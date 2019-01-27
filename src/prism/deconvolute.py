from prism.mixture import BetaBinomialMixture
from collections import Counter

import numpy as np
import prism.util as util
import cleanlog

logger = cleanlog.ColoredLogger('deconvolute')

def postfiltered(pattern_counter_generator, full_pattern_proportion):
    """"""
    for header, pattern_counter in pattern_counter_generator:
        if len(pattern_counter) == 1:
            continue

        num_CpGs = len(header.split(';'))
        counts = pattern_counter.values()
        pattern1, pattern2 = pattern_counter.most_common(2)[0][0], pattern_counter.most_common(2)[1][0]

        total = sum(counts)
        if (
            util.is_fully_methylated_or_unmethylated(pattern1) and
            util.is_fully_methylated_or_unmethylated(pattern2) and
            pattern_counter[pattern1] + pattern_counter[pattern2] > total * full_pattern_proportion
        ):
            yield header, pattern_counter


def methylated_pattern(p1, p2):
    return [p1, p2][sum(c == '1' for c in str(p1)) < sum(c == '1' for c in str(p2))]


def parse_met_file(fp, full_pattern_proportion=0.8):
    depths, counts, headers = [], [], []
    for header, pattern_counter in postfiltered(util.pattern_counters_from_met(fp), full_pattern_proportion):
        p1, p2 = pattern_counter.most_common(2)[0][0], pattern_counter.most_common(2)[1][0]
        depth = pattern_counter[p1] + pattern_counter[p2]
        count = pattern_counter[methylated_pattern(p1, p2)]

        depths.append(depth)
        counts.append(count)
        headers.append(header)
    
    return np.array(depths), np.array(counts), np.array(headers)


def common_intersection(headers_list):
    common_headers = set(headers_list[0])
    for i in range(1, len(headers_list)):
        common_headers = common_headers & set(headers_list[i])
    common_headers = list(common_headers)

    return common_headers

def merge_met_files(met_files, full_pattern_proportion, intersection_method):
    n_samples = len(met_files)

    depths_list, counts_list, headers_list = [], [], []
    for met_file in met_files:
        logger.debug('Parsing %s.' % met_file)
        depths, counts, headers = parse_met_file(met_file, full_pattern_proportion=full_pattern_proportion)
        depths_list.append(depths)
        counts_list.append(counts)
        headers_list.append(headers)

    header_depth_dicts, header_count_dicts = [], []
    for depths, counts, headers in zip(depths_list, counts_list, headers_list):
        header_depth_dict, header_count_dict = dict(), dict()

        for depth, count, header in zip(depths, counts, headers):
            header_depth_dict[header] = depth
            header_count_dict[header] = count

        header_depth_dicts.append(header_depth_dict)
        header_count_dicts.append(header_count_dict)

    if intersection_method == 'common':
        common_headers = common_intersection(headers_list)
        merged_depths = np.array([
            [header_depth_dicts[i][h] for i in range(n_samples)] for h in common_headers
        ])

        merged_counts = np.array([
            [header_count_dicts[i][h] for i in range(n_samples)] for h in common_headers
        ])

    elif intersection_method == 'jaccard':
        assert len(headers_list) == 2, 'Extracting common headers by jaccard similarity is applicable only for two samples.'
        common_h1, common_h2, common_headers = util.get_common_headers_by_jaccard_similarity(headers_list[0], headers_list[1])
        merged_depths = np.array([
            [header_depth_dicts[0][h1], header_depth_dicts[1][h2]] for h1, h2 in zip(common_h1, common_h2)
        ])

        merged_counts = np.array([
            [header_count_dicts[0][h1], header_count_dicts[1][h2]] for h1, h2 in zip(common_h1, common_h2)
        ])


    return merged_depths, merged_counts, common_headers

def mark_outlier_clusters(model, outlier_dispersion_cutoff):
    return np.array([any(d > 0.2) for d in model.get_dispersions()])

def merge_subclones(subclones, cluster_a, cluster_b):
    """Merge two clusters containing cluster a and cluster b.
    If a and b are already in the same subclone, just return the subclones unchanged."""
    merged_subclones = []

    for subclone in subclones:
        # Cluster a and b is already in the same subclone.
        if cluster_a in subclone and cluster_b in subclone:
            return subclones

        if cluster_a in subclone:
            subclone_containing_a = subclone
        elif cluster_b in subclone:
            subclone_containing_b = subclone
        else:
            merged_subclones.append(subclone)

    # Append merged subclone to the result.
    merged_subclones.append(subclone_containing_a | subclone_containing_b)

    return merged_subclones

def identify_subclone(model, merge_cutoff, outlier_cluster_mask):
    midpoint_distance = lambda a, b: np.sqrt(np.square((a + b) / 2 - 0.5).sum())

    cluster_means = model.get_means() 
    num_clusters = model.get_n_components()
    subclones = [{i} for i in range(num_clusters)]

    # For 1-dimensional analysis, skip merging clusters.
    if model.get_n_dimensions() == 1:
        return subclones, outlier_cluster_mask[:]
    
    # Examine all pairs of clusters, and merge if their midpoint is close to (0.5, ..., 0.5).
    for cluster_a in range(num_clusters):
        for cluster_b in range(cluster_a + 1, num_clusters):
            if midpoint_distance(cluster_means[cluster_a], cluster_means[cluster_b]) < merge_cutoff:
                subclones = merge_subclones(subclones, cluster_a, cluster_b)
    
    # Mark a subclone as outlier subclone if it contains outlier cluster.
    outlier_subclone_mask = []
    for subclone in subclones:
        is_outlier_subclone = any([outlier_cluster_mask[cluster] for cluster in subclone])
        outlier_subclone_mask.append(is_outlier_subclone)
    
    # Reorder subclone so that non-outlier subclones have indices starting from 0.
    final_subclones = []
    final_outlier_subclone_mask = []
    for subclone_index, subclone in enumerate(subclones):
        if not outlier_subclone_mask[subclone_index]:
            final_subclones.append(subclone)
            final_outlier_subclone_mask.append(False)

    for subclone_index, subclone in enumerate(subclones):
        if outlier_subclone_mask[subclone_index]:
            final_subclones.append(subclone)
            final_outlier_subclone_mask.append(True)

    return final_subclones, final_outlier_subclone_mask

def posthoc_process(model, merge_cutoff, outlier_dispersion_cutoff):
    """"""
    outlier_cluster_mask = mark_outlier_clusters(model, outlier_dispersion_cutoff)
    subclones, outlier_subclone_mask = identify_subclone(model, merge_cutoff, outlier_cluster_mask)
    return subclones, outlier_subclone_mask

def get_subclone_assignment(subclones, assignment, outlier_subclone_mask):
    for subclone_index, subclone in enumerate(subclones):
        if assignment in subclone:

            # If the assigned subclone is outlier, return -1.
            if outlier_subclone_mask[subclone_index]:
                return -1
            # Else, return the index of subclone.
            else:
                return subclone_index
    
def run(input_fps, full_pattern_proportion=0.8, merge_cutoff=0.05, outlier_dispersion_cutoff=0.2, num_max_cluster=15, seed=12345, intersection_method='common', verbose=False, output_fp=None):
    if verbose:
        logger.setLevel(cleanlog.DEBUG)
    
    merged_depths, merged_counts, common_headers = merge_met_files(input_fps, full_pattern_proportion=full_pattern_proportion, intersection_method=intersection_method)
    logger.debug('Total %d fingerprint epiloci will be used for deconvolution.' % len(common_headers))
    
    models = []
    for n_components in range(1, num_max_cluster + 1):
        logger.debug('Fitting beta-binomial mixture model with %d clusters.' % n_components)
        bbmm = BetaBinomialMixture(n_components=n_components, seed=seed)
        bbmm.fit(merged_depths, merged_counts, common_headers)

        models.append(bbmm)

    selected_model = min(models, key=lambda m: m.bic())
    logger.debug('The best model had %d clusters.' % selected_model.get_n_components())

    subclones, outlier_subclone_mask = posthoc_process(selected_model, merge_cutoff, outlier_dispersion_cutoff)
    cluster_assignments = selected_model.predict_proba(merged_depths, merged_counts).argmax(axis=0)
    subclone_assignments = [get_subclone_assignment(subclones, assignment, outlier_subclone_mask) for assignment in cluster_assignments]

    logger.debug('The clusters were processed, and resulted in %d subclones.' % len(subclones))

    if output_fp is None:
        return selected_model
    
    headers = [
        'epilocus',
        'cluster',
        'subclone',
        'depths',
        'fingerprint_counts'
        'fingerprint_fractions',
    ]
    with open(output_fp, 'w') as outFile:
        print('\t'.join(headers), file=outFile)
        for header, cluster, subclone, depths, counts in zip(common_headers, cluster_assignments, subclone_assignments, merged_depths, merged_counts):
            ffs = ",".join(map(lambda t: str(t[1] / t[0]), zip(depths, counts)))
            print(f'{header}\t{cluster}\t{subclone}\t{",".join(map(str, depths))}\t{",".join(map(str, counts))}\t{ffs}', file=outFile)
    return selected_model