import cycler
import sys
import cleanlog

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pybedtools as pb
import prism.util as util

logger = cleanlog.ColoredLogger('annotate')

def scatter_1d(subclone_assignments, num_subclones, depths, fingerprint_fractions, annotation_names, annotation_mask, width, height, dpi):
    """Generate annotated scatterplot for one-dimensional PRISM analysis.

    :param list subclone_assignments: List containing subclone assignment status for each fingerprint epilocus.
    :param int num_subclones: Number of subclones.
    :param list depths: List of depths (pattern counts).
    :param list fingerprint_fractions: List of fingerprint fractions.
    :param list annotation_names: List of the names of annotations.
    :param list annotation_mask: List containing annotation status for each fingerprint epilocus.
    :param float width: Figure width.
    :param float height: Figure height.
    :param int dpi: Figure DPI.
    """
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel('Fraction of fingerprint')
    ax.set_ylabel('Depth')

    for subclone_index in range(num_subclones):
        mask = (subclone_assignments == subclone_index)
        ax.scatter(
            fingerprint_fractions[mask, 0],
            depths[mask, 0],
            s=60, alpha=0.2, linewidth=0,
        )
    
    plt.gca().set_prop_cycle(None)

    for name, mask in zip(annotation_names, annotation_mask):
        ax.scatter(fingerprint_fractions[mask, 0], depths[mask, 0], label=name, marker='x', s=60, lw=1.33)
    

def scatter_2d(subclone_assignments, num_subclones, fingerprint_fractions, annotation_names, annotation_mask, width, height, dpi):
    """Generate annotated scatterplot for two-dimensional PRISM analysis.

    :param list subclone_assignments: List containing subclone assignment status for each fingerprint epilocus.
    :param int num_subclones: Number of subclones.
    :param list fingerprint_fractions: List of fingerprint fractions.
    :param list annotation_names: List of the names of annotations.
    :param list annotation_mask: List containing annotation status for each fingerprint epilocus.
    :param float width: Figure width.
    :param float height: Figure height.
    :param int dpi: Figure DPI.
    """
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel('FF1')
    ax.set_ylabel('FF2')

    for subclone_index in range(num_subclones):
        mask = (subclone_assignments == subclone_index)
        ax.scatter(
            fingerprint_fractions[mask, 0],
            fingerprint_fractions[mask, 1],
            s=60, alpha=0.2, linewidth=0,
        )

    plt.gca().set_prop_cycle(None)
    
    for name, mask in zip(annotation_names, annotation_mask):
        ax.scatter(fingerprint_fractions[mask, 0], fingerprint_fractions[mask, 1], label=name, marker='x', s=60, lw=1.33)

def run(input_fp, output_fp, bed_fps, annotation_names, output_figure_fp=None, dpi=400, width=4, height=4, scale=1, font_family=None):
    if len(bed_fps) != len(annotation_names):
        logger.error('The number of bed files and their names should match. (Given %d bed files and %d names.)' % (len(bed_fps), len(annotation_names)))
        
    util.preset_rc(scale=scale, font_family=font_family)

    headers = []
    subclone_assignments = []
    depths = []
    fingerprint_fractions = []
    with open(input_fp) as inFile:
        inFile.readline()
        for line in inFile.readlines():
            header, cluster, subclone, d, c, ffs = util.parse_result_line(line)

            headers.append(header)
            subclone_assignments.append(subclone)
            depths.append(d)
            fingerprint_fractions.append(ffs)
    
    header_bed = util.prepare_header_bed(headers)

    num_annotations = len(annotation_names)
    beds = [pb.BedTool(fp) for fp in bed_fps]
    annotation_mask = np.array([[False] * len(headers) for _ in range(num_annotations)])

    for annotation_index, bed in enumerate(beds):
        for epiloci_index, interval in enumerate(header_bed.intersect(bed, c=True)):
            overlap_counts = int(interval.fields[-1])
            # If there's overlap with current bed file with annotations, mark as annotated.
            annotation_mask[annotation_index][epiloci_index] = (overlap_counts != 0)
    
    with open(output_fp, 'w') as outFile:
        with open(input_fp) as inFile:
            inFile.readline()
            for epiloci_index, line in enumerate(inFile.readlines()):
                for annotation_index in range(num_annotations):
                    annotation = ','.join([annotation_names[i] for i in range(num_annotations) if annotation_mask[i][epiloci_index] == True])
                    print(line.strip() + '\t%s' % annotation, file=outFile)
    
    if output_figure_fp is not None:
        subclone_assignments = np.array(subclone_assignments)
        fingerprint_fractions = np.array(fingerprint_fractions)
        depths = np.array(depths)
        num_subclones = len(set(subclone_assignments))
        n_dim = len(fingerprint_fractions[0])
        
        if n_dim == 1:
            scatter_1d(subclone_assignments, num_subclones, depths, fingerprint_fractions, annotation_names, annotation_mask, width, height, dpi)
        elif n_dim == 2:
            scatter_2d(subclone_assignments, num_subclones, fingerprint_fractions, annotation_names, annotation_mask, width, height, dpi)

        plt.legend()
        plt.tight_layout()
        plt.savefig(output_figure_fp, dpi=dpi)