import cleanlog
import prism.util as util
import prism.proofreading as proofreading
import multiprocessing as mp
import numpy as np

from collections import Counter
from io import StringIO

logger = cleanlog.ColoredLogger('preprocess')

def prefiltered(pattern_counter_generator, full_pattern_proportion, no_prefilter):
    """"""
    if no_prefilter:
        for header, pattern_counter in pattern_counter_generator:
            yield header, pattern_counter
    
    else:
        for header, pattern_counter in pattern_counter_generator:
            counts = pattern_counter.values()
            pattern1, pattern2 = pattern_counter.most_common(2)[0][0], pattern_counter.most_common(2)[1][0]

            total = sum(counts)
            if (
                util.is_fully_methylated_or_unmethylated(pattern1) and
                util.is_fully_methylated_or_unmethylated(pattern2) and
                pattern_counter[pattern1] + pattern_counter[pattern2] > total * full_pattern_proportion
            ):
                yield header, pattern_counter

def param_generator(fp, no_prefilter, full_pattern_proportion, error, bisulfite_conversion_rate, processivity, recruitment_efficiency):
    """"""
    hmm_params = {
        'e_m': error,
        'e_b': 1 - bisulfite_conversion_rate,
        'p': processivity,
        'q': recruitment_efficiency,
    }

    for header, pattern_counter in prefiltered(util.pattern_counters_from_met(fp), full_pattern_proportion, no_prefilter):
        binary_patterns = []
        for pattern, count in pattern_counter.items():
            binary_patterns += [pattern for _ in range(count)]
        binary_patterns = np.genfromtxt(StringIO('\n'.join(binary_patterns)), dtype=np.int8, delimiter=[1] * len(binary_patterns[0]))

        yield header, binary_patterns, hmm_params


def run(
    input_fp,
    output_fp,
    no_prefilter,
    full_pattern_proportion,
    error,
    bisulfite_conversion_rate,
    processivity,
    recruitment_efficiency,
    threads,
    seed,
    verbose,
):
    if verbose:
        logger.setLevel(cleanlog.DEBUG)

    with mp.Pool(processes=threads) as p:
        results = list(p.imap(  
            proofreading.proofread,
            param_generator(input_fp, no_prefilter, full_pattern_proportion, error, bisulfite_conversion_rate, processivity, recruitment_efficiency),
            chunksize=100,
        ))

    with open(output_fp, 'w') as outFile:
        for header, proofread_patterns in results:
            stringified_patterns = [''.join(map(str, p)) for p in proofread_patterns]

            print('>' + header, file=outFile)
            print('\n'.join(stringified_patterns), file=outFile)