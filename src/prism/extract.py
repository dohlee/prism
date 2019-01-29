import argparse
import pysam
import cleanlog
import re
import sys
import os

import prism.util as util

from collections import namedtuple, defaultdict, Counter
from prism.util import Region

logger = cleanlog.ColoredLogger('extract')

def extend_region(region, read):
    """Returns the union of genomic regions covered by given region and read.

    :param Region region: Genomic region to be extended.
    :param AlignedSegment read: Read used for the extension.

    :returns: Extended genomic region, which is the union of region and read.
    """
    assert region.reference_name == read.reference_name, \
           'Unable to extend region. Contigs are incompatible: %s vs %s' % (region.reference_name, read_reference_name)

    new_start = min(region.reference_start, read.reference_start)
    new_end = max(region.reference_end, read.reference_end)

    return Region(region.reference_name, new_start, new_end)


def has_overlap(region, read):
    """Returns True if given region and read overlap, otherwise False.

    :param Region region: Genomic region to test overlap. 
    :param AlignedSegment read: Read to test overlap.

    :returns: True if given region and read overlap, otherwise False.
    """
    if region.reference_name != read.reference_name:
        return False

    a0, a1 = region.reference_start, region.reference_end
    b0, b1 = read.reference_start, read.reference_end

    return a0 <= b1 and b0 <= a1

def character_indices(string, characters):
    """Returns a sorted list of indices of characters in the string.

    :param string string: A string where the characters are expected to appear.
    :param list characters: List of characters to return indices.

    :returns: A list of indices of string where the characters appear.
    """
    i = []
    for character in characters:
        i.extend([m.start() for m in re.finditer(character, string)])

    return list(sorted(i))

def get_cpg_coordinates(read, paired=False):
    """Returns genomic coordinates of cytosines of CpGs in the read.

    :param AlignedSegment read: Read to examine for the location of CpGs.
    :param bool paired: True if the sequencing library is paired, otherwise False.

    :returns: Absolute genomic coordinates of CpGs that appear in the read.
    """
    start = read.reference_start
    methylation_string = read.get_tag('XM')

    if not paired:
        # We should care whether the read is mapped to forward strand or reverse strand
        # because the coordinates of two cytosines in complementary CpG pair will differ by
        # a single base pair.
        # Also note that we use 1-based, inclusive coordinate system in MET files.
        if not read.is_reverse:
            # +1 to convert to 1-based system.
            return tuple(start + i + 1 for i in character_indices(methylation_string, ['z', 'Z']))
        else:
            # -1 then +1, so apparently no differences are made in the coordinate.
            return tuple(start + i for i in character_indices(methylation_string, ['z', 'Z']))

    else:
        # Flag 99 and 147 indicates reverse strand. 
        if read.flag in [99, 147]:
            return tuple(start + i + 1 for i in character_indices(methylation_string, ['z', 'Z']))
        # Forward strand.
        else:
            return tuple(start + i for i in character_indices(methylation_string, ['z', 'Z']))


def discard_noninformative_reads(cpgs_reads_dict, depth_cutoff, num_cpg_cutoff):
    """Discards noninformative reads with low depth or few CpGs.

    :param dict cpgs_reads_dict: Dictionary mapping a set of CpG coordinates to a group of reads.
    :param int depth_cutoff: Minimum depth of a group of reads to be retained.
    :param int num_cpg_cutoff: Mininum number of CpGs appearing in the reads to be retained.

    :returns: Retained (informative) CpG-Read group mapping dictionary.
    """
    retained = defaultdict(list)

    for cpgs, reads in cpgs_reads_dict.items():
        # Discard if (1) depth or (2) the number of CpG is not sufficient.
        if len(reads) < depth_cutoff or len(cpgs) < num_cpg_cutoff:
            continue

        pattern_counter = Counter([util.get_binary_string_representation(read) for read in reads])
        # If there is only one kind of pattern, discard it.
        if len(pattern_counter) == 1:
            continue
        
        retained[cpgs] = reads
        # p1, p2 = pattern_counter.most_common(2)[0][0], pattern_counter.most_common(2)[1][0]

        # # (3) the two most frequent pattern should be fully methylated/unmethylated.
        # if len(set(str(p1))) != 1 or len(set(str(p2))) != 1:
        #     continue

        # # (4) The two most frequent patterns should account for more than 80% of reads (by default).
        # if pattern_counter[p1] + pattern_counter[p2] > sum(pattern_counter.values()) * full_pattern_proportion:
        #     retained[cpgs] = reads

    return retained


def save_met_file(handle, output_path, depth_cutoff, num_cpg_cutoff, paired, contigs):
    """Read through the bam file and identify group of reads mapped to the same genomic region.
    Group of reads with insufficient depth and insufficient number of CpGs will be discarded.
    NOTE: BAM file should be sorted to guarantee proper generation of MET file.

    :param AlignmentFile handle: pysam AlignmentFile object for BAM file.
    :param string output_path: Output MET file path.
    :param int depth_cutoff: Minimum depth for a group of reads to be retained.
    :param num_cpg_cutoff: Minimum number of CpGs appearing in the reads to be retained.
    :param bool paired: True if the sequencing library is paired, otherwise False.
    :param list contigs: A list of contigs to examine.
    """
    out = open(output_path, 'w')

    # We are going to sweep through the genome.
    current_region = None

    # To gather reads harboring the same set of CpGs,
    # tuples of CpG coordinates will be mapped to a list of corresponding reads.
    reads = []
    cpgs_reads_dict = defaultdict(list)

    for contig in contigs:
        logger.debug('Processing contig %s.' % contig)
        count = 0

        for read in handle.fetch(contig):
            # The first iteration.
            if current_region is None:
                reads.append(read)
                current_region = Region(read.reference_name, read.reference_start, read.reference_end)
                cpgs_reads_dict[get_cpg_coordinates(read, paired)].append(read)
                continue

            # If examining current region is not finished,
            # extend the current region state with current read,
            # and map read information to dictionary with the set of coordinates of CpGs it harbors as a key.
            if has_overlap(current_region, read):
                reads.append(read)
                current_region = extend_region(current_region, read)
                cpgs_reads_dict[get_cpg_coordinates(read, paired)].append(read)

            # If examining current region is finished (i.e., met reads mapped at distant region),
            # discard noninformative read groups from the dictionary and
            # write read groups to the output file.
            else:
                retained = discard_noninformative_reads(cpgs_reads_dict, depth_cutoff, num_cpg_cutoff)
                
                for cpgs, reads in retained.items():
                    print('>' + ';'.join([current_region.reference_name + ':' + str(cpg) for cpg in cpgs]), file=out)
                    print('\n'.join([util.get_binary_string_representation(read) for read in reads]), file=out)
                    count += 1

                # Refresh the current information to examine the next region.
                reads = [read]
                current_region = Region(read.reference_name, read.reference_start, read.reference_end)
                cpgs_reads_dict = defaultdict(list)
                cpgs_reads_dict[get_cpg_coordinates(read, paired)].append(read)

        logger.debug('Found %d loci in contig %s' % (count, contig))
    out.close()

def run(
    input_fp,
    output_fp,
    depth_cutoff,
    num_cpg_cutoff,
    prepend_chr,
    paired,
    verbose,
):
    if verbose:
        logger.setLevel(cleanlog.DEBUG)

    if prepend_chr:
        contigs = list(map(lambda x: 'chr' + str(x), range(1, 23)))
    else:
        contigs = list(map(str, range(1, 23)))
    
    handle = util.get_bam_handle(input_fp)

    save_met_file(
        handle=handle,
        output_path=output_fp,
        depth_cutoff=depth_cutoff,
        num_cpg_cutoff=num_cpg_cutoff,
        paired=paired,
        contigs=contigs,
    )
