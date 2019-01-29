import pysam
import cycler

import numpy as np
import pybedtools as pb
import matplotlib as mpl
import matplotlib.pyplot as plt

from collections import namedtuple, Counter, defaultdict

class Region(namedtuple('Region', ['reference_name', 'reference_start', 'reference_end'])):
    def __new__(cls, reference_name, reference_start, reference_end):
        assert reference_start <= reference_end, 'Invalid genomic region %s:%d-%d.' % (reference_name, reference_start, reference_end)
        return super(Region, cls).__new__(cls, reference_name, reference_start, reference_end)

    def __str__(self):
        return '%s:%d-%d' % (self.reference_name, self.reference_start, self.reference_end)

def get_bam_handle(bam_path):
    """Returns `pysam.AlignmentFile` object for BAM file."""
    return pysam.AlignmentFile(bam_path)


def get_binary_string_representation(read):
    """Returns a binary representation string of methylation state of a read.
    Unmethylated state = 0, methylated state = 1.
    """
    return ''.join([['0', '1'][c == 'Z'] for c in read.get_tag('XM') if c in ['z', 'Z']])

def is_fully_methylated_or_unmethylated(pattern):
    return len(set(str(pattern))) == 1

def pattern_counters_from_met(fp):
    """Parse MET file and generate header, and counter of patterns."""
    header, binary_patterns = None, []
    with open(fp) as inFile:
        for line in inFile.readlines():
            if line.startswith('>'):
                if header is not None:
                    pattern_counter = Counter(binary_patterns)
                    yield header[1:], pattern_counter
                    header, binary_patterns = line.strip(), []

                else:
                    header = line.strip()

            else:
                binary_patterns.append(line.strip())

        pattern_counter = Counter(binary_patterns)
        yield header[1:], pattern_counter

def random_pattern(l, n):
    """Generate `n` random methylation pattern with length `l`."""
    s = set()
    ps = []
    p = np.random.randint(2, size=l)
    s.add(''.join(map(str, p)))
    ps.append(p)
    while True:
        if len(ps) == n:
            break
        p = list(np.random.randint(2, size=l))
        if ''.join(map(str, p)) not in s:
            s.add(''.join(map(str, p)))
            ps.append(p)

    return np.array(list(ps)).flatten()

def jaccard_similarity(s1, s2):
    """Compute Jaccard similarity of two sets."""
    s1 = set(s1)
    s2 = set(s2)
    return len(s1 & s2) / len(s1 | s2)

def extract_chromosome(header):
    """Extract chromosome information for epilocus header."""
    return header.split(';')[0].split(':')[0]

def merge_two_headers(h1, h2):
    """Union two epiloci headers."""
    chromosome = h1.split(';')[0].split(':')[0]

    start1 = int(h1.split(';')[0].split(':')[1])
    start2 = int(h2.split(';')[0].split(':')[1])

    end1 = int(h1.split(';')[-1].split(':')[1])
    end2 = int(h2.split(';')[-1].split(':')[1])

    return '%s:%d;%s:%d' % (chromosome, min(start1, start2), chromosome, max(end1, end2))

def get_common_headers_by_jaccard_similarity(headers1, headers2, cutoff=0.5):
    """Given two lists of epiloci headers, determine which epiloci are common."""
    common_h1, common_h2, merged_headers = [], [], []
    chr_header1_dict = defaultdict(list)
    chr_header2_dict = defaultdict(list)
    
    for h1 in headers1:
        chr_header1_dict[extract_chromosome(h1)].append(h1)
        
    for h2 in headers2:
        chr_header2_dict[extract_chromosome(h2)].append(h2)
        
    for chromosome in chr_header1_dict.keys():
        for h1 in chr_header1_dict[chromosome]:
            for h2 in chr_header2_dict[chromosome]:
                cpgs1 = h1.split(';')
                cpgs2 = h2.split(';')

                if jaccard_similarity(cpgs1, cpgs2) >= cutoff:
                    common_h1.append(h1)
                    common_h2.append(h2)
                    merged_headers.append(merge_two_headers(h1, h2))
    
    return common_h1, common_h2, merged_headers

def preset_rc(scale=1, font_family=None):
    """Set visually plausible matplotlib rc file."""
    plt.rc('axes', linewidth=1.33, labelsize=14)
    plt.rc('xtick', labelsize=8 * scale)
    plt.rc('ytick', labelsize=8 * scale)

    plt.rc('xtick', bottom=True)
    plt.rc('xtick.major', size=5 * scale, width=1.33)
    plt.rc('xtick.minor', size=5 * scale, width=1.33)

    plt.rc('ytick', left=True)
    plt.rc('ytick.major', size=5 * scale, width=1.33)
    plt.rc('ytick.minor', size=5 * scale, width=1.33)

    plt.rc('legend', fontsize=7 * scale)
    plt.rc('grid', color='grey', linewidth=0.5, alpha=0.33)

    if font_family:
        plt.rc('font', family=font_family)

    color_palette = [
        '#005AC8',
        '#AA0A3C',
        '#0AB45A',
        '#FA7850',
        '#8214A0',
        '#FA78FA',
        '#A0FA82',
        '#006E82',
        '#00A0FA',
        '#14D2DC',
        '#F0F032',
        '#FAE6BE',
    ]

    mpl.rcParams['axes.prop_cycle'] = cycler.cycler(color=color_palette)

def parse_result_line(line):
    """Parse each line in the PRISM result file."""
    fields = line.strip().split()

    header = fields[0]
    cluster = int(fields[1])
    subclone = int(fields[2])
    depths = np.array(list(map(int, fields[3].split(','))))
    counts = np.array(list(map(int, fields[4].split(','))))
    fingerprint_fractions = list(map(float, fields[5].split(',')))

    return header, cluster, subclone, depths, counts, fingerprint_fractions

def prepare_header_bed(headers):
    """Converts a list of epiloci headers to a BedTool object."""
    bed_strings = []
    for header in headers:
        chrom = header.split(':')[0]
        # BED file uses 0-based exclusive coordinate system.
        start = int(header.split(';')[0].split(':')[1]) - 1  
        end = int(header.split(';')[-1].split(':')[1])

        bed_strings.append('%s\t%d\t%d' % (chrom, start, end))
    
    return pb.BedTool('\n'.join(bed_strings), from_string=True)