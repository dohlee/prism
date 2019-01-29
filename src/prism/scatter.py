import cycler

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import prism.util as util

def scatter_1d(subclone_assignments, num_subclones, depths, fingerprint_fractions, width, height, dpi):
    """Generate annotated scatterplot for one-dimensional PRISM analysis.

    :param list subclone_assignments: List containing subclone assignment status for each fingerprint epilocus.
    :param int num_subclones: Number of subclones.
    :param list depths: List of depths (pattern counts).
    :param list fingerprint_fractions: List of fingerprint fractions.
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
            label='S%d' % (subclone_index + 1),
            s=60, edgecolors='black', linewidth=0.3,
        )
    

def scatter_2d(subclone_assignments, num_subclones, fingerprint_fractions, width, height, dpi):
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
            label='S%d' % (subclone_index + 1),
            s=60, edgecolors='black', linewidth=0.3,
        )

def run(input_fp, output_fp, dpi=400, width=4, height=4, scale=1, font_family=None):
    util.preset_rc(scale=scale, font_family=font_family)

    subclone_assignments = []
    depths = []
    fingerprint_fractions = []

    with open(input_fp) as inFile:
        inFile.readline()
        for line in inFile.readlines():
            header, cluster, subclone, d, c, ffs = util.parse_result_line(line)

            subclone_assignments.append(subclone)
            depths.append(d)
            fingerprint_fractions.append(ffs)
        
    subclone_assignments = np.array(subclone_assignments)
    fingerprint_fractions = np.array(fingerprint_fractions)
    depths = np.array(depths)
    num_subclones = len(set(subclone_assignments))
    n_dim = len(fingerprint_fractions[0])
    
    if n_dim == 1:
        scatter_1d(subclone_assignments, num_subclones, depths, fingerprint_fractions, width, height, dpi)
    elif n_dim == 2:
        scatter_2d(subclone_assignments, num_subclones, fingerprint_fractions, width, height, dpi)

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_fp, dpi=dpi)