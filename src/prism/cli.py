"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mprism` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``prism.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``prism.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse

import prism.extract
import prism.deconvolute
import prism.preprocess
import prism.scatter
import prism.annotate

# Create the top-level parser.
parser = argparse.ArgumentParser(description='Decompose a bulk sample into putative epigenetic subclones.')
subparsers = parser.add_subparsers(description='Commands', dest='subcommand')

def subcommand(args=[], parent=subparsers):
    def decorator(func):
        parser = parent.add_parser(func.__name__, description=func.__doc__)
        for arg in args:
            parser.add_argument(*arg[0], **arg[1])
        parser.set_defaults(func=func)
    return decorator

def argument(*name_or_flags, **kwargs):
    return ([*name_or_flags], kwargs)

def argument_list(short, long, default=None, required=False, **kwargs):
    return argument(short, long, nargs='+', default=default, required=required, **kwargs)

def argument_string(short, long, default=None, required=False, **kwargs):
    return argument(short, long, default=default, required=required, **kwargs)

def argument_int(short, long, default=None, required=False, **kwargs):
    return argument(short, long, type=int, default=default, required=required, **kwargs)

def argument_bool(short, long, default=None, required=False, **kwargs):
    return argument(short, long, action='store_true', default=default, required=required, **kwargs)

def argument_float(short, long, default=None, required=False, **kwargs):
    return argument(short, long, type=float, default=default, required=required, **kwargs)

@subcommand([
    argument_string('-i', '--input', help='Input BAM file.', required=True),
    argument_string('-o', '--output', help='Output MET file, which stores the epiloci information.', required=True),
    argument_int(
        '-d', 
        '--depth-cutoff',
        default=20,
        help='Depth cutoff for epiloci. Every epiloci should have at least this much reads mapped to it.'
        'Note that filtering by depth and CpG counts is always done even --no-prefilter option is specified.'
    ),
    argument_int(
        '-c',
        '--num-cpg-cutoff',
        default=4,
        help='CpG count cutoff for epiloci. Reads mapped to every epiloci should harbor at least this much CpGs.'
        'Note that filtering by depth and CpG counts is always done even --no-prefilter option is specified.'
    ),
    argument_bool('-u', '--prepend-chr', default=False, help='Prepend "chr" to the contigs. Useful when using UCSC reference genome.'),
    argument_bool('-x', '--paired', default=False, help='Specify if the reads are paired.'),
    argument_bool('-v', '--verbose', default=False, help='Increase verbosity.'),
])
def extract(args):
    """Run epiloci extraction step of PRISM workflow.
    """
    prism.extract.run(
        input_fp=args.input,
        output_fp=args.output,
        depth_cutoff=args.depth_cutoff,
        num_cpg_cutoff=args.num_cpg_cutoff,
        prepend_chr=args.prepend_chr,
        paired=args.paired,
        verbose=args.verbose,
    )

@subcommand([
    argument_string('-i', '--input', help='Input MET file, which stores the raw epiloci information.', required=True),
    argument_string('-o', '--output', help='Output corrected MET file, of which patterns are in sillico proofread.', required=True),
    argument_bool(
        '-f',
        '--no-prefilter',
        default=False, 
        help='Do not apply pre-filter to reduce the number of epiloci passed to in sillico proofreading.'
        'This may cause a significant increase in running time of processbam step.'
    ),
    argument_float(
        '-r',
        '--full-pattern-proportion',
        default=0.5,
        help='Cutoff for the proportion of fully methylated/unmethylated patterns in the pre-filtering step.'
        'This is ignored when --no-prefilter option is specified. (Default: 0.5)'
    ),
    argument_float('-e', '--error', default=0.01, help='Expected error rate of methylation maintenance.'),
    argument_float('-b', '--bisulfite-conversion-rate', default=0.99, help='Bisulfite conversion rate of the sample.'),
    argument_float('-p', '--processivity', default=0.96, help='Expected processivity of DNMT1 enzyme.'),
    argument_float('-q', '--recruitment-efficiency', default=0.96, help='Expected recruitment efficiency of DNMT1 enzyme.'),
    argument_int('-t', '--threads', default=1, help='Number of threads to use.'),
    argument_int('-s', '--seed', default=12345, help='It is recommended to set random seed for reproducibility.'),
    argument_bool('-v', '--verbose', default=False, help='Increase verbosity.'),
])
def preprocess(args):
    """Run preprocessing step of PRISM workflow.
    """
    prism.preprocess.run(
        input_fp=args.input,
        output_fp=args.output,
        no_prefilter=args.no_prefilter,
        full_pattern_proportion=args.full_pattern_proportion,
        error=args.error,
        bisulfite_conversion_rate=args.bisulfite_conversion_rate,
        processivity=args.processivity,
        recruitment_efficiency=args.recruitment_efficiency,
        threads=args.threads,
        seed=args.seed,
        verbose=args.verbose,
    )

@subcommand([
    argument_list('-i', '--input', help='A list of input MET files.', required=True),
    argument_string('-o', '--output', help='Output PRISM profile result.', required=True),
    argument_int('-m', '--num-max-cluster', default=15, help='Maximum number of clusters to examine.'),
    argument_float(
        '-r', 
        '--full-pattern-proportion',
        default=0.8,
        help='Cutoff for the proportion of fully methylated/unmethylated patterns in the filtering step.'
    ),
    argument_float(
        '-a',
        '--merge-cutoff',
        default=0.05,
        help='Cutoff for the distance from the midpoint of two clusters and (0.5, ..., 0.5) to treat clusters as single subclone.'
    ),
    argument_float(
        '-d',
        '--outlier-dispersion-cutoff',
        default=0.2,
        help='Cutoff for the dispersion to call outlier cluster.'
    ),
    argument_string(
        '-t',
        '--intersection-method',
        default='common',
        help='Method for calling common epiloci when analyzing more than one samples jointly.'
    ),
    argument_string('-c', '--copynumber', help='Called copy number file.', default=None),
    argument_string('-p', '--cn-prior', help='Prior for probability of methylation patterns in copy-number-gained segments.', default='random'),
    argument_int('-s', '--seed', default=12345, help='It is recommended to set random seed for reproducibility.'),
    argument_bool('-v', '--verbose', default=False, help='Increase verbosity.'),
])
def deconvolute(args):
    """Solve beta-binomial mixture of FF to detect subclones.
    """
    prism.deconvolute.run(
        input_fps=args.input,
        output_fp=args.output,
        full_pattern_proportion=args.full_pattern_proportion,
        merge_cutoff=args.merge_cutoff,
        outlier_dispersion_cutoff=args.outlier_dispersion_cutoff,
        intersection_method=args.intersection_method,
        copynumber=args.copynumber,
        cn_prior=args.cn_prior,
        num_max_cluster=args.num_max_cluster,
        seed=args.seed,
        verbose=args.verbose,
    )

@subcommand([
    argument_string('-i', '--input', help='Input PRISM deconvolution result file.', required=True),
    argument_string('-o', '--output', help='Output scatterplot result file.', required=True),
    argument_int('-x', '--width', help='Figure width in inches.', default=4),
    argument_int('-y', '--height', help='Figure height in inches.', default=4),
    argument_float('-s', '--scale', help='Figure scale.', default=1),
    argument_string('-f', '--font-family', help='Font faily used for the plot.', default=None),
])
def scatter(args):
    """Generate scatterplot showing the fractions of fingerprint patterns.
    """
    prism.scatter.run(
        input_fp=args.input,
        output_fp=args.output,
        width=args.width,
        height=args.height,
        scale=args.scale,
        font_family=args.font_family,
    )
    
@subcommand([
    argument_string('-i', '--input', help='Input PRISM deconvolution result file.', required=True),
    argument_string('-o', '--output', help='Output annotated result file.', required=True),
    argument_list('-b', '--beds', help='BED file with intervals to annotate.', required=True),
    argument_list('-n', '--annotation-names', help='Names describing the annotations.', required=True),
    argument_string('-g', '--figure', help='Output figure name.'),
    argument_int('-x', '--width', help='Figure width in inches.', default=4),
    argument_int('-y', '--height', help='Figure height in inches.', default=4),
    argument_float('-s', '--scale', help='Figure scale.', default=1),
    argument_string('-f', '--font-family', help='Font faily used for the plot.', default=None)
])
def annotate(args):
    """Annotate epiloci with given annotation bed files.
    Optionally, a scatterplot with annotation marks can be generated.
    """
    prism.annotate.run(
        input_fp=args.input,
        output_fp=args.output,
        bed_fps=args.beds,
        annotation_names=args.annotation_names,
        output_figure_fp=args.figure,
        width=args.width,
        height=args.height,
        scale=args.scale,
        font_family=args.font_family,
    )

def main(args=None):
    args = parser.parse_args(args=args)
    args.func(args)
