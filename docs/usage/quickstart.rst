Quickstart
==========

Install PRISM with pip.

.. code-block:: console

    $ pip install subclone-prism

Run analysis.

.. code-block:: console

    ------------------------------
    Extract epiloci from BAM file.
    ------------------------------
    $ prism extract -i sample.bam -o sample.met

    -------------------------------------------------------------------
    Preprocess epiloci to get finer estimates of epigenetic subclones
    and to rescue more fingerprint epiloci from noisy methylation data.
    -------------------------------------------------------------------
    $ prism preprocess -i sample.met -o sample.corrected.met

    ----------------------------------------------
    Infer the subclonal composition of the sample.
    ----------------------------------------------
    If you want 1-sample deconvolution, run:
    $ prism deconvolute -i sample.corrected.met -o sample.prism.result
    or if you want 2-sample deconvolution, run:
    $ prism deconvolute -i sample1.corrected.met sample2.corrected.met -o sample.prism.result

    --------------------------------------------
    Scatterplot for visualization of the result.
    --------------------------------------------
    $ prism scatter -i sample.prism.result -o sample.png

    --------------------------------------------------------------------
    Annotation of fingerprint epiloci for functional characterization of
    discovered epigenetic subclones.
    --------------------------------------------------------------------
    $ prism annotate -i sample.prism.result -o sample.prism.annotated.result \
    --beds annotation_a.bed annotation_b.bed \
    --annotation-names ANNOTATION-A ANNOTATION-B

The PRISM analysis is done in three steps: ``extract`` - ``preprocess`` - ``deconvolute``.
