Deconvolute
===========

.. note::
    
    This section is under construction.

``prism deconvolute`` command infers the subclonal composition of the sample. Simply give methylation pattern-corrected epiloci file.

.. code-block:: console

    $ prism deconvolute -i sample.corrected.met -o sample.prism.result

Another feature of PRISM is that is can utilize two or more samples from a single tumor at the same time, and jointly infer subclonal composition. To provoke joint-analysis, specify a list of ``corrected.met`` files.

.. code-block:: console

    $ prism deconvolute -i sample1.corrected.met sample2.corrected.met -o sample.prism.result

For a more detailed description about all options, run ``prism deconvolute -h``.

