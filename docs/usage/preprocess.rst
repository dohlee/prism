Preprocess
==========

.. note::

    This section is under construction.

``prism preprocess`` command corrects for the errors in methylation patterns in order to amplify the number of fingerprint epiloci and calibrate for the subclone size estimates.

.. code-block:: console

    $ prism preprocess -i sample.met -o sample.corrected.met

You may benefit from multithreading with ``-t/--threads`` option.

.. code-block:: console

    $ prism preprocess -i sample.met -o sample.corrected.met -t 30

This step is resource intensive, so by default PRISM tries to pre-filter the epilocus that is not likely to be a fingerprint epilocus. This pre-filtering of course can be turned off by ``-f/--no-prefilter`` option and this indeed gives a better estimates of subclones. However, please be warned, depending on your data size, this step will last long. To help you deciding whether or not to apply prefiltering, with 30 threads (``-t 30``) ~200M met file took about 5 hours to be preprocessed without prefiltering.

.. code-block:: console

    $ prism preprocess -i sample.met -o sample.corrected.met --no-prefilter -t 30

For a more detailed description about all options, run ``prism preprocess -h``.
