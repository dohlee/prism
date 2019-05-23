Extract
=======

.. note::
    
    This section is under construction.

``prism extract`` command extracts viable epiloci from a BAM file. In other words, it extracts genomic regions harboring a sufficient number of mapped reads (>= d) with a sufficient number of CpGs (>= c). A resulting file with those epiloci information is generated, whose file extension will be ``.met`` afterwards. To extract epiloci with default settings (d = 20, c = 4), simply run the command below:

.. code-block:: console

    $ prism extract -i sample.bam -o sample.met

If you want to specify your own cutoffs for the required sequencing depth and the number of CpGs, say, d = 15 and c = 3, run as follows:

.. code-block:: console

    $ prism extract -i sample.bam -o sample.met -d 15 -c 3
