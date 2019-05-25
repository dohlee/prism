Scatter
=======

.. note::
    
    This section is under construction.

``prism scatter`` command generates a scatterplot of the PRISM analysis result. You need a result of ``prism deconvolute``. The dimension of analysis (i.e., the number of samples you gave to ``prism deconvolute`` command) should not be more than three to visualize it. Note that the file extension of output file should be a general one for image files such as png, jpg, or pdf.

.. code-block:: console

    $ prism scatter -i sample.prism.result -o sample.png

.. image:: _images/scatter.png
    :align: center
