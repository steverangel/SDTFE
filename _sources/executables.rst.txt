Executables
===========

Example Data
------------

Download and extract example data.

.. code-block:: bash

   curl -O http://users.eecs.northwestern.edu/~emr126/sdtfe_examples.tar.gz
   tar xvfz sdtfe_examples.tar.gz



dtfe
----

.. code-block::

   Usage: dtfe [ path\_to\_file n\_particles grid\_dim center\_x center\_y center\_z field\_width field\_depth particle\_mass mc\_box\_width n\_mc\_samples sample\_factor ]

Run the examples with all the parameters specified:

.. code-block:: bash

   ./bin/dtfe ../data/913571938961.bin 216683 768 2556.9 1510.4 1986.6 6.0 4.0 1 0.01 5 0.25


Or by terminal prompts:

.. code-block:: bash

   ./bin/dtfe