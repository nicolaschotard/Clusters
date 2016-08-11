Data format
-----------

The data format used in all scripts is based on the `Astropy Table
<http://docs.astropy.org/en/stable/table/>`_ format. The main table is
built from the LSST DM butler as shown in the diagram [1]_ below:

.. [1] Diagram created using https://www.jgraph.com/. Use the draw,io application and the last xml file to update the diagram if needed. 

.. image:: https://cdn.rawgit.com/nicolaschotard/Clusters/master/lib/data-table-1.1.svg
   :scale: 100 %
   :alt: Data table construction
   :align: center

For each filter `f`, an Astropy Table is create for all available
patch `p`. Since we have the same amount of patch for all filter,
which contain the exact same amount of sources, all table created from
a patch will be of the same size for all filter. Each table (1,p), (2,p), etc. 
