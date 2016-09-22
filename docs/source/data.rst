Data format
===========

Overview
--------

The data format used in all scripts is based on the `Astropy Table
<http://docs.astropy.org/en/stable/table/>`_ format. In the `Build the
table`_ section, we show how these Astropy Tables are create from the
DM bulter. This work is automaticaly done for the `meas` and `forced`
catalogs, and then saved together in one single `hdf5
<http://www.h5py.org/>`_ file. The procedures to write and read these
files, and work with the loaded tables, are describe in the `Work with
the table`_ section.


Build the table
---------------

The main table is built from the LSST DM butler as shown in the
diagram [#]_ below:

.. image:: https://cdn.rawgit.com/nicolaschotard/Clusters/master/docs/source/_static/data-table-1.1.svg
   :scale: 100 %
   :alt: Data table construction
   :align: center

For each filter `f`, an Astropy Table is create for all available
patches `p1`, `p2`, etc. Since we have the same amount of patches for
all filters, which contain the exact same amount of sources, all table
(1,p), (2,p), etc., created from a patch will be of the same size for
all filter. The Astropy Tables created from all individual
filters/pathes set will then be verticaly stacked. This means that if
we respecively have `N1`, `N2` and `N3` sources for the patches 1, 2
and 3, the output table will contains `N1` + `N2` + `N3`
sources. After the first stack, we end up with one table per filter,
all containing the same number of sources. These per-filter tables are
finaly stacked together to form the final table shown on the
right-hand side of the diagram.

--------

.. [#] Diagram created using https://www.jgraph.com/. Use the
       https://www.draw.io/?demo=1 application and the last xml file
       from this repository to update the diagram if needed.


Work with the table
-------------------

.. include:: data_tuto.rst
