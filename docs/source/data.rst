Data format
===========

Overview
--------

The data format used in all scripts is based on the `Astropy Table
<http://docs.astropy.org/en/stable/table/>`_ format. In the `Build the
table`_ section, we show how these Astropy Tables are created from the
DM butler. This work is automatically done for the ``deepCoadd_meas``,
``deepCoadd_forced_src``, and ``forced_src`` catalogs, and then saved
together in one single `hdf5 <http://www.h5py.org/>`_ file. The
procedures to write and read these files, and work with the loaded
tables, are described in the `Working with the table`_ section.


Build the table
---------------

The main table is built from the LSST DM butler as shown in the
diagram [#]_ below:

.. image:: https://cdn.rawgit.com/nicolaschotard/Clusters/master/docs/source/_static/data-table-1.1.svg
   :scale: 100 %
   :alt: Data table construction
   :align: center

For each filter ``f``, an Astropy Table is created for all available
patches ``p1``, ``p2``, etc. Since we have the same number of patches
for all filters, which contain the exact same number of sources, all
table (1,p), (2,p), etc., created from a patch will be of the same
size for all filter. The Astropy Tables created from all individual
filters/patches set will then be vertically stacked. This means that if
we respectively have ``N1``, ``N2`` and ``N3`` sources for the patches
1, 2 and 3, the output table will contains ``N1`` + ``N2`` + ``N3``
sources. After the first stack, we end up with one table per filter,
all containing the same number of sources. These per-filter tables are
finally stacked together to form the final table shown on the
right-hand side of the diagram.

--------

.. [#] Diagram created using https://www.jgraph.com/. Use the
       https://www.draw.io/?demo=1 application and the last xml file
       from this repository to update the diagram if needed.

WCS
---

The WCS computed during the data processing is also stored in the
``hdf5`` file in the ``wcs`` path. If you load the data using the
`read_hdf5 <clusters.html#clusters.data.read_hdf5>`_ function, the
output will be a dictionary containing a ``wcs`` key, which refers to
an `astropy.wcs.WCS
<http://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS>`_
object. The `skycoord_to_pixel
<clusters.html#clusters.data.skycoord_to_pixel>`_ and
`pixel_to_skycoord <clusters.html#clusters.data.pixel_to_skycoord>`_
functions take this ``wcs`` object to convert the input coordinates
into an output format (sky <-> pixel).

All the tables (correspondig to the catalogs listed previously)
already contain three coordinates columns:

- ``coord_ra`` and ``coord_dec``: they are the (ra, dec) coordinates in radian;
- ``coord_ra_deg`` and ``coord_dec``: they are the (ra, dec) coordinates in degree;
- ``x_Src`` and ``y_Src``: they are the (x, y) position in pixel.


Working with the table
----------------------

.. include:: data_tuto.rst
