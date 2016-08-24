Data format
===========

Overview
--------

The data format used in all scripts is based on the `Astropy Table
<http://docs.astropy.org/en/stable/table/>`_ format. In the `Building
the table`_ section, we show how these Astropy Tables are create from
the DM bulter. This work is automaticaly done for the `meas` and
`forced` catalgue, and then saved together in one single `hdf5`
file. The procedure to write and read these files is describe in the
`Saving / reading the table`_ section. In the `Working with the
table`_, we show how you can use these tables to run your analysis or
make plots.


Building the table
------------------

The main table is
built from the LSST DM butler as shown in the diagram [#]_ below:

.. image:: https://cdn.rawgit.com/nicolaschotard/Clusters/master/docs/source/_static/data-table-1.1.svg
   :scale: 100 %
   :alt: Data table construction
   :align: center

For each filter `f`, an Astropy Table is create for all available
patch `p`. Since we have the same amount of patch for all filter,
which contain the exact same amount of sources, all table (1,p),
(2,p), etc., created from a patch will be of the same size for all
filter. The Astropy Tables created from all individual filters/pathes
set will then be verticaly stacked. This means that if we respecively
have `N1`, `N2` and `N3` sources for the patches 1, 2 and 3, the
output table will contains `N1` + `N2` + `N3` sources. After the first
stack, we end up with one table per filter, all containing the same
number of sources. These per-filter table are finaly stacked together
to form the final table show on the right-hand side of the diagram.

--------

.. [#] Diagram created using https://www.jgraph.com/. Use the
       https://www.draw.io/?demo=1 application and the last xml file
       from this repository to update the diagram if needed.

Saving / reading the table
--------------------------

Working with the table
----------------------

Such an Asptroy table is great to work with, and can be used for all
kind of analysis in the context of our cluster study. You can apply
filters, group by column, concatenate them, etc.
