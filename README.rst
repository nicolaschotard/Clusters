.. image:: https://readthedocs.org/projects/clusters/badge/?version=latest
   :target: http://clusters.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://landscape.io/github/nicolaschotard/Clusters/master/landscape.svg?style=flat
   :target: https://landscape.io/github/nicolaschotard/Clusters/master
   :alt: Code Health

.. image:: https://travis-ci.org/nicolaschotard/Clusters.svg?branch=master
   :target: https://travis-ci.org/nicolaschotard/Clusters
   :alt: Travis CI build status (Linux)

.. image:: https://codecov.io/gh/nicolaschotard/Clusters/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/nicolaschotard/Clusters

____

**WARNING**: Package under development

____

.. inclusion-marker-do-not-remove

Clusters
--------

Python package wrapping up the ongoing cluster analysis of the french
LSST/DESC group. For more info, see the two following github
repositories:

- A collection of `notebooks <https://github.com/lsst-france/LSST_notebooks>`_ for LSST
- The `ReprocessingTaskForce <https://github.com/DarkEnergyScienceCollaboration/ReprocessingTaskForce>`_ repository

See also the private `Trello board
<https://trello.com/b/Lhg6VAq2/clusters>`_ that we use to share our
work.

Installation
------------

To install::

  git clone https://github.com/nicolaschotard/Clusters.git
  pip install Clusters/

To install in a local directory ``mypath``, use::

  pip install --prefix='mypath' Clusters/

and do not forget to add it to your PYTHONPATH.

To upgrade to a new version (after a ``git pull`` or a local modification), use::

  pip install --upgrade (--prefix='mypath') Clusters/

To install a release version (no release version available yet)::

  pip install http://github.com/nicolaschotard/Cluster/archive/v0.1.tar.gz

Also works with the master::

  pip install (--upgrade) https://github.com/nicolaschotard/Clusters/archive/master.zip

In the future, release versions will be listed at this `location
<http://github.com/nicolaschotard/Clusters/releases>`_.

To install dependencies, simply do::

  pip install -r requirements.txt


Dependencies
------------

``Clusters`` has for now the following dependencies:

- The LSST DM `stack <https://developer.lsst.io/build-ci/lsstsw.html>`_
- Python 2.7 and libraries listed in the `requirements <requirements.txt>`_ file
- `LEPHARE <http://cesam.lam.fr/lephare/lephare.html>`_


Configuration file
------------------

All the scripts will take the same input YAML file, which contains
necessary informations for the analysis or simply for plotting purpose,
such as the name of the studied cluster. Keys are listed below and are
case-sensitive. Additional keys are simply ignored. You can find
examples of these configuration files in the `config
<https://github.com/nicolaschotard/Clusters/blob/master/configs>`_
directory, or clicking `here
<https://github.com/nicolaschotard/Clusters/blob/master/configs/MACSJ2243.3-0935.yaml>`_
for MACSJ2243.3-0935.

+--------------------+--------+-------------------------------------------------------------------+
| General keys       | Type   | Description [units]                                               |
+====================+========+===================================================================+
| ``"cluster"``      | string | Name of the cluster                                               |
+--------------------+--------+-------------------------------------------------------------------+
| ``"ra"``           | float  | RA coordinate of the cluster **[deg]**                            |
+--------------------+--------+-------------------------------------------------------------------+
| ``"dec"``          | float  | DEC coordinate of the cluster **[deg]**                           |
+--------------------+--------+-------------------------------------------------------------------+
| ``"redshift"``     | float  | Cluster redshift                                                  |
+--------------------+--------+-------------------------------------------------------------------+
| ``"butler"``       | string | Absolute path to the intput data (butler)                         |
+--------------------+--------+-------------------------------------------------------------------+
| ``"filter"``       | list   | List of filters to be considered, e.g., 'ugriz' (Megacam filters) |
+--------------------+--------+-------------------------------------------------------------------+
| ``"patch"``        | list   | List of patches to study                                          |
+--------------------+--------+-------------------------------------------------------------------+

The following list of optional keys can also be added to the
configuration file. They correspond to specific configurations of the
different steps of the analysis. While the previous list will most
likely stay unchanged, the following one will be completed with new
keys as this analysis will progress.

+----------------------+--------+------------------------------------------------------------------+
| Optional keys        | Type   | Description [units]                                              |
+======================+========+==================================================================+
| ``"keys"``           | dict   | Dictionnary containing list of keys for the catalogs (see below) |
+----------------------+--------+------------------------------------------------------------------+
| ``"zpara"``          | list   | List of paths to ``zphota`` configuration files (see below)      |
+----------------------+--------+------------------------------------------------------------------+
| ``"zspectro_file"``  | string | File containing spectroz sample for LePhare training             |
+----------------------+--------+------------------------------------------------------------------+

- ``keys`` is a dictionary having the name of the different catalogs
  like **deepCoadd_meas**, **deepCoadd_forced_src** and
  **forced_src**. The list of keys for a given catalog can include:

  - "the_full_name_of_a_key";
  - "\*_a_part_of_a_key_name" or "an_other_part_of_a_key_name\*"
    preceded or followed by a \*;
  - a combination of all the above: ["key1", "ke\*", "\*ey"];
  - or a "*" to get all keys available in a catalog, which is the
    default value for all catalogs.




General usage
-------------

``Clusters`` consists in several command-line executables that you
have to run in the right order.

- Get the input data and dump them in a hdf5 file containing astropy
  tables (see the `data format section
  <http://clusters.readthedocs.io/en/latest/data.html>`_ of the
  documentation for detail)::

    clusters_data config.yaml (--output data.hdf5)

You can adapt the content of the output file using the ``keys``
parameter of the config.yaml file.

- Correct the data for Milky Way extinction::

    clusters_extinction.py config.yaml data.hdf5 (--output extinction.hdf5)

- Get the photometric redshift using LEPHARE::

    clusters_zphot.py config.yaml data.hdf5 (--extinction extinction.hdf) (--output zphot.hdf5)

The configuration file(s) used in LEPHARE can be given with the option
``--zpara``. The code will loop over the different files and run
LEPHARE for each of them. All results are saved in the same ``hdf5``
file. This list of configuration files can also be given in the
CONFIG.yaml file (see above). ``--zpara`` will overwrite what is given
in the configuration file.

- Extract background galaxies from the whole sample: remove the
  cluster galaxies (red sequence) and other foreground galaxies using
  the photometric redshifts::

    clusters_getbackground config.yaml input.hdf5 output.hdf5

- Compute the shear::

    clusters_shear config.yaml input.hdf5 output.hdf5

- A pipeline script which run all the above step in a raw with standard options::

    clusters_pipeline config.yaml

With any command, you can run with ``-h`` or ``--help`` to see all the
optional arguments, e.g., ``clusters_data.py -h``.


An example
----------

If you have installed ``Clusters`` but do not have any data to run it
on, you can use one of our re-processing outputs for
MACSJ2243.3-0935. The corresponding configuration file is stored
`there <configs/MACSJ2243.3-0935.yaml>`_. To use it, you either need
to be connected at CC-IN2P3, or change the path to the butler inside
the config file (if you have your own data). You could also mount sps
on your personal computer (see this `how to
<http://lsstnotes.readthedocs.io/en/latest/sshfs.html>`_).

The first step of the ``Clusters`` package is ``clusters_data.py``,
which will get the data from the DM butler, convert them into a
``astropy`` tables and save them in a single ``hdf5`` file. To do so,
you need the LSST DM stack to be installed. If you want to skip this
part and try the code whithout having to install the DM stack, you
could also use the output of this first step for MACSJ2243.3-0935 that
we have stored under::

  /sps/lsst/data/clusters/MACSJ2243.3-0935/analysis/output_v1/MACSJ2243.3-0935_data.hdf5

A `short tutorial
<http://clusters.readthedocs.io/en/latest/data.html#work-with-the-table>`_
explains how to use this ``hdf5`` file to start an analysis.
