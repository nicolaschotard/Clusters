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

Python package wrapping up the ongoing cluster analysis of the cluster
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


Dependencies
------------

``Clusters`` has for now the following dependencies (see the quick
installs below):

- Python 2.7 and libraries listed in the `requirements <requirements.txt>`_ file
- The LSST DM `stack <https://developer.lsst.io/build-ci/lsstsw.html>`_. 
- `LEPHARE <http://cesam.lam.fr/lephare/lephare.html>`_ (Photometric
  Analysis for Redshift Estimate)

Python
``````

To install the python dependencies, simply do::

  pip install -r requirements.txt


DM stack quick install
``````````````````````

This four-step procedure should allow you to install and configure a
light version of the DM stack, but complete enough to use the
``Clusters`` package. It should take ~10 minutes.

- Get and install miniconda, if you do not have it already::

    .. warning:: If you install the stack in an 'lsst' conda
                 environment, you must also install the cluster
                 package in this environment.
  
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda

- Install the needed part of the DM stack (we do not need the entire
  stack)::
    
    conda config --add channels http://conda.lsst.codes/stack/0.12.0
    conda create -q -n lsst python=2.7
    source activate lsst
    conda install -q gcc lsst-daf-persistence lsst-log lsst-afw lsst-skypix lsst-meas-algorithms lsst-pipe-tasks

- Install a compatible version of the ``obs_cfht`` package, which is
  not yet avalaible in the conda repository::

    setup daf_persistence
    git clone https://github.com/lsst/obs_cfht.git
    cd obs_cfht
    git checkout b7ab2c4
    setup -k -r .
    scons opt=3
    eups declare -r . -t yourname
  
- To use this install of the DM stack, do not forget these following
  setups::
  
    export PATH="$HOME/miniconda/bin:$PATH"
    source activate lsst
    source eups-setups.sh
    setup daf_persistence
    setup obs_cfht -t yourname

If these steps went well, you should be able to use
``clusters_data.py`` on one of the outputs of the DM stack (see below
to get some data).

LEPHARE quick install
`````````````````````

You can download and install a pre-configured version of LEPHARE as
followed:

- for linux system::
    
    wget https://lapp-owncloud.in2p3.fr/index.php/s/MDaXObLSD9IVQ1B/download -O lephare.tar.gz
    tar zxf lephare.tar.gz

- for mac::

    wget http://lpsc.in2p3.fr/upload/doc/fdd5e8/lephare_macosx.tar.gz -O lephare.tar.gz
    tar zxf lephare.tar.gz
    
When the download is complete, exctract the ``lephare`` directory where it
suits you (``mypath`` in this example), and set the following
environment variables (use setenv if needed)::

    export LEPHAREWORK="mypath/lephare/lephare_work"
    export LEPHAREDIR="mypath/lephare/lephare_dev"
    export PATH="$PATH:mypath/lephare/lephare_dev/source"

You should now be able to run ``clusters_zphot.py`` (only tested on
linux systems).


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


Test the code
-------------

If you have installed all the dependencies previoulsy mentionned,
download the following test data set::

  wget https://lapp-owncloud.in2p3.fr/index.php/s/xG2AoS2jggbmP0k/download -O testdata.tar.gz
  tar zxf testdata.tar.gz

The ``testdata`` directory contains a subset of the reprocessing data
available for MACSJ2243.3-0935. It can be used as a test set of the
code, but is not complete enough to run the full analysis. Here is the
full structure and content of the directory, which has the exact same
structure as a regulare DM stack output directory::

  testdata/
  ├── input
  │   ├── _mapper
  │   └── registry.sqlite3
  ├── output
  │   ├── coadd_dir
  │   │   ├── deepCoadd
  │   │   │   ├── g
  │   │   │   │   └── 0
  │   │   │   │       ├── 1,5
  │   │   │   │       └── 1,5.fits
  │   │   │   └── skyMap.pickle
  │   │   ├── deepCoadd-results
  │   │   │   └── g
  │   │   │       └── 0
  │   │   │           └── 1,5
  │   │   │               ├── bkgd-g-0-1,5.fits
  │   │   │               ├── calexp-g-0-1,5.fits
  │   │   │               ├── detectMD-g-0-1,5.boost
  │   │   │               ├── det-g-0-1,5.fits
  │   │   │               ├── forced_src-g-0-1,5.fits
  │   │   │               ├── meas-g-0-1,5.fits
  │   │   │               ├── measMD-g-0-1,5.boost
  │   │   │               └── srcMatch-g-0-1,5.fits
  │   │   ├── forced
  │   │   │   └── 08BO01
  │   │   │       └── SCL-2241_P1
  │   │   │           └── 2008-09-03
  │   │   │               └── g
  │   │   │                   └── 0
  │   │   │                       ├── FORCEDSRC-1022175-00.fits
  │   │   │                       ├── FORCEDSRC-1022175-09.fits
  │   │   │                       ├── FORCEDSRC-1022176-00.fits
  │   │   │                       ├── FORCEDSRC-1022176-09.fits
  │   │   │                       ├── FORCEDSRC-1022177-00.fits
  │   │   │                       ├── FORCEDSRC-1022177-09.fits
  │   │   │                       ├── FORCEDSRC-1022178-00.fits
  │   │   │                       ├── FORCEDSRC-1022178-09.fits
  │   │   │                       ├── FORCEDSRC-1022179-00.fits
  │   │   │                       ├── FORCEDSRC-1022179-09.fits
  │   │   │                       ├── FORCEDSRC-1022180-00.fits
  │   │   │                       └── FORCEDSRC-1022180-09.fits
  │   │   └── _parent -> ../
  │   └── _parent -> ../input/
  └── travis_test.yaml

With this data set, you should be able to test most of the
``Clusters`` parts, starting with the ``clusters_data.py`` script.


Get the data
------------

Raw DM stack outputs
`````````````````````

If you have installed ``Clusters`` but do not have any data to run it
on, you can use one of our re-processing outputs for
MACSJ2243.3-0935. The corresponding configuration file is stored
`there <configs/MACSJ2243.3-0935.yaml>`_. To use it, you either need
to be connected at CC-IN2P3, or change the path to the butler inside
the config file (if you already have a copy of this data). You could
also mount sps on your personal computer (see this `how to
<http://lsstnotes.readthedocs.io/en/latest/sshfs.html>`_).


``clusters_data.py`` output
```````````````````````````

The first step of the ``Clusters`` package is ``clusters_data.py``,
which will get the data from the DM butler, convert them into
``astropy`` tables and save them in a single ``hdf5`` file. To do so,
you need the LSST DM stack to be installed. If you want to skip this
part and try the code whithout having to install the DM stack, you
could also use the outputs of this first step that you can download
from `this repository
<https://lsst-web.ncsa.illinois.edu/~nchotard/data/clusters/>`_, which
contains the following files::

  |-- CL0016
  |   |-- [4.4G]  CL0016_data.hdf5                     # full data set
  |   |-- [334M]  CL0016_filtered_data.hdf5            # only quality-filtered galaxies
  |   `-- [ 312]  CL0016.yaml                          # configuration file
  |-- MACSJ224330935
  |   |-- [5.6G]  MACSJ2243.3-0935_data.hdf5           # full data set
  |   |-- [367M]  MACSJ2243.3-0935_filtered_data.hdf5  # only quality-filtered galaxies
  |   |-- [ 329]  MACSJ2243.3-0935.yaml                # configuration file


This `short tutorial
<http://clusters.readthedocs.io/en/latest/data.html#work-with-the-table>`_
explains how to use these ``hdf5`` files to start an analysis.
