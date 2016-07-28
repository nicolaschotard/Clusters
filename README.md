# Clusters

***Warning***: This package is under development!

Python package wrapping up the ongoing cluster analysis of the french LSST/DESC group. For more info, see the two following github repos:

- https://github.com/lsst-france/LSST_notebooks
- https://github.com/DarkEnergyScienceCollaboration/ReprocessingTaskForce

## Installation

To install (use --prefix to install it locally)

```
git clone https://github.com/nicolaschotard/Clusters.git
pip install Clusters/
```

To install a release version:

```
pip install http://github.com/nicolaschotard/Cluster/archive/v0.1.tar.gz
```

Release versions are listed
[here](http://github.com/nicolaschotard/Clusters/releases).

`Clusters` has for now (too) many dependencies:

- The lsst DM stack! (see [here](https://developer.lsst.io/build-ci/lsstsw.html))
- Python 2.7
- numpy
- matplotlib
- LEPHARE (see [here](http://cesam.lam.fr/lephare/lephare.html))
- and many other packages that will be replaced at some point


Usage
-----

`Clusters` consists of several command-line executables that you have
to run in the right order.

Get the input data and dump them in a pickle file (will change soon).

```
clusters_data config.yaml output.pkl
```

Correct the data for Milky Way extinction (output.pkl is the output of the previous step)

```
clusters_extinction.py config.yaml output.pkl --plot
```

**Next parts will come soon**

Get the photometric redshift using LEPHARE

```
clusters_zphot config.yaml input.pkl output.pkl
```

Exctract background galaxies from the whole sample (remove the red sequence)

```
clusters_background config.yaml input.pkl output.pkl
```

Etc.

With any command, you can run with `-h` or `--help` options to see all the
optional arguments, e.g., `cubefit --help`.

Input format
------------

All the scripts will take the same input YAML file. Keys are
case-sensitive. Additional keys are simply ignored. Each script will
aslo have a set of option that you can use.

| Parameter        | Type     | Description [units]                   |
| ---------------- | ------   | ------------------------------------- |
| `"cluster"`      | *string* | Name of the cluster |
| `"ra"`           | *float*  | RA coordinate of the cluster **[deg]** |
| `"dec"`          | *float*  | DEC coordinate of the cluster **[deg]** |
| `"redshift"`     | *float*  | Redshift the cluster |
| `"filters"`     | *string*  | Filter list to study, e.g., 'ugriz' (Megacam filters) |
| `"butler"`     | *string*  | Absolute path to the intput data (butler) |
| `"patches"`     | *list*  | List of patches to study |
