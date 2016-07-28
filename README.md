# Clusters

Cluster package analysis

## Installation

To install a release version:

```
pip install http://github.com/nicolaschotard/Cluster/archive/v0.1.tar.gz
```

Release versions are listed
[here](http://github.com/nicolaschotard/Cluster/releases). `Clusters`
has for now (too) many dependencies:

- The lsststack ! (see [here](https://developer.lsst.io/build-ci/lsstsw.html))
- Python 2.7
- numpy
- matplotlib
- LEPHARE (see [here](http://cesam.lam.fr/lephare/lephare.html))
- and many other packages that will be replaced at some point


Usage
-----

`Clusters` consists of several command-line executables that you have
to run in the right order.

Get the input data and put them in a nice format.

```
clusters_data config.json output.pkl
```

Get the photometric redshift using LEPHARE

```
clusters_zphot config.json input.pkl output.pkl
```

Exctract background galaxies from the whole sample (remove the red sequence)

```
clusters_background config.json input.pkl output.pkl
```

Etc.

With any command, you can run with `-h` or `--help` options to see all the
optional arguments, e.g., `cubefit --help`.

Input format
------------

All the scripts will take the same input JSON file. Keys are
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
