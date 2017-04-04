"""Test the reddening module."""

import os
from clusters.mains import data, extinction, zphot

CONFIG = "testdata/travis_test.yaml"
DATAFILE = "travis_test_data.hdf5"

# Test the data module

def test_load_config():
    """Try to load a real config (yaml) file."""
    data.cdata.load_config(CONFIG)


def test_catalogs_class(config=CONFIG, datafile=DATAFILE):
    """Test the clusters.data.Catalogs class."""
    if not os.path.exists('testdata'):
        get_testdata = """
        - wget https://lapp-owncloud.in2p3.fr/index.php/s/xG2AoS2jggbmP0k/download
        - tar zxvf download
        - rm -f download
        """
        raise IOError("No test data. Try: %s" % get_testdata)
    catalogs = ['forced_src', 'deepCoadd_meas', 'deepCoadd_forced_src']
    config = data.cdata.load_config(config)
    cats = data.cdata.Catalogs(config['butler'])
    cats.load_catalogs(catalogs, matchid=True, **config)
    cats.load_catalogs(None, show=True)
    cats.save_catalogs(datafile.split('.')[0])


def test_data_functions(datafile=DATAFILE):
    """Test functions of data.py."""
    # Read the hdf5 file and load the catalogs
    catalogs = data.cdata.read_hdf5(datafile)

    # Apply filters to the catalogs (keep galaxies)
    fcatalogs = data.cdata.filter_table(catalogs)

    # Make sure we can get data from the filtered catalogs
    ra = fcatalogs['deepCoadd_forced_src']['coord_ra'].tolist()
    dec = fcatalogs['deepCoadd_forced_src']['coord_dec'].tolist()
    wcs = data.cdata.load_wcs(fcatalogs['wcs'])

    # Transformations: coordinates <-> pixel
    x, y = data.cdata.skycoord_to_pixel([ra, dec], wcs, unit='rad')
    data.cdata.pixel_to_skycoord(x, y, wcs)


# Test the pipeline


def test_main(config=CONFIG, datafile=DATAFILE):
    """Test the pipeline."""

    data.load_data([config, "--output", datafile, "--overwrite"])
    data.load_data([config, "--show", "--overwrite"])
    filtered_data = datafile.replace('.hdf5', '_filtered_data.hdf5')
    extinction.extinction([config, filtered_data, "--overwrite"])
    zphot.photometric_redshift([config, filtered_data, "--overwrite"])
#    main.getbackground([config, filtered_data, "--overwrite"])    
