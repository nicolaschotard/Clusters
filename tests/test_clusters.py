"""Test the reddening module."""

import os
from Clusters import main
from Clusters import data


CONFIG = "configs/MACSJ2243.3-0935.yaml"

# Test the data module

def test_load_config():
    """Try to load a real config (yaml) file."""
    data.load_config(CONFIG)


def test_catalogs_class(config="testdata/travis_test.yaml", datafile="travis_data"):
    """Test the Clusters.data.Catalogs class."""
    if not os.path.exists('testdata'):
        get_testdata = """
        - wget https://lapp-owncloud.in2p3.fr/index.php/s/xG2AoS2jggbmP0k/download
        - tar zxvf download
        - rm -f download
        """
        raise IOError("No test data. Try: %s" % get_testdata)
    catalogs = ['forced_src', 'deepCoadd_meas', 'deepCoadd_forced_src']
    config = data.load_config(config)
    cats = data.Catalogs(config['butler'])
    cats.load_catalogs(catalogs, matchid=True, **config)
    cats.load_catalogs(None, show=True)
    cats.save_catalogs(datafile)


def test_data_functions(datafile="travis_data.hdf5"):
    """Test functions of data.py."""
    # Read the hdf5 file and load the catalogs
    catalogs = data.read_hdf5(datafile)

    # Apply filters to the catalogs (keep galaxies)
    fcatalogs = data.filter_table(catalogs)

    # Make sure we can get data from the filtered catalogs
    ra = fcatalogs['deepCoadd_forced_src']['coord_ra'].tolist()
    dec = fcatalogs['deepCoadd_forced_src']['coord_dec'].tolist()
    wcs = data.load_wcs(fcatalogs['wcs'])

    # Transformations: coordinates <-> pixel
    x, y = data.skycoord_to_pixel([ra, dec], wcs, unit='rad')
    data.pixel_to_skycoord(x, y, wcs)


# Test the pipeline


#def test_main(config="testdata/travis_test.yaml", datafile="travis_test_data.hdf5"):
#    """Test the pipeline."""
#    main.load_data([config, "--output", datafile, "--overwrite"])
#    main.load_data([config, "--show", "--overwrite"])
#    filtered_data = datafile.replace('.hdf5', '_filtered_data.hdf5')
#    main.extinction([config, filtered_data, "--plot", "--overwrite"])
#    extinction_data = filtered_data.replace('.hdf5', '_extinction.hdf5')
#    main.photometric_redshift([config, filtered_data, "--extinction",
#                               extinction_data, "--overwrite"])
