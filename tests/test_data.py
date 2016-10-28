"""Test the reddening module."""

import os
from Clusters import data


CONFIG = "configs/MACSJ2243.3-0935.yaml"

# Test the data module

def test_load_config():
    """Try to load a real config (yaml) file."""
    data.load_config(CONFIG)


def test_catalogs_class():
    """Test the Clusters.data.Catalogs class."""
    if not os.path.exists('testdata'):
        get_testdata = """
        - wget https://lapp-owncloud.in2p3.fr/index.php/s/xG2AoS2jggbmP0k/download
        - tar zxvf download
        - rm -f download
        """
        raise IOError("No test data. Try: %s" % get_testdata)
    catalogs = ['forced_src', 'deepCoadd_meas', 'deepCoadd_forced_src']
    config = data.load_config("testdata/travis_test.yaml")
    cats = data.Catalogs(config['butler'])
    cats.show_keys(catalogs)
    cats.load_catalogs(catalogs, matchid=True, **config)
    cats.save_catalogs("travis_data")
    catalogs = data.read_hdf5("travis_data.hdf5")
