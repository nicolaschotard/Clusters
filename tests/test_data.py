"""Test the reddening module."""

from Clusters import data


CONFIG = "configs/MACSJ2243.3-0935.yaml"

# Test the data module

def test_load_config():
    """Try to load a real config (yaml) file."""
    data.load_config(CONFIG)
