#!/usr/bin/env python

"""Setup script."""

import os
import glob
import yaml

from setuptools import setup, find_packages


README = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/README.rst'
VERSION = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/version.yaml'

# Get __version__ from version.py without importing package itself.
__version__ = yaml.load(open(VERSION))['version']

# Package name
name = 'Clusters'

# Packages (subdirectories in clusters/)
packages = find_packages()

# Scripts (in scripts/)
scripts = glob.glob("scripts/*.py")

package_data = {}

setup(name=name,
      version=__version__,
      description=("Cluster analysis tools and scripts"),
      license="MIT",
      classifiers=["Topic :: Scientific :: Astronomy",
                   "Intended Audience :: Science/Research"],
      url="https://github.com/nicolaschotard/Clusters",
      author="Nicolas Chotard, Dominique Boutigny, Celine Combet",
      author_email="nchotard@in2p3.fr",
      packages=packages,
      scripts=scripts,
      package_data=package_data,
      long_description=open(README).read(),
      setup_requires=['pytest-runner'],
      tests_require=['pytest']
)
