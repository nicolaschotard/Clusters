#!/usr/bin/env python
"""Setup script."""
    
import os
import glob
import yaml

from setuptools import setup

README = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/README.rst'
VERSION = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/version.yaml'

# Get __version__ from version.py without importing package itself.
__version__ = yaml.load(open(VERSION))['version']

# Package name
name = 'Clusters'

# Packages (subdirectories in clusters/)
packages = ["clusters"]

# Modules (all python files in clusters/)
modules = [m.replace("clusters/", "%s." % name).replace('.py', '') for m in glob.glob("clusters/*.py")]

# Scripts (in scripts/)
scripts = glob.glob("scripts/*.py")

cmdclass = {}
command_options = {}
package_data = {}

print __version__

setup(name=name,
      version=__version__,
      description=("Cluster analysis tools and scripts"),
      license="MIT",
      classifiers=["Topic :: Scientific :: Astronomy",
                   "Intended Audience :: Science/Research"],
      url="https://github.com/nicolaschotard/Clusters",
      author="Nicolas Chotard, Dominique Boutigny, Celine Combet",
      author_email="nchotard@in2p3.fr",
      package_dir={name: 'clusters'},
      packages=packages,
      py_modules=modules,
      scripts=scripts,
      package_data=package_data,
      cmdclass=cmdclass,
      command_options=command_options,
      long_description=open(README).read(),
)
