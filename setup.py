#!/usr/bin/env python

import os
import sys
import glob

from setuptools import setup
from setuptools.command.test import test as TestCommand

# Get __version__ from version.py without importing package itself.
with open('/'.join(os.path.realpath(__file__).split('/')[:-1]) + \
     '/version.txt') as f:
    exec(f.read())

# Package name
name = 'Clusters'

# Packages (subdirectories in lib/)
packages = ["lib"]

# Modules (all python files in lib/)
modules = [m.replace("lib/", "%s." % name).replace('.py', '') for m in glob.glob("lib/*.py")]

# Scripts (in scripts/)
scripts = glob.glob("scripts/*.py")

cmdclass = {}
command_options = {}
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
      package_dir={name: 'lib'},
      packages=packages,
      py_modules=modules,
      scripts=scripts,
      package_data=package_data,
      cmdclass=cmdclass,
      command_options=command_options,
      long_description=open('README.rst').read(),
)
