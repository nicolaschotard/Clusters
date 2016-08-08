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

def readme():
    with open('README.rst') as f:
        return f.read()
    
# Package name
name = 'Clusters'

# Packages (subdirectories in lib/)
packages = []

# Modules (all python files in lib/)
modules = [m.replace("lib/", "%s." % name).replace('.py', '') for m in glob.glob("lib/*.py")]

# Scripts (in scripts/)
scripts = glob.glob("scripts/*.py")

cmdclass = {}
command_options = {}
package_data = {}

# build_sphinx command in case sphinx is installed
try:
    from sphinx.setup_command import BuildDoc
    cmdclass.update({'build_sphinx': BuildDoc})
    command_options.update({'build_sphinx': {
        'version': ('setup.py', __version__),
        'release': ('setup.py', __version__)}
    })
except ImportError:
    pass

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
      command_options=command_options
)
