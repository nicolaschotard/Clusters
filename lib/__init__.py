#!/usr/bin/env python

import os, glob

# Automatically import all modules (python files)
__all__ = [os.path.basename(m).replace('.py', '') for m in glob.glob("lib/*.py")
           if '__init__' not in m]

# Set to True if you want to import all previous modules directly
importAll = True

if importAll:
    for pkg in __all__:
        __import__(__name__ + '.' + pkg)
