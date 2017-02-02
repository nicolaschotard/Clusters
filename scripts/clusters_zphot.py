#!/usr/bin/env python
"""Comput photometric redshift using LEPHARE."""

import sys
from clusters.mains import zphot

sys.exit(zphot.photometric_redshift())
