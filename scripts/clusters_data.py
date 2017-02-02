#!/usr/bin/env python
"""Load the cluster data."""

import sys
from clusters.mains import data

sys.exit(data.load_data())
