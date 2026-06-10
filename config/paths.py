#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
config/paths.py
------------------------------------------------------------
Central path configuration for the Python preprocessing stage.

All scripts resolve input/output locations from environment
variables instead of hardcoded absolute paths. Set these in your
shell before running, or rely on the local-relative defaults below.

    DATA_DIR     Xenium spatial input datasets   (was: Datasets_April)
    EXTDATA_DIR  External validation datasets     (was: External_Data)
    RESULTS_DIR  Pipeline outputs / results        (was: Results_New)

Usage inside a script:
    from config.paths import DATA_DIR, RESULTS_DIR
    # or, if config is not a package on sys.path:
    import os
    DATA_DIR = os.environ.get("DATA_DIR", "./data")

Bash:        export DATA_DIR=/path/to/Xenium_datasets
PowerShell:  $env:DATA_DIR = "D:\\path\\to\\Xenium_datasets"
============================================================
"""

import os

DATA_DIR    = os.environ.get("DATA_DIR",    "./data")
EXTDATA_DIR = os.environ.get("EXTDATA_DIR", "./external_data")
RESULTS_DIR = os.environ.get("RESULTS_DIR", "./results")
