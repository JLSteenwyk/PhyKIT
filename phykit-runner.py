#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convenience wrapper for running phykit directly from source tree."""
import sys

from phykit.phykit import main

if __name__ == "__main__":
    main(sys.argv[1:])
