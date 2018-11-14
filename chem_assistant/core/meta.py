#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json

__all__ = ['Meta']

class Meta(object):
    """Class for parsing *meta.json* files, present in every directory where computational calculations are desired.

Usage:
    >>> meta = Meta('meta.json')
    >>> meta.using # returns coordinates i.e. 'initial.xyz'
"""
    
    def __init__(self, filename):
        with open(filename, "r") as f:
            object.__setattr__(self, '_meta', json.load(f))       
            # https://stackoverflow.com/questions/16171289/python-how-to-use-setattr-on-dictionary-object-that-is-part-of-the-class-th           
 
    def __getattr__(self, key):
        try:
            return self._meta[key]
        except KeyError:
            raise AttributeError(f'No {key} in this class') 

    def __setattr__(self, key, value):
        self._meta[key] = value

# Example meta.json
# {
#   "runtype": "optimisation",
#   "is_fmo": true, 
#   "frags": 5,
#   "supercomp": "raijin",
#   "package": "gamess", 
#   "using": "initial.xyz"
# }
