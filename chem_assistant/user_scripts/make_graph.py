#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['plot_opt']

from ..interfaces.gamess_results import GamessResults
from ..interfaces.psi_results import PsiResults
from ..interfaces.gaussian_results import GaussianResults
from .grep_results import get_results_class

import os
import re

"""
File: make_graph.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Graph of energy vs iteration of geom opt in GAMESS 
"""

def get_data(file):
    r = get_results_class(file)
    iters = []
    energies = []
    if isinstance(r, GamessResults):
        with open(file, "r") as f:
            for line in f.readlines():
                if 'NSERCH:' in line:
                    _, iteration, _, energy, *_ = line.split()
                    iters.append(int(iteration))
                    energies.append(float(energy))
    elif isinstance(r, GaussianResults):
        energies_per_iter = []
        with open(file, "r") as f:
            for line in f.readlines():
                if 'Step number' in line:
                    iters.append(int(line.split()[2]))
                if re.search('^\sE=\s*-?[0-9]*.[0-9]*', line):
                    energies_per_iter.append(float(line.split()[1]))
                if 'Step number' in line: # find next iteration
                    energies.append(energies_per_iter[-1])
                    energies_per_iter.clear()
    return {'iters': iters, 'energies': energies} 


def plot_opt(file):

    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    import seaborn as sns 
    # large import- only import if needed!

    data = get_data(file)
    sns.set_style('darkgrid')
    plt.figure(figsize=(8,6))
    sns.lineplot(x = "iters", y = "energies", data = data)
    plt.ticklabel_format(useOffset = False)
    plt.xlabel('Iterations')
    plt.ylabel('Energy, E$_h$')
    plt.show()
