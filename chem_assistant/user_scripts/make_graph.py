#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['plot_opt']

"""
File: make_graph.py 
Author: Tom Mason
Email: tommason14@gmail.com
Github: https:github.com/tommason14
Description: Graph of energy vs iteration of geom opt in GAMESS 
"""
# Plotting

import pandas as pd
import matplotlib.pyplot as plt

def get_data(file):
    iters = []
    energies = []
    with open(file, "r") as f:
        for line in f.readlines():
            if 'NSERCH:' in line:
                _, iteration, _, energy, *_ = line.split()
                iters.append(int(iteration))
                energies.append(float(energy))
    return iters, energies 

def plot_opt(file):
    iters, energies = get_data(file)
    df = pd.DataFrame({'iters': iters, 'energies': energies})
    plot = df.plot.line(x = 'iters', y = 'energies', legend = False)
    plot.set_xlabel('Iterations')
    plot.set_ylabel('Energy, E$_\mathrm{h}$')
    plot.yaxis.major.formatter._useMathText = True
    plt.show()
