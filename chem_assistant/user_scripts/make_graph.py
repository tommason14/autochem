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
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import seaborn as sns

def get_data(file):
    iters = []
    energies = []
    with open(file, "r") as f:
        for line in f.readlines():
            if 'NSERCH:' in line:
                _, iteration, _, energy, *_ = line.split()
                iters.append(int(iteration))
                energies.append(float(energy))
    return {'iters': iters, 'energies': energies} 

def plot_opt(file):
    data = get_data(file)
    sns.set_style('darkgrid')
    plt.figure(figsize=(8,6))
    sns.lineplot(x = "iters", y = "energies", data = data)
    plt.ticklabel_format(useOffset = False)
    plt.xlabel('Iterations')
    plt.ylabel('Energy, E$_h$')
    plt.show()
