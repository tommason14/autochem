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

import pandas as pd
from ..core.utils import module_exists
if module_exists('plotnine'):
    from plotnine import *
else:
    import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

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
    if module_exists('plotnine'):
        plot = ggplot(df) + geom_line(aes(x = 'iters', y = 'energies'))+\
               labs(x = 'Iterations',  y = 'Energy, E$_\mathrm{h}$')
        print(plot)
    else:
        plot = df.plot.line(x = 'iters', y = 'energies', legend = False)
        plot.set_xlabel('Iterations')
        plot.set_ylabel('Energy, E$_\mathrm{h}$')
        plot.yaxis.major.formatter._useMathText = True
        plt.show()
