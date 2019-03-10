from dfply import *
from plotnine import *
import sys

# Pretty
import matplotlib.font_manager as fm
fpath = '/System/Library/Fonts/Avenir.ttc'
font = fm.FontProperties(fname=fpath)

# Unnecessary maplotlib warnings
import warnings
warnings.catch_warnings()
warnings.simplefilter('ignore')


def bar(df, x, y, fill = None, position = None):

    plot = None
    
    if fill is not None and position is not None:
        plot = (
        ggplot(df) +
        geom_bar(aes(x, y, fill = fill), stat = 'identity', position = position) +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_fill_brewer(palette = 'Set2', type = 'qual'))
    elif fill is None and position is not None:
        plot = (
        ggplot(df) +
        geom_bar(aes(x, y), stat = 'identity', position = position) +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_fill_brewer(palette = 'Set2', type = 'qual'))
    elif fill is not None and position is None:
        plot = (
        ggplot(df) +
        geom_bar(aes(x, y, fill = fill), stat = 'identity', position = 'dodge') +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_fill_brewer(palette = 'Set2', type = 'qual'))
    else:
        plot = (
        ggplot(df) +
        geom_bar(aes(x, y), stat = 'identity') +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_fill_brewer(palette = 'Set2', type = 'qual'))
    return plot

def bar_with_errors(df, x, y, ymin = None, ymax = None, fill = None, group = None):
    if any(thing is None for thing in (ymin,  ymax, fill, group)):
        print('Error: Incorrect arguments passed')
        print('Syntax: bar_with_errors(df, x, y, ymin, ymax, fill = value, group = value)')
        sys.exit()
    else:
        plot = (
        ggplot(df) +
        geom_bar(aes(x, y, fill = fill), stat = 'identity', position = 'dodge') +
        geom_errorbar(aes(x, ymin = ymin, ymax = ymax, group = group), position = position_dodge(0.9)) +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_fill_brewer(palette = 'Set2', type = 'qual'))

        return plot

def scatter(df, x, y, color):
    return (
        ggplot(df) +
        geom_point(aes(x, y, color = color)) +
        theme_seaborn() +
        theme(text = element_text(fontproperties = font)) +
        scale_color_gradient(low = '#1b9e77', mid = '#d95f02', high = '#7570b3'))






