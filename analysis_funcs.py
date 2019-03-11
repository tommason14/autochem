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

def lollipop(df, x, xend, y, yend, color):
    return (
    ggplot(df) + 
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = color)) +
    geom_point(aes(x = xend, y = yend, color = color)) +
    theme_seaborn() +
    theme(text = element_text(fontproperties = font)) +
    scale_color_brewer(type = 'qual', palette = 'Dark2'))


##### DFPLY FUNCTIONS ####

# --- Image rendering for high res screens ---
%config InlineBackend.figure_format = 'retina'


# --- Helper Functions ---
def rename_cations(series, d = None):
    if d is not None:
        return series.replace(d)
    return series.replace({'c4mim': 'C$_{\\mathrm{4}}$mim',
                           'choline': 'Choline'})

def rename_anions(series, d = None):
    if d is not None:
        return series.replace(d)
    return series.replace({'acetate': 'Acetate',
                           'dhp': 'DHP',
                           'mes': 'Mesylate'})

@pipe
def rename_ions(df):
    vals = {'Cation': {'ch': 'Choline', 'c4mim': 'C$_{\\mathrm{4}}$mim'},
            'Anion': {'ac': 'Acetate', 'dhp': 'DHP', 'mes': 'Mesylate'}
           }

    def conditional_mutate(df, value, col, new_col, new_value):
        df.loc[df[col].str.contains(value), new_col] = new_value
        return df

    for new_col, d in vals.items():
        for string, new_value in d.items():
            df = conditional_mutate(df, string, new_col, new_col, new_value)
    return df


@pipe
def make_IL_column(df, lst):
    """
    Uses a list of dictionaries to define the ionic liquids being used.
    i.e.
      lst = [{'Cation': 'ch', 'Anion': 'ac', 'IL': '[ch][ac]'},
             {'Cation': 'ch', 'Anion': 'dhp', 'IL': '[ch][dhp]'},
             {'Cation': 'ch', 'Anion': 'mes', 'IL': '[ch][mes]'},
             {'Cation': 'c$_{\\mathrm{4}}$mim', 'Anion': 'ac', 'IL': '[C$_{\\mathrm{4}}$mim][ac]'},
             {'Cation': 'c$_{\\mathrm{4}}$mim', 'Anion': 'dhp', 'IL': '[C$_{\\mathrm{4}}$mim][dhp]'},
             {'Cation': 'c$_{\\mathrm{4}}$mim', 'Anion': 'mes', 'IL': '[C$_{\\mathrm{4}}$mim][mes]'}]
    """
    
    def add_il(row):
      
        for item in lst:
            if item['Cation'] in row['Cation'].lower() and item['Anion'] in row['Anion'].lower():
                return item['IL']

    df = df.assign(IL = df.apply(add_il, axis = 1))
    return df


@make_symbolic
def confidence(column):
    """
    95% confidence intervals defined as:
        1.96 * standard deviation from the mean / sqrt(number of items)
        
    1.96 assumes a normal distribution.
    
    https://www.itl.nist.gov/div898/handbook/prc/section1/prc14.htm
    http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html
    """
    return 1.96 * sd(column) * (n(column) ** -0.5) 
    
@make_symbolic
def order_column(series):
    return pd.Categorical(series, categories = series.unique(), ordered = True)

@make_symbolic
def bp(series, as_percent = False):
    """
    Takes in energies in kJ/mol, produces 
    probabilities according to a Boltzmann distribution
    """
    R = 8.3145
    T = 298.15
    exponent = np.exp((-1 * series * 1000) / (R * T))
    summed = exponent.sum()
    if as_percent:
        return (exponent / summed) * 100
    return exponent / summed


clear_x_axis = theme(axis_text_x = element_blank(),
                     axis_ticks_major_x = element_blank(),
                     axis_title_x = element_blank())
clear_y_axis = theme(axis_text_y = element_blank(),
                     axis_ticks_major_y = element_blank(),
                     axis_title_y = element_blank())

# Useful regex!


# df = df >> mutate(Cation = X.Config.str.extract('(^[^_]+(?=_))'), 
#                   Anion = X.Config.str.extract('((?<=_)[^_]+(?=_))'))

@pipe
def ions_from_filename(df, series):
    return df >> mutate(Cation = series.str.extract('(^[^_]+(?=_))'),
                         Anion = series.str.extract('((?<=_)[^_]+(?=_))'))
                    
