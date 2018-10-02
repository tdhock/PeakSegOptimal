import numpy as np
import matplotlib.pyplot as plt

def arfpop_stats(fit):
    df = {'penalty': [fit['penalty']],
          'changepoints_n': [len(fit['spikes'])],
          'cost': [fit['cost'][len(fit['cost']) - 1] - len(fit['spikes']) * fit['penalty']]}

    return df

def update_path_stats(df, new_fit):
    new_values = arfpop_stats(new_fit)
    if len(df) == 0:
        df = new_values
    else:
        df['penalty'].append(new_values['penalty'][0])
        df['changepoints_n'].append(new_values['changepoints_n'][0])
        df['cost'].append(new_values['cost'][0])
    return df

def get_num_changepts(penalty, df):
    return df['changepoints_n'][df['penalty'].index(penalty)]


def get_cost(penalty, df):
    return df['cost'][df['penalty'].index(penalty)]

def plot_penalty_path(fits):
    stats = fits['path_stats']
    inds = np.argsort(stats['penalty'])
    x = [stats['penalty'][i] for i in inds]
    y = [stats['changepoints_n'][i] for i in inds]
    plt.figure()
    plt.plot(x, y)