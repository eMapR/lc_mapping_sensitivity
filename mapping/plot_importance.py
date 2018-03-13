# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 09:14:48 2018

@author: shooper
"""

import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def split_name(name, index):
    split = name.split('_')
    try:
        return split[index]
    except:
        return None

def main(var_info_txt):#, sample_txt):
    
    df = pd.read_csv(var_info_txt, sep='\t', 
                     usecols=['var_name', 'importance', 'rank'], 
                     index_col='rank')
    df = df[df.var_name.apply(lambda x: len(x.split('_')) >= 3)]
    
    # Calc sum of each variable type
    df['run_index'] = df.var_name.apply(lambda x: split_name(x, index=0))
    df['date_range'] = df.var_name.apply(lambda x: split_name(x, index=1))
    df['out_index'] = df.var_name.apply(lambda x: split_name(x, index=-1))
    
    sum_run = {}
    for ind in df.run_index.unique(): 
        sum_run[ind] = df.loc[df.run_index == ind, 'importance'].sum()
    sum_run = pd.Series(sum_run).sort_values()
    
    sum_date = {}
    for ind in df.date_range.unique(): 
        sum_date[ind] = df.loc[df.date_range == ind, 'importance'].sum()
    sum_date = pd.Series(sum_date).sort_values()
    
    sum_out = {}
    for ind in df.out_index.unique(): 
        sum_out[ind] = df.loc[df.out_index == ind, 'importance'].sum()
    sum_out = pd.Series(sum_out).sort_values()
    
    # plot 
    #   make a bar chart of all 21 index/date combos
    run_indices = df.run_index.unique()
    date_ranges = df.date_range.unique()
    out_indices = [i for i in df.out_index.unique() if i != 'ts']
    
    n_indices = len(out_indices)
    bar_width = 1
    var_buffer = n_indices * bar_width/3. # distance between groups of bars
    var_width = n_indices * bar_width # width of ech group
    
    colors = sns.color_palette()[:n_indices]
    df.sort_values('var_name', inplace=True) # sort so we can match labels to groups
    c = 0 # count of date/index combos
    # plot each run_index/date combo
    for run_index in run_indices:
        for date_range in date_ranges:
            # get x axis locations for this set of bars
            x_loc = np.array([(c * (var_width + var_buffer) + i) for i in range(n_indices)])
            # select the variables from this index/date combo
            ind_vars = df.loc[
                            (df.run_index == run_index) &
                            (df.date_range == date_range) &
                            (df.var_name.apply(
                                lambda x: 'delta' not in x and 'ts' not in x))
                            ]
            bars = plt.bar(x_loc, 
                           ind_vars.importance, 
                           bar_width, 
                           color=colors,
                           edgecolor='none')
            c += 1 # increment count of date/index combos 
    
    # Plot lables
    #   get labels in format {run_index}_{date_range}
    var_labels = sorted([(r + '_' + d) for d in date_ranges for r in run_indices])
    #   center each label on each respective group of bars
    x_loc = [i * (var_buffer + var_width) - (var_buffer + var_width)/2. for i in range(len(var_labels))]
    plt.xticks(x_loc, var_labels, size='small', rotation=45, multialignment='right')
    plt.subplots_adjust(bottom=0.2)
    
    # plot legend
    plt.legend(bars, out_indices, bbox_to_anchor=(1.13, 1)) # Place legend outside
    
    # plot title
    region_id = re.search('region_\d\d', var_info_txt).group().split('_')[1]
    title = 'Variable importance for region %s' % region_id
    plt.title(title)
    
    out_dir = os.path.join(os.path.dirname(var_info_txt), 'evaluation')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    region = os.path.basename(os.path.dirname(os.path.dirname(var_info_txt))).split('_')[1]
    out_png = os.path.join(out_dir, 'importance_plot_%s.png' % region)
    plt.savefig(out_png, dpi=300)
    print'\nImportance plot saved to', out_png
    
    out_txt = out_png.replace('_plot', '')
    df.sort_index().to_csv(out_txt, sep='\t')
    
    print '\nSum of run indices:\n', sum_run.sort_values()
    print '\nSum of date indices:\n', sum_date.sort_values()
    print '\nSum of out indices:\n', sum_out.sort_values()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
    # 
    