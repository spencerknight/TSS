import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm, tqdm_pandas
import time
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from matplotlib import rc, font_manager
import subprocess
import os

def trapezoid(time_series, value_series):
    '''
    Finds the trapezoidal area across a time/value series (e.g., Trapezoid rule)
    '''
    baseline = value_series[0]
    answer = 0
    for i in range(len(time_series)-1):
        answer += (time_series[i+1]-time_series[i])*((value_series[i] + value_series[i+1])/2 - baseline)
    return answer

def sloper(time_series, value_series):
    '''
    Finds the slope from the lefthand direction for a time/value series (e.g., lefthand derivative approximation)
    '''
    answer = []
    for i in range(len(time_series)-1):
        answer.append((value_series[i+1]-value_series[i])/(time_series[i+1]-time_series[i]))

def normalizer(arr):
    '''
    Normalizes a read array from 0 to 1. arr is a list or numpy array
    '''
    return [(ix-min(arr))/(max(arr)-min(arr)) for ix in arr]

def plotter3():
    r = pd.read_pickle('marko_new_years.pkl')
    r = r[r['anticorr_trapezoid?'] == True]
    #r = r.sort_values(by = 'slope_slope_r', ascending = False)
    r.index = range(len(r))
    t_rna = [0.0,1.0,2.0,4.0,6.0,9.0,12.0]
    t_ribo = [0.0,0.5,1.0,2.0,4.0,6.0,8.0,9.0,12.0]
    rnas = r['0_1_rna_array'].tolist()
    ribos = r['0_1_ribo_array'].tolist()
    genes = r['Gene'].tolist()
    translation_type = r['Translation_annotation'].tolist()
    for i,rna in enumerate(tqdm(rnas)):
        times = []
        values = []
        hues = []
        times += t_ribo + t_rna
        values += list(ribos[i]) + list(rna)
        hues += ['{} Footprint Read'.format(genes[i]) for _ in range(len(t_ribo))] + ['{} RNA-Seq Read'.format(genes[i]) for _ in range(len(t_rna))]
        df = pd.DataFrame()
        df['Time'] = times
        df['Normalized_Read'] = values
        df[''] = hues
        df0 = df[df[''] == '{} Footprint Read'.format(genes[i])]
        df0['Footprint Read'] = df0['Normalized_Read']
        df1 = df[df[''] == '{} RNA-Seq Read'.format(genes[i])]
        df1['RNA-Seq Read'] = df1['Normalized_Read']
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_ylim([-0.04, 1.04])
        plt.plot(df1['Time'], df1['RNA-Seq Read'])
        plt.plot(df0['Time'], df0['Footprint Read'])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        ax.set_xlabel('Time, h')
        ax.set_ylabel('Normalized Read')
        timestr = time.strftime("%Y%m%d-%H%M%S")
        plt.title('RNA-Seq v. Footprinting Reads for {}'.format(genes[i] + '_' + translation_type[i]))
        plt.savefig('new_ts_plots/regression/Time_series_{}.pdf'.format(timestr + '_' + genes[i] + '_' + translation_type[i]))

def venner(col_a, col_b, col_c, data):
    '''
    Overlap values for three columns for making Venn diagram
    '''
    a_only = len(data[data[col_a] == True][data[col_b] == False][data[col_c] == False])
    b_only = len(data[data[col_b] == True][data[col_a] == False][data[col_c] == False])
    c_only = len(data[data[col_c] == True][data[col_b] == False][data[col_a] == False])
    a_b_only = len(data[data[col_a] == True][data[col_b] == True][data[col_c] == False])
    a_c_only = len(data[data[col_a] == True][data[col_c] == True][data[col_b] == False])
    b_c_only = len(data[data[col_b] == True][data[col_c] == True][data[col_a] == False])
    a_b_c = len(data[data[col_a] == True][data[col_b] == True][data[col_c] == True])

    title_font = {'fontname':'Letter Gothic Std', 'size':'20', 'color':'black', 'weight':'bold',
    'verticalalignment':'bottom'}
    axis_font = {'fontname':'Letter Gothic Std', 'size':'18', 'weight': 'bold'}

    v = venn3(subsets=(a_only, b_only, a_b_only, c_only, a_c_only, b_c_only, a_b_c), set_labels = (col_a, col_b, col_c))
    c = venn3_circles(subsets=(a_only, b_only, a_b_only, c_only, a_c_only, b_c_only, a_b_c), linestyle = 'solid')
    num_labels = ['100', '101', '110', '111', '010', '011', '001']
    lett_labels = ['A', 'B', 'C']
    for label in num_labels:
        label = v.get_label_by_id(label)          # Those are subset labels (i.e. numbers)
        label.set_family('Letter Gothic Std')

    for label in lett_labels:
        label = v.get_label_by_id(label)          # Those are subset labels (i.e. numbers)
        label.set_fontsize(16)
        label.set_family('Letter Gothic Std')
        if label == 'C':
            x,y = label.get_position()
            label.set_position((x + 10, y - 0.5))

    plt.title('Anticorrelation Metric Overlap', **title_font)
    plt.savefig('hits_overlap.png', dpi = 600)

def venner2():

    title_font = {'fontname':'Letter Gothic Std', 'size':'20', 'color':'black', 'weight':'bold',
    'verticalalignment':'bottom'}
    axis_font = {'fontname':'Letter Gothic Std', 'size':'18', 'weight': 'bold'}

    v = venn2(subsets=(607, 797, 789), set_labels = ('Pre-removal', 'Post-removal'))
    c = venn2_circles(subsets=(607, 797, 789), linestyle = 'solid')
    num_labels = ['11', '10', '01']
    lett_labels = ['A', 'B']

    for label in num_labels:
        label = v.get_label_by_id(label)          # Those are subset labels (i.e. numbers)
        label.set_family('Letter Gothic Std')

    for label in lett_labels:
        label = v.get_label_by_id(label)          # Those are subset labels (i.e. numbers)
        label.set_fontsize(16)
        label.set_family('Letter Gothic Std')

    plt.title('Overlap pre- v. post-UCXXX removal', **title_font)
    plt.savefig('derivative_uc_overlap.png', dpi = 600)

def merger():
    translation_types = ['Annotated', 'Downstream', 'Extension', 'Internal', 'Isoform', 'New', 'Truncation', 'Upstream']
    os.chdir('new_ts_plots/trapezoid')
    for t in tqdm(translation_types):
        subprocess.call('pdftk *{}.pdf'.format(t) + ' cat output merged_{}.pdf'.format(t), shell=True)

if __name__ == "__main__":
    #df = pd.read_pickle('marko_new_years.pkl')
    #df['Regression Slope'] = df['anticorr_slope?']
    #df['Trapezoid Rule'] = df['anticorr_trapezoid?']
    #df['Derivative Anticorrelation'] = df['anticorr_slope_slope?']
    #venner('Regression Slope', 'Trapezoid Rule', 'Derivative Anticorrelation', data=df)
    #plotter3()
    merger()
    #venner2()
    print 'DONE!'
