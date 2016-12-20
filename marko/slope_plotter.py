import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def opposites(a, b):
    if a < 0:
        if b > 0:
            return True
    elif a > 0:
        if b < 0:
            return True
    else:
        return False

def plotter1():
    df = pd.read_pickle('rna_ribo_comp.pkl')
    df = df[df['nn_ribo_r'] > 0.5]
    df = df[df['nn_rna_r'] > 0.5]
    fig, ax = plt.subplots()
    plt.title('Normalized Slope for Ribosome Footprinting v. RNA-Seq')
    sns.jointplot(x = 'nn_rna_slope', y = 'nn_ribo_slope', data = df)
    ax.set_xlabel('RNA-Seq slope')
    ax.set_ylabel('Ribosome Footprint slope')
    plt.savefig('RNARibo_slope_nn_r05.png', dpi = 300)

def plotter2():
    df = pd.read_pickle('rna_ribo_comp.pkl')
    tp = [1,2,4,6,9,12]
    t_len = len(tp)
    for i,t in enumerate(tp):
        dx = pd.DataFrame()
        index = t_len-i
        coi = 'rna_delta_{}h_nn'.format(t)
        x,y, hues = [], [], []
        for j in range (index-1):
            x += df[coi].tolist()
        print 'X_length: {}'.format(len(x))
        dx['RNA-Seq Delta 0-{}h'.format(t)] = x
        for k in range(i+1, t_len):
            y += df['ribo_delta_{}h_nn'.format(tp[k])].tolist()
            hues += ['delta_footprint_{}h'.format(tp[k]) for _ in range(len(df))]
        dx['Ribosome Footprint Delta'] = y
        dx[''] = hues
        fig,ax = plt.subplots()
        plt.title('Delta Ribosome Footprint v. RNA-Seq')
        ax = sns.lmplot(x = 'RNA-Seq Delta 0-{}h'.format(t), y = 'Ribosome Footprint Delta', data = dx, hue = '', fit_reg = False, legend = False)
        plt.legend(loc = 'upper left')
        plt.savefig('Delta_plot_{}_nn'.format(t), dpi = 300)

def plotter3():
    df = pd.read_pickle('rna_ribo_comp.pkl')
    df['truth'] = df.apply(lambda row: opposites(row['nn_rna_slope'], row['nn_ribo_slope']), axis = 1)
    df = df[df['truth'] == True]
    fig, ax = plt.subplots()
    plt.title('Normalized Slope for Ribosome Footprinting v. RNA-Seq')
    sns.jointplot(x = 'nn_rna_slope', y = 'nn_ribo_slope', data = df)
    ax.set_xlabel('RNA-Seq slope')
    ax.set_ylabel('Ribosome Footprint slope')
    plt.savefig('RNARibo_slope_nn_opps.png', dpi = 300)

if __name__ == "__main__":
    df = pd.read_pickle('rna_ribo_comp.pkl')
    plotter3()
