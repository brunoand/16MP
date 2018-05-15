#!/usr/bin/python3

import skbio
import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import warnings
import sys

Input = sys.argv[1]
Meta = sys.argv[2]
Output = sys.argv[3]
warnings.filterwarnings("ignore")
metadata = pd.read_csv(Meta, sep = '\t', index_col = 0)
my_obj =  DistanceMatrix.read(Input, 'lsmat')
PC = pcoa(my_obj)



def plot_PCoA(matrix, ID_column):
    figure = matrix.plot(metadata, ID_column, axis_labels=('PC 1', 'PC 2', 'PC 3'), cmap='jet', s=50)
    figure.set_size_inches(12.5, 8.5)
    figure.text(0,0.9, r'Samples colored by {}'.format(ID_column), fontsize=16)
    figure.savefig(Output + 'PCOA_{}.pdf'.format(ID_column), bbox_inches='tight')
    
for x in metadata.columns:
    plot_PCoA(PC, x)
