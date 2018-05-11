#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg') 
import sys



def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j, regex=True)
    return text


def parse(Dataframe, Output, Percentage):
    taxa_all = Dataframe.reset_index()
    taxa_all = taxa_all.drop(['Score','OTU'], axis=1)
    taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= 's__.{1,}', value='', regex=True)
    taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= '(.__|\s.{1,}|\[|\]|\-.{1,})', value='', regex=True)
    taxa_all2 = pd.DataFrame(taxa_all['Taxonomy'].str.split(';').values.tolist(), columns = ['Domain', 'Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'], index = taxa_all.index)

    taxa_all2 = taxa_all2.drop('Species', axis=1)
    taxa_all2 = taxa_all2.merge(taxa_all, left_index=True, right_index = True).drop('Taxonomy', axis=1)
    for x in ['Phylum', 'Order', 'Class', 'Family', 'Genus']:
        rel_abundance(taxa_all2, x, Output, Percentage)


def rel_abundance(dataframe, tax_level, Output, Percentage):
    frame = dataframe.groupby(tax_level).sum()
    Df_Rel_ab = (frame*100)/frame.sum(axis=0)
    Df_Rel_ab['Mean'] = Df_Rel_ab.mean(axis = 1)
    Df_Rel_ab = Df_Rel_ab[Df_Rel_ab.Mean >=Percentage]
    Df_Rel_ab = Df_Rel_ab.sort_values(by=['Mean'], ascending = False).drop(['Mean'], axis=1)
    sns.set(rc={'figure.figsize':(15,10)},font_scale = 2)
    sns.set_style("whitegrid", {'axes.grid' : False})
    barplot = Df_Rel_ab.transpose().plot(kind='bar', stacked=True, cmap="tab20b", edgecolor='black', linewidth=1)
    barplot.legend(loc=9, bbox_to_anchor=(1.12, 1), prop={'size': 15})
    barplot.axes.set_title("Relative abundance at the {} level(>0.5% of abundance)".format(tax_level),fontsize=22)
    barplot.set_xlabel("Percentage",fontsize=20)
    barplot.set_ylabel("Samples",fontsize=20)
    plt.xticks(fontsize=18, rotation = 45)
    plt.yticks(fontsize=14)
    barplot.figure.savefig(Output + 'Relative_abundance_{}.png'.format(tax_level), bbox_inches='tight')
    Df_Rel_ab.to_csv('Relative abundance' + tax_level + '.tsv', sep = '\t')
    frame.to_csv('Counts' + tax_level + '.tsv', sep = '\t')

Input_OTU = sys.argv[1]
Taxonomy = sys.argv[2]
Output = sys.argv[3]
Output_plots = sys.argv[4]
Percentage = sys.argv[5]


OTU_table = pd.read_csv(Input_OTU, sep = '\t')
header_taxonomy = ['OTU', 'Taxonomy', 'Score', 'Size']
Taxonomy = pd.read_csv(Taxonomy, sep = '\t', names = header_taxonomy)
Pattern = {"D_0":"d", "D_1":"p", "D_2":"c", "D_3":"o", "D_4":"f", "D_5":"g", "D_6":"s" }

Taxonomy['Taxonomy']= replace_all(Taxonomy['Taxonomy'], Pattern)
Taxonomy = Taxonomy.merge(OTU_table, right_on='OTU', left_on = 'OTU').drop(['Size'], axis=1).set_index('OTU')
Taxonomy.to_csv(Output, sep = '\t')
plot = parse(Taxonomy, Output_plots)
