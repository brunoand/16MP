import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

Input = sys.argv[1]
Output = sys.argv[2]

taxa_all = pd.read_csv('OTU_table_taxonomy.txt', sep = '\t', index_col=0)
taxa_all = taxa_all.drop(['Score','OTU ID'], axis=1)
taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= 's__.{1,}', value='', regex=True)
taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= '(.__|\s.{1,}|\[|\]|\-.{1,})', value='', regex=True)
taxa_all2 = pd.DataFrame(taxa_all['Taxonomy'].str.split(';').values.tolist(), columns = ['Domain', 'Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'], index = taxa_all.index)

taxa_all2 = taxa_all2.drop('Species', axis=1)
taxa_all2 = taxa_all2.merge(taxa_all, left_index=True, right_index = True).drop('Taxonomy', axis=1)
for x in ['Phylum', 'Order', 'Class', 'Family', 'Genus']:
   rel_abundance(taxa_all2, x, Output)

def rel_abundance(dataframe, tax_level, Output):
    frame = dataframe.groupby(tax_level).sum()
    Df_Rel_ab = (frame*100)/frame.sum(axis=0)
    Df_Rel_ab['Mean'] = Df_Rel_ab.mean(axis = 1)
    Df_Rel_ab = Df_Rel_ab[Df_Rel_ab.Mean >=0.5]
    Df_Rel_ab = Df_Rel_ab.sort_values(by=['Mean'], ascending = False).drop(['Mean'], axis=1)
    sns.set(rc={'figure.figsize':(15,10)},font_scale = 2)
    sns.set_style("whitegrid", {'axes.grid' : False})
    barplot = Df_Rel_ab.transpose().plot(kind='bar', stacked=True, cmap="tab20b", edgecolor='black', linewidth=1)
   # plt.gca().invert_yaxis()
    barplot.legend(loc=9, bbox_to_anchor=(1.12, 1), prop={'size': 15})
    barplot.axes.set_title("Relative abundance at the {} level(>0.5% of abundance)".format(tax_level),fontsize=22)
    barplot.set_xlabel("Percentage",fontsize=20)
    barplot.set_ylabel("Samples",fontsize=20)
    plt.xticks(fontsize=14, rotation = 45)
    plt.yticks(fontsize=14)
    barplot.figure.savefig(Output + 'Relative_abundance_{}.png'.format(tax_level), bbox_inches='tight')