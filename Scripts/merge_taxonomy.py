#!/usr/bin/env python

import pandas as pd
import sys

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j, regex=True)
    return text


Input_OTU = sys.argv[1]
Taxonomy = sys.argv[2]
Output = sys.argv[3]

OTU_table = pd.read_csv(Input_OTU, sep = '\t')
header_taxonomy = ['OTU', 'Taxonomy', 'Score', 'Size']
Taxonomy = pd.read_csv(Taxonomy, sep = '\t', names = header_taxonomy)
Pattern = {"D_0":"d", "D_1":"p", "D_2":"c", "D_3":"o", "D_4":"f", "D_5":"g", "D_6":"s" }

Taxonomy['Taxonomy']= replace_all(Taxonomy['Taxonomy'], Pattern)
Taxonomy.merge(OTU_table, right_on='OTU', left_on = 'OTU').drop(['Size'], axis=1).set_index('OTU').to_csv(Output, sep = '\t')
