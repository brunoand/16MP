#!/usr/bin/env python
import pandas as pd
import re
import sys
Input = sys.argv[1]
Output = sys.argv[2]
Fasta = sys.argv[3]
header = list(pd.read_csv(Input, sep = '\t', nrows=0))
header.insert(0, "Sequences")
table = pd.read_csv(Input, sep = '\t', skiprows=1, names = header).reset_index().rename(columns={'index':'OTU'})
table['OTU'] = table['OTU'].astype(str).replace('^', 'OTU_', regex = True)
table2 = table
table2['Header'] = table['OTU'].replace('^', '>', regex = True)
OTU_fasta = table2[['Header','Sequences']]
OTU_fasta.set_index('Header').to_csv(Fasta, sep = '\n', header = False)
table.drop(['Sequences', 'Header'], axis = 1).set_index('OTU').to_csv(Output, sep = '\t')
