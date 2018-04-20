import pandas as pd
import re
import sys
Input = sys.argv[1]
header = list(pd.read_csv(Input, sep = '\t', nrows=0))
header.insert(0, "Sequences")
table = pd.read_csv(Input, sep = '\t', skiprows=1, names = header).reset_index().rename(columns={'index':'OTU'})
table['OTU'] = table['OTU'].astype(str).replace('^', 'OTU_', regex = True)
OTU_fasta = table[['OTU','Sequences']]
OTU_fasta['OTU'] = OTU_fasta['OTU'].replace('^', '>', regex = True)
OTU_fasta.set_index('OTU').to_csv('OTU.fasta', sep = '\n', header = False)
table.drop(['Sequences'], axis = 1).set_index('OTU').to_csv('Dada2_OTU.tsv', sep = '\t')
