# coding: utf-8
import pandas as pd
import numpy as np
import sys
def compute(cell, frac):
    print(cell)
    a = pd.read_csv(cell + '/output-epw/predictions-gbm.csv', header=0)
    a = a.set_index(['promoter_name', 'enhancer_name'])
    print(sum(a['label']), sum(a['prediction']),
            sum(a.ix[a['label']==1,'prediction'] == a.ix[a['label']==1, 'label']),
            sum(a['prediction']) - sum(a.ix[a['label']==1,'prediction'] == a.ix[a['label']==1, 'label']))
    b = pd.read_table(cell+'/output-epw/temp_pos_0.%d.bed' % frac, header=None)
    b.columns= ['chrom', 'start', 'end', 'promoter_name', 'enhancer_name']
    b = b.set_index(['promoter_name', 'enhancer_name'])
    c = a.index.intersection(b.index)
    print(len(c), np.sum(a.ix[c, :]))
    b = pd.read_table(cell+'/output-epw/temp_0.%d.bed' % frac, header=None)
    b.columns= ['chrom', 'start', 'end', 'promoter_name', 'enhancer_name']
    b = b.set_index(['promoter_name', 'enhancer_name'])
    c = a.index.intersection(b.index)
    print(len(c), np.sum(a.ix[c, :]))
    print('='*20)
    print()

for cell in ['K562', 'GM12878', 'HeLa-S3', 'IMR90', 'NHEK', 'HUVEC']:
    compute(cell, int(sys.argv[1]))
