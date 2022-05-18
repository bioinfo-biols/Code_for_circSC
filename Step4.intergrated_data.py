#####4.Intergrated circRNAs and gene expression to matrix####
import pandas as pd
import glob
import os
import numpy as np
import scipy.stats as stats

"""
description='Perform the circRNAs and genes related-analysis, but we only provided the basic analysis'
"""
datadir = ''  # PATH of pipeline directory
files = []
for fn in glob.glob(datadir + '/*_bsj.txt'):
    files.append(fn)

test_factor = pd.read_table(datadir + 'brain_factor.csv', sep=',', header=0, index_col=0)
test_meta = pd.read_table(datadir + 'brain_metadata.csv', sep=',', header=0, index_col=0)
test_meta = pd.concat([test_meta, test_factor], axis=1).fillna(0)
test_meta['cell_cluster'] = [test_cell[i] + str(i) for i in test_meta['seurat_clusters']]

test_GSE = {}
tese_gene_counts = pd.DataFrame()
for fn in glob.glob(datadir + '/*_count.csv'):
    tmps = pd.read_table(fn, sep=',', header=0, index_col=0)
    tese_gene_counts = pd.concat([tese_gene_counts, tmps], axis=1).fillna(0)
    prefix = os.path.split(fn)[-1].split('_')[0]
    if prefix != 'merge':
        with open(fn, 'r') as al:
            head = al.readline().rstrip().split(',')
            for ele in head:
                if ele:
                    test_GSE[ele] = prefix
for ele in test_cluster.index:
    test_cluster.loc[ele, 'Dataset'] = test_GSE[ele]
test_bsjs = pd.DataFrame()
for ele in files:
    tmp = pd.read_table(ele, sep=',', header=0, index_col=0)
    test_bsjs = pd.concat([test_bsjs, tmp], axis=1).fillna(0)
test_bsjs = test_bsjs[list(set(test_bsjs.columns) & set(test_cluster.index))]
test_bsjs['sum'] = test_bsjs.apply(lambda x: x.sum(), axis=1)
test_bsjs = test_bsjs[test_bsjs['sum'] > 1]
del test_bsjs['sum']

test_cluster = test_cluster.loc[list(set(test_bsjs.columns) & set(test_cluster.index))]
test_cexps = pd.DataFrame()
for e in test_bsjs.columns:
    test_cexps[e] = np.log2(test_bsjs[e] * 1.0 / test_factor.loc[e, 'sizeFactor'] + 1)
test_gexps = pd.DataFrame()
for e in tese_gene_counts.columns:
    test_gexps[e] = np.log2(test_gexps[e] * 1.0 / test_factor.loc[e, 'sizeFactor'] + 1)

test_gexps.to_csv(datadir + 'brain_gexp.csv')
test_cexps.to_csv(datadir + 'brain_cexp.csv')
test_cluster.to_csv(datadir + 'brain_cmeta.csv')
test_bsjs['CIRI'] = test_bsjs.index
test_bsjs = pd.melt(test_bsjs, id_vars='CIRI')
test_bsjs = test_bsjs[test_bsjs['value'] > 0]
test_bsjs['cell'] = test_bsjs['variable'].apply(lambda x: test_cluster.loc[x, 'ident'])
cell_num = {}
for i in set(test_cluster['ident']): jn
tmp = test_cluster[test_cluster['ident'] == i]
cell_num[i] = len(tmp)
circ_por = pd.DataFrame(index=set(test_bsjs['CIRI']), columns=set(test_bsjs['cell']))
for i in set(test_bsjs['CIRI']):
    tmp = test_bsjs[test_bsjs['CIRI'] == i]
    for j in set(tmp['cell']):
        tpmT = tmp[tmp['cell'] == j]
        circ_por.loc[i, j] = len(tpmT) / cell_num[j]
circ_por = circ_por.fillna(0)
circ_por.to_csv(datadir + '/brain_circ_por.csv')
circ_por['sum'] = circ_por.apply(lambda x: x.sum(), axis=1)

circ_mean = pd.DataFrame(index=set(test_bsjs['CIRI']), columns=set(test_bsjs['cell']))
for i in set(test_bsjs['CIRI']):
    tmp = test_bsjs[test_bsjs['CIRI'] == i]
    for j in set(tmp['cell']):
        tpmT = tmp[tmp['cell'] == j]
        circ_mean.loc[i, j] = np.mean(tpmT['exp'])
circ_mean = circ_mean.fillna(0)
circ_mean['sum'] = circ_mean.apply(lambda x: x.sum(), axis=1)
circ_mean = circ_mean.sort_values('sum', ascending=False)
circ_mean.to_csv(datadir + '/brain_circ_mean.csv')


