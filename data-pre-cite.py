import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import numpy as np
import muon

#read datasets
cite = sc.read("/home/ubuntu/data/cite_data.h5ad")

#classification
rna = cite[:, cite.var['feature_types'] == 'GEX'].copy()
adt = cite[:, cite.var['feature_types'] == 'ADT'].copy()

#rna preprocessing
rna.X = rna.layers['counts'].copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=4000, batch_key='batch')
rna_hvg = rna[:, rna.var.highly_variable].copy()

#adt preprocessing
adt.X = adt.layers['counts'].copy()
muon.prot.pp.clr(adt)
adt.layers['clr'] = adt.X.copy()

#save files
rna.write('/home/ubuntu/multimodal-metrics/multigrate/data/cite/rna_cite.h5ad')
rna_hvg.write('/home/ubuntu/multimodal-metrics/multigrate/data/cite/rna_cite_hvg.h5ad')
adt.write('/home/ubuntu/multimodal-metrics/multigrate/data/cite/adt.h5ad')