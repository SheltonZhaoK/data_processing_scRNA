#Build-in function
import os.path, sys,random
from datetime import date

#Utility
import numpy as np
import pandas as pd
import scanpy as sc
'''
#FFNN
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dense, Dropout, Conv1D, Flatten, MaxPooling1D
from tensorflow.keras import optimizers
from tensorflow.keras.regularizers import l2
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OrdinalEncoder
'''
def check_file(file1):
   if not(os.path.exists(file1)):
      error("%s doesn't exist." % file1)

def error(message):
   print('There is a problem: %s' % message)
   sys.exit(1)

def read_data(inputFile):
   check_file(inputFile)
   adata = sc.read_h5ad(inputFile)   
   #adata = sc.read_10x_mtx(inputFile, var_names='gene_symbols',cache=False)
   return adata

def prepare_data(adata, min_cells, min_features, max_genes_counts):
   sc.pp.filter_cells(adata, min_genes = min_features)
   sc.pp.filter_genes(adata, min_cells = min_cells)
   
   adata.var['mt'] = adata.var_names.str.startswith('MT-')
   sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

   adata = adata[adata.obs.n_genes_by_counts < max_genes_counts, :]
   adata = adata[adata.obs.pct_counts_mt < 5, :]
   return adata

def linear_transformation(adata):
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)
   sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
   adata.raw = adata
   
   adata = adata[:, adata.var.highly_variable]
   sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
   sc.pp.scale(adata, max_value=10)

   sc.tl.pca(adata, svd_solver='arpack')
   return adata

def none_linear_transformation(adata):
   sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) 
   sc.tl.umap(adata)
   return adata

def clustering(adata):
   sc.tl.leiden(adata)
   return adata

def find_maker_genes(adata):
   sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
   marker_rank = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(15)
   print(marker_rank)
   
def format_ouput():
   print('')
   
def main(inputFile, min_cells, min_features, max_genes_counts):
   
   adata = read_data(inputFile)
   adata = prepare_data(adata, min_cells, min_features, max_genes_counts)
   adata = linear_transformation(adata)
   adata = none_linear_transformation(adata)
   adata = clustering(adata)
   find_maker_genes(adata)

   adata.write('/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/output/azizi_scanpy_data_processed.h5ad')
   
   
if __name__ == "__main__":
   dataDir = '/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/output/'
   inputFile = dataDir + 'azizi_scanpy_data_processed.h5ad'
   #inputFile = '/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/data/filtered_gene_bc_matrices/hg19/'   
   
   min_cells = 3
   min_features = 200
   max_genes_counts = 2500
   main(inputFile, min_cells, min_features, max_genes_counts)
   
   print(date.today())
   
