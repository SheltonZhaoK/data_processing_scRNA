import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#file that will store the analysis results
results_file = 'python_result.h5ad' 

adata = sc.read_10x_mtx(
    '/deac/csc/khuriGrp/zhaok220/data_processing_scRNA/data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)

adata.var_names_make_unique()

#genes that yield the highest fraction of counts in each single cell
#across all cells
sc.pl.highest_expr_genes(adata, n_top=20,save='python_plots.pdf')

#filtered out genes that are detected in less than 3 cells
#and cells that dedected less than 200 cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#quality control
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#violin plot of the number of genes expressed in the count matrix per cell
#the total counts per cell
#the percentage of counts in mitochondrial genes
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='python_plots.pdf' )

#filter cells that have unique feature counts over 2,500 or less than 200
#filter cells that have >5% mitochondrial counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]


#Total-count normalize (library-size correct) the data matrix 
#to 10,000 reads per cell and logrithmize to prepare for pca
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='python_plots.pdf')

adata.raw = adata

#filter highly-variable genes, regress out total counts and mt cells, and scale
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

#pca
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save='python_plots.pdf')

adata.obsm.to_df()[['X_pca1', 'X_pca2','X_pca3','X_pca4','X_pca5','X_pca6','X_pca7'
                    ,'X_pca8','X_pca9','X_pca10']].to_csv('pbmc3k_corrected_X_pca.csv')

#compute and embed the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'],save='python_plots.pdf', use_raw=False)

#clustering the neighborhood graph/non-linear dimensinality reduction
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])

#find marker genee in each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='python_plots.pdf')

new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False,save='python_plots.pdf')



