import re
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy

'''
Parse arguments.
Input file is required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Cluster cell types using mcmicro marker expression data.')
    parser.add_argument('-i', '--input', help="Input CSV of mcmicro marker expression data for cells", type=str, required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved', type=str, required=False)
    parser.add_argument('-m', '--markers', help='A text file with a marker on each line to specify which markers to use for clustering', type=str, required=False)
    parser.add_argument('-k', '--neighbors', help='the number of nearest neighbors to use when clustering. The default is 30.', default=30, type=int, required=False)
    parser.add_argument('-c', '--method', help='Include a column with the method name in the output files.', action="store_true", required=False)
    parser.add_argument('-y', '--config', help='A yaml config file that states whether the input data should be log/logicle transformed.', type=str, required=False)
    parser.add_argument('--force-transform', help='Log transform the input data. If omitted, and --no-transform is omitted, log transform is only performed if the max value in the input data is >1000.', action='store_true', required=False)
    parser.add_argument('--no-transform', help='Do not perform Log transformation on the input data. If omitted, and --force-transform is omitted, log transform is only performed if the max value in the input data is >1000.', action='store_true', required=False)
    args = parser.parse_args()
    return args


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Get markers to use for clustering from a text file where each marker is on a line and corresponds exactly to the column name in the input data file.
Returns a list of markers to use for clustering.
'''
def get_markers(markers_file):
    markers = [] # list of markers in file

    # read markers from file
    f = open(markers_file, 'r')
    for line in f:
        markers.append(line.strip())

    return markers


'''
Clean data in input file.
NOTE: Currently we are doing this with pandas however, using csv might be faster, or numpy.

Exclude the following data from clustering:
    - X_centroid, …, Extent, Orientation - morphological features
    - Any of the following DNA stains
        - DNA0, DNA1, DNA2, …
        - Hoechst0, Hoechst1, ....
        - DAPI0, DAPI1, …
    - AF* (e.g., AF488, AF555, etc.) - autofluorescence
    - A* (e.g., A488, A555, etc.) - secondary antibody staining only

To include any of these markers in the clustering, provide their exact names in a file passed in with the '-m' flag
'''
def clean(input_file):

    # constants

    # a default list of features to exclude from clustering
    # we do want to also save these data separately as well. modify below to save 
    FEATURES_TO_REMOVE = ['X_centroid', 'Y_centroid', # morphological features
                        'column_centroid', 'row_centroid', 
                        'Area', 'MajorAxisLength', 
                        'MinorAxisLength', 'Eccentricity', 
                        'Solidity', 'Extent', 'Orientation', 
                        'DNA.*', 'Hoechst.*', 'DAP.*', # DNA stain
                        'AF.*', # autofluorescence
                        'A\d{3}.*'] # secondary antibody staining only (iy has to have 3 digist after)

    # load csv
    data = pd.read_csv(input_file)

    # if markers provided, keep only those features and the Cell IDs. It is important that the CellID column is first.
    if args.markers:
        if CELL_ID not in markers: # add cell ID to list of columns to keep
            markers.insert(0, CELL_ID)
        elif markers.index(CELL_ID) != 0: # if cell ID column is included but not first, move it to the front
            markers.insert(0, markers.pop(markers.index(CELL_ID)))
        meta_data = data.drop(columns = markers, axis = 1).copy() # keep meta data separate
        data = data[markers]
    else:
        # find any columns in the input csv that should be excluded from clustering be default
        # NOTE: may want to replace this with regex, it might be faster.
        col_to_remove = []
        cols = data.columns
        for feature in FEATURES_TO_REMOVE:
            r = re.compile(feature)
            col_to_remove.extend(list(filter(r.match, cols)))
        
        # save meta data that we need.
        meta_data = data[col_to_remove].copy()
        # drop all columns that should be excluded
        data = data.drop(columns=col_to_remove, axis=1)

    # save cleaned data to csv
    data.to_csv(f'{output}/{clean_data_file}', index=False)
    # save meta data separately
    meta_data.to_csv(f'{output}/{meta_data_file}', index=False)


'''
Write CELLS_FILE from leidenCluster() adata
'''
def writeCells(adata):
    cells = pd.DataFrame(adata.obs[CELL_ID].astype(int)) # extract cell IDs to dataframe
    cells[CLUSTER] = adata.obs[LEIDEN] # extract and add cluster assignments to cells dataframe

    # add in method column if requested
    if args.method:
        cells[METHOD] = SCANPY

    cells.to_csv(f'{output}/{cells_file}', index=False)


'''
Write CLUSTERS_FILE from leidenCluster() adata
'''
def writeClusters(adata):
    clusters = pd.DataFrame(columns=adata.var_names, index=adata.obs[LEIDEN].cat.categories)   
    clusters.index.name = CLUSTER # name indices as cluster column
    for cluster in adata.obs.leiden.cat.categories: # this assumes that LEIDEN = 'leiden' if the name is changed, replace it for 'leiden' in this line
        clusters.loc[cluster] = adata[adata.obs[LEIDEN].isin([cluster]),:].X.mean(0)
    
    # add in method column if requested
    if args.method:
        clusters[METHOD] = SCANPY

    clusters.to_csv(f'{output}/{clusters_file}')

'''
Write adata into h5ad
'''
def writeH5ad(adata):
    adata.write_h5ad(f'{output}/{h5ad_file}')

'''
Get max value in dataframe.
'''
def getMax(df):
    return max([n for n in df.max(axis = 0)])

'''
normalization function taken from xyz_utils
'''
def clr_normalize_each_cell(adata, inplace=True):
    """
    Normalize count vector for each cell, i.e. for each row of .X
    """

    def seurat_clr(x):
        # TODO: support sparseness
        s = np.sum(np.log1p(x[x > 0]))
        exp = np.exp(s / len(x))
        return np.log1p(x / exp)

    if not inplace:
        adata = adata.copy()

    # apply to dense or sparse matrix, along axis. returns dense matrix
    adata.X = np.apply_along_axis(
        seurat_clr, 0, (adata.X.A if scipy.sparse.issparse(adata.X) else adata.X) # seurat tutorial suggested normalization across cells, so I think the axis should be 0 instead of 1
    )
    return adata

'''
Cluster data using the Leiden algorithm via scanpy
'''
def leidenCluster():

    sc.settings.verbosity = 3 # print out information
    adata_init = sc.read(f'{output}/{clean_data_file}', cache=True) # load in clean data

    # move CellID info into .obs
    # this assumes that 'CELL_ID' is the first column in the csv
    adata_init.obs[CELL_ID] = adata_init.X[:,0]
    adata = ad.AnnData(np.delete(adata_init.X, 0, 1), obs=adata_init.obs, var=adata_init.var.drop([CELL_ID]))

    # log transform the data according to parameter. If 'auto,' transform only if the max value >1000. Don't do anything if transform == 'false'. Write transform decision to yaml file.
    if transform == 'true':
        sc.pp.log1p(adata, base=10)
        writeConfig(True)
    elif transform == 'auto' and getMax(adata.X) > 1000:
        sc.pp.log1p(adata, base=10)
        writeConfig(True)
    else:
        writeConfig(False)

    # compute neighbors and cluster
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=10) # compute neighbors, using the first 10 principle components and the number of neighbors provided in the command line. Default is 30.
    sc.tl.leiden(adata, key_added = LEIDEN, resolution=1.0) # run leidan clustering. default resolution is 1.0

    # write cell/cluster information to 'CELLS_FILE'
    writeCells(adata)

    # write cluster mean feature expression to 'CLUSTERS_FILE'
    writeClusters(adata)



'''
Cluster data with our preferred way of clustering.
'''

def louvainCluster():

    sc.settings.verbosity = 3 # print out information
    adata_init = sc.read(f'{output}/{clean_data_file}', cache=True) # load in clean data
    meta_data = pd.read_csv(f'{output}/{meta_data_file}') # read in meta data

    # make coordinates ready for addition later
    coordinates = meta_data.loc[:,['X_centroid','Y_centroid']]
    coordinates.columns = ['X','Y']

    # move CellID info into .obs
    # this assumes that 'CELL_ID' is the first column in the csv
    # I think this way of creating anndata is giving us trouble. 
    adata_init.obs[CELL_ID] = adata_init.X[:,0]
    adata = ad.AnnData(np.delete(adata_init.X, 0, 1), obs=adata_init.obs, var=adata_init.var.drop([CELL_ID]))

    # add meta data and spatial coordinates back to adata
    adata.obs = pd.concat([adata.obs, meta_data.set_index(adata.obs.index)], axis=1)
    adata.obsm = {"spatial": coordinates.loc[:,['X','Y']].to_numpy()}

    # save untransformed data in the backup slots
    adata.raw = adata # store raw data to use for scaling for multiple different methods.
    adata.layers['raw'] = adata.X

    #transform the data. We're just going to assume tranformation is true
    adata = clr_normalize_each_cell(adata)
    adata.layers['clr'] = adata.X
    writeConfig(True)

    # compute PCA
    sc.pp.pca(adata)

    # compute neighbors and cluster
    # using the first 15 components since we have generally have more markers
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs = 15, use_rep  = 'X_pca') 
    sc.tl.umap(adata) # compute UMAP
    sc.tl.louvain(adata, resolution=0.6) # run louvain clustering. default resolution is 1.0

    # write cell/cluster information to 'CELLS_FILE'
    # writeCells(adata)

    # write cluster mean feature expression to 'CLUSTERS_FILE'
    # writeClusters(adata)

    # write the h5ad file
    writeH5ad(adata)


'''
Write to a yaml file whether the data was transformed or not.
'''
def writeConfig(transformed):
    os.mkdir('qc')
    with open('qc/config.yml', 'a') as f:
        f.write('---\n')
        if transformed:
            f.write('transform: true')
        else:
            f.write('transform: false')


'''
Read config.yml file contents.
'''
def readConfig(file):
    f = open(file, 'r')
    lines = f.readlines()

    # find line with 'transform:' in it
    for l in lines:
        if 'transform:' in l.strip():
            transform = l.split(':')[-1].strip() # get last value after colon

    return transform


'''
Main.
'''
if __name__ == '__main__':
    args = parseArgs() # parse arguments

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # get list of markers if provided
    if args.markers is not None:
        markers = get_markers(args.markers)

    # assess log transform parameter
    if args.force_transform and not args.no_transform:
        transform = 'true'
    elif not args.force_transform and args.no_transform:
        transform = 'false'
    elif args.config is not None:
        transform = readConfig(args.config)
    else:
        transform = 'auto'

    # constants
    CELL_ID = 'CellID' # column name holding cell IDs
    CLUSTER = 'Cluster' # column name holding cluster number
    LEIDEN = 'leiden' # obs name for cluster assignment
    METHOD = 'Method' # name of column containing the method for clustering
    SCANPY = 'Scanpy' # the name of this method
    
    # output file names
    data_prefix = getDataName(args.input) # get the name of the input data file to add as a prefix to the output file names
    clean_data_file = f'{data_prefix}-clean.csv' # name of output cleaned data CSV file
    meta_data_file = f'{data_prefix}-meta_data.csv' # name of output meta data csv file
    clusters_file = f'{data_prefix}-clusters.csv' # name of output CSV file that contains the mean expression of each feaute, for each cluster
    cells_file = f'{data_prefix}-cells.csv' # name of output CSV file that contains each cell ID and it's cluster assignation
    h5ad_file = f'{data_prefix}-scanpy.h5ad'

    # clean input data file
    clean(args.input)

    # cluster using scanpy implementation of Leiden algorithm
    # leidenCluster()

    # cluster using Louvain, then output the h5ad file
    louvainCluster()

