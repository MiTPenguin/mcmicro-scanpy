# mcmicro-scanpy
An MCMICRO module implementation of scanpy for clustering cells using the Leiden algorithm.

Example usage:
```
docker run --rm -v "$PWD":/data labsyspharm/mc-scanpy:1.0.1 python3 /app/cluster.py -i /data/unmicst-exemplar-001.csv -o /data/ -c
```

## TODO:
1. ~~Change the output to a scanpy h5ad object instead of cell lists~~
2. ~~Include method for transform using CLR instead of log~~
3. ~~Publish to dockerhub (mitpenguin/mc_scanpy:latest)~~
4. Fix writecell and writecluster to allow different cluster inputs
5. Specify version of scanpy and pandas package via requirements txt
6. Add sample tagging into the meta data for easy h5ad combination
7. Make graphical output for quick QC (umap cluster, etc.)

~~We probably need to test this just by the image first. so we can get a whole run through of the other methods first.~~

### modification to cluster.py

1. ~~Add methodology for writing h5ad~~
2. ~~Add methodology for custom clustering~~
3. ~~When it makes h5ad files, xy coordinates need to be added in.~~
~~follow example in immunai-product/research/dgic/01_parsed_data.ipynb~~
 
## Output Files
- `cells.csv` contains the cluster assignment for each cell
- `clusters.csv` contains each clusters' mean values for every feature 
(if the max feature value is >1000 then the values will be log transformed for clustering and remain transformed in this output file)
- `scanpy.h5ad` contains the scanpy anndata objects for single-cell analysis

## Paramenter Reference
```
usage: cluster.py [-h] -i INPUT [-o OUTPUT] [-m MARKERS] [-k NEIGHBORS] [-c] [-y CONFIG] [--force-transform] [--no-transform]

Cluster cell types using mcmicro marker expression data.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV of mcmicro marker expression data for cells
  -o OUTPUT, --output OUTPUT
                        The directory to which output files will be saved
  -m MARKERS, --markers MARKERS
                        A text file with a marker on each line to specify which markers to use for clustering
  -k NEIGHBORS, --neighbors NEIGHBORS
                        the number of nearest neighbors to use when clustering. The default is 30.
  -c, --method          Include a column with the method name in the output files.
  -y CONFIG, --config CONFIG
                        A yaml config file that states whether the input data should be log/logicle transformed.
  --force-transform     Log transform the input data. If omitted, and --no-transform is omitted, log transform is only performed if the max value in the input data is >1000.
  --no-transform        Do not perform Log transformation on the input data. If omitted, and --force-transform is omitted, log transform is only performed if the max value in the input data is >1000.
```

## Default input from mcmicro

```
  downstream:
    -
      name: scanpy
      container: labsyspharm/mc-scanpy
      version: 1.0.1
      cmd: python3 /app/cluster.py -c
      input: -i
```

Taken from the `mcmicro/config` repo. output is defaulted to just the current directory of the script execution if
nothing is provided. so far it doesn't seem ilke by default the command provides a lot of alternative output.

To use the customized version of this module, add the following to `params.yml`

```
downstream:
    -
      name: scanpy
      container: mitpenguin/mc-scanpy
      version: 1.2
```

change the `stop-at` value to `downstream` in the yml file as well.

## Build

Build from this repo as you would any other docker git repo

```
docker build -t mitpenguin/mc-scanpy:1.2 .
```
