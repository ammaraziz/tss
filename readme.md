## Pipeline Overview

`tss` is an Snakemake-based Nextstrain/Augur phylogenetic reconstruction pipeline for seasional Influenza. It is based on the [Nextstrain seasonal-influ build](https://github.com/nextstrain/seasonal-flu). 

Steps:

1. Align to annotated reference (`.gbk`) and trim:
	- H1N1pdm A/California/07/2009
	- H3N2 A/Beijing/32/1992
	- Bvic B/Hong Kong/02/1993
	- Byam B/Singapore/11/1994
2. Construct tree using `IQTree`
3. Refine tree using sequencing data (`.nwk`)
4. Reconstruct ancestral sequences using `TreeTime`
5. Translate sequences/mutations to amino acid (.json)
6. Assign clades based on nt/aa mutations (`.json`)
7. Calculate recency based on date (`.json`)
8. Add mutations (aa, nt), clades and receny data to tree (`.nhx`)
9. Plot tree using `ggtree` (`.pdf`)

## Usage

### Input

1. Multi sequence fasta file
	- Do not include any spaces in the file names(replace with `_`)
	- Replace diacritic letters with their standard equvalent. eg `ô` with `o`
	- Replace or delete apostrophes `'`
	- Most standard unicode characters are accepted (`*?_-`) 
2. Metadata in `.tsv` format
Required fields:
	- strain: virus name, must match exactly the fasta headers
	- date: sample date in YYYY-MM-DD. If unknown leave blank.
	- month_abbr(default): Virus names on tree will be colored according to this column
	- type(default): Tip shape and color will be set according to this column 
Other columns/data can be included. No empty rows are allowed.

See examples of each in `data` folder.

### Running the Pipeline

1. Place your formatted `.fasta` and `.tsv` files into the `data` folder
2. Modify the `snakefile` to point to these files.
3. Select lineage/segment. eg `H3N2` and `ha`

For example:

```
input_fasta = Path("data/YOURDATA.fasta")
input_metadata = Path("data/YOURDATA.tsv")
lineage = 'h3n2'
segment = 'ha'
```

Finally, in the terminal navigate to the `tss` directory. Run the pipeline:

```
snakemake -j 15
```

Where `-j` is the number of cpus to use. 

### Output

Output directory is `tss/results`. In this folder you will find the below.

Of interest:
- YOURDATA.fasta - aligned, trimmed fasta file
- YOURDATA.fasta.insertions.csv - insertions that were trimmed from the alignment
- YOURDATA.pdf - plot of the tree
- YOURDATA.tree - `.nhx` tree file

Augur specific files:
- aa_muts.json 
- branch_lengths.json
- nt_muts.json
- recency.json
- tree_raw.nwk
- tree_refined.nwk - refined tree
- data/ha_clades.json - clade data, used for assigning clades to NA.

Log files. Useful for debugging errors otherwise ignore
- YOURDATA-delim.fasta.log
- YOURDATA-delim.iqtree.log
- YOURDATA.fasta.log
- refine.log
- tree.log
- align.log

### Common Issues and Questions

1. Hitting the recussion limit

This happens with very large trees.

If you see an error about 'recursion limit' from augur, then increase the limit. In bash run:

```
export AUGUR_RECURSION_LIMIT=10000
```

For increasing python recurrsion limit, add this to the top of the `snakefile`

```
import sys

sys.setrecursionlimit(5000)
```

2. Colors do not match/taxa labels not colored

When running the pipeline, a list of missing taxa names will be printed if the `strain` in meta data does not match the fasta headers. The `strain` column needs to match *exactly* the fasta headers.

3. How are the colors/shapes/categories set?

There are two `.tsv` files which control the colors and shapes used in the pdf plot. They can be found in `scripts` folder. The numbers for 'shapes' are ggplot2/r specific. [See this link for more information](http://www.sthda.com/english/wiki/ggplot2-point-shapes)

4. I get a weird error, a strange tree or random data.

Delete the output from `results`, rerun the pipeline. 

Look at the `.log` files, it will indicate what/where the pipeline failed.

5. Can you help with X?

Create an issue but please be nice.

## Install

### Install conda

https://docs.conda.io/en/main/miniconda.html

Python version > 3.8 is required

### Create a conda environemnt and install packages:

```
conda env create -n tss --file tss.yml
```

### Install R packages

```
RPACKAGES=$(Rscript <(echo 'chooseCRANmirror(ind=2)
install.packages("BiocManager")
install.packages("optparse")
install.packages("ggplot2")
BiocManager::install("ggtree")'))
```

### Download the pipeline:

Either download the pipeline from github or run:

```
git clone github.com/ammaraziz/tss
```

Place this somewhere smart like your `bin` folder in your home directory.

### Dry run

In the `tss` directory perform a dry run:

```
snakemake -nq
```

You should see:

```
Building DAG of jobs...
Job stats:
job            count    min threads    max threads
-----------  -------  -------------  -------------
align              1              1              1
all                1              1              1
ancestral          1              1              1
clades             1              1              1
export             1              1              1
plot_ggtree        1              1              1
recency            1              1              1
refine             1              1              1
translate          1              1              1
tree               1              1              1
total             10              1              1
```

Or similar.

