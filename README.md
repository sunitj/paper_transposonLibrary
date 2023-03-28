# Transposon Library Paper

Authors: Carlos Voogdt, Surya Tripathi, K.C Huang et al.

This repository contains the code and data for the community comparison tree figure in the paper *"[WORKING TITLE] Transposon Library Paper"* by Carlos Voogdt, Surya Tripathi, K.C Huang etal. Here we compile a list of species across communities described in Zimmermann et.al, Maier et.al and Cheng et.al., generate a unique superset (see: [genomes_superset.csv](data/imported/genomes_superset.csv)), select a representative genome for each species, place them on a tree, and then prune the tree to show only the species of interest.

## Environment

1. Create a new virtual environment. Here is an example command using conda:

```bash
conda create -n trees python=3.11
```

2. Activate the environment:

```bash
conda activate trees
```

3. Install the required packages:

```bash
git clone REPO SOURCE
cd TransposonLibrary_Voodgt_Tripathi_etal
pip install -U .
conda install -c conda-forge ncbi-datasets-cli
```

## How to generate the tree?

Once you have a list of species names that you'd like to visualize in a tree (example: [genomes.list](data/imported/genomes.list)). Upload them to the [NCBI Tax Identifier](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) tool to get representative taxon ids. Then, create a file with one taxon id per line (see: [taxids.list](data/imported/taxids.list)). Finally, run the notebooks in the following order:

### Download the genomes

This notebook first identifies a representative genome for each taxon id provided as input, then downloads the genomes from NCBI using the `ncbi-datasets` command line tool.

Notebook: [download_genomes.ipynb](place_genome_on_tree/download_genomes.ipynb)

Data generated: `data/generated/download_genomes/`

```bash
data/generated/download_genomes/
├── genome_accessions.list
├── genomes.zip (not uploaded to github)
└── genomes_metadata.csv
```

### Place genomes on the GTDB tree

We use gtdb-tk to place these genomes on a tree. This notebook first extracts the representative genomes from the zip file generated in the previous step, then processes the genomes using gtdb-tk's `classify_wf`.

Notebook: [process_genomes.ipynb](place_genome_on_tree/process_genomes.ipynb)

Data generated: `data/generated/process_genomes/`

```bash
data/generated/process_genomes
├── gtdb.TransposonLibrary_20210331.bac120.classify.tree
└── gtdb.TransposonLibrary_20210331.bac120.summary.tsv
```

### Prune the tree

This notebook uses the `treeViz.py` python script to prune the tree and generate a figure as a pdf file. It also procduces additional files (`.tree` - newick format tree) that can be used with a tree viewer such as [iTOL](https://itol.embl.de/) to make customizations.

Notebook: [prune_tree.ipynb](place_genome_on_tree/prune_tree.ipynb)
Dependency: [treeViz.py](place_genome_on_tree/prune_tree.ipynb) python script
Data generated: `data/generated/prune_tree/`

```bash
data/generated/prune_tree/
├── TransposonLibrary_20210331.circular_w_bgcolor.pruned.pdf
├── TransposonLibrary_20210331.circular_w_bgcolor.pruned.tree
├── TransposonLibrary_20210331.rect_no_color.pruned.pdf
├── TransposonLibrary_20210331.rect_no_color.pruned.tree
├── TransposonLibrary_20210331.rect_w_bgcolor.pruned.pdf
├── TransposonLibrary_20210331.rect_w_bgcolor.pruned.tree
└── TransposonLibrary_20210331.taxa_levels.csv
```
