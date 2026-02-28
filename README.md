# MeganServer Access

A simple Python package for downloading and analysing datasets from a MeganServer. It provides:

- command-line utilities for listing and fetching classification blocks from a megan-server
- functions to merge, clean and convert datasets into OTU/taxonomy tables
- a `MicrobiomeDataAnalyzer` class with statistics and plotting helpers for microbiome studies

This repository contains helper modules for working with megan-server outputs and performing downstream reporting and visualization.

For more information about the megan-server refer to the following repo: https://github.com/husonlab/megan-ce/tree/master/src/megan/ms
---

## Features

- Authenticate and list available datasets from the megan-server
- Select multiple datasets interactively and merge them into a single comparison table
- Convert raw MEGAN classification blocks into OTU and taxonomy matrices
- Request taxonomic lineages from NCBI using `ete3`
- Perform common microbiome analyses: PCoA, ANOVA, t‑tests, Wilcoxon, PERMANOVA, beta diversity
- Plot rank‑level abundance and top taxa charts

## Requirements

- Python 3.8+
- [`pandas`](https://pandas.pydata.org/)
- [`requests`](https://pypi.org/project/requests/)
- [`ete3`](http://etetoolkit.org/) (requires `ncbi_taxonomy` database)
- [`numpy`, `scipy`, `matplotlib`, `scikit-bio`]

You can install the dependencies with pip:

```bash
pip install pandas requests ete3 numpy scipy matplotlib scikit-bio
```

> **Note:** `ete3` may require additional system dependencies (e.g., `==libsqlite3-dev` on Linux) and a one‑time database update via `ncbi.update_taxonomy_database()`.

## Installation

There is no published package yet, so clone the repository and use the modules directly:

```bash
git clone https://github.com/negorov15/megan_access.git
cd megan_access
# optionally create a virtual environment
python3 -m venv venv && source venv/bin/activate
pip install pandas requests ete3 numpy scipy matplotlib scikit-bio
```

## Usage

### Command‑line access

Run `megan_accessor.py` to interactively query the megan-server:

```bash
python3 megan_accessor.py
```

You will be prompted for your `login` and `password` (both can be set to 'guest'), shown a list of available datasets and asked to select one or more by name. After choosing the classification (e.g. `taxonomy`, `SEED`, `GTDB`), the script downloads each file, cleans it and merges the results into one csv file.

### Examples

Read merged comparison files and derive OTU and taxonomy tables:

```python
import pandas as pd
from modificator import otu_table, tax_table
from microbiome_class import MicrobiomeDataAnalyzer

# Suppose we have two different merged datasets: Alice and Bob
otumat_alice = otu_table('alice.csv')
taxmat_alice = tax_table('alice.csv')

otumat_bob = otu_table('bob.csv')
taxmat_bob = tax_table('bob.csv')
```

Define sample metadata:

```python
metadata_dict_alice = {
    'SampleID': ['Alice00-1mio.daa', 'Alice01-1mio.daa',
                 'Alice03-1mio.daa', 'Alice06-1mio.daa',
                 'Alice08-1mio.daa', 'Alice34-1mio.daa'],
    'Group': ['Without treatment', 'With treatment', 'With treatment',
              'With treatment', 'Without treatment', 'Without treatment'],
    'Property': ['0-', '1+', '3+', '6+', '8-', '34-']
}
metadata_alice = pd.DataFrame(metadata_dict_alice)

metadata_dict_bob = {
    'SampleID': ['Bob00-1mio.daa', 'Bob01-1mio.daa',
                 'Bob03-1mio.daa', 'Bob06-1mio.daa',
                 'Bob08-1mio.daa', 'Bob34-1mio.daa'],
    'Group': ['Without treatment', 'With treatment', 'With treatment',
              'With treatment', 'Without treatment', 'Without treatment'],
    'Property': ['0-', '1+', '3+', '6+', '8-', '34-']
}
metadata_bob = pd.DataFrame(metadata_dict_bob)

sample_alice = MicrobiomeDataAnalyzer(otumat_alice, taxmat_alice, metadata_alice)
sample_bob = MicrobiomeDataAnalyzer(otumat_bob, taxmat_bob, metadata_bob)
```

Perform statistical analysis and visualization:

```python
print(sample_alice.permanova())
print(sample_bob.permanova())
print(sample_alice.t_test())
print(sample_bob.t_test())
sample_bob.plot_rank('Phylum')
sample_bob.plot_top(5)
```

## Data formats

### MEGAN Classification Blocks
Raw output from the server. Tab-delimited text with multiple header lines:
- **Lines 1–2**: metadata/comments (skipped by `dataset_modifier()`)
- **Line 3+**: data rows with columns `readID`, `taxid`, `weight`, etc.

Example:
```
#Files	Datasets
#CreationDate	2025-01-15
readID	taxid	weight
read_001	562	100.5
read_002	1280	89.3
```

### Merged Comparison Table
Produced by `merge_data()`. Tab-delimited CSV with a `Taxa` column (NCBI taxon IDs) and one column per sample (read counts):

Example (2 samples merged):
```
Taxa	alice.csv	bob.csv
562	150	200
1280	89	120
2157	45	78
```

### OTU Table
Output of `otu_table()`. Rows = OTUs (named `OTU1`, `OTU2`, ...), columns = samples with abundance counts. Used directly by `MicrobiomeDataAnalyzer`:

|        | alice.csv | bob.csv |
|--------|-----------|---------|
| OTU1   | 150       | 200     |
| OTU2   | 89        | 120     |
| OTU3   | 45        | 78      |

### Taxonomy Table
Output of `tax_table()`. Rows = OTUs, columns = taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species). Missing ranks are `NaN`:

|      | Kingdom   | Phylum        | Class            | Order       | Family      | Genus     | Species           |
|------|-----------|---------------|------------------|-------------|-------------|-----------|-------------------|
| OTU1 | Bacteria  | Proteobacteria| Gammaproteobacteria | Enterobacterales | Enterobacteriaceae | Escherichia | Escherichia coli  |
| OTU2 | Bacteria  | Firmicutes    | Bacilli          | Bacillales  | Bacillaceae | Bacillus  | NaN               |

### Metadata Table
User-defined pandas DataFrame describing samples. Must contain at least:
- `SampleID`: identifier matching column names in OTU table
- `Group`: experimental grouping (e.g., "Treatment", "Control")
- Optional: additional columns like `Property` for time points or other covariates

Example:
```python
metadata = pd.DataFrame({
    'SampleID': ['alice.csv', 'bob.csv'],
    'Group': ['Control', 'Treatment'],
    'Property': ['Day0', 'Day7']
})
```

---