# omataxonomy

omataxonomy combines the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
and the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org),
including every GTDB genome, into a single taxonomy that you can query from
Python. GTDB taxa and genomes get numeric taxids, reproducible for a given
GTDB release, so they can be used together with NCBI taxids. It is
developed for the [OMA](https://omabrowser.org) project and can be used
independently of it.

```python
from omataxonomy import Taxonomy

tax = Taxonomy()
tax.get_name_lineage(["Homo sapiens", "f__Leptotrichiaceae", "RS_GCF_000012805.1"])
```

```{toctree}
:maxdepth: 2

background
installation
database
tutorial
api
```
