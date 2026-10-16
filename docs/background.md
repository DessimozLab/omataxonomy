# Why combine NCBI and GTDB?

There are two taxonomies that matter in prokaryotic genomics, and neither one
covers everything.

## The NCBI Taxonomy

The [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) is the reference
taxonomy of the INSDC sequence databases (GenBank, ENA, DDBJ). Every organism
with a deposited sequence has an NCBI taxid, so it covers all of life,
including eukaryotes and viruses. It is assembled by curators from the
published literature. For Bacteria and Archaea this has well-known drawbacks:

- Many groups are not monophyletic. Historical names such as *Clostridium*
  have been used for organisms that, according to genome phylogenies, belong
  to many unrelated lineages.
- Ranks are applied unevenly. A family in one part of the tree can be much
  older or much younger than a family elsewhere.
- Uncultured organisms and metagenome-assembled genomes are often only
  placed as "unclassified" or "environmental samples".

## The Genome Taxonomy Database (GTDB)

[GTDB](https://gtdb.ecogenomic.org) is built from a phylogeny of conserved
marker genes (`bac120` for Bacteria and `ar53` for Archaea). Its taxa are
monophyletic, and their ranks are normalised by relative evolutionary
divergence, so a GTDB family means roughly the same thing everywhere in the
tree. Each species is defined around a representative genome. GTDB has its
own conventions:

- Every name carries a rank prefix: `d__` (domain), `p__`, `c__`, `o__`,
  `f__`, `g__` and `s__`.
- When a traditional group turns out to be polyphyletic, it is split and the
  parts get alphabetical suffixes. In GTDB release 232 there are dozens of
  `g__Clostridium_*` genera (`g__Clostridium_A`, `g__Clostridium_B`, ...).
- Species without a valid name get placeholders derived from an accession,
  e.g. `s__Desulfitibacter sp001515945`.

However, GTDB only covers Bacteria and Archaea, and it has **no numeric
taxon identifiers**. Taxa are identified only by their names, and genomes by
their assembly accessions.

## The problem

Many applications need both taxonomies at once: eukaryotes and viruses from
NCBI, and prokaryotes classified with GTDB. With two separate taxonomies,
they run into the following problems:

- The two classifications of the same organism disagree. Here is the
  ammonia-oxidising bacterium *Nitrosococcus oceani* (genome
  `GCA_000012805.1`) in both of them. NCBI places it in the order
  Chromatiales; GTDB places it in its own order, Nitrosococcales:

  ```pycon
  >>> from omataxonomy import Taxonomy
  >>> tax = Taxonomy()
  >>> def named_lineage(name):
  ...     taxid = tax.get_name_translator([name])[name][0]
  ...     lineage = tax.get_lineage(taxid)
  ...     return tax.translate_to_names(lineage)
  >>> named_lineage("Nitrosococcus oceani")
  ['root', 'cellular organisms', 'Bacteria', 'Pseudomonadati', 'Pseudomonadota', 'Gammaproteobacteria', 'Chromatiales', 'Chromatiaceae', 'Nitrosococcus', 'Nitrosococcus oceani']
  >>> named_lineage("GCA_000012805.1")
  ['root', 'd__Bacteria', 'p__Pseudomonadota', 'c__Gammaproteobacteria', 'o__Nitrosococcales', 'f__Nitrosococcaceae', 'g__Nitrosococcus', 's__Nitrosococcus oceani', 'RS_GCF_000012805.1']

  ```

- Without numeric identifiers, GTDB taxa cannot be stored where tools expect
  an integer taxid, and they cannot be mixed with NCBI taxids.
- Many names exist in both taxonomies (`Gammaproteobacteria`,
  `Nitrosococcus`, ...) but stand for different sets of organisms.

## What omataxonomy does

omataxonomy merges both taxonomies into a single tree, stored in one SQLite
database:

- The full NCBI Taxonomy is kept unchanged, with its positive taxids.
- The GTDB tree (Bacteria and Archaea) is added as a second subtree directly
  below the root. The prefixed GTDB names keep the two apart: `Bacteria` is
  NCBI's domain and `d__Bacteria` is GTDB's.
- **Every GTDB genome becomes a leaf** of rank `subspecies` below its GTDB
  species. Its name is the GTDB accession (e.g. `RS_GCF_000012805.1`).
- GTDB taxa and genomes get **stable negative taxids**. These are computed by
  hashing the name (for genomes, the accession number), so a taxon keeps its
  taxid in every release in which its name stays the same. Hash collisions
  are resolved deterministically.
- Genomes can also be found by their GenBank accession, by their NCBI
  organism name, and by the NCBI taxid recorded for them in the GTDB
  metadata (as the synonym `ncbi_taxid:<taxid>`):

  ```pycon
  >>> tax.get_name_translator(["d__Bacteria", "Bacteria"])
  {'Bacteria': [2], 'd__Bacteria': [-2086509484]}
  >>> tax.get_name_translator(["RS_GCF_000012805.1", "GCA_000012805.1", "ncbi_taxid:323261"])
  {'RS_GCF_000012805.1': [-540457981], 'GCA_000012805.1': [-540457981], 'ncbi_taxid:323261': [-540457981]}

  ```

- [UniProt species mnemonics](https://www.uniprot.org/docs/speclist) such as
  `HUMAN` or `ECOLI` are attached to the taxa, and you can query by them.

:::{note}
The NCBI and GTDB subtrees are not linked to each other. The last common
ancestor of a GTDB genome and an NCBI taxon is therefore always `root`, even
if both are bacteria. Use GTDB names for prokaryotes and NCBI names for
everything else.
:::

The query interface follows ete's
[`NCBITaxa`](https://etetoolkit.github.io/ete/tutorial/tutorial_taxonomy.html)
module. omataxonomy reads the data from the database as it needs it instead
of loading it into memory, so it uses little memory even though the
combined taxonomy has more than four million nodes.
