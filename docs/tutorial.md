# Tutorial

This tutorial walks through the most common queries. All outputs were
produced with a database built from GTDB release 232 (see
{doc}`database`). With another release, the GTDB parts of the outputs can
differ.

## Getting started

Everything starts from a {class}`~omataxonomy.Taxonomy` object. It opens the
database, or builds it the first time. We also import
{func}`pprint.pprint` to print the larger results in a readable way.

```pycon
>>> from pprint import pprint
>>> from omataxonomy import Taxonomy
>>> tax = Taxonomy()

```

## Names and taxids

{meth}`~omataxonomy.Taxonomy.get_name_translator` maps names to taxids. The
values are lists, because a name can be ambiguous. NCBI taxa have positive
taxids and GTDB taxa and genomes have negative ones. The lookup ignores
case. If a name is not a scientific name, its synonyms are searched, so a
GTDB genome is also found by its GenBank accession:

```pycon
>>> pprint(tax.get_name_translator(["Homo sapiens", "f__Leptotrichiaceae", "GB_GCA_001515945.1", "GCA_000012805.1"]))
{'GB_GCA_001515945.1': [-1326969492],
 'GCA_000012805.1': [-540457981],
 'Homo sapiens': [9606],
 'f__Leptotrichiaceae': [-769982259]}

```

Names that cannot be found are left out of the result instead of raising an
error. Check for them if your input can contain typos.

The other direction is {meth}`~omataxonomy.Taxonomy.get_taxid_translator`.
{meth}`~omataxonomy.Taxonomy.translate_to_names` does the same for a list
and keeps its order:

```pycon
>>> pprint(tax.get_taxid_translator([9606, -540457981]))
{-540457981: 'RS_GCF_000012805.1', 9606: 'Homo sapiens'}
>>> tax.translate_to_names([9606, -540457981])
['Homo sapiens', 'RS_GCF_000012805.1']

```

Note that a genome is named by its GTDB accession, which includes the `RS_`
(RefSeq) or `GB_` (GenBank) prefix.

## Lineages

{meth}`~omataxonomy.Taxonomy.get_lineage` returns the taxids from the root
down to the given taxon. Here is the house mouse:

```pycon
>>> lineage = tax.get_lineage(10090)
>>> lineage[:5]
[1, 131567, 2759, 33154, 33208]
>>> tax.translate_to_names(lineage)[-4:]
['Murinae', 'Mus', 'Mus', 'Mus musculus']

```

(`Mus` appears twice because NCBI has both a genus and a subgenus *Mus*.)
{meth}`~omataxonomy.Taxonomy.get_rank` returns the rank of each taxid:

```pycon
>>> ranks = tax.get_rank(lineage)
>>> [ranks[t] for t in lineage[-4:]]
['subfamily', 'genus', 'subgenus', 'species']

```

If you start from names, {meth}`~omataxonomy.Taxonomy.get_name_lineage` does
the translation for you. It returns the lineage of each name as a list of
names. NCBI and GTDB names can be mixed:

```pycon
>>> lineages = tax.get_name_lineage(["RS_GCF_000012805.1", "RS_GCF_020885915.1", "Gallus"])
>>> pprint(lineages, width=100)
{'Gallus': ['root',
            'cellular organisms',
            'Eukaryota',
            'Opisthokonta',
            'Metazoa',
            'Eumetazoa',
            'Bilateria',
            'Deuterostomia',
            'Chordata',
            'Craniata',
            'Vertebrata',
            'Gnathostomata',
            'Teleostomi',
            'Euteleostomi',
            'Sarcopterygii',
            'Dipnotetrapodomorpha',
            'Tetrapoda',
            'Amniota',
            'Sauropsida',
            'Sauria',
            'Archelosauria',
            'Archosauria',
            'Dinosauria',
            'Saurischia',
            'Theropoda',
            'Coelurosauria',
            'Aves',
            'Neognathae',
            'Galloanserae',
            'Galliformes',
            'Phasianidae',
            'Phasianinae',
            'Gallus'],
 'RS_GCF_000012805.1': ['root',
                        'd__Bacteria',
                        'p__Pseudomonadota',
                        'c__Gammaproteobacteria',
                        'o__Nitrosococcales',
                        'f__Nitrosococcaceae',
                        'g__Nitrosococcus',
                        's__Nitrosococcus oceani',
                        'RS_GCF_000012805.1'],
 'RS_GCF_020885915.1': ['root',
                        'd__Archaea',
                        'p__Halobacteriota',
                        'c__Methanoliparia',
                        'o__Methanoliparales',
                        'f__Methanoliparaceae',
                        'g__Methanoliparum',
                        's__Methanoliparum thermophilum',
                        'RS_GCF_020885915.1']}

```

The GTDB lineages go directly from `root` to `d__Bacteria` or `d__Archaea`.
They never pass through NCBI's `cellular organisms` (see
{doc}`background`).

## UniProt species mnemonics

UniProt identifies species by short mnemonic codes, such as `HUMAN` in the
protein entry name `INS_HUMAN`. omataxonomy translates between taxids and
mnemonics in both directions. Taxa that have no mnemonic are left out of
the result:

```pycon
>>> pprint(tax.get_mnemonic_names([9606, 9823, 43715, 658031, 73382]))
{9606: 'HUMAN', 9823: 'PIG', 43715: 'ACRAC', 658031: 'AEDFL'}
>>> pprint(tax.get_mnemonic_translator(["MOUSE", "YEAST", "EREGS", "CAPSP"]))
{'CAPSP': [-50085752], 'EREGS': [284811], 'MOUSE': [10090], 'YEAST': [559292]}

```

Mnemonics that UniProt gives to an NCBI taxid are moved to the matching
GTDB genome when there is one. Above, `CAPSP` resolves to a negative taxid,
i.e. to a GTDB genome. If several GTDB genomes share that NCBI taxid, the
GTDB species representative is chosen.

## Descendants of a taxon

{meth}`~omataxonomy.Taxonomy.get_descendant_taxa` accepts a taxid or a name.
By default it returns the taxids of all leaves below that taxon. For a GTDB
taxon, those leaves are the individual genomes, which can be a lot. With
`collapse_subspecies=True`, the genomes are collapsed into their species:

```pycon
>>> species = list(tax.get_descendant_taxa("f__Nitrosococcaceae", collapse_subspecies=True))
>>> tax.translate_to_names(species)
['s__Nitrosococcus oceani', 's__Nitrosococcus watsonii', 's__Nitrosococcus halophilus', 's__Nitrosococcus wardiae', 's__Nitrosacidococcus sp945836005', 's__Nitrosacidococcus tergens', 's__Nitrosoglobus terrae']

```

With `return_tree=True`, you get an ete
[`PhyloTree`](https://etetoolkit.github.io/ete/reference/reference_phylo.html)
instead. Its nodes are named by taxid and carry the properties `sci_name`,
`taxid`, `rank` and `lineage`:

```pycon
>>> subtree = tax.get_descendant_taxa("f__Nitrosococcaceae", collapse_subspecies=True, return_tree=True)
>>> print(subtree.write())
((-1292655344,-1504426611,-744880077,-1412537045),(-1911723685,-1206911645),(-1293032416));
>>> print(subtree.to_str(props=["sci_name", "rank"], compact=True))
                                                     ╭╴s__Nitrosococcus oceani,species
                            ╭╴g__Nitrosococcus,genus╶┼╴s__Nitrosococcus watsonii,species
                            │                        ├╴s__Nitrosococcus halophilus,species
╴f__Nitrosococcaceae,family╶┤                        ╰╴s__Nitrosococcus wardiae,species
                            ├╴g__Nitrosacidococcus,genus╶┬╴s__Nitrosacidococcus sp945836005,species
                            │                            ╰╴s__Nitrosacidococcus tergens,species
                            ╰╴g__Nitrosoglobus,genus╶╌╴s__Nitrosoglobus terrae,species

```

## The tree connecting a set of taxa

{meth}`~omataxonomy.Taxonomy.get_topology` returns the smallest subtree of
the taxonomy that connects the given taxids. Here we combine two archaeal
groups with the Nitrosococcus genome from above. The names are translated
to taxids first:

```pycon
>>> name2taxids = tax.get_name_translator(["p__Hydrothermarchaeota", "c__Methanoliparia", "RS_GCF_000012805.1"])
>>> taxids = [taxid for ids in name2taxids.values() for taxid in ids]

```

`intermediate_nodes=True` keeps the full lineage of each taxon, and
`collapse_subspecies=True` ends it at the species instead of the genome:

```pycon
>>> tree = tax.get_topology(taxids, intermediate_nodes=True, collapse_subspecies=True)
>>> print(tree.to_str(props=["sci_name", "rank"], compact=True))
╴root,no rank╶┬╴d__Bacteria,superkingdom╶╌╴p__Pseudomonadota,phylum╶╌╴c__Gammaproteobacteria,class╶╌╴o__Nitrosococcales,order╶╌╴f__Nitrosococcaceae,family╶╌╴g__Nitrosococcus,genus╶╌╴s__Nitrosococcus oceani,species
              ╰╴d__Archaea,superkingdom╶┬╴p__Halobacteriota,phylum╶╌╴c__Methanoliparia,class
                                        ╰╴p__Hydrothermarchaeota,phylum

```

Without `intermediate_nodes`, nodes with a single child are removed. Only
the requested taxa and the points where lineages split remain:

```pycon
>>> tree = tax.get_topology(taxids, collapse_subspecies=True)
>>> print(tree.to_str(props=["sci_name", "rank"], compact=True))
              ╭╴d__Archaea,superkingdom╶┬╴p__Hydrothermarchaeota,phylum
╴root,no rank╶┤                         ╰╴c__Methanoliparia,class
              ╰╴RS_GCF_000012805.1,subspecies

```

## Annotating your own tree

If you have a tree whose leaves are taxa, e.g. a species tree,
{meth}`~omataxonomy.Taxonomy.annotate_tree` adds the taxonomic information
to it. The method reads the taxid of each node from the node property named
by `taxid_attr`.

In this example, the leaves are named by GTDB taxa. We tell
{class}`ete4.PhyloTree` to set the `species` property of each leaf to the
taxid of its name, and then pass `taxid_attr="species"`:

```pycon
>>> from ete4 import PhyloTree
>>> def name_to_taxid(name):
...     return tax.get_name_translator([name])[name][0]
>>> tree = PhyloTree("((c__Thorarchaeia, c__Lokiarchaeia), s__Caballeronia udeis);",
...                  sp_naming_function=name_to_taxid)
>>> tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")

```

Every node, including the internal ones, now has `taxid`, `sci_name`,
`rank` and `lineage` properties. Internal nodes are assigned the last common
ancestor of their leaves. Here, Thorarchaeia and Lokiarchaeia are joined in
the phylum Asgardarchaeota. The `name` column shows the original node names
(`⊗` means the node has no name):

```pycon
>>> print(tree.to_str(props=["name", "sci_name", "rank"], compact=True))
                ╭╴⊗,p__Asgardarchaeota,phylum╶┬╴c__Thorarchaeia,c__Thorarchaeia,class
╴⊗,root,no rank╶┤                             ╰╴c__Lokiarchaeia,c__Lokiarchaeia,class
                ╰╴s__Caballeronia udeis,s__Caballeronia udeis,species

```

Genomes and NCBI taxa can be mixed in one tree. As explained in
{doc}`background`, a GTDB genome and an NCBI taxon only meet at the root:

```pycon
>>> tree = PhyloTree("(RS_GCF_000012805.1, (Homo sapiens, Gallus gallus));",
...                  sp_naming_function=name_to_taxid)
>>> _ = tax.annotate_tree(tree, taxid_attr="species")
>>> print(tree.to_str(props=["sci_name", "rank"], compact=True))
╴root,no rank╶┬╴RS_GCF_000012805.1,subspecies
              ╰╴Amniota,clade╶┬╴Homo sapiens,species
                              ╰╴Gallus gallus,species

```

`annotate_tree` also returns three dictionaries that map each taxid in the
tree to its name, its lineage and its rank. They are handy for further
processing:

```pycon
>>> tax2rank[-72785574], tax2name[-72785574]
('class', 'c__Lokiarchaeia')

```
