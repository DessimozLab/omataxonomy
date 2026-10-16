

"""
>>> import itertools

>>> from omataxonomy import Taxonomy
>>> from ete4 import PhyloTree
>>> tax = Taxonomy()

>>> descendants = tax.get_descendant_taxa('f__Nitrosococcaceae', collapse_subspecies=True, return_tree=True)
>>> print(descendants.write())
((-1292655344,-1504426611,-744880077,-1412537045),(-1911723685,-1206911645),(-1293032416));

>>> print(descendants.to_str(props=['sci_name', 'taxid', 'rank'], compact=True))
                                                                            ╭╴s__Nitrosococcus oceani,-1292655344,species
                                        ╭╴g__Nitrosococcus,-348385529,genus╶┼╴s__Nitrosococcus watsonii,-1504426611,species
                                        │                                   ├╴s__Nitrosococcus halophilus,-744880077,species
╴f__Nitrosococcaceae,-1150330694,family╶┤                                   ╰╴s__Nitrosococcus wardiae,-1412537045,species
                                        ├╴g__Nitrosacidococcus,-255272695,genus╶┬╴s__Nitrosacidococcus sp945836005,-1911723685,species
                                        │                                       ╰╴s__Nitrosacidococcus tergens,-1206911645,species
                                        ╰╴g__Nitrosoglobus,-969350119,genus╶╌╴s__Nitrosoglobus terrae,-1293032416,species


>>> taxids = list(itertools.chain.from_iterable(
...     tax.get_name_translator(["p__Hydrothermarchaeota", "c__Methanoliparia", "RS_GCF_000012805.1"])
...    .values()
... ))
>>> tree = tax.get_topology(taxids, intermediate_nodes=True, collapse_subspecies=True, annotate=True)
>>> print(tree.to_str(props=["taxid", "sci_name", "rank"], compact=True))
╴1,root,no rank╶┬╴-2086509484,d__Bacteria,superkingdom╶╌╴-2022353309,p__Pseudomonadota,phylum╶╌╴-1336395690,c__Gammaproteobacteria,class╶╌╴-1490952751,o__Nitrosococcales,order╶╌╴-1150330694,f__Nitrosococcaceae,family╶╌╴-348385529,g__Nitrosococcus,genus╶╌╴-1292655344,s__Nitrosococcus oceani,species
                ╰╴-1206591166,d__Archaea,superkingdom╶┬╴-2104852476,p__Halobacteriota,phylum╶╌╴-1789775371,c__Methanoliparia,class
                                                      ╰╴-995994921,p__Hydrothermarchaeota,phylum


>>> tree = PhyloTree('((c__Thorarchaeia, c__Lokiarchaeia), s__Caballeronia udeis);',
...                    sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
>>> tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
>>> print(tree.to_str(props=["taxid", "name", "sci_name", "rank"], compact=True))
                  ╭╴-2035759834,⊗,p__Asgardarchaeota,phylum╶┬╴-825485653,c__Thorarchaeia,c__Thorarchaeia,class
╴1,⊗,root,no rank╶┤                                         ╰╴-72785574,c__Lokiarchaeia,c__Lokiarchaeia,class
                  ╰╴-1479973048,s__Caballeronia udeis,s__Caballeronia udeis,species

>>> tree = PhyloTree('(RS_GCF_000012805.1, (Homo sapiens, Gallus gallus));',
...                  sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
>>> tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
>>> print(tree.to_str(props=["taxid", "name", "sci_name", "rank"], compact=True))
╴1,⊗,root,no rank╶┬╴-540457981,RS_GCF_000012805.1,RS_GCF_000012805.1,subspecies
                  ╰╴32524,⊗,Amniota,clade╶┬╴9606,Homo sapiens,Homo sapiens,species
                                          ╰╴9031,Gallus gallus,Gallus gallus,species

>>> print(tax.get_name_lineage(['RS_GCF_000012805.1', 'RS_GCF_020885915.1', "Homo sapiens", "f__Leptotrichiaceae", "Gallus"]))
{'f__Leptotrichiaceae': ['root', 'd__Bacteria', 'p__Fusobacteriota', 'c__Fusobacteriia', 'o__Fusobacteriales', 'f__Leptotrichiaceae'], 'Gallus': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Sauropsida', 'Sauria', 'Archelosauria', 'Archosauria', 'Dinosauria', 'Saurischia', 'Theropoda', 'Coelurosauria', 'Aves', 'Neognathae', 'Galloanserae', 'Galliformes', 'Phasianidae', 'Phasianinae', 'Gallus'], 'Homo sapiens': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Homo', 'Homo sapiens'], 'RS_GCF_000012805.1': ['root', 'd__Bacteria', 'p__Pseudomonadota', 'c__Gammaproteobacteria', 'o__Nitrosococcales', 'f__Nitrosococcaceae', 'g__Nitrosococcus', 's__Nitrosococcus oceani', 'RS_GCF_000012805.1'], 'RS_GCF_020885915.1': ['root', 'd__Archaea', 'p__Halobacteriota', 'c__Methanoliparia', 'o__Methanoliparales', 'f__Methanoliparaceae', 'g__Methanoliparum', 's__Methanoliparum thermophilum', 'RS_GCF_020885915.1']}
>>> print(tax.get_mnemonic_names([9606, 43715, 658031, 73382, 9823]))
{9606: 'HUMAN', 9823: 'PIG', 43715: 'ACRAC', 658031: 'AEDFL'}

>>> lin = tax.get_lineage(10090)
>>> t2n = tax.get_taxid_translator(lin)
>>> print([t2n[x] for x in lin])
['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Glires', 'Rodentia', 'Myomorpha', 'Muroidea', 'Muridae', 'Murinae', 'Mus', 'Mus', 'Mus musculus']

>>> print(tax.get_mnemonic_translator(['MOUSE', 'YEAST', 'EREGS', 'CAPSP', ]))
{'CAPSP': [-50085752], 'EREGS': [284811], 'MOUSE': [10090], 'YEAST': [559292]}

"""