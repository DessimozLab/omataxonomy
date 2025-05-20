

"""
>>> import itertools

>>> from omataxonomy import Taxonomy
>>> from ete4 import PhyloTree
>>> tax = Taxonomy()

>>> descendants = tax.get_descendant_taxa('c__Thorarchaeia', collapse_subspecies=True, return_tree=True)
>>> print(descendants.write())
(((-2036468532,-2077252572,-34509770,-159862293,-1725536879,-1970875077,-1152849540,-1977275412),(-134830589,-2068828550),(-1110741460,-611357532,-254183316,-1296315779,-598729758,-1402255420,-2078812169),(-1837152244,-1829728714,-1098403373),(-1578072040),(-1232904002,-2112901554),(-1726229641),(-1795683282,-1192980426,-180605867,-1703035231,-1394921732),(-49029718,-1705185043,-325256247),(-2142362126,-1638457287,-1884908377,-431893286,-935030981,-2056141946,-1469378617),(-929230146,-801809351,-1293299626,-1418265015),(-1712588471,-1117767839),(-1269442283),(-1070083721),(-750786656),(-2131970766)));

>>> print(descendants.to_str(props=['sci_name', 'taxid', 'rank'], compact=True))
                                                                                                         ╭╴s__MP8T-1 sp029858305,-2036468532,species
                                                                                                         ├╴s__MP8T-1 sp016839275,-2077252572,species
                                                                                                         ├╴s__MP8T-1 sp030612005,-34509770,species
                                                                           ╭╴g__MP8T-1,-1740241678,genus╶┼╴s__MP8T-1 sp018238205,-159862293,species
                                                                           │                             ├╴s__MP8T-1 sp004524565,-1725536879,species
                                                                           │                             ├╴s__MP8T-1 sp003345545,-1970875077,species
                                                                           │                             ├╴s__MP8T-1 sp018238825,-1152849540,species
                                                                           │                             ╰╴s__MP8T-1 sp016933575,-1977275412,species
                                                                           ├╴g__WJIL01,-1033566068,genus╶┬╴s__WJIL01 sp018335335,-134830589,species
                                                                           │                             ╰╴s__WJIL01 sp014730275,-2068828550,species
                                                                           │                            ╭╴s__WTCK01 sp018238375,-1110741460,species
                                                                           │                            ├╴s__WTCK01 sp023134215,-611357532,species
                                                                           │                            ├╴s__WTCK01 sp019057475,-254183316,species
                                                                           ├╴g__WTCK01,-440441776,genus╶┼╴s__WTCK01 sp030612365,-1296315779,species
                                                                           │                            ├╴s__WTCK01 sp013138615,-598729758,species
                                                                           │                            ├╴s__WTCK01 sp030669595,-1402255420,species
                                                                           │                            ╰╴s__WTCK01 sp030642465,-2078812169,species
                                                                           │                              ╭╴s__MP5-2226 sp024275415,-1837152244,species
                                                                           ├╴g__MP5-2226,-887393441,genus╶┼╴s__MP5-2226 sp029210805,-1829728714,species
                                                                           │                              ╰╴s__MP5-2226 sp030638785,-1098403373,species
                                                                           ├╴g__JAWAQK01,-1920664692,genus╶╌╴s__JAWAQK01 sp038893995,-1578072040,species
                                                                           ├╴g__OWC5,-1946171831,genus╶┬╴s__OWC5 sp003345555,-1232904002,species
                                                                           │                           ╰╴s__OWC5 sp003345595,-2112901554,species
                                                                           ├╴g__JALVPN01,-1556843369,genus╶╌╴s__JALVPN01 sp027031765,-1726229641,species
                                                                           │                               ╭╴s__SMTZ1-83 sp016839605,-1795683282,species
╴o__Thorarchaeales,-634304458,order╶╌╴f__Thorarchaeaceae,-362013849,family╶┤                               ├╴s__SMTZ1-83 sp020355105,-1192980426,species
                                                                           ├╴g__SMTZ1-83,-1143256210,genus╶┼╴s__SMTZ1-83 sp016840795,-180605867,species
                                                                           │                               ├╴s__SMTZ1-83 sp011364985,-1703035231,species
                                                                           │                               ╰╴s__SMTZ1-83 sp029856425,-1394921732,species
                                                                           │                               ╭╴s__JAJRWI01 sp030593165,-49029718,species
                                                                           ├╴g__JAJRWI01,-1773114923,genus╶┼╴s__JAJRWI01 sp024276825,-1705185043,species
                                                                           │                               ╰╴s__JAJRWI01 sp020348985,-325256247,species
                                                                           │                              ╭╴s__SMTZ1-45 sp004376265,-2142362126,species
                                                                           │                              ├╴s__SMTZ1-45 sp001940705,-1638457287,species
                                                                           │                              ├╴s__SMTZ1-45 sp019894665,-1884908377,species
                                                                           ├╴g__SMTZ1-45,-708346970,genus╶┼╴s__SMTZ1-45 sp016840685,-431893286,species
                                                                           │                              ├╴s__SMTZ1-45 sp016840855,-935030981,species
                                                                           │                              ├╴s__SMTZ1-45 sp001563335,-2056141946,species
                                                                           │                              ╰╴s__SMTZ1-45 sp011364905,-1469378617,species
                                                                           │                             ╭╴s__B65-G9 sp029855935,-929230146,species
                                                                           ├╴g__B65-G9,-1352786013,genus╶┼╴s__B65-G9 sp029855925,-801809351,species
                                                                           │                             ├╴s__B65-G9 sp003662765,-1293299626,species
                                                                           │                             ╰╴s__B65-G9 sp021498125,-1418265015,species
                                                                           ├╴g__SHMX01,-919784897,genus╶┬╴s__SHMX01 sp014728035,-1712588471,species
                                                                           │                            ╰╴s__SHMX01 sp008080745,-1117767839,species
                                                                           ├╴g__FT1-004,-1823890864,genus╶╌╴s__FT1-004 sp016839445,-1269442283,species
                                                                           ├╴g__JACAEL01,-1993894532,genus╶╌╴s__JACAEL01 sp013388835,-1070083721,species
                                                                           ├╴g__TEKIR-14,-842377461,genus╶╌╴s__TEKIR-14 sp004524445,-750786656,species
                                                                           ╰╴g__TEKIR-12S,-1684714758,genus╶╌╴s__TEKIR-12S sp004524435,-2131970766,species


>>> taxids = list(itertools.chain.from_iterable(
...     tax.get_name_translator(["p__Hydrothermarchaeota", "o__Peptococcales", "f__Korarchaeaceae", "GB_GCA_011358815.1"])
...    .values()
... ))
>>> tree = tax.get_topology(taxids, intermediate_nodes=True, collapse_subspecies=True, annotate=True)
>>> print(tree.to_str(props=["taxid", "sci_name", "rank"], compact=True))
╴1,root,no rank╶┬╴-2086509484,d__Bacteria,superkingdom╶╌╴-1778889416,p__Bacillota,phylum╶╌╴-1552719505,c__Peptococcia,class╶╌╴-854041149,o__Peptococcales,order
                ╰╴-1206591166,d__Archaea,superkingdom╶┬╴-1857155222,p__Korarchaeota,phylum╶╌╴-1303715095,c__Korarchaeia,class╶╌╴-759496122,o__Korarchaeales,order╶╌╴-1431518377,f__Korarchaeaceae,family
                                                      ╰╴-995994921,p__Hydrothermarchaeota,phylum


>>> tree = PhyloTree('((c__Thorarchaeia, c__Lokiarchaeia), s__Caballeronia udeis);',
...                    sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
>>> tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
>>> print(tree.to_str(props=["taxid", "name", "sci_name", "rank"], compact=True))
                  ╭╴-2035759834,⊗,p__Asgardarchaeota,phylum╶┬╴-825485653,c__Thorarchaeia,c__Thorarchaeia,class
╴1,⊗,root,no rank╶┤                                         ╰╴-72785574,c__Lokiarchaeia,c__Lokiarchaeia,class
                  ╰╴-1479973048,s__Caballeronia udeis,s__Caballeronia udeis,species

>>> tree = PhyloTree('(RS_GCF_006228565.1, (Homo sapiens, Gallus gallus));',
...                  sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
>>> tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
>>> print(tree.to_str(props=["taxid", "name", "sci_name", "rank"], compact=True))
╴1,⊗,root,no rank╶┬╴-1948361583,RS_GCF_006228565.1,RS_GCF_006228565.1,subspecies
                  ╰╴32524,⊗,Amniota,clade╶┬╴9606,Homo sapiens,Homo sapiens,species
                                          ╰╴9031,Gallus gallus,Gallus gallus,species

>>> print(tax.get_name_lineage(['RS_GCF_006228565.1', 'GB_GCA_001515945.1', "Homo sapiens", "f__Leptotrichiaceae", "Gallus"]))
{'f__Leptotrichiaceae': ['root', 'd__Bacteria', 'p__Fusobacteriota', 'c__Fusobacteriia', 'o__Fusobacteriales', 'f__Leptotrichiaceae'], 'Gallus': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Sauropsida', 'Sauria', 'Archelosauria', 'Archosauria', 'Dinosauria', 'Saurischia', 'Theropoda', 'Coelurosauria', 'Aves', 'Neognathae', 'Galloanserae', 'Galliformes', 'Phasianidae', 'Phasianinae', 'Gallus'], 'GB_GCA_001515945.1': ['root', 'd__Bacteria', 'p__Bacillota', 'c__DSM-16504', 'o__Desulfitibacterales', 'f__Desulfitibacteraceae', 'g__Desulfitibacter', 's__Desulfitibacter sp001515945', 'GB_GCA_001515945.1'], 'Homo sapiens': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Homo', 'Homo sapiens'], 'RS_GCF_006228565.1': ['root', 'd__Bacteria', 'p__Bacillota', 'c__Moorellia', 'o__Moorellales', 'f__Moorellaceae', 'g__Moorella', 's__Moorella thermoacetica', 'RS_GCF_006228565.1']}
>>> print(tax.get_mnemonic_names([9606, 43715, 658031, 73382, 9823]))
{9606: 'HUMAN', 9823: 'PIG', 43715: 'ACRAC', 658031: 'AEDFL'}

>>> lin = tax.get_lineage(10090)
>>> t2n = tax.get_taxid_translator(lin)
>>> print([t2n[x] for x in lin])
['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Glires', 'Rodentia', 'Myomorpha', 'Muroidea', 'Muridae', 'Murinae', 'Mus', 'Mus', 'Mus musculus']

>>> print(tax.get_mnemonic_translator(['MOUSE', 'YEAST', 'EREGS', 'CAPSP', ]))
{'CAPSP': [-50085752], 'EREGS': [284811], 'MOUSE': [10090], 'YEAST': [559292]}

"""