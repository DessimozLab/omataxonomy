# OMA Taxonomy 

omataxonomy is a library based on 
[ete's ncbi_taxonomy module](https://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html) 
that is used internally for the [OMA project](https://omabrowser.org), but that can also 
be used in different contexts.
Essentially, it combines the NCBI Taxonomy and the GTDB taxonomy (including all their genome as subspecies). 
For GTDB we generate stable taxon ids by hashing the scientific names. 

omataxonomy stores the data in a sqlite database under `${HOME}/.config/omataxonomy/` 
and therefor uses little resources when being used as not all the data will be loaded into memory.

## Install

OMA Taxonomy can be installed directly from pip:

    pip install omataxonomy

## Usage

    from omataxonomy import Taxonomy
    tax = Taxonomy()
    print(tax.get_name_lineage(['RS_GCF_006228565.1', 'GB_GCA_001515945.1', "f__Leptotrichiaceae", "Homo sapiens", "Gallus"]))

{'f__Leptotrichiaceae': ['root', 'd__Bacteria', 'p__Fusobacteriota', 'o__Fusobacteriales', 'f__Leptotrichiaceae'], 'Gallus': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Sauropsida', 'Sauria', 'Archelosauria', 'Archosauria', 'Dinosauria', 'Saurischia', 'Theropoda', 'Coelurosauria', 'Aves', 'Neognathae', 'Galloanserae', 'Galliformes', 'Phasianidae', 'Phasianinae', 'Gallus'], 'GB_GCA_001515945.1': ['root', 'd__Bacteria', 'p__Firmicutes_B', 'c__Moorellia', 'o__Desulfitibacterales', 's__Desulfitibacter sp001515945', 'GB_GCA_001515945.1'], 'Homo sapiens': ['root', 'cellular organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Homo', 'Homo sapiens'], 'RS_GCF_006228565.1': ['root', 'd__Bacteria', 'p__Firmicutes_B', 'c__Moorellia', 'o__Moorellales', 'f__Moorellaceae', 'g__Moorella', 's__Moorella thermoacetica', 'RS_GCF_006228565.1']}
