import collections
import tempfile
import os
import dendropy
import numpy
import csv
import re
import tarfile
from urllib.request import urlretrieve
import logging
from hashlib import md5

logger = logging.getLogger(__name__)

GTDB_DOWNLOAD_BASE_URL = "https://data.gtdb.ecogenomic.org/releases"
NCBI_DOMAIN_ROOT = {"Bacteria": 2, "Archaea": 2157}
PREFIX_2_RANK = {"d__": "superkingdom", "p__": "phylum", "c__": "class", "o__": "order",
                 "f__": "family", "g__": "genus", "s__": "species"}


def _load_meta(fn):
    with open(fn, 'rt', newline="") as fh:
        reader = csv.DictReader(fh, dialect="excel-tab")
        meta = {row['accession']: row for row in reader}
    return meta


def build_tree_dump(tree_fn, meta_fn, offset):
    tree = dendropy.Tree.get_from_path(tree_fn, schema="newick", preserve_underscores=True)
    meta = _load_meta(meta_fn)
    tree = extend_tree_with_all_genomes(tree, meta)
    tree = _fix_tree_single_nodes_per_rank(tree)
    with open("gtdbnodes.dmp", "at") as nodes_dump, open("gtdbnames.dmp", "at") as names_dump:
        write_dump_files(tree, offset, nodes_dump, names_dump)


def hash_to_int(val):
    h = md5(str(val).encode('utf-8'))
    return -abs(int(numpy.frombuffer(h.digest()[:4], dtype='i4')[0]))


def write_dump_files(tree, offset, node_dump, names_dump):
    re_label = re.compile(r"(?P<bootstrap>[\d.]+)?:?(?P<name>[\w -]+)?")
    for node in tree.preorder_node_iter():

        offset -= 1
        if node.label is not None:
            m = re_label.match(node.label)
            if m.group("name") is not None:
                node.annotations.add_new("scientific_name", m.group("name"))
                node.annotations.add_new("rank", PREFIX_2_RANK[m.group("name")[:3]])
                node.annotations.add_new("taxid", hash_to_int(m.group('name')))
            else:
                taxid = offset
                offset -= 1
                node.annotations.add_new("scientific_name", f"x__tax{taxid}")
                node.annotations.add_new("taxid", taxid)

        rank = node.annotations.get_value('rank')
        taxid = node.annotations.require_value('taxid')
        try:
            p = node.parent_node
            if p is None:
                raise AttributeError()
        except AttributeError:
            parent_taxid = 1
            lineage_parent_taxid = 1
            lineage_parent_dist = 0
        else:
            parent_taxid = p.annotations.get_value('taxid')
            lineage_parent_taxid = 0
            lineage_parent_dist = node.edge_length
            # get the lineage parent for nodes that are part of the lineage
            while rank is not None:
                if p.annotations.get_value('rank') is not None:
                    lineage_parent_taxid = p.annotations.get_value('taxid')
                    break
                lineage_parent_dist += p.edge_length
                p = p.parent_node
        # tax_id					-- node id in GenBank omataxonomy database
        # parent tax_id				-- parent node id in GenBank omataxonomy database
        # rank					-- rank of this node (superkingdom, kingdom, ...)
        # embl code				-- locus-name prefix; not unique
        # division id				-- see division.dmp file
        # inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
        # genetic code id				-- see gencode.dmp file
        # inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
        # mitochondrial genetic code id		-- see gencode.dmp file
        # inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
        # GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        # hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        # comments				-- free-text comments and citations
        ### additional fields
        #  length of branch in GTDB tree to parent node
        #  parent node of lineage (skipping intermediate nodes)
        #  length of branch in GTDB tree to next lineage node == sum of branches to next lineage node
        #  is_reference_genome (1 or 0)    -- 1 if it is the reference genome of a species
        is_ref = 1 if node.annotations.get_value('is_reference', False) else 0
        nodes_buffer = [taxid, parent_taxid, rank, "XX", 0, 0, 11, 1, 1, 0, (0 if rank is not None else 1),
                        0, "", node.edge_length, lineage_parent_taxid, lineage_parent_dist, is_ref]
        node_dump.write("\t|\t".join(map(str, nodes_buffer)) + "\n")
        for key in ("scientific_name", "ncbi_taxid", "mnemonic_code", "ncbi_organism_name"):
            val = node.annotations.get_value(key)
            if val is None:
                continue
            buf = [taxid, val, "", key]
            names_dump.write("\t|\t".join(map(str, buf)) + "\n")


def _fix_tree_single_nodes_per_rank(tree):
    for node in tree.preorder_internal_node_iter():
        try:
            bootstrap, name = node.label.split(":")
        except ValueError:
            continue
        if len(name.split(':')) > 1:
            tok = name.split(":")
            p = node.parent_node
            p.remove_child(node)
            for i, lev in enumerate(tok):
                n = dendropy.Node(label=f"{bootstrap}:{lev}")
                n.edge_length = node.edge_length if i == 0 else 0
                p.add_child(n)
                p = n
            for c in node.child_node_iter():
                n.add_child(c)
    return tree


def extend_tree_with_all_genomes(tree, meta):
    same_rep = collections.defaultdict(list)
    for g in meta.values():
        same_rep[g['gtdb_genome_representative']].append(g)

    for n in tree.leaf_node_iter():
        parent = n.parent_node
        parent.remove_child(n)

        s_node = dendropy.Node(label=meta[n.taxon.label]['gtdb_taxonomy'].split(';')[-1])
        s_node.edge_length = n.edge_length
        for res in same_rep[n.taxon.label]:
            c = dendropy.Node(taxon=dendropy.Taxon(res['accession']))
            c.edge_length = 0.001
            c.annotations.add_new("is_reference", res['accession'] == n.taxon.label)
            c.annotations.add_new("rank", "subspecies")
            c.annotations.add_new("scientific_name", res['accession'])
            c.annotations.add_new("taxid", hash_to_int(res['accession']))
            for tag in ("ncbi_taxid", "ncbi_organism_name", "ncbi_strain_identifiers", "ncbi_genbank_assembly_accession"):
                c.annotations.add_new(tag, res[tag])
            s_node.add_child(c)
        parent.add_child(s_node)
    return tree


def download_gtdb_release(rel=None, dom=None, target_folder=None):
    from . import cwd
    if dom is None:
        dom = ["ar53", "bac120"]
    if isinstance(dom, str):
        dom = [dom]
    if rel is None:
        rel = "latest"
    if target_folder is not None:
        target_folder = "./"
    with cwd(target_folder):
        for d in dom:
            tree_fn, meta_fn = download_dom(d, rel)
            build_tree_dump(tree_fn, meta_fn, offset=-10000 if d.startswith("ar") else -50000)


def download_dom(dom, rel):
    logger.info(f"download tax release '{rel}' from {GTDB_DOWNLOAD_BASE_URL}")
    files = [f"{dom}{fn}" for fn in (".tree.tar.gz", "_metadata.tar.gz")]
    extracted_files = []
    for file in files:
        url = f"{GTDB_DOWNLOAD_BASE_URL}/{rel}/{file}"
        urlretrieve(url, file)
        tar_fh = tarfile.open(file)
        tar_files = tar_fh.getmembers()
        tar_fh.extractall(members=tar_files)
        extracted_files.extend([z.name for z in tar_files])
        tar_fh.close()
    return extracted_files


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    download_gtdb_release()