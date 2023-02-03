#!/usr/bin/env python

import sys
import os

import pickle
from collections import defaultdict, Counter
import sqlite3
import warnings
from ete3 import PhyloTree
from .build_db import update_db, update_mnemonic_codes

__all__ = ["Taxonomy", "is_taxadb_up_to_date"]


DB_VERSION = 3
DEFAULT_DB = os.path.expanduser(os.path.join('~', '.config', 'omataxonomy', 'taxonomy.sqlite'))
DEFAULT_DUMP_NAME = "omatax.tar.gz"
DEFAULT_RELEASE_ENV = "DARWIN_GENOMES_PATH"


def is_taxadb_up_to_date(dbfile=DEFAULT_DB):
    """Check if a valid and up-to-date taxa.sqlite database exists
    If dbfile= is not specified, DEFAULT_TAXADB is assumed
    """
    db = sqlite3.connect(dbfile)
    try:
        r = db.execute('SELECT version FROM stats;')
        version = r.fetchone()[0]
    except (sqlite3.OperationalError, ValueError, IndexError, TypeError):
        version = None

    db.close()

    if version != DB_VERSION:
        return False
    return True


class Taxonomy:

    def __init__(self, dbfile=None, taxdump_file=None, memory=False):

        if not dbfile:
            self.dbfile = DEFAULT_DB
        else:
            self.dbfile = dbfile

        if taxdump_file:
            self.update_taxonomy_database(taxdump_file)

        if dbfile != DEFAULT_DB and not os.path.exists(self.dbfile):
            print('omataxonomy database not present yet (first time used?)', file=sys.stderr)
            self.update_taxonomy_database(taxdump_file=taxdump_file)

        if not os.path.exists(self.dbfile):
            raise ValueError("Cannot open omataxonomy database: %s" % self.dbfile)

        self.db = None
        self._connect()

        if not is_taxadb_up_to_date(self.dbfile):
            print('omataxonomy database format is outdated. Upgrading', file=sys.stderr)
            self.update_taxonomy_database(taxdump_file)

        if memory:
            filedb = self.db
            self.db = sqlite3.connect(':memory:')
            filedb.backup(self.db)

    def update_taxonomy_database(self, taxdump_file=None):
        """Updates the omataxonomy database by downloading and parsing the latest
        taxonomies from NCBI and GTDB.
        :param None taxdump_file: an alternative location of the taxdump.tar.gz file.
        """
        if taxdump_file is None:
            taxdump_file = os.path.join(os.path.dirname(self.dbfile), DEFAULT_DUMP_NAME)
        update_db(self.dbfile, targz_file=taxdump_file)

    def _connect(self):
        self.db = sqlite3.connect(self.dbfile)

    def _translate_merged(self, all_taxids):
        conv_all_taxids = set(map(int, all_taxids))
        cmd = 'select taxid_old, taxid_new FROM merged WHERE taxid_old IN (%s)' % ','.join(map(str, all_taxids))

        result = self.db.execute(cmd)
        conversion = {}
        for old, new in result.fetchall():
            conv_all_taxids.discard(int(old))
            conv_all_taxids.add(int(new))
            conversion[int(old)] = int(new)
        return conv_all_taxids, conversion

    def get_rank(self, taxids):
        """return a dictionary converting a list of taxids into their
        corresponding omataxonomy rank

        :param taxids: an iterable of taxids to be converted"""

        all_ids = set(taxids)
        all_ids.discard(None)
        all_ids.discard("")
        query = ','.join(['"%s"' % v for v in all_ids])
        cmd = "select taxid, rank FROM species WHERE taxid IN (%s);" % query
        result = self.db.execute(cmd)
        id2rank = {}
        for tax, rank in result.fetchall():
            id2rank[tax] = rank
        return id2rank

    def get_lineage_translator(self, taxids):
        """Given a valid taxid number, return its corresponding lineage track as a
        hierarchically sorted list of parent taxids.
        """
        all_ids = set(taxids)
        all_ids.discard(None)
        all_ids.discard("")
        query = ','.join([f'"{v}"' for v in all_ids])
        result = self.db.execute(f'SELECT taxid, track FROM species WHERE taxid IN ({query});')
        id2lineages = {}
        for tax, track in result.fetchall():
            id2lineages[tax] = list(map(int, reversed(track.split(","))))
        return id2lineages

    def get_name_lineage(self, taxnames):
        """Given a valid taxname, return its corresponding lineage track as a
        hierarchically sorted list of parent taxnames.
        """
        name_lineages = {}
        name2taxid = self.get_name_translator(taxnames)
        for key, value in name2taxid.items():
            lineage = self.get_lineage(value[0])
            names = self.get_taxid_translator(lineage)
            name_lineages.update({key: [names[taxid] for taxid in lineage]})

        return name_lineages

    def get_lineage(self, taxid):
        """Given a valid taxid number, return its corresponding lineage track as a
        hierarchically sorted list of parent taxids.
        """
        if not taxid:
            return None
        taxid = int(taxid)
        result = self.db.execute('SELECT track FROM species WHERE taxid=%s' % taxid)
        raw_track = result.fetchone()
        if not raw_track:
            # perhaps is an obsolete taxid
            _, merged_conversion = self._translate_merged([taxid])
            if taxid in merged_conversion:
                result = self.db.execute('SELECT track FROM species WHERE taxid=%s' % merged_conversion[taxid])
                raw_track = result.fetchone()
            if not raw_track:
                raise ValueError("%s taxid not found" % taxid)
            else:
                warnings.warn(f"taxid {taxid} was translated into {merged_conversion[taxid]}")
        track = list(map(int, raw_track[0].split(",")))
        return list(reversed(track))

    def get_common_names(self, taxids):
        query = ','.join([f'"{v}"' for v in taxids])
        cmd = f"select taxid, common FROM species WHERE taxid IN ({query});"
        result = self.db.execute(cmd)
        id2name = {}
        for tax, common_name in result.fetchall():
            if common_name:
                id2name[tax] = common_name
        return id2name

    def get_mnemonic_names(self, taxids):
        query = ','.join([f"{v}" for v in taxids])
        cmd = f"select taxid, mnemonic FROM species WHERE taxid IN ({query});"
        result = self.db.execute(cmd)
        id2name = {}
        for tax, os_code in result.fetchall():
            if os_code:
                id2name[tax] = os_code
        return id2name

    def get_mnemonic_translator(self, codes):
        """Given a list of mnemonic codes, returns a dictionary
        with thier corresponding taxids."""
        code2id = {}
        code2origcode = {c.upper(): c for c in codes if c not in (None, '')}
        query = ','.join([f'"{c}"' for c in code2origcode.keys()])
        cmd = f'select mnemonic, taxid from species where mnemonic IN ({query})'
        result = self.db.execute(cmd)
        for code, taxid in result.fetchall():
            ocode = code2origcode[code.upper()]
            code2id.setdefault(ocode, []).append(taxid)
        return code2id

    def get_taxid_translator(self, taxids, try_synonyms=True):
        """Given a list of taxids, returns a dictionary with their corresponding
        scientific names.
        """

        all_ids = set(map(int, taxids))
        all_ids.discard(None)
        all_ids.discard("")
        query = ','.join([f'"{v}"' for v in all_ids])
        cmd = "select taxid, spname FROM species WHERE taxid IN (%s);" % query
        result = self.db.execute(cmd)
        id2name = {}
        for tax, spname in result.fetchall():
            id2name[tax] = spname

        return id2name

    def get_name_translator(self, names):
        """
        Given a list of taxid scientific names, returns a dictionary translating them into their corresponding taxids.
        Exact name match is required for translation.
        """

        name2id = {}
        name2origname = {}
        for n in names:
            name2origname[n.lower()] = n

        names = set(name2origname.keys())

        query = ','.join(['"%s"' % n for n in name2origname.keys()])
        cmd = 'select spname, taxid from species where spname IN (%s)' % query
        result = self.db.execute(cmd)
        for sp, taxid in result.fetchall():
            oname = name2origname[sp.lower()]
            name2id.setdefault(oname, []).append(taxid)
        missing = names - set([n.lower() for n in name2id.keys()])
        if missing:
            query = ','.join(['"%s"' % n for n in missing])
            result = self.db.execute('select spname, taxid from synonym where spname IN (%s)' % query)
            for sp, taxid in result.fetchall():
                oname = name2origname[sp.lower()]
                name2id.setdefault(oname, []).append(taxid)
        return name2id

    def translate_to_names(self, taxids):
        """
        Given a list of taxid numbers, returns another list with their corresponding scientific names.
        """
        id2name = self.get_taxid_translator(taxids)
        names = []
        for sp in taxids:
            names.append(id2name.get(sp, sp))
        return names

    def get_descendant_taxa(self, parent, intermediate_nodes=False, rank_limit=None, collapse_subspecies=False,
                            return_tree=False):
        """
        given a parent taxid or scientific species name, returns a list of all its descendants taxids.
        If intermediate_nodes is set to True, internal nodes will also be dumped.
        """
        try:
            taxid = int(parent)
        except ValueError:
            try:
                taxid = self.get_name_translator([parent])[parent][0]
            except KeyError:
                raise ValueError('%s not found!' % parent)

        # checks if taxid is a deprecated one, and converts into the right one.
        _, conversion = self._translate_merged([taxid])  # try to find taxid in synonyms table
        if conversion:
            taxid = conversion[taxid]

        with open(self.dbfile + ".traverse.pkl", "rb") as CACHED_TRAVERSE:
            prepostorder = pickle.load(CACHED_TRAVERSE)
        descendants = {}
        found = 0
        for tid in prepostorder:
            if tid == taxid:
                found += 1
            elif found == 1:
                descendants[tid] = descendants.get(tid, 0) + 1
            elif found == 2:
                break

        if not found:
            raise ValueError("taxid not found:%s" % taxid)
        elif found == 1:
            return [taxid]
        if rank_limit or collapse_subspecies or return_tree:
            tree = self.get_topology(list(descendants.keys()), intermediate_nodes=intermediate_nodes,
                                     collapse_subspecies=collapse_subspecies, rank_limit=rank_limit)
            if return_tree:
                return tree
            elif intermediate_nodes:
                return list(map(int, [n.name for n in tree.get_descendants()]))
            else:
                return map(int, [n.name for n in tree])

        elif intermediate_nodes:
            return [tid for tid, count in descendants.items()]
        else:
            return [tid for tid, count in descendants.items() if count == 1]

    def get_topology(self, taxids, intermediate_nodes=False, rank_limit=None, collapse_subspecies=False, annotate=True):
        """Given a list of taxid numbers, return the minimal pruned GTDB omataxonomy tree
        containing all of them.

        :param False intermediate_nodes: If True, single child nodes
            representing the complete lineage of leaf nodes are kept.
            Otherwise, the tree is pruned to contain the first common
            ancestor of each group.

        :param None rank_limit: If valid NCBI rank name is provided,
            the tree is pruned at that given level. For instance, use
            rank="species" to get rid of sub-species or strain leaf
            nodes.

        :param False collapse_subspecies: If True, any item under the
            species rank will be collapsed into the species upper
            node.
        """
        taxids, merged_conversion = self._translate_merged(taxids)
        if len(taxids) == 1:
            root_taxid = int(list(taxids)[0])
            with open(self.dbfile + ".traverse.pkl", "rb") as CACHED_TRAVERSE:
                prepostorder = pickle.load(CACHED_TRAVERSE)
            descendants = {}
            found = 0
            nodes = {}
            hit = 0
            visited = set()
            start = prepostorder.index(root_taxid)
            try:
                end = prepostorder.index(root_taxid, start + 1)
                subtree = prepostorder[start:end + 1]
            except ValueError:
                # If root taxid is not found in postorder, must be a tip node
                subtree = [root_taxid]
            leaves = set([v for v, count in Counter(subtree).items() if count == 1])
            nodes[root_taxid] = PhyloTree(name=str(root_taxid))
            current_parent = nodes[root_taxid]
            for tid in subtree:
                if tid in visited:
                    current_parent = nodes[tid].up
                else:
                    visited.add(tid)
                    nodes[tid] = PhyloTree(name=str(tid))
                    current_parent.add_child(nodes[tid])
                    if tid not in leaves:
                        current_parent = nodes[tid]
            root = nodes[root_taxid]
        else:
            taxids = set(map(int, taxids))
            sp2track = {}
            elem2node = {}
            id2lineage = self.get_lineage_translator(taxids)
            all_taxids = set()
            for lineage in id2lineage.values():
                all_taxids.update(lineage)
            id2rank = self.get_rank(all_taxids)
            for sp in taxids:
                track = []
                lineage = id2lineage[sp]

                for elem in lineage:
                    if elem not in elem2node:
                        node = elem2node.setdefault(elem, PhyloTree())
                        node.name = str(elem)
                        node.taxid = elem
                        node.add_feature("rank", str(id2rank.get(int(elem), "no rank")))
                    else:
                        node = elem2node[elem]
                    track.append(node)
                sp2track[sp] = track
            # generate parent child relationships
            for sp, track in sp2track.items():
                parent = None
                for elem in track:
                    if parent and elem not in parent.children:
                        parent.add_child(elem)
                    if rank_limit and elem.rank == rank_limit:
                        break
                    parent = elem
            root = elem2node[1]

        # remove onechild-nodes
        if not intermediate_nodes:
            for n in root.get_descendants():
                if len(n.children) == 1 and int(n.name) not in taxids:
                    n.delete(prevent_nondicotomic=False)

        if len(root.children) == 1:
            tree = root.children[0].detach()
        else:
            tree = root

        if collapse_subspecies:
            to_detach = []
            for node in tree.traverse():
                if node.rank == "species":
                    to_detach.extend(node.children)
            for n in to_detach:
                n.detach()

        if annotate:
            self.annotate_tree(tree)

        return tree

    def annotate_tree(self, t, taxid_attr="name", tax2name=None, tax2track=None, tax2rank=None):
        """Annotate a tree containing taxids as leaf names by adding the  'taxid',
        'sci_name', 'lineage', 'named_lineage' and 'rank' additional attributes.
        :param t: a Tree (or Tree derived) instance.
        :param name taxid_attr: Allows to set a custom node attribute
            containing the taxid number associated to each node (i.e.
            species in PhyloTree instances).
        :param tax2name,tax2track,tax2rank: Use these arguments to
            provide pre-calculated dictionaries providing translation
            from taxid number and names,track lineages and ranks.
        """

        taxids = set()
        for n in t.traverse():
            try:
                tid = int(getattr(n, taxid_attr))
            except (ValueError,AttributeError):
                pass
            else:
                taxids.add(tid)
        merged_conversion = {}

        taxids, merged_conversion = self._translate_merged(taxids)

        if not tax2name or taxids - set(map(int, list(tax2name.keys()))):
            tax2name = self.get_taxid_translator(taxids)
        if not tax2track or taxids - set(map(int, list(tax2track.keys()))):
            tax2track = self.get_lineage_translator(taxids)

        all_taxid_codes = set([_tax for _lin in list(tax2track.values()) for _tax in _lin])
        extra_tax2name = self.get_taxid_translator(list(all_taxid_codes - set(tax2name.keys())))
        tax2name.update(extra_tax2name)
        tax2common_name = self.get_common_names(tax2name.keys())

        if not tax2rank:
            tax2rank = self.get_rank(list(tax2name.keys()))

        n2leaves = t.get_cached_content()

        for n in t.traverse('postorder'):
            try:
                node_taxid = int(getattr(n, taxid_attr))
            except (ValueError, AttributeError):
                node_taxid = None

            n.add_feature('taxid', node_taxid)
            if node_taxid:
                if node_taxid in merged_conversion:
                    node_taxid = merged_conversion[node_taxid]
                n.add_features(sci_name=tax2name.get(node_taxid, getattr(n, taxid_attr, '')),
                            common_name=tax2common_name.get(node_taxid, ''),
                            lineage=tax2track.get(node_taxid, []),
                            rank=tax2rank.get(node_taxid, 'Unknown'),
                            named_lineage=[tax2name.get(tax, str(tax)) for tax in tax2track.get(node_taxid, [])])
            elif n.is_leaf():
                n.add_features(sci_name=getattr(n, taxid_attr, 'NA'),
                            common_name='',
                            lineage=[],
                            rank='Unknown',
                            named_lineage=[])
            else:
                lineage = self._common_lineage([lf.lineage for lf in n2leaves[n]])
                ancestor = lineage[-1]
                n.add_features(sci_name=tax2name.get(ancestor, str(ancestor)),
                            common_name=tax2common_name.get(ancestor, ''),
                            taxid=ancestor,
                            lineage=lineage,
                            rank=tax2rank.get(ancestor, 'Unknown'),
                            named_lineage=[tax2name.get(tax, str(tax)) for tax in lineage])

        return tax2name, tax2track, tax2rank

    def _common_lineage(self, vectors):
        occurrence = defaultdict(int)
        pos = defaultdict(set)
        for v in vectors:
            for i, taxid in enumerate(v):
                occurrence[taxid] += 1
                pos[taxid].add(i)

        common = [taxid for taxid, ocu in occurrence.items() if ocu == len(vectors)]
        if not common:
            return [""]
        else:
            sorted_lineage = sorted(common, key=lambda x: min(pos[x]))
            return sorted_lineage

        # OLD APPROACH:

        # visited = defaultdict(int)
        # for index, name in [(ei, e) for v in vectors for ei, e in enumerate(v)]:
        #     visited[(name, index)] += 1

        # def _sort(a, b):
        #     if a[1] > b[1]:
        #         return 1
        #     elif a[1] < b[1]:
        #         return -1
        #     else:
        #         if a[0][1] > b[0][1]:
        #             return 1
        #         elif a[0][1] < b[0][1]:
        #             return -1
        #     return 0

        # matches = sorted(visited.items(), _sort)

        # if matches:
        #     best_match = matches[-1]
        # else:
        #     return "", set()

        # if best_match[1] != len(vectors):
        #     return "", set()
        # else:
        #     return best_match[0][0], [m[0][0] for m in matches if m[1] == len(vectors)]

    def get_broken_branches(self, t, taxa_lineages, n2content=None):
        """Returns a list of GTDB lineage names that are not monophyletic in the
        provided tree, as well as the list of affected branches and their size.
        CURRENTLY EXPERIMENTAL
        """
        if not n2content:
            n2content = t.get_cached_content()

        tax2node = defaultdict(set)

        unknown = set()
        for leaf in t.iter_leaves():
            if leaf.sci_name.lower() != "unknown":
                lineage = taxa_lineages[leaf.taxid]
                for index, tax in enumerate(lineage):
                    tax2node[tax].add(leaf)
            else:
                unknown.add(leaf)

        broken_branches = defaultdict(set)
        broken_clades = set()
        for tax, leaves in tax2node.items():
            if len(leaves) > 1:
                common = t.get_common_ancestor(leaves)
            else:
                common = list(leaves)[0]
            if (leaves ^ set(n2content[common])) - unknown:
                broken_branches[common].add(tax)
                broken_clades.add(tax)

        broken_clade_sizes = [len(tax2node[tax]) for tax in broken_clades]
        return broken_branches, broken_clades, broken_clade_sizes

    # def annotate_tree_with_taxa(self, t, name2taxa_file, tax2name=None, tax2track=None, attr_name="name"):
    #     if name2taxa_file:
    #         names2taxid = dict([map(strip, line.split("\t"))
    #                             for line in open(name2taxa_file)])
    #     else:
    #         names2taxid = dict([(n.name, getattr(n, attr_name)) for n in t.iter_leaves()])

    #     not_found = 0
    #     for n in t.iter_leaves():
    #         n.add_features(taxid=names2taxid.get(n.name, 0))
    #         n.add_features(species=n.taxid)
    #         if n.taxid == 0:
    #             not_found += 1
    #     if not_found:
    #         print >>sys.stderr, "WARNING: %s nodes where not found within NCBI omataxonomy!!" %not_found

    #     return self.annotate_tree(t, tax2name, tax2track, attr_name="taxid")

    def update_mnemonic_codes(self, speclist=None, extra_file=None):
        """updates the mnemonic species codes.

        :param speclist: uniprot species code file. If not provided, the latest version
                         from UniProtKB is downloaded and used.

        :param extra_file: a path to a file that contains extra mappings in TSV format,
                           with two columns (CODE and taxid/sciname).
        """
        update_mnemonic_codes(self.db, speclist, extra_file)


class EnvReleaseTaxonomy(Taxonomy):
    def __init__(self, env_var=None):
        if env_var is None:
            env_var = DEFAULT_RELEASE_ENV
        try:
            db_file = os.path.join(os.environ[env_var], os.path.basename(DEFAULT_DB))
        except KeyError:
            raise KeyError(f"Environment variable '{env_var}' is not defined. "
                           f"Needed for {__class__.__name__} initialization")
        super().__init__(db_file)


if __name__ == "__main__":
    # from .. import PhyloTree
    tax = Taxonomy()
    #tax.update_taxonomy_database(DEFAULT_DUMP)

    descendants = tax.get_descendant_taxa('c__Thorarchaeia', collapse_subspecies=True, return_tree=True)
    print(descendants.write(features=None))
    print(descendants.get_ascii(attributes=['sci_name', 'taxid', 'rank']))

    import itertools
    taxids = list(itertools.chain.from_iterable(
        tax.get_name_translator(["p__Hydrothermarchaeota", "o__Peptococcales", "f__Korarchaeaceae", "GB_GCA_011358815.1"])
        .values()
    ))
    tree = tax.get_topology(taxids, intermediate_nodes=True, collapse_subspecies=True, annotate=True)
    print(tree.get_ascii(attributes=["taxid", "sci_name", "rank"]))

    tree = PhyloTree('((c__Thorarchaeia, c__Lokiarchaeia), s__Caballeronia udeis);',
                     sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
    tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
    print(tree.get_ascii(attributes=["taxid", "name", "sci_name", "rank"]))

    tree = PhyloTree('(RS_GCF_006228565.1, (Homo sapiens, Gallus gallus));',
                     sp_naming_function=lambda name: tax.get_name_translator([name])[name][0])
    tax2name, tax2track, tax2rank = tax.annotate_tree(tree, taxid_attr="species")
    print(tree.get_ascii(attributes=["taxid", "name", "sci_name", "rank"]))

    print(tax.get_name_lineage(['RS_GCF_006228565.1', 'GB_GCA_001515945.1', "Homo sapiens", "f__Leptotrichiaceae", "Gallus"]))

    print(tax.get_mnemonic_names([9606, 43715, 658031, 73382, 9823]))

    lin = tax.get_lineage(10090)
    t2n = tax.get_taxid_translator(lin)
    print([t2n[x] for x in lin])

    #tax.update_mnemonic_codes()
    print(tax.get_mnemonic_translator(['MOUSE', 'YEAST', 'ASHGO', 'CAPSP', ]))
