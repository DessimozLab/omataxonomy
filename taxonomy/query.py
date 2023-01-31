#!/usr/bin/env python

import sys
import os

import pickle
from collections import defaultdict, Counter
from itertools import chain
from tqdm import tqdm
import sqlite3
import tarfile
import warnings
from ete4 import PhyloTree

__all__ = ["Taxonomy", "is_taxadb_up_to_date"]

DB_VERSION = 2
DEFAULT_DB = os.path.join(os.path.dirname(os.path.realpath(__file__)), '.data', 'taxa.sqlite')
DEFAULT_DUMP = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'omatax.tar.gz')


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
            print('taxonomy database not present yet (first time used?)', file=sys.stderr)
            self.update_taxonomy_database(taxdump_file=DEFAULT_DUMP)

        if not os.path.exists(self.dbfile):
            raise ValueError("Cannot open taxonomy database: %s" % self.dbfile)

        self.db = None
        self._connect()

        if not is_taxadb_up_to_date(self.dbfile):
            print('taxonomy database format is outdated. Upgrading', file=sys.stderr)
            self.update_taxonomy_database(taxdump_file)

        if memory:
            filedb = self.db
            self.db = sqlite3.connect(':memory:')
            filedb.backup(self.db)

    def update_taxonomy_database(self, taxdump_file=None):
        """Updates the GTDB taxonomy database by downloading and parsing the latest
        gtdbtaxdump.tar.gz file from gtdbdump folder.
        :param None taxdump_file: an alternative location of the gtdbtaxdump.tar.gz file.
        """
        if not taxdump_file:
            update_db(self.dbfile)
        else:
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
        'return a dictionary converting a list of taxids into their corresponding GTDB taxonomy rank'

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
        name_lineages = []
        name2taxid = self.get_name_translator(taxnames)
        for key, value in name2taxid.items():
            lineage = self.get_lineage(value[0])
            names = self.get_taxid_translator(lineage)
            name_lineages.append({key: [names[taxid] for taxid in lineage]})

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
        code2origcode = {c.upper(): c for c in codes}
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
        """Given a list of taxid numbers, return the minimal pruned GTDB taxonomy tree
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
                        node.add_prop("rank", str(id2rank.get(int(elem), "no rank")))
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
                    if rank_limit and elem.props.get('rank') == rank_limit:
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
                if node.props.get('rank') == "species":
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

            n.add_prop('taxid', node_taxid)
            if node_taxid:
                if node_taxid in merged_conversion:
                    node_taxid = merged_conversion[node_taxid]
                n.add_props(sci_name=tax2name.get(node_taxid, getattr(n, taxid_attr, '')),
                            common_name=tax2common_name.get(node_taxid, ''),
                            lineage=tax2track.get(node_taxid, []),
                            rank=tax2rank.get(node_taxid, 'Unknown'),
                            named_lineage=[tax2name.get(tax, str(tax)) for tax in tax2track.get(node_taxid, [])])
            elif n.is_leaf():
                n.add_props(sci_name=getattr(n, taxid_attr, 'NA'),
                            common_name='',
                            lineage=[],
                            rank='Unknown',
                            named_lineage=[])
            else:
                lineage = self._common_lineage([lf.props.get('lineage') for lf in n2leaves[n]])
                ancestor = lineage[-1]
                n.add_props(sci_name=tax2name.get(ancestor, str(ancestor)),
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
    #         print >>sys.stderr, "WARNING: %s nodes where not found within NCBI taxonomy!!" %not_found

    #     return self.annotate_tree(t, tax2name, tax2track, attr_name="taxid")

    def update_mnemonic_codes(self, speclist=None, extra_file=None):
        """updates the mnemonic species codes.

        :param speclist: uniprot species code file. If not provided, the latest version
                         from UniProtKB is downloaded and used.

        :param extra_file: a path to a file that contains extra mappings in TSV format,
                           with two columns (CODE and taxid/sciname).
        """
        update_mnemonic_codes(self.db, speclist, extra_file)


def load_gtdb_tree_from_dump(tar):
    from ete4 import Tree
    # Download: gtdbdump/gtdbr202dump.tar.z
    parent2child = {}
    name2node = {}
    node2taxname = {}
    synonyms = set()
    name2rank = {}
    node2common = {}
    print("Loading node names...")
    unique_nocase_synonyms = set()
    for line in tar.extractfile("names.dmp"):
        line = str(line.decode())
        fields = [_f.strip() for _f in line.split("|")]
        nodename = fields[0]
        name_type = fields[3].lower()
        taxname = fields[1]

        # Clean up tax names so we make sure the don't include quotes. See https://github.com/etetoolkit/ete/issues/469
        taxname = taxname.rstrip('"').lstrip('"')

        if name_type in ("scientific name", "scientific_name"):
            node2taxname[nodename] = taxname
        if name_type == "genbank common name":
            node2common[nodename] = taxname
        elif name_type in set(["mnemonic_code", "ncbi_taxid",  "ncbi_organism_name", "genbank equivalent name",
                               "anamorph", "genbank synonym", "genbank anamorph", "teleomorph"]):
            if name_type == "ncbi_taxid":
                taxname = f"ncbi_taxid:{taxname}"

            # Keep track synonyms, but ignore duplicate case-insensitive names. See https://github.com/etetoolkit/ete/issues/469
            synonym_key = (nodename, taxname.lower())
            if synonym_key not in unique_nocase_synonyms:
                unique_nocase_synonyms.add(synonym_key)
                synonyms.add((nodename, taxname))

    print(len(node2taxname), "names loaded.")
    print(len(synonyms), "synonyms loaded.")

    print("Loading nodes...")
    for line in tar.extractfile("nodes.dmp"):
        line = str(line.decode())
        fields = line.split("|")
        nodename = fields[0].strip()
        parentname = fields[1].strip()
        n = Tree()
        n.name = nodename
        # n.taxname = node2taxname[nodename]
        n.add_prop('taxname', node2taxname[nodename])
        if nodename in node2common:
            n.add_prop('common_name', node2taxname[nodename])
        n.add_prop('rank', fields[2].strip())
        parent2child[nodename] = parentname
        name2node[nodename] = n
    print(len(name2node), "nodes loaded.")

    print("Linking nodes...")
    for node in name2node:
        if node == "1":
            t = name2node[node]
        else:
            parent = parent2child[node]
            parent_node = name2node[parent]
            parent_node.add_child(name2node[node])
    print("Tree is loaded.")
    return t, synonyms


def generate_table(t):
    OUT = open("taxa.tab", "w")
    for j, n in enumerate(t.traverse()):
        if j % 1000 == 0:
            print("\r", j, "generating entries...", end=' ')
        temp_node = n
        track = []
        while temp_node:
            temp_rank = temp_node.props.get("rank")
            if temp_rank not in (None, "None"):
                track.append(temp_node.name)
            temp_node = temp_node.up
        if n.up:
            print('\t'.join(
                [n.name, n.up.name, n.props.get('taxname'), n.props.get("common_name", ''), n.props.get("rank"),
                 ','.join(track)]), file=OUT)
        else:
            print('\t'.join([n.name, "", n.props.get('taxname'), n.props.get("common_name", ''), n.props.get("rank"),
                             ','.join(track)]), file=OUT)
    OUT.close()


def update_db(dbfile, targz_file=None):
    basepath = os.path.split(dbfile)[0]
    if basepath and not os.path.exists(basepath):
        os.mkdir(basepath)

    try:
        tar = tarfile.open(targz_file, 'r')
    except:
        raise ValueError("Please provide taxa dump tar.gz file")

    t, synonyms = load_gtdb_tree_from_dump(tar)

    prepostorder_full = [int(node.name) for post, node in t.iter_prepostorder()]
    with open(dbfile + ".full.traverse.pkl", "wb") as fh:
        pickle.dump(prepostorder_full, fh, 5)
    prepostorder_lineage = [int(node.name) for post, node in t.iter_prepostorder() if node.props.get('rank') not in (None, "None")]
    with open(dbfile + ".traverse.pkl", "wb") as fh:
        pickle.dump(prepostorder_lineage, fh, 5)

    print("Updating database: %s ..." % dbfile)
    generate_table(t)

    with open("syn.tab", "w") as SYN:
        SYN.write('\n'.join(["%s\t%s" %(v[0],v[1]) for v in synonyms]))

    with open("merged.tab", "w") as merged:
        for line in tar.extractfile("merged.dmp"):
            line = str(line.decode())
            out_line = '\t'.join([_f.strip() for _f in line.split('|')[:2]])
            merged.write(out_line+'\n')
    try:
        upload_data(dbfile)
    except:
        raise
    else:
        os.system("rm taxa.tab")
        # remove only downloaded taxdump file
        if not targz_file:
            os.system("rm gtdbtaxdump.tar.gz")


def upload_data(dbfile):
    print()
    print('Uploading to', dbfile)
    basepath = os.path.split(dbfile)[0]
    if basepath and not os.path.exists(basepath):
        os.mkdir(basepath)

    db = sqlite3.connect(dbfile)

    create_cmd = """
    DROP TABLE IF EXISTS stats;
    DROP TABLE IF EXISTS species;
    DROP TABLE IF EXISTS synonym;
    DROP TABLE IF EXISTS merged;
    CREATE TABLE stats (version INT PRIMARY KEY);
    CREATE TABLE species (taxid INT PRIMARY KEY, 
                          parent INT, 
                          spname VARCHAR(50) COLLATE NOCASE, 
                          common VARCHAR(50) COLLATE NOCASE, 
                          rank VARCHAR(50),
                          mnemonic VARCHAR(5), 
                          track TEXT);
    CREATE TABLE synonym (taxid INT, 
                          spname VARCHAR(50) COLLATE NOCASE, 
                          PRIMARY KEY (spname, taxid));
    CREATE TABLE merged (taxid_old INT, taxid_new INT);
    CREATE INDEX spname1 ON species (spname COLLATE NOCASE);
    CREATE INDEX spname2 ON synonym (spname COLLATE NOCASE);
    CREATE INDEX mnemonic ON species (mnemonic);
    """
    for cmd in create_cmd.split(';'):
        db.execute(cmd)
    print()

    db.execute("INSERT INTO stats (version) VALUES (%d);" % DB_VERSION)
    db.commit()

    try:
        with open("syn.tab", 'rt') as fh:
            for i, line in tqdm(enumerate(fh), desc="inserting synonymes"):
                taxid, spname = line.strip('\n').split('\t')
                db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))
        db.commit()
    except FileNotFoundError:
        print("no synonym table found. skipping", file=sys.stderr)

    try:
        with open("merged.tab", "rt") as fh:
            for i, line in tqdm(enumerate(fh), desc="inserting taxid merges"):
                taxid_old, taxid_new = line.strip('\n').split('\t')
                db.execute("INSERT INTO merged (taxid_old, taxid_new) VALUES (?, ?);", (taxid_old, taxid_new))
        db.commit()
    except FileNotFoundError:
        print("no merged.tab found. skipping", file=sys.stderr)

    with open("taxa.tab", "rt") as fh:
        for i, line in tqdm(enumerate(fh), desc="inserting taxids"):
            taxid, parentid, spname, common, rank, lineage = line.strip('\n').split('\t')
            db.execute("INSERT INTO species (taxid, parent, spname, common, rank, mnemonic, track) VALUES (?, ?, ?, ?, ?, ?, ?);",
                       (taxid, parentid, spname, common, rank, "", lineage))
    db.commit()
    update_mnemonic_codes(db)
    print("\rdatabase created", file=sys.stderr)


def update_mnemonic_codes(db, speclist=None, extra_file=None):
    from taxonomy.mnemonic import iter_mnemonic_species_codes, iter_extra_mnemonic_species_codes
    db.execute("UPDATE species set mnemonic = ''")
    try:
        for i, (os_code, taxid) in tqdm(enumerate(chain(iter_mnemonic_species_codes(speclist),
                                                  iter_extra_mnemonic_species_codes(extra_file))),
                                        desc="inserting mnemonic codes"):
            try:
                taxid = int(taxid)
                result = db.execute("SELECT taxid FROM synonym WHERE spname=? ORDER BY taxid DESC", (f"ncbi_taxid:{taxid}",))
                e = result.fetchone()
                if e is not None:
                    taxid = e[0]
                db.execute("UPDATE species SET mnemonic = ? WHERE taxid = ?", (os_code, taxid))
            except ValueError:
                db.execute("UPDATE species SET mnemonic = ? where spname = ?", (os_code, taxid))
        db.commit()
    except Exception as e:
        print(f"update mnemoinc codes failed: {e}")
        db.rollback()


if __name__ == "__main__":
    # from .. import PhyloTree
    tax = Taxonomy()
    #tax.update_taxonomy_database(DEFAULT_DUMP)

    descendants = tax.get_descendant_taxa('c__Thorarchaeia', collapse_subspecies=True, return_tree=True)
    print(descendants.write(properties=None))
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

    print(tax.get_name_lineage(['RS_GCF_006228565.1', 'GB_GCA_001515945.1', "Homo sapiens", "Gallus"]))

    print(tax.get_mnemonic_names([9606, 43715, 658031, 73382, 9823]))

    lin = tax.get_lineage(10090)
    t2n = tax.get_taxid_translator(lin)
    print([t2n[x] for x in lin])

    # tax.update_mnemonic_codes()
    print(tax.get_mnemonic_translator(['MOUSE', 'YEAST', 'ASHGO', 'CAPSP', ]))
