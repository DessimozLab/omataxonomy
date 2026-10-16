# Building the taxonomy database

All queries use a local SQLite database. This page explains how it is built,
where it is stored and how to choose the GTDB release.

## The first run

The first time you create a {class}`~omataxonomy.Taxonomy` object, the
database is built:

```python
from omataxonomy import Taxonomy

tax = Taxonomy()
```

This

1. downloads the current NCBI taxdump (`taxdump.tar.gz`) and checks it
   against its md5 checksum,
2. downloads the GTDB metadata tables for Bacteria (`bac120_metadata`) and
   Archaea (`ar53_metadata`) of the **latest** GTDB release,
3. assigns the negative taxids to the GTDB taxa and genomes and merges both
   taxonomies into one tree,
4. writes the database and fills in the UniProt species mnemonics from
   UniProt's `speclist.txt`.

This takes a while, mostly for the downloads and for loading more than four
million nodes. The result is stored as two files that must stay together:

```text
~/.config/omataxonomy/taxonomy.sqlite               # the database
~/.config/omataxonomy/taxonomy.sqlite.traverse.pkl  # tree traversal used by get_descendant_taxa()
```

Later calls to `Taxonomy()` just open the existing database.

## Storing the database elsewhere

Pass `dbfile` to use a different location. If the file does not exist yet,
it is built there:

```python
tax = Taxonomy(dbfile="/data/taxonomy/taxonomy.sqlite")
```

## Choosing a GTDB release

By default, the latest GTDB release is used. Use `gtdb_release` to choose a
specific one. A release can be given as `232`, `"232"`, `"r232"`,
`"release232"`, or with a minor version such as `"214.1"`:

```python
tax = Taxonomy(dbfile="/data/taxonomy/gtdb-r232/taxonomy.sqlite", gtdb_release=232)
```

Choosing a release makes your results reproducible: the GTDB part of the
tree, and thus all GTDB taxids, stays the same. omataxonomy's own test suite
pins a release for this reason.

:::{important}
`gtdb_release` is only used **when a database is built**. If `dbfile`
already exists, it is opened as it is, whatever release it was built from.
The database does not record which GTDB release it contains. Put the release
in the path, as in the example above, to keep track of it.
:::

Only the GTDB release can be chosen. The NCBI Taxonomy is always the
current one at the time the database is built.

### Several releases side by side

Databases are independent files, so you can keep one per release and open
whichever you need:

```python
tax_r220 = Taxonomy(dbfile="/data/taxonomy/gtdb-r220/taxonomy.sqlite", gtdb_release=220)
tax_r232 = Taxonomy(dbfile="/data/taxonomy/gtdb-r232/taxonomy.sqlite", gtdb_release=232)
```

The taxid of a GTDB taxon is computed from its name, so a taxon has the same
taxid in both databases as long as GTDB did not rename it.

## Updating an existing database

To rebuild a database in place with fresh data, call
{meth}`~omataxonomy.Taxonomy.update_taxonomy_database`. This downloads the
current NCBI Taxonomy and the requested (or latest) GTDB release again:

```python
tax = Taxonomy()
tax.update_taxonomy_database()                   # latest GTDB release
tax.update_taxonomy_database(gtdb_release=232)   # a specific GTDB release
```

When a new version of omataxonomy changes the database format, the database
is rebuilt automatically the next time it is opened.

## Building from a prepared archive

The downloads can also be separated from building the database. This is
useful when you want to build several databases from the same data, or
build the database on a machine without internet access.
{func}`~omataxonomy.build_combined_tarball` writes the merged NCBI and GTDB
dump files to a single archive:

```python
from omataxonomy import build_combined_tarball

build_combined_tarball("omatax-gtdb-r232.tar.gz", gtdb_release=232)
```

Then pass the archive as `taxdump_file` to build a database from it:

```python
tax = Taxonomy(dbfile="/data/taxonomy/gtdb-r232/taxonomy.sqlite",
               taxdump_file="omatax-gtdb-r232.tar.gz")
```

:::{note}
Whenever `taxdump_file` is given, the database is rebuilt from that archive.
Pass it only once, to build the database, and open the database later with
`Taxonomy(dbfile=...)` alone.
:::

Archives made by older versions of omataxonomy can contain duplicated GTDB
taxids. Such archives are rejected with a `ValueError` and have to be
rebuilt.

## Keeping the database in memory

`Taxonomy(memory=True)` copies the whole database into an in-memory SQLite
database. This makes queries faster when you run many of them, but it needs
as much RAM as the size of the database file.

## Database from an environment variable

{class}`~omataxonomy.EnvReleaseTaxonomy` opens `taxonomy.sqlite` in the
directory named by an environment variable (`DARWIN_GENOMES_PATH` by
default). This is how the OMA pipeline uses one database for a whole OMA
release:

```python
from omataxonomy import EnvReleaseTaxonomy

tax = EnvReleaseTaxonomy()                 # $DARWIN_GENOMES_PATH/taxonomy.sqlite
tax = EnvReleaseTaxonomy("MY_TAXDB_DIR")   # $MY_TAXDB_DIR/taxonomy.sqlite
```

## Updating mnemonic codes

The UniProt species mnemonics are loaded while the database is built. To
reload them, e.g. after a new UniProt release, call
{meth}`~omataxonomy.Taxonomy.update_mnemonic_codes`. You can pass a local
copy of `speclist.txt`, and a tab-separated file with extra codes. Its first
column is the code and its second column a taxid or a scientific name:

```python
tax.update_mnemonic_codes(speclist="speclist.txt", extra_file="extra_codes.tsv")
```

If no `extra_file` is given, omataxonomy reads
`$DARWIN_GENOMES_PATH/extra_mnemonic_codes.tsv` if that file exists.
