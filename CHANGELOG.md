# Changelog

All notable changes to omataxonomy are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2026-10-16

### Upgrade notes

- Existing databases are rebuilt automatically. `DB_VERSION` went from 3 to 4,
  so the first `Taxonomy()` call after upgrading rebuilds the database. This
  downloads the NCBI and GTDB data again.
- `omatax.tar.gz` archives built before 0.5.0 can contain duplicate taxids, and
  loading one now raises a `ValueError`. Delete the old archive, or let
  `update_db` generate a new one.
- GTDB taxa and genomes that were involved in a taxid collision get new
  taxids. Taxids are reproducible for a given GTDB release, but a taxon can get
  a different taxid in another release if the set of collisions changes.

### Added

- `Taxonomy`, `update_db` and `build_combined_tarball` take a `gtdb_release`
  argument, for example `232`, `"r232"` or `"214.1"`. The default is still the
  latest release. Before this, only `latest` produced a valid download URL.
- Documentation on Read the Docs: <https://omataxonomy.readthedocs.io>, with
  background, installation, database building and a tutorial. The examples are
  tested as doctests.
- Support for Python 3.13, 3.14 and 3.15.

### Changed

- The project moved from Poetry to `uv` (build backend `uv_build`), and
  `uv.lock` is now committed. CI tests Python 3.9–3.14, each with locked,
  highest and lowest-direct dependency resolution.
- The license metadata is now the SPDX identifier `GPL-3.0-only`, and project
  URLs for the docs, the repository and this changelog were added.
- Releases are published to PyPI with trusted publishing.

### Fixed

- GTDB taxid hash collisions. GTDB taxids are a 31-bit hash of the taxon name
  or genome accession. In GTDB r232, 318 taxids were shared by two keys. The
  entry written last silently replaced the other one, so the replaced taxon
  could not be found, and its genomes could end up under an unrelated node.
  Collisions are now resolved deterministically. The alphabetically first key
  keeps the hash, and the others get a new hash. Taxids that never collided
  keep their values.
- `build_combined_tarball` handles a relative target path. Before, the archive
  was written into a temporary directory and deleted along with it.
- `update_db` no longer fails when given a bare archive file name.
- The md5 checksum of the NCBI taxdump is also checked right after downloading,
  and a mismatch raises an `IOError`.
- `update_mnemonic_codes(speclist=...)` reads a local `speclist.txt`. Before, a
  local file was never loaded.
- The inverted `target_folder` check in `download_gtdb_release`.
- A malformed logging call in `upload_data`.

## [0.4.0] - 2025-05-20

### Changed

- Upgraded from ete3 to ete4. Several functions that returned `PhyloTree`
  objects now return ete4 `Tree` objects. The APIs of ete3 and ete4 differ, so
  this may break your code.
- Upgraded numpy from 1.26 to 2.0.

### Added

- Several new tests.

[0.5.0]: https://github.com/DessimozLab/omataxonomy/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/DessimozLab/omataxonomy/compare/v0.3.2...v0.4.0
