# Installation

omataxonomy requires Python 3.9 or newer. Install it from PyPI:

```bash
pip install omataxonomy
```

or add it to a project managed by [uv](https://docs.astral.sh/uv/):

```bash
uv add omataxonomy
```

To install the development version from GitHub:

```bash
pip install git+https://github.com/DessimozLab/omataxonomy.git
```

That is all for the package itself. The taxonomy data is not part of the
package. The first time you use it, omataxonomy downloads the data and builds
a local database (see {doc}`database`). Make sure that you have:

- internet access to `ftp.ncbi.nih.gov`, `data.ace.uq.edu.au` (GTDB) and
  `ftp.uniprot.org`,
- a few GB of free disk space. The database for GTDB release 232 takes about
  1.3 GB, and the downloads are stored temporarily while it is being built.
