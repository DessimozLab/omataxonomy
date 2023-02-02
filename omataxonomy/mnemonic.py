import re
import tempfile
import urllib.request
import os
import csv
from pathlib import Path
DEFAULT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt"


def iter_mnemonic_species_codes(url=None):
    if url is None:
        url = DEFAULT_URL
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfile = Path(tmpdir) / "speclist.txt"
        urllib.request.urlretrieve(url, tmpfile)
        end_header_seen = False
        spec_re = re.compile(r"(?P<code>[A-Z][A-Z0-9]{2,4})\s+[EBA]\s+(?P<taxid>\d+):")
        with open(tmpfile, "rt") as fh:
            for line in fh:
                if not end_header_seen:
                    if line.startswith("_____ _ _______  ___"):
                        end_header_seen = True
                    continue
                m = spec_re.match(line)
                if m is not None:
                    yield m.group('code'), m.group('taxid')


def iter_extra_mnemonic_species_codes(fname):
    if fname is None:
        fname = Path(os.getenv("DARWIN_GENOMES_PATH", ".")) / "extra_mnemonic_codes.tsv"
    fname = Path(fname)
    if fname.is_file():
        with open(fname, "rt") as fh:
            dialect = csv.Sniffer().sniff(fh.read(1000))
            fh.seek(0)
            reader = csv.reader(fh, dialect=dialect)
            yield from reader