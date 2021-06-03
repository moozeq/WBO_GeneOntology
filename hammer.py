import json
import subprocess
from pathlib import Path

from Bio.SearchIO.HmmerIO import Hmmer3TabParser, Hmmer3TextParser

from utils import Term, download_and_unpack


class Hammer:
    EXTENSIONS = ['h3f', 'h3i', 'h3m', 'h3p']

    def __init__(self, db_name: str):
        self.db_name = db_name

    def hmmscann(self, input_file: str, out: str, evalue: float = 1e-20):
        if Path(out).exists():
            return out
        Term.info(f'Running `hmmscan` for = {input_file}')
        # subprocess.run(['hmmscan', '--tblout', out, '-E', str(evalue), self.db_name, input_file])
        subprocess.run(['hmmscan', '-o', out, '-E', str(evalue), self.db_name, input_file])
        if Path(out).exists():
            Term.ok(f'Success, results saved at = {out}')
            return out
        return ''

    @classmethod
    def hmmer_db_from_url(cls, url: str = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'):
        hmmerdb = 'hmmerdb/Pfam-A.hmm'
        if not Path(hmmerdb).exists():
            Term.info(f'HMMER Pfam database not found, downloading...')
            if download_and_unpack(url, hmmerdb):
                Term.ok(f'HMMER Pfam database stored at = {hmmerdb}')

        if any(Path(f'{hmmerdb}.{ext}').exists() for ext in Hammer.EXTENSIONS):
            Term.info(f'Loading BLAST database from: {hmmerdb}')
            return cls(hmmerdb)
        Term.info(f'HMMER Pfam database found, but need to be prepared by `hmmpress` first')
        subprocess.run(['hmmpress', hmmerdb])
        if any(Path(f'{hmmerdb}.{ext}').exists() for ext in Hammer.EXTENSIONS):
            Term.info(f'Loading BLAST database from: {hmmerdb}')
            return cls(hmmerdb)

    @staticmethod
    def parse_table(table_file: str):
        with open(table_file) as f:
            parser = Hmmer3TabParser(f)
            return {entry.id: entry.hit_keys for entry in parser}

    @staticmethod
    def parse_text(text_file: str):
        with open(text_file) as f:
            parser = Hmmer3TextParser(f)
            return {
                entry.id: {
                    'hmmer': {
                        hit.id: hit.description
                        for hit in entry.hits
                    }
                }
                for entry in parser
                if entry.hit_keys
            }


if __name__ == '__main__':
    hmmer = Hammer.hmmer_db_from_url()
    hmmer.hmmscann('data/blasted.fasta', 'data/hmmscanned.txt')
    entries = Hammer.parse_text('data/hmmscanned.txt')
    with open('data/study_annotated.json', 'w') as f:
        json.dump(entries, f, indent=4)
    with open('data/studyPfam.txt', 'w') as f:
        f.write('\n'.join(entries.keys()))
    Term.info(f'Got {len(entries)} genes for study')
