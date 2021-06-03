import gzip
import shutil
from pathlib import Path

import requests


class Term:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    @staticmethod
    def error(msg: str):
        print(f'{Term.FAIL}[ERROR] {msg}{Term.ENDC}')

    @staticmethod
    def ok(msg: str):
        print(f'{Term.OKGREEN}[OK] {msg}{Term.ENDC}')

    @staticmethod
    def info(msg: str):
        print(f'[INFO] {msg}')


def read_seqs(filename: str):
    from Bio import SeqIO
    seqs = [record for record in SeqIO.parse(filename, 'fasta')]
    return seqs


def write_seqs(records: list, filename: str):
    from Bio import SeqIO
    SeqIO.write(records, filename, 'fasta')
    return filename


def download_and_unpack(url: str, out: str):
    unpack = True if url.endswith('.gz') else False
    out_target = f'{out}.gz' if unpack else out
    if not Path(out).exists():
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(out_target, 'wb') as f:
                for chunk in r.iter_content(chunk_size=16384):
                    f.write(chunk)
        if unpack:
            with gzip.open(out_target, 'rb') as gf:
                with open(out, 'wb') as f:
                    shutil.copyfileobj(gf, f)
            Path(out_target).unlink()

    return Path(out).exists()
