import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List

from utils import read_seqs, Term, write_seqs


@dataclass
class Entry:
    query: str
    subject: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    q_start: int
    q_end: int
    s_start: int
    s_end: int
    evalue: float
    score: float

    def __post_init__(self):
        self.identity = float(self.identity)
        self.alignment_length = int(self.alignment_length)
        self.mismatches = int(self.mismatches)
        self.gap_opens = int(self.gap_opens)
        self.q_start = int(self.q_start)
        self.q_end = int(self.q_end)
        self.s_start = int(self.s_start)
        self.s_end = int(self.s_end)
        self.evalue = float(self.evalue)
        self.score = float(self.score)


class Blaster:
    BLAST_BINS_PATH = '/usr/local/ncbi/blast/bin/'
    EXTENSIONS = ['ndb', 'nhr', 'nin', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto']
    ENV = {'PATH': BLAST_BINS_PATH}

    def __init__(self, db_name: str, input_file: str, dbtype: str, gencode: int):
        self.db_name = db_name
        self.input_file = input_file
        self.dbtype = dbtype
        self.gencode = gencode

    def search_tblastn(self, input_file: str, out: str, evalue: float = 1e-20):
        """
        E-value by default is 1e-20 so very similar sequences will be found.

        Output format is simple .tsv
        """
        from Bio.Blast.Applications import NcbitblastnCommandline

        tblastn_cline = NcbitblastnCommandline(
            query=input_file,
            db=self.db_name,
            db_gencode=self.gencode,
            out=out,
            evalue=evalue,
            outfmt=6
        )
        stdout, stderr = tblastn_cline(env=Blaster.ENV)
        if stderr:
            Term.error(f'Searching using `tblastn` failed:\n{stderr}')
        else:
            Term.ok(f'Searching using `tblastn` succeed, results at: {out}')

    def gather_hits_seqs(self, results_tsv: str, out: str, identity_threshold: float = 90.0, aa: bool = False) -> str:
        hits = Blaster.get_hits(results_tsv, identity_threshold)
        hits_genes = set(hit.subject for hit in hits)

        seqs = read_seqs(self.input_file)
        seqs = [
            seq
            for seq in seqs
            if seq.id in hits_genes
        ]
        if aa:
            aa_seqs = []
            for seq in seqs:
                aa_seq = seq.translate(table=self.gencode)
                aa_seq.id = seq.id
                aa_seq.description = seq.description
                aa_seqs.append(aa_seq)
            seqs = aa_seqs

        return write_seqs(seqs, out)

    @classmethod
    def blastdb_from_file(cls, input_file: str, dbtype: str = 'nucl', gencode: int = 1):
        """
        Pick right gencode:
            - 1:    standard
            - 11:   bacterial
        """
        from Bio.Blast.Applications import NcbimakeblastdbCommandline
        title = Path(input_file).stem
        out = f'blastdb/{title}'

        if any(Path(f'{out}.{ext}').exists() for ext in Blaster.EXTENSIONS):
            Term.info(f'Loading BLAST database from: {out}')
            return cls(out, input_file, dbtype, gencode)

        Path(out).parent.mkdir(parents=True, exist_ok=True)

        build_db_cline = NcbimakeblastdbCommandline(
            dbtype=dbtype,
            input_file=input_file,
            parse_seqids=True,
            out=out,
            title=title
        )
        stdout, stderr = build_db_cline(env=Blaster.ENV)
        if stderr:
            Term.error(f'Building BLAST db failed:\n{stderr}')
        else:
            Term.ok(f'Building BLAST db succeed, database at: {out}')
        return cls(out, input_file, dbtype, gencode)

    @staticmethod
    def get_hits(results_tsv: str, identity_threshold: float = 90.0) -> List[Entry]:
        with open(results_tsv, newline='') as res:
            reader = csv.reader(res, delimiter='\t')
            hits = [
                entry
                for row in reader
                if (entry := Entry(*row)).identity >= identity_threshold
            ]
            return hits


if __name__ == '__main__':
    db = Blaster.blastdb_from_file('data/genes_e_coli.fasta', gencode=11)
    db.search_tblastn('data/protein_fragments.fasta', 'data/blasted.tsv')
    db.gather_hits_seqs('data/blasted.tsv', 'data/blasted.fasta', aa=True)
