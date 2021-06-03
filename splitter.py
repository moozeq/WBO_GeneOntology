from utils import read_seqs

if __name__ == '__main__':
    genes = read_seqs('data/genes_e_coli.fasta')
    genes = [gene.id for gene in genes]
    with open('data/population.txt', 'w') as f:
        f.write('\n'.join(genes))
    with open('data/studyA.txt', 'w') as f:
        f.write('\n'.join(gene for gene in genes if gene[-1] == 'A'))
    with open('data/studyB.txt', 'w') as f:
        f.write('\n'.join(gene for gene in genes if gene[-1] == 'B'))
