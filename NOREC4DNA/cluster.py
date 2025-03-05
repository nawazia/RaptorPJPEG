import argparse
import os

from clustering_dna_storage.utils import get_fastq_records, get_badread_strand_id, create_fasta_file
from clustering_dna_storage.clustering import Clustering
from clustering_dna_storage.checksum import CheckSum4

def cluster(path, oligo_len, num_of_oligos):
    records = get_fastq_records(fastq_filepath=path)
    reads = [str(i.seq) for i in records]
    ids = [get_badread_strand_id(i) for i in records]

    clustering = Clustering(strand_pool=reads, reference_length=oligo_len, strand_pool_ids=ids, n_reference_strands=num_of_oligos, distance_threshold=40)
    clustering.run_pipeline()

    checksum = CheckSum4(reference_length=oligo_len)
    decoded_strands = checksum.decode(candidates=clustering.candidates, n_reference_strands=num_of_oligos, clustered_seqs=clustering.clustered_seqs, n_guesses=5, guesses=True)
    ids = ["" for _ in range(len(decoded_strands))]
    print(path)
    # create_fasta_file(ids, decoded_strands, output_filepath='decoded_raptor.fasta')
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", metavar="file", type=str, help="the Badread fastq to cluster")
    parser.add_argument("oligo_len", type=int, help="Length of oligos")
    parser.add_argument("num_of_oligos", type=int, help="Number of oligos")
    opt = parser.parse_args()

    cluster(opt.path, opt.oligo_len, opt.num_of_oligos)
