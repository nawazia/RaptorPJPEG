from utils import (
    create_fasta_file, read_synthesized_strands_from_file,
    get_fastq_records, get_badread_strand_id,
    run_badread_simulation
    )
from checksum import CheckSum4
from clustering import Clustering
import gzip
import shutil

def main():
    checksum = CheckSum4(reference_length=204)
    original_strands, original_strand_ids = read_synthesized_strands_from_file(file_path="/Users/i/Downloads/DNA_STorage/code/libjpeg/cat_FFDX.jpg_RU10.fasta")

    encoded_strands = checksum.encode(original_strands)
    create_fasta_file(ids=original_strand_ids, strands=encoded_strands, output_filepath="/Users/i/Downloads/DNA_STorage/code/libjpeg/cat_FFDX_checksum.fasta")

    with gzip.open("/Users/i/Downloads/DNA_STorage/code/libjpeg/reads.fastq.gz", 'rb') as f_in:
        with open("/Users/i/Downloads/DNA_STorage/code/libjpeg/reads.fastq", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    records = get_fastq_records(fastq_filepath="/Users/i/Downloads/DNA_STorage/code/libjpeg/reads.fastq")

    reads = [str(i.seq) for i in records]
    ids = [get_badread_strand_id(i) for i in records]

    clustering = Clustering(strand_pool=reads, reference_length=208, original_strands=encoded_strands, strand_pool_ids=ids, distance_threshold=40)

    clustering.run_pipeline()

    # We want to recover the original strands after checksum
    decoded_strands = checksum.decode(candidates=clustering.candidates, n_reference_strands=len(original_strands), clustered_seqs=clustering.clustered_seqs, n_guesses=5, guesses=True)
    ids = ["" for i in range(len(decoded_strands))]
    create_fasta_file(ids, decoded_strands, output_filepath="/Users/i/Downloads/DNA_STorage/code/libjpeg/cat_FFDX_BR.fasta")
    return 0

if __name__ == "__main__":
    main()
