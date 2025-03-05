
from tqdm import tqdm


def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def remove_adapter(subseq, adapter):
    """
    Given the strand, remove the largest direct match of the adapter.
    """
    
    kmer_length = len(subseq)

    while kmer_length >= 3:
        kmers = get_kmers(subseq, kmer_length)

        if any([i for i in kmers if i in adapter]):
            for ind, i in enumerate(kmers):
                if i in adapter:
                    return ind+kmer_length

        kmer_length -= 1

    return 0

def remove_adapters_from_strands(strands, original_strand_length, ids=None,
                                starting_adapter="AATGTACTTCGTTCAGTTACGTATTGCT"):
    """
    Removes the start adapters from the strand
    """

    # Add some adapter logging here
    cleaned_strands = []

    if ids:
        cleaned_ids = []

    for ind in tqdm(range(len(strands))):

        strand = strands[ind]
        overhang = len(strand) - original_strand_length

        if overhang > len(starting_adapter)*2:
            continue

        if overhang > 5:
            start_index = remove_adapter(strand[:len(starting_adapter) + 5], starting_adapter)
            strand = strand[start_index:]

        cleaned_strands.append(strand)

        if ids:
            cleaned_ids.append(ids[ind])
        
    return cleaned_strands, cleaned_ids