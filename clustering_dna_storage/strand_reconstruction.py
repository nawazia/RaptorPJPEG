import random
import operator
import subprocess
import os
from Bio import AlignIO, Align
from utils import reverse_complement
from Levenshtein import ratio
from typing import List
from tqdm import tqdm


def align(seqA: str, seqB: str, identity: bool = True) -> any:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    
    return aligner.align(seqA, seqB)[0]

def multiple_alignment_muscle(
        cluster: List[str], out: bool = False, running_on_hpc: bool = False) -> List[any]:
    
    # write cluster to file
    file = open("clm.fasta","w") 
    
    for i,c in enumerate(cluster):
        file.write(">S%d\n" % i)
        file.write(c)
        file.write("\n")
    file.close()
    
    if running_on_hpc:
        # Will need to install muscle on new device and change the path
        muscle_exe = os.path.join(os.environ['HOME'], 'muscle.exe')
    else:
        # Change the muscle path here
        muscle_exe = "/Users/i/Downloads/muscle-osx-arm64.v5.3"
        if not muscle_exe:
            print("Please add the muscle path!")
            exit()

    output_alignment = "clmout.fasta"

    try:
        output_message = subprocess.run(
            args=[
                f"{muscle_exe}", "-align", "clm.fasta", "-output", "clmout.fasta"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            timeout=120
        )
    except Exception as e:
        print(e)
        return None

    msa = AlignIO.read(output_alignment, "fasta")

    if out:
        print(msa)
    
    alignedcluster = []
    
    for i in msa:
        alignedcluster += [i.seq]
    
    return alignedcluster


def majority_merge(reads: List[str], weight: float = 0.4) -> str:
    
    assert [len(i) == len(reads[0]) for i in reads]

    res = ""
    for i in range(len(reads[0])):
        counts = {'A':0,'C':0,'G':0,'T':0,'-':0,'N':0}
        for j in range(len(reads)):
            if i >= len(reads[j]):
                continue
            counts[reads[j][i]] +=1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res

def make_prediction(cluster: List[str], sample_size: int=5) -> str:
    #assert sample_size <= len(cluster)
    cluster = random.sample(cluster, min(sample_size, len(cluster)))
    return majority_merge(multiple_alignment_muscle(cluster))

def get_clustered_seqs(clusters: List[List[int]], reversed_markers: List[bool], strand_pool: List[str]) -> List[List[str]]:
    assert len(reversed_markers) == len(strand_pool)
    return [
        [reverse_complement(strand_pool[index]) if reversed_markers[index] else strand_pool[index] for index in cluster] for cluster in clusters]

def get_candidates(clustered_seqs: List[List[str]], n_samples: int) -> List[str]:
    candidates = []
    for cluster_seqs in tqdm(clustered_seqs):
        candidates.append(make_prediction(cluster_seqs, n_samples))
    return candidates

def get_candidate_orientation(original_strands: List[str], candidates: List[str]) -> List[bool]:
    """Loop through all the candidates, find their best matching reference and reverse the orientation if needbe"""

    reversed_markers = []
    for candidate in candidates:

        rev = reverse_complement(candidate)

        if candidate in original_strands:
            reversed_markers.append(False)
            continue

        if rev in original_strands:
            reversed_markers.append(True)
            continue

        flag = False
        best_ratio = 0.0
        for strand in original_strands:

            rec_straight = ratio(candidate, strand)
            rec_rev = ratio(rev, strand)

            if rec_straight > best_ratio:
                best_ratio = rec_straight
                flag = False

            if rec_rev > best_ratio:
                best_ratio = rec_rev
                flag = True

        reversed_markers.append(flag)

    return reversed_markers      
