
from typing import List, Tuple
from math import floor
from utils import reverse_complement
from strand_reconstruction import make_prediction
from tqdm import tqdm


class CheckSum4():

    def __init__(self, reference_length: int):
        self.reference_length = reference_length
        self.base_maps = {"A": 0, "C": 1, "G": 2, "T": 3}
        self.number_maps = ["A", "C", "G", "T"]

    def get_quaternary_sum(self, seq: str) -> int:
        return sum([self.base_maps[base] + 1 for base in seq]) % 256

    def encode_to_quaternary(self, num: int) -> str:
        seq = ""
        ind = 3
        while not num == 0 or len(seq) < 4:
            freq = floor(num / (4 ** ind))
            seq = seq + self.number_maps[freq]
            num = num % (4 ** ind)
            ind -= 1
        return seq
    
    def decode_quateranary(self, checksum_seq: str) -> int:
        checksum = 0
        for ind, base in enumerate(checksum_seq):
            checksum += 4**(3 - ind) * self.base_maps[base]
        return checksum


    def encode(self, strands: List[str]) -> List[str]:
        """Add 4 nt seq to the end of each strand corresponding to the checksum"""
        encoded_strands = []
        for strand in strands:
            checksum = self.get_quaternary_sum(seq=strand)
            checksum_seq = self.encode_to_quaternary(checksum)
            encoded_strand = strand + checksum_seq
            encoded_strands.append(encoded_strand)
        
        return encoded_strands
    
    def make_checksum_guesses_from_cluster(
            self, cluster: List[str], n_guesses: int) -> Tuple[str, bool, bool]:
        
        for i in range(n_guesses):
            candidate = make_prediction(cluster, sample_size=10)
            rev_candidate = reverse_complement(candidate)

            if self.verify_checksum(candidate):
                return candidate, True, False
            if self.verify_checksum(rev_candidate):
                return rev_candidate, True, True        
        return candidate, False, False
    
    def verify_checksum(self, candidate: str):
        return self.decode_quateranary(
            candidate[-4:]) == self.get_quaternary_sum(
                candidate[:-4]) and len(candidate) == self.reference_length + 4

    def decode(
            self, candidates: List[str], n_reference_strands: int,
            clustered_seqs: List[List[str]], n_guesses: int,
            guesses: bool = False) ->  List[str]:
        
        decoded_strands = set()
        n_reversed = 0
        found_indices = []
        for ind, candidate in tqdm(enumerate(candidates), total=len(candidates)):
            
            rev_candidate = reverse_complement(candidate)
            if self.verify_checksum(candidate):
                decoded_strands.add(candidate[:-4])
                found_indices.append(ind)
                continue
            elif self.verify_checksum(rev_candidate):
                decoded_strands.add(rev_candidate[:-4])
                n_reversed += 1
                found_indices.append(ind)
                continue
        
        print(f"{len(found_indices)} direct checksum matches found")
        guess_indices = [i for i in range(len(candidates)) if i not in found_indices]

        if guesses:
            print("Making guesses")
            for ind in tqdm(guess_indices):
                candidate = candidates[ind]
                candidate, decoded, reversed = self.make_checksum_guesses_from_cluster(
                    cluster=clustered_seqs[ind], n_guesses=n_guesses)        
                if decoded:
                    decoded_strands.add(candidate[:-4])
                    if reversed:
                        n_reversed += 1

        decoded_strands = list(decoded_strands)
        print(f"{len(decoded_strands)  - len(found_indices)} extra found after guessing")
        print(f"{n_reversed / len(decoded_strands)} were reversed")
        print(f"{len(decoded_strands)} Valid checksum strands found")

        return decoded_strands
