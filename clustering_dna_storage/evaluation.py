
from Levenshtein import editops, ratio, distance
from tqdm import tqdm
from utils import reverse_complement
import numpy as np
from typing import Dict


def get_recovery_rate(consensus_strand: str, original_strand: str) -> float:
    """Gets recovery rate based on two strands. Chooses the length of the original strand to evaluate identity"""

    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
                1 for i in range(min_length)
                if consensus_strand[i] == original_strand[i]]
                ) / len(original_strand)


def count_ids_errors(str1: str, str2: str):
    edit_operations = editops(str1, str2)
    
    insertions = sum(1 for op in edit_operations if op[0] == 'insert')
    deletions = sum(1 for op in edit_operations if op[0] == 'delete')
    substitutions = sum(1 for op in edit_operations if op[0] == 'replace')

    return {'Insertions': insertions, 'Deletions': deletions, 'Substitutions': substitutions}

def evaluate_candidates(original_strands: list[str], candidates: list[str],
                        metric='edit_distance_ratio') -> Dict[str, np.ndarray]:
    """
    Calculate and return reference recoveries, matching reference indices, and recovery rates per candidate

    Args:
        original_strands: list of reference strands to compare to
        candidates: list of strand guesses
        metric: edit_distance_ratio or identity
    Returns:f
        {
            "reference_recoveries": best recovery percentage of each reference strand
            "reference_strand_indices": matching reference indices for each candidate
            "recovery_rates": recovery rate for each candidate with its best matching reference
        } 
    """

    assert metric == 'edit_distance_ratio' or metric == 'identity'

    reference_strand_indices = np.zeros(len(candidates), dtype=int)
    recovery_rates = np.zeros(len(candidates))
    reference_recoveries = np.zeros(len(original_strands))

    for ind, candidate in tqdm(enumerate(candidates)):

        if candidate in original_strands:
            original_strand_index = original_strands.index(candidate)
            reference_strand_indices[ind] = original_strand_index
            recovery_rates[ind] = 1.0
            reference_recoveries[original_strand_index] = 1.0
            continue
        
        else:
            best_recovery = 0.0
            best_matching_index = 0
            for ind_, original_strand in enumerate(original_strands):
                
                if metric == 'edit_distance_ratio':
                    recovery = ratio(original_strand, candidate)
                elif metric == 'identity':
                    recovery = get_recovery_rate(candidate, original_strand)

                if recovery > best_recovery:                 
                    best_recovery, best_matching_index = recovery, ind_                                          

            reference_strand_indices[ind] = best_matching_index
            reference_recoveries[best_matching_index] = best_recovery
            recovery_rates[ind] = best_recovery

    return {
        'reference_recoveries': reference_recoveries,
        'reference_strand_indices': reference_strand_indices,
        'recovery_rates': recovery_rates
    }

    
