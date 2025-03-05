
from Levenshtein import distance
from tqdm import tqdm
from strand_reconstruction import multiple_alignment_muscle, majority_merge
from utils import reverse_complement
import random
import numpy as np
from typing import List, Dict


def get_edit_distance_matrix(strands):
    """
    Returns the edit distance matrix for the strands
    O(n^2)
    """
    n_strands = len(strands)
    edit_distance_matrix = np.zeros([n_strands, n_strands])
    for i in range(n_strands - 1):
        for j in range(i + 1, n_strands):
            edit_distance_matrix[i,j] = edit_distance_matrix[j, i] = distance(strands[i], strands[j])

    return edit_distance_matrix

def calculate_centroid(strands: list[str]):
    edit_distance_matrix = get_edit_distance_matrix(strands)

    distances = [sum(edit_distance_matrix[i, :]) for i in range(len(edit_distance_matrix))]
    return strands[distances.index(min(distances))]

def calculate_centroid(strands: list[str], edit_distance_matrix: np.ndarray):
    distances = [sum(edit_distance_matrix[i, :]) for i in range(len(strands))]
    return strands[distances.index(min(distances))]

def update_distance_matrix(
        added_strand: str, cluster_strands: list[str], distance_matrix: np.ndarray):
    """Adds the distances of the added strand to the cluster distance matrix"""

    # List of lists. keep going till next. if sum < main centroids, start going to next hmm

    new_strand_index = len(cluster_strands)
    for ind, cluster_strand in enumerate(cluster_strands):
        distance_matrix[ind, new_strand_index] = distance_matrix[new_strand_index, ind] = distance(added_strand, cluster_strand)

    return distance_matrix


def cluster_strands(strand_pool: List[str], distance_threshold: int = 40, use_centroids: bool = False,
                    sort_order: bool = True) -> Dict[str, List]:
    """
    Agglomeratively cluster strands using the Levenshtien distance.

    Args:
        strand_pool: List of the strands in the pool to be clustered
        distance_threshold: Threshold for edit distance to cluster similar strands together
        use_centroids: Calculate and update centroids while clustering
        sort_order: Sort from largest to smallest cluster every 100 iterations
    Returns:
        clusters: List of clusters with indices of the strand in each cluster
        reversed_markers: boolean list of whether a strand is reversed or not
        cluster_heads: The representative strand for each cluster. If using centroids, this becomes the centroid
    """

    if use_centroids:
        print("This feature is not implemented yet!")
        return

    clusters = []
    cluster_heads = []
    n_strands_pool = len(strand_pool)
    reversed_markers = np.zeros(n_strands_pool, dtype=bool)
    within_clusters = False

    print(f"Total strands {len(strand_pool)}")

    for ind, strand in enumerate(tqdm(strand_pool)):

        within_clusters = False
        reversed_flag = False
        rev_strand = reverse_complement(strand)
        
        for cluster_ind, cluster_head in enumerate(cluster_heads):
            
            if distance(cluster_head, strand) <= distance_threshold:
                within_clusters = True
                break

            elif distance(cluster_head, rev_strand) <= distance_threshold:
                within_clusters = True
                reversed_flag = True
                break

        if within_clusters:            
            clusters[cluster_ind].append(ind)
            
        if not within_clusters:
            clusters.append([ind])
            cluster_heads.append(strand)

        if reversed_flag:
            reversed_markers[ind] = reversed_flag

        if sort_order:
            if ind % 100 == 0 or ind == (n_strands_pool - 1):
                clusters, cluster_heads = zip(*sorted(
                    zip(clusters, cluster_heads), key=lambda x: len(x[0]), reverse=True))

                clusters = list(clusters)
                cluster_heads = list(cluster_heads)
                
    
    print(f"Number of clusters = {len(cluster_heads)}")

    return {
        "clusters": clusters,
        "reversed_markers": reversed_markers,
        "cluster_heads": cluster_heads
    }
           
def make_prediction(cluster, sample_size=5):
    cluster = random.sample(cluster, sample_size)
    return majority_merge(multiple_alignment_muscle(cluster))

