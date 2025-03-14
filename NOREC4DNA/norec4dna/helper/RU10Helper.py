#!/usr/bin/python
# -*- coding: latin-1 -*-
import typing
import numpy as np
from functools import lru_cache
from math import ceil, floor, pow, sqrt, log

from norec4dna.distributions import RaptorDistribution

int63 = int(pow(2, 63) - 1)
int31 = int(pow(2, 31) - 1)


def choose_packet_numbers(number_of_chunks: int, code_block_index: int, dist: RaptorDistribution,
                          systematic: bool = False, max_l: int = None) -> typing.List[int]:
    """
    Parameters
    ----------
    number_of_chunks : int
        k
    code_block_index : int
        seed for random generator
    dist : RaptorDistribution
        Distribution used to generate random values and determine the degree distribution

    Returns
    -------
    indices : list
        Return the sorted list of indices of intermediate (k + s + h) chunks used.
        Size of list = degree d.
    """
    if systematic:
        d, a, b = systematic_ru10_triple_generator(number_of_chunks, code_block_index, dist)
    else:
        # always
        d, a, b = ru10_triple_generator(number_of_chunks, code_block_index, dist, max_l)
    
    # l = k + s + h
    if max_l is None:
        l, _, _ = intermediate_symbols(number_of_chunks, dist)
    else:
        l = max_l
    lprime: np.uint32 = np.uint32(dist.smallestPrimeGreaterOrEqual(l))

    # upper bound of degree is # of intemediate symbols
    if d > l:
        d = l
    l = np.uint32(l)
    # initialise indices to size d.
    indices: typing.List[int] = [0] * d
    # Adjust b to ensure it is within the range [0, l - 1] by repeatedly adding a and taking modulo lprime.
    while b >= l:
        b = (b + a) % lprime

    # Store b as first index
    indices[0] = b

    # for all d
    for idx in range(1, d):
        # add a to b and mod lprime
        b = (b + a) % lprime
        # ensure b in range [0, l - 1]
        while b >= l:
            b = (b + a) % lprime
        # store b
        indices[idx] = b
    # print(sorted(indices))
    return sorted(indices)


@lru_cache(maxsize=None)
def intermediate_symbols(k, dist) -> typing.Tuple[int, int, int]:
    """
    Parameters
    ----------
    k : int
        Number of chunks.
    dist : Distribution
        Class containing prime & central binomial utlity functions.

    Returns
    -------
    (k + s + h, s, h) : tuple
        k + s + h = the total number of intermediate symbols.
        s = the number of LDPC symbols.
        h = the number of Half symbols.
    """
    # x is the smallest positive integer such that X*(X-1) >= 2*K
    x = int(floor(sqrt(2 * np.float64(k))))
    if x < 1:
        x = 1

    while (x * (x - 1)) < (2 * k):
        x += 1

    # find smallest prime ≥ ceil(0.01k) + x
    s = int(ceil(0.01 * np.float64(k))) + x
    s = dist.smallestPrimeGreaterOrEqual(s)

    # initially set h to floor(log4(s+k))
    h = int(floor(log(np.float64(s) + np.float64(k)) / log(4)))
    # increase h until h = number of unique combinations of Half symbols ≥ k + s
    # dist.centerBinomial = central binomial coefficient (2h, h)
    while dist.centerBinomial(h) < k + s:
        h += 1
    return k + s + h, s, h


def ru10_triple_generator(k: int, x: int, dist: RaptorDistribution, max_l: typing.Optional[int] = None) \
        -> typing.Tuple[int, np.uint32, np.uint32]:
    """
    Parameters
    ----------
    k : int
        Number of source chunks
    x : int
        Seed value for random number generator
    dist : RaptorDistribution
        Raptor distribution, used for generating random numbers and determining degree
    max_l : int
        max number of intermediate symbols, if None l is calculated with intermediate_symbols()

    Returns
    -------
    d : int
        Degree - number of intermediate blocks to be combined.
    a : int
        Random int
    b : int
        Random int
    """
    # l = k + s + h, total number of intermediate symbols
    if max_l is None:
        l, _, _ = intermediate_symbols(k, dist)
    else:
        l = max_l
    lprime = dist.smallestPrimeGreaterOrEqual(l)

    # set the random generator seed, the same x will produce the same (d, a, b) output
    rng: np.random = np.random
    rng.seed(x)

    # v = a random 32-bit integer in the range [0, 1048575] (i.e., 2^20 − 1). 
    # This value is used to determine the degree d.
    v = np.uint32(r_int63(rng) % 1048576)
    a = np.uint32(1 + (r_int63(rng) % (lprime - 1)))
    b = np.uint32(r_int63(rng) % lprime)

    # calculate degree from random value v.
    d = dist.deg(v)
    return d, a, b


def systematic_ru10_triple_generator(k: int, x: int, dist: RaptorDistribution) -> typing.Tuple[int, int, int]:
    l, _, _ = intermediate_symbols(k, dist)
    lprime = dist.smallestPrimeGreaterOrEqual(l)
    q = 65521  # largest prime < 2 ^ 16
    jk = dist.systematicIndextable[k]

    a = int(53591 + jk * 997) % int(q)
    b = 10267 * (jk + 1) % q
    y = int(b + (x * a)) % q
    v: int = dist.raptor_rand(y, 0, 1048576)
    d: int = dist.deg(v)
    a: int = 1 + dist.raptor_rand(y, 1, int(lprime - 1))
    b: int = dist.raptor_rand(y, 2, int(lprime))
    return d, a, b


def r_int63(rng: np.random) -> int:
    return rng.randint(0, int31)


def from_true_false_list(tf_list: typing.List[bool]) -> typing.List[int]:
    return [i for i, x in enumerate(tf_list) if x]


if __name__ == '__main__':
    print(choose_packet_numbers(500, 123, RaptorDistribution.RaptorDistribution(500)), False, None)
