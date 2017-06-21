import numpy as np
from scipy.stats.mstats import gmean
from scipy.stats import rankdata


def generate_permutations(p, k, n):
    """
    ------------------------------------------------------------
    Parameters:
    ------------------------------------------------------------
    p : int
        number of permutations
    k : int
        number of replicates
    n : int
        length of ranked list
    ------------------------------------------------------------
    Returns:
    ------------------------------------------------------------
    np.array of shape [p, n] containing rank products
    """
    rank_products = np.empty([p, n], dtype=np.float32)
    # generate p permutations of k rank lists of length n
    array = np.array([range(1, n+1) for i in range(k)])
    # shuffle in place p times
    for i in range(p):
        for row in array:
            np.random.shuffle(row)
        # calculate rank product of each col
        rank_products[i, :] = gmean(array, axis=0)
    return rank_products


def count_smaller_than(observed, permutations):
    """
    count how many times the rank products of the genes in the permutations
    are smaller or equal to the observed rank product
    ------------------------------------------------------------
    Parameters:
    ------------------------------------------------------------
    observed: np.array
        observed rank products
    permutations: np.array
        permuted rank products
    ------------------------------------------------------------
    Returns:
    ------------------------------------------------------------
    np.array c
    """
    assert len(observed) == permutations.shape[1]
    return (permutations <= observed).sum(axis=0)


def get_average_expected_RP(c, p):
    """
    calculate the average expected value for the rank product
    """
    return c / p


def signif_rank_product(ranks_array, n_permutations=10000):
    # calculate rank products of array:
    rank_products = gmean(ranks_array, axis=1)
    # generate permutations of rank_products
    n_samples, n_replicates = ranks_array.shape
    permuted = generate_permutations(n_permutations,
                                     k=n_replicates,
                                     n=n_samples)
    c = count_smaller_than(observed=rank_products, permutations=permuted)
    av_expected = get_average_expected_RP(c, n_permutations)
    # calculate ranks
    ranks = rankdata(rank_products)
    return av_expected / ranks


