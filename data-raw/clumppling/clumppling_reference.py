# Reference results from the clumppling python package (v2.0.10), used to
# test the native R/C++ port in tidygenclust. This module is imported with
# reticulate by regenerate_python_fixtures.R, in the conda environment created
# by setup_clumppling_env.sh.
#
# Indices returned to R are 1-based.

import numpy as np


def _relabel(communities):
    """Relabel communities by decreasing size, as detect_communities() does."""
    from collections import Counter
    counts = Counter(communities)
    order = sorted(counts.keys(), key=lambda x: (-counts[x], communities.index(x)))
    mapping = {old: new for new, old in enumerate(order)}
    return [int(mapping[c]) for c in communities]


def ilp_cases(cases):
    """Optimal alignments (align_ILP, align_ILP_diff1) and alignQ_wrtP without
    merging, for a list of dicts with matrices P and Q."""
    from clumppling.core import align_ILP, align_ILP_diff1, alignQ_wrtP
    out = []
    for case in cases:
        P = np.array(case["P"], dtype=float, ndmin=2)
        Q = np.array(case["Q"], dtype=float, ndmin=2)
        obj, idx = align_ILP(P, Q)
        idx = [int(i) for i in idx]
        res = {"obj": float(obj), "idx": [i + 1 for i in idx]}
        aligned_q, _ = alignQ_wrtP(P, Q, idx, merge=True)
        res["merge_q"] = aligned_q
        if Q.shape[1] - P.shape[1] == 1:
            obj1, idx1 = align_ILP_diff1(P, Q)
            res["obj_diff1"] = float(obj1)
            res["idx_diff1"] = [int(i) + 1 for i in idx1]
        if Q.shape[1] > P.shape[1]:
            aligned_q, pattern = alignQ_wrtP(P, Q, idx, merge=False)
            res["nomerge_q"] = aligned_q
            res["nomerge_pattern"] = [int(i) + 1 for i in pattern]
        out.append(res)
    return out


def graph_cases(adjs, seeds):
    """Community detection on adjacency matrices (with ones on the diagonal):
    louvain (one result per seed), markov clustering, the 'custom' example
    method (markov clustering with pruning threshold 0.2, here run with
    cdlib's copy of markov_clustering), and the test for community
    structure."""
    from clumppling.core import cd_default
    from clumppling.detectMode import test_comm_struc
    from cdlib.algorithms.internal import markov_clustering as mc
    out = []
    for A in adjs:
        A = np.array(A, dtype=float)
        res = {"louvain": []}
        for seed in seeds:
            np.random.seed(seed)
            res["louvain"].append(_relabel(cd_default(A.copy(), method="louvain")))
        res["markov_clustering"] = _relabel(cd_default(A.copy(), method="markov_clustering"))
        result = mc.run_mcl(A.copy(), pruning_threshold=0.2)
        communities = [-1] * A.shape[0]
        for cluster_idx, nodes in enumerate(mc.get_clusters(result)):
            for node in nodes:
                communities[node] = cluster_idx
        res["custom"] = _relabel(communities)
        res["test_comm"] = bool(test_comm_struc(A.copy(), alpha=0.01))
        out.append(res)
    return out


def within_k_cases(Q_list, prepare_q_list, methods, seed):
    """Alignment within K and mode detection for each K, for a list of Q
    matrices prepared as in tidygenclust (prepare_q_list is
    py_rclumppling._prepare_q_list)."""
    from clumppling.core import align_within_k, detect_communities
    from clumppling.alignWithinK import aligned_res_to_dicts
    from clumppling.detectMode import construct_cost_mat, cost_to_adj, test_comm_struc
    Q_names, Q_mats, K_range, K2IDs = prepare_q_list(Q_list)
    out = {}
    for K in K_range:
        prefixes = tuple(f"{i}_" for i in K2IDs[K])
        idx_K = [i for i, name in enumerate(Q_names) if name.startswith(prefixes)]
        names_K = [Q_names[i] for i in idx_K]
        Q_K = [Q_mats[i] for i in idx_K]
        aligned_res = align_within_k(Q_K)
        alignment, cost = aligned_res_to_dicts(aligned_res, names_K)
        cost_mat = construct_cost_mat(cost, names_K)
        n = len(names_K)
        pattern = np.zeros((K, n, n), dtype=int)
        for i in range(n):
            pattern[:, i, i] = np.arange(K) + 1
            for j in range(n):
                if i != j:
                    pattern[:, i, j] = np.array(alignment[(names_K[i], names_K[j])]) + 1
        communities = {}
        for method in methods:
            for test_comm in (False, True):
                np.random.seed(seed)
                try:
                    comm = detect_communities(cost_mat, test_comm=test_comm,
                                              method=method, res=1.0,
                                              min_threshold=1e-6,
                                              max_threshold=1e-2)
                    comm = [int(c) for c in comm]
                except Exception:
                    comm = None
                communities[f"{method}_{test_comm}"] = comm
        out[f"K{K}"] = {
            "q_names": names_K,
            "cost_mat": cost_mat,
            "pattern": pattern,
            "test_comm": bool(test_comm_struc(cost_to_adj(cost_mat, norm=True), alpha=0.01)),
            "communities": communities,
        }
    return out


def detect_communities_cases(cost_mats):
    """detect_communities() on small cost matrices, to check the special
    cases (one or two replicates, costs below or above the thresholds, all
    costs equal)."""
    from clumppling.core import detect_communities
    out = []
    for C in cost_mats:
        np.random.seed(1)
        C = np.array(C, dtype=float, ndmin=2)
        out.append([int(c) for c in detect_communities(C, test_comm=False)])
    return out
