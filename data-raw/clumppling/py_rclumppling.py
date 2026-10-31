# this is derived from https://github.com/PopGenClustering/Clumppling/blob/main/clumppling/__main__.py
# @author of the original file Xiran Liu
# modified by Andrea Manica for use with tidygenclust
#
# Previous versions of tidygenclust used this module to run clumppling; it is
# now only used to compute the reference results for the tests of the native
# port (see data-raw/clumppling/README.md).
#
# Rather than running the clumppling command line tool (which reads all its
# inputs from, and writes all its outputs to, files), we call the clumppling
# functions directly, so that Q matrices can be passed from R as numpy arrays,
# and the results returned to R as python objects. Files are only used when
# the input is a directory of Q files, in which case we rely on clumppling's
# own parser to read them.

import os
import numpy as np


def _prepare_q_list(Q_list, tolerance=1e-6):
	"""Validate a list of Q matrices passed from R, and give them the same names
	and ids that clumppling's process_files() would give to the processed files.

	Returns Q_names, Q_mats (both ordered by increasing K), K_range and K2IDs.
	"""
	from clumppling.parseInput.utils import validate_membership_matrix, rescale_membership_matrix

	grouped = {}
	expected_nrows = None
	for i, Q in enumerate(Q_list):
		# copy, as rescale_membership_matrix modifies the matrix in place
		Q = np.array(Q, dtype=float, ndmin=2)
		validate_membership_matrix(Q, tolerance)
		Q = rescale_membership_matrix(Q)
		n_rows, n_cols = Q.shape
		if n_cols == 1:
			print(f"Ignored Q matrix {i+1} with K=1.")
			continue
		if expected_nrows is None:
			expected_nrows = n_rows
		elif n_rows != expected_nrows:
			raise ValueError(f"Q matrix {i+1} has {n_rows} individuals, expected {expected_nrows}.")
		grouped.setdefault(n_cols, []).append(Q)

	Q_names = []
	Q_mats = []
	K2IDs = {}
	global_id = 1
	for k in sorted(grouped.keys()):
		K2IDs[k] = []
		for i, Q in enumerate(grouped[k]):
			Q_names.append(f"{global_id}_K{k}R{i+1}")
			Q_mats.append(Q)
			K2IDs[k].append(global_id)
			global_id += 1
	return Q_names, Q_mats, sorted(grouped.keys()), K2IDs


def _read_q_files(input_dir, processed_dir, fmt, extension, remove_missing):
	"""Parse a directory of Q files with clumppling's own parser (which writes
	the processed matrices to processed_dir), and read them back in."""
	from clumppling.parseInput import process_files, extract_meta_input

	os.makedirs(processed_dir, exist_ok=True)
	process_files(input_dir=input_dir, output_dir=processed_dir,
				  fmt=fmt, extension=extension,
				  skip_missing=remove_missing,
				  delimiter=" ", skip_rows=0, label_cols=[0, 1, 3],
				  mat_start_col=5, reorder_indices=None)
	Q_names, K_range, K2IDs = extract_meta_input(processed_dir)
	Q_mats = [np.loadtxt(os.path.join(processed_dir, f"{name}.Q"), ndmin=2)
			  for name in Q_names]
	return Q_names, Q_mats, K_range, K2IDs


def clumppling_run(Q_list=None, input_dir=None, processed_dir=None,
				   fmt="admixture", extension=".Q", remove_missing=True,
				   cd_method="louvain", cd_res=1.0, test_comm=False,
				   comm_min=1e-6, comm_max=1e-2, merge=True, use_rep=True,
				   use_best_pair=True, random_seed=None):
	"""Run the clumppling pipeline (align within K, detect modes, align across
	K) either on a list of Q matrices (Q_list) or on a directory of Q files
	(input_dir, which is parsed into processed_dir)."""
	from clumppling.core import (align_within_k, detect_communities,
		extract_modes_and_stats, align_across_k, reorderQ_across_k)
	from clumppling.alignWithinK import aligned_res_to_dicts
	from clumppling.detectMode import construct_cost_mat, community_labels_to_modes
	from clumppling.utils import get_modes_all_K, unnest_list
	import pandas as pd

	if random_seed is not None:
		np.random.seed(random_seed)

	if Q_list is not None:
		Q_names, Q_mats, K_range, K2IDs = _prepare_q_list(Q_list)
	else:
		Q_names, Q_mats, K_range, K2IDs = _read_q_files(input_dir, processed_dir,
			fmt, extension, remove_missing)
	if len(K_range) < 2:
		raise ValueError(f"At least two K values are required for alignment across K: K={K_range}")

	# align within K, and detect and extract modes
	cd_res_list = []
	for K in K_range:
		prefixes = tuple(f"{i}_" for i in K2IDs[K])
		idx_K = [i for i, name in enumerate(Q_names) if name.startswith(prefixes)]
		Q_names_K = [Q_names[i] for i in idx_K]
		Q_list_K = [Q_mats[i] for i in idx_K]

		aligned_res = align_within_k(Q_list_K)
		alignment_withinK, cost_withinK = aligned_res_to_dicts(aligned_res, Q_names_K)
		cost_mat = construct_cost_mat(cost_withinK, Q_names_K)
		communities = detect_communities(cost_mat, test_comm=test_comm,
			method=cd_method, res=cd_res,
			min_threshold=comm_min, max_threshold=comm_max)
		modes = community_labels_to_modes(communities)
		cd_res_list.append(extract_modes_and_stats(modes, cost_mat, Q_names_K,
			Q_list_K, alignment_withinK, label_prefix=f"K{K}M"))

	# align modes across K
	mode_names_list, Q_rep_modes_list, Q_avg_modes_list = get_modes_all_K(K_range, cd_res_list)
	Q_modes_list = Q_rep_modes_list if use_rep else Q_avg_modes_list
	alignment_acrossK, cost_acrossK, best_acrossK, major_acrossK = align_across_k(
		K_range, Q_modes_list, mode_names_list, merge=merge)

	# reorder the Q matrices of the modes according to the alignment
	if use_best_pair:
		anchor_pairs = best_acrossK["Best Pair"].tolist()
	else:
		anchor_pairs = major_acrossK["Major Pair"].tolist()
	aligned_Qs_allK, _ = reorderQ_across_k(K_range, Q_modes_list,
		mode_names_list, alignment_acrossK, anchor_pairs)

	mode_names = unnest_list(mode_names_list)
	return {
		"N": int(Q_mats[0].shape[0]),
		"K_range": [int(K) for K in K_range],
		"mode_names": mode_names,
		"aligned_modes": [aligned_Qs_allK[name] for name in mode_names],
		"cost_acrossK": {pair: float(cost) for pair, cost in cost_acrossK.items()},
		"mode_stats": pd.concat([res["mode_stats"] for res in cd_res_list], ignore_index=True),
		"mode_alignment": pd.concat([res["align_info"] for res in cd_res_list], ignore_index=True),
	}
