import os
import sys
import json
import time
import shutil
import logging
import builtins   
import argparse
import datetime
from pkg_resources import resource_stream
import numpy as np

from clumppling.funcs import *
import warnings



def clumppling_run(args):

    input_path = args.input_path
    output_path = args.output_path
    input_format = args.input_format
    
    if os.path.exists(input_path+".zip"):
        shutil.unpack_archive(input_path+".zip",input_path)
    if not os.path.exists(input_path):
        sys.exit("ERROR: Input file {} doesn't exist.".format(input_path))
    if input_path==output_path:
        sys.exit("ERROR: The input and output file paths are the same: {}. Please use different paths.".format(input_path))
    if not input_format in ["structure","fastStructure","admixture","generalQ"]:
        sys.exit("ERROR: Input data format is not supported. \nPlease specify input_format as one of the following: structure, admixture, fastStructure, and generalQ.")  

    overwrite = False
    if os.path.exists(output_path):
        overwrite = True
        shutil.rmtree(output_path)

    if os.path.exists(output_path+".zip"):
        os.remove(output_path+".zip")
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if args.vis is not None:
        visualization = bool(args.vis)
    else:
        visualization = True

    # load default parameters and process
    parameters = vars(args)
    parameters = process_parameters(parameters)
    plot_params = ['plot_modes','plot_modes_withinK','plot_major_modes','plot_all_modes']
    if visualization==False:
        for param in plot_params:
            parameters[param] = False
    disp = display_parameters(input_path,input_format,output_path,parameters)
    
    #%% Set-up 
    # log outputs
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=os.path.join(output_path,'output.log'), level=logging.INFO, format='')
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger('matplotlib.pyplot').disabled = True
    logging.getLogger('matplotlib.pyplot').disabled = True

    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#    logging.info("{} Program starts.".format(current_time))

    if overwrite:
        logging.info('\033[93m'+"Overwriting existing output directory {}.".format(output_path)+'\033[0m') 
    
#    logging.info("==================================")
#    logging.info(disp)
    logging.info("======= Running Clumppling =======")

    tot_tic = time.time()   
    tic = time.time()
    logging.info(">>> Processing input data files and checking arguments")
    
    
    #%% Loading membership data    
    # load input data
    Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs = load_inputs(data_path=input_path, 
                                                                   output_path=output_path, 
                                                                   input_format=input_format)
    
    distruct_cmap = np.loadtxt(resource_stream('clumppling', 'files/default_colors.txt'),dtype=str,comments=None)

    # set colormap
    if parameters['custom_cmap']!='':
        custom_cmap = [s.strip() for s in parameters['custom_cmap'].split(",")]
        while len(custom_cmap) < K_max:
            logging.info("The provided colormap does not have enough colors for all clusters. Colors are recycled.")
            custom_cmap.extend(custom_cmap)
        cmap = cm.colors.ListedColormap(custom_cmap)
    else:
        cmap = cm.colors.ListedColormap(distruct_cmap[:K_max])

    if visualization:
        # visualize results
        fig_path = os.path.join(output_path,"visualization")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        plot_colorbar(cmap,K_max,fig_path)


    toc = time.time()
    logging.info("Time: %.3fs",toc-tic)
    
    
    #%% Alignment within-K and mode detection
    
    tic = time.time()
    logging.info(">>> Aligning replicates within K and detecting modes")

    # align within-K
    alignment_withinK, cost_withinK = align_withinK(output_path,Q_list,Q_files,K_range,K2IDs)

    # detect mode
    modes_allK, cost_matrices, msg = detect_modes(cost_withinK,Q_files,K_range,K2IDs,default_cd=parameters['cd_default'],cd_param=parameters['cd_param'])

    # extract representative/consensus replicates
    mode_labels, rep_modes, repQ_modes, avgQ_modes, alignment_to_modes, stats, costs = extract_modes(Q_list,Q_files,modes_allK,alignment_withinK,cost_matrices,K_range,K2IDs, output_path)
    
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)
    
    #%% Alignment across-K
        
    # use representative membership as the consensus
    logging.info(">>>  Aligning modes across K")
    tic = time.time()
    # align across-K
    acrossK_path = os.path.join(output_path,"alignment_acrossK")
    if parameters['use_rep']:
        alignment_acrossK_rep, cost_acrossK_rep, best_acrossK_rep = align_ILP_modes_acrossK(repQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="rep",merge=parameters['merge_cls'])
    else:
        alignment_acrossK_avg, cost_acrossK_avg, best_acrossK_avg = align_ILP_modes_acrossK(avgQ_modes,mode_labels,K_range,acrossK_path,cons_suffix="avg",merge=parameters['merge_cls'])
    toc = time.time()
    logging.info("Time: %.3fs", toc-tic)

    # set up summaries for visualisation
    if parameters['use_rep']:
        plot_structure_on_multipartite(K_range,mode_labels,stats,repQ_modes,alignment_acrossK_rep,cost_acrossK_rep,best_acrossK_rep,"rep",output_path,False,cmap)
    else:
        plot_structure_on_multipartite(K_range,mode_labels,stats,avgQ_modes,alignment_acrossK_avg,cost_acrossK_avg,best_acrossK_avg,"avg",output_path,False,cmap)

    # zip all files
#    logging.info(">>>  Zipping files")
#    shutil.make_archive(output_path, "zip", output_path)

    tot_toc = time.time()
    logging.info("======== Total Time: %.3fs ========", tot_toc-tot_tic)

    if parameters['use_rep']:
      return args, cmap, Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs, alignment_withinK, cost_withinK, modes_allK, cost_matrices, msg, mode_labels, rep_modes, repQ_modes, avgQ_modes, alignment_to_modes, stats, costs, alignment_acrossK_rep, cost_acrossK_rep, best_acrossK_rep
    else:
      return args, cmap, Q_list, K_list, Q_files, R, N, K_range, K_max, K2IDs, alignment_withinK, cost_withinK, modes_allK, cost_matrices, msg, mode_labels, rep_modes, repQ_modes, avgQ_modes, alignment_to_modes, stats, costs, alignment_acrossK_avg, cost_acrossK_avg, best_acrossK_avg
      
