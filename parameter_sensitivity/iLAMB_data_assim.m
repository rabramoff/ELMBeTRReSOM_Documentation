%% workflow
% 0. use the right branch for sparse grid tuning

% =============== this matlab code will do the following =================
% 1. extract representative grids to tune
% 2. generate parameters ensemble
% 3. generate sparse grid domain/surface dataset 
%    and aggregate multiple ensemble into one big-ensemble
% ========================================================================

% 4. run simulation on edison

% 6. post-process 1: split big-ensemble into multiple ensembles, run ilamb for multiple ensemble
%Step 1 Edison: ilamb.global_sparse_grid_post.sh, concatenate monthly outputs
%Step 2 Edison command line: 
%./exensvar_QZ.bash /global/cscratch1/sd/qzhu/acme_scratch/edison/bgccomparison/ xx1 1
%pick benchmark variables for ilamb
%Generate sparse grid DATA
%Step 3 MAC:
%cd /Users/QZhu/Desktop/ilamb_iext/bscripts_QZ
% ./run_glbflow_QZ.bash (for cru and gswp seperately)
%generate sparse grid DATA, split output into ensembles and pfts, and run ilamb, aggregate results
% 7. pick best ensembles (go back to step 3, iterate)

addpath('../');
%% step 1 extract representative grids to tune
nsites = step1_representative_sites(1); % extract n * hundred representative grids

%% step 2 generate parameters ensemble
step2_generate_parameter_ensemble('params_range.txt',12); % sample 2^8 parameters

%% step 3 gerate sparse grid domain/surface dataset and aggregate multiple ensemble into one big-ensemble
step3_generate_surface_domain('global_sparse_grid.txt','ALM_param.txt', 5); % tune 5 parameter in total, generate two domains and surfdata
