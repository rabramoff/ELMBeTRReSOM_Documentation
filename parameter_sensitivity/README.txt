These Matlab scripts were written by Qing Zhu to run ensembles of parameter values for parameter sensitivity analysis. I have not modified these scripts to run parameter ensembles on ELM-BeTR-ReSOM yet. 

Step 1 identifies major representative eco-regions by sampling climatology (not needed for single site testing)
Step 2 generates parameter ensembles from user-defined ranges in params_range.txt
Step 3 generates a grid where each grid cell is a different model run - this allows for many model runs to be executed in parallel

Usually parameters for ELM are stored in a dedicated parameter file like clm_param.nc, but for this analysis the parameters are integrated into the surface dataset so that each grid cell's surface data contains all the parameter information needed in the run. Because this is a structural change to the model, there is a particular branch that is currently set up to run these ensembles. It is qzhu-lbl/lnd/bgc-param-est-v7-qzhu-before-alloc-bf (https://github.com/E3SM-Project/E3SM/tree/qzhu-lbl/lnd/bgc-param-est-v7-qzhu-before-alloc-bf). Therefore, some code development or carefully-done merging has to happen in my branch in order to support this analysis.

The parameters that I would test first are maximum microbial growth rate (gmax_mic), mineral surface area (minsite), maximum enzyme production rate (pmax_enz), the activation energy parameters (ea_*) (see SummsParaType.F90) and the binding affinity for DOM/DOC (Kd; see TracerParamSetMod.F90).