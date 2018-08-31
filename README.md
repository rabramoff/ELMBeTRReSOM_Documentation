README
--------------------------------------------------------------------------------
The soil decomposition model embedded in ELM-BeTR is called the Reaction-network-based model of Soil Organic Matter and microbes (ReSOM) model. 

1. The simplest version can be found [here](https://github.com/jinyun1tang/one_bug_model)

2. A version with somewhat more options and some scripts for simple global change experiments can be found [here](https://github.com/rabramoff/ReSOM_Interface)

3. The fortran version of the model integrated into the reactive transport model BeTR is [here]( https://github.com/rabramoff/sbetr/tree/rzacplsbetr) <- you will find the code for the model itself in src/Applications/soil-farm/SUMMS (an older name for the model). There are three folders in there - summsNlayer which contains scripts that couple the model to the reactive transport model, summs1layer which contains the main model scripts (the whole structure of the model has changed substantially from the matlab version to comply with the common format for moving CNP around in a large-scale model), and summspara which has the parameter values. See Model_Description.docx and Model_organization.pptx for more details.

4. The fortran version of the Earth System Model E3SM can be found [here](https://github.com/rabramoff/ACME/tree/rzacplsumms) (this calls repository #3 above as a submodule): 
