# Root and Treated Root causal Effects (R-TRE)

Treatments ideally mitigate pathogenesis, or the detrimental effects of the root causes of disease. However, existing definitions of treatment effect fail to account for pathogenic mechanism. We therefore introduce the Treated Root causal Effects (TRE) metric which measures the ability of a treatment to modify root causal effects. We leverage TREs to automatically identify treatment targets and cluster patients who respond similarly to treatment. The R-TRE algorithm learns a partially linear causal model to extract the root causal effects of each variable and then estimates TREs for target discovery and downstream subtyping. R-TRE does not require an invertible structural equation model.

# Installation

> library(devtools)

> install_github("ericstrobl/RTRE")

> library(RTRE)

# Run R-TRE

> DAG = generate_DAG(p,2) # generate a DAG and linear SEM with p=10 vertices

> dataT = sample_DAG_start(10000,DAG) # draw 10000 samples from the SEM + discretize some variables to ensure non-invertibility

> out = RTRE(dataT$data,DAG$graph,10) # recover the root causal effects and TREs; can also replace `DAG$graph` with a learned graph
