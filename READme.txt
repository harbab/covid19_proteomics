CODE AND DATASETS FOR PEER REVIEW

The information here is organised in three segments:
1. Code 
- contains the code used for making the figures

2. Data 
- contains the datasets needed to reproduce the findings

3. Results 
- contains the figure outputs as generated in R

CODE
The code is organised according to the figures in the manuscript.
Make sure that you set the working directory to correspond to the location of the code_submission folder before using it.

The code should be reproducible if you run each R file separately. However, if you get errors, try running the files in this order:
figure_1_2_4_covid_hirief_serum_analyses.Rmd
figure_3_sars_cov_2_infection_experiment.Rmd
figure_5_systematic_review_summary.R
figure_6_meta_analysis_summary.Rmd

The folder covimapp contains all the functions used to calculate the estimates for covimapp, excluding the app.R file for building the user interface.

DATA
In the meta-analysis per dataset, we only included an example of our study because we do not own the datasets shared with us by others.
However, we provided datasets with summary values that you can load in R and reproduce the code.

We are legally not allowed to share the datasets with clinical information because they contain sensitive data that in theory could allow tracing of the patients.
However, we kept the code used for cleaning the clinical variables and analyses that include such variables, so that the reviewers can have a look at it. If you run those chunks that use clinical variables, you will get an error because the code will not be able to use the clinical information.
