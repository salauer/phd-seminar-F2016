Title: Assessing Incremental Value of Biomarkers with Multi-phase Nested Case-Control Studies

We propose robust statistical procedures of constructing prediction models and evaluating the incremental values (IncV) of new markers for three-phase nested case control (NCC) studies. In Phase I, clinical markers, are collected at baseline on the full cohort. In Phase II, for each case, m1 controls are randomly selected from the risk set of the case matched on several matching variables. Biomarkers are available for only the Phase II sub-cohort. In Phase III, m2 controls are further randomly selected from the m1 controls chosen in the Phase II sub-cohort. Genetic markers are available for only the Phase III sub-cohort.

This .zip package includes two data files (data0.csv and NCC.dat.csv), an R file including all the R functions needed in the analysis, and an example R file (example.R) to run the analysis.

The R function execute.fun is to construct several prediction models: (1) clinical model with clinical markers alone, denoted by model "Z"; (2) biomarker model with clinical and biological markers, denoted by model "B"; (3) genetic model with clinical and genetic markers, denoted by model "G"; and (4) full model with all the available clinical, biological and genetic markers, denoted by "full" model. The output of this function includes

- $beta: provide, for each prediction model, point estimates of relative risk parameters and their naive and adjusted standard error estimates. 
- $acc: provide, for each prediction model, point estimates of four accuracy measures: area under ROC curve (AUC) and true positive rate (TPR), positive predictive values (PPV) and negative predictive values (NPV) at the cut-off value achieving certain level (u0) of false positive rate (FPR). In addition, their naive and adjusted standard error estimates are also included.
- $cv: if choosing yes.cv=TRUE, provide, for each prediction model, cross-validated estimates of the four accuracy measures described above as well as the corresponding ROC curves.
- $IncV: if choosing yes.IncV=TRUE, provide point estimates as well as naive and adjusted standard error estimates for IncV parameters of new markers by comparing different pairs of prediction models.

The arguments needed in the execute.fun include

- two data files
(1) data0: this file lists the following information
	- Column 1: ID of participants
	- Column 2: Censored event time, i.e., the minimum between the event time and censoring time
	- Column 3: Censoring indicator variable; = 1, if the event time is censored; = 0, otherwise
	- The remaining columns: the markers including clinical markers, biomarkers, and genetic markers.
(2) NCC.dat: this file lists the following information:
	- Column 1: = 1, if the subject is selected in the sub-cohort as a case; = 0, otherwise
	- Column 2: = 1, if the subject is selected in the Phase II nested case control sub-cohort as a control; = 0, otherwise
	- Column 3: = 1, if the subject is selected in the Phase III nested case control sub-cohort as a control; = 0, otherwise
	- The remaining columns: the matching variables which might be needed in the NCC sampling
- cov.list: list of marker names (“nm”) used for different prediction models and the cohort or sub-cohort (“ncc”) where all the involved markers are available. For example, in Model "Z" where the clinical markers are available on Phase I full cohort, “ncc=S0”; in Model "B" where the clinical and biological markers are both available on Phase II sub-cohort, “ncc=S1”; in Model "G" where the clinical and genetic markers are both available on Phase III sub-cohort, “ncc=S2”; in the full model where all the markers are available on Phase III sub-cohort, “ncc=S2”.
- model: either time-specific GLM model (which you can specify the link function) or Cox's PH model
- u0: the level of FPR to achieve for obtaining TPR, PPV and NPV
- t0: the prediction time
- execute.rtn: =“ ALL" (default), if both point estimates and standard error estimates are required; =“ EST", if only point estimates are required.
- yes.DepCen: = FALSE (default), if the distribution of censoring time is estimated by the Kaplan-Meier estimates assuming independent censoring; = TRUE, if allowing for marker-dependent censoring
- cov.DepCen: if yes.DepCen=FALSE, cov.DepCen = list(Cov=NULL,Stage=NULL) (default); if yes.DepCen=TRUE, list of the names of markers on which the distribution of censoring time depends as well as which phase these markers are available.
- control:
	yes.cv: = TRUE, if cross-validated estimates are required; = FALSE (default), otherwise
	yes.IncV: = TRUE (default), if IncVs are required; = FALSE (default), otherwise








