# viral-dynamics-research
Instructions:
1. Damped Gomperts - ... and Two Stage Decay - ... files are Monolix scripts that will run the respective model fitting scheme for the stated age group based on our data in Monolix.
  1a. data folder contains all data from the project.
  1b. model contains the model files used by Monolix.
2. loocv_monolix.R runs the leave one out cross validation in R and REQUIRES an active Monolix suite to run, as it runs Monolix at each iteration.
3. Trajectories and viralLoadPlots will reproduce plots from the manuscript in Rstudio.

"A simple model for viral decay dynamics and the distribution of
infected cell life spans in SHIV-infected infant rhesus macaques" abstract:

The viral decay following the initiation of antiretroviral therapy in a nonhuman primate
model of HIV infection and infant rhesus macaques infected with SHIV.C.CH505 does not
follow simple exponential dynamics. Several mathematical models have been proposed to
describe its more complex behavior, the most popular of which is two-phase exponential
decay. The underlying assumption in two-phase exponential decay is that there are two
classes of infected cells with different lifespans. However, excepting of CD4+ T cells, there is
not a consensus on all of the cell types that can become productively infected, and the fit of
the two-phase exponential decay to the SHIV.C.CH505 data was relatively poor. Therefore,
we propose a new model for viral decay, inspired by the Gompertz model where the decay
rate itself is a dynamic variable. We modify the Gompertz model to include a linear term
that modulates the decay rate. We show that this simple model performs as well as the two-
phase exponential decay model on HIV and SIV data sets, and outperforms it for the infant
rhesus macaque SHIV.C.CH505 infection data set. We also show that by using a stochastic
differential equation formulation, the modified Gompertz model can be interpreted as being
driven by a population of infected cells with a continuous distribution of cell lifespans, and
estimate this distribution for the SHIV.C.CH505-infected infant rhesus macaques. Thus, we
find that the dynamics of viral decay in this model of infant HIV infection and treatment
may be explained by a distribution of cell lifespans, rather than two distinct cell types.
