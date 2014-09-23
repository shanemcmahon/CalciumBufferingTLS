# This script executes a total least squares fit of a calcium buffering model using a stochastic gradient descent algorithm. It is divided into three sections. The first section contains user defined parameters, and must be specified by the user. The second section defines functions used for the fitting process, and should not be modified. The third section executes the fitting procedure, and should only be modified with caution. After modifying section 1, run this source file. A series of prompts will ask you to choose the files containing the fluorescence measurements, and standard errors for each measurement, and estimates for the change in total calcium concentration, and their standard errors.

#############################################################################################################################
#############################################################################################################################

# Begin Section 1: user defined parameters

#############################################################################################################################
#############################################################################################################################

# The fitting algorithm estimates values for the unknown fluorescence measurement errors at each point by maximizing the likelihood function evaluated at f.initial.i + f.initial.error.i and f.final.i + f.final.error.i for each i of n total measurements. We must place some restrictions on the algorithm so that the resulting test points are physically reasonable. That is, the normalized fluorescence must be positive and must be less than one. Here we specifiy the values F.MIN and F.MAX that limit the permited values of the error estimates so that F.MIN < f.initial.i + f.initial.error < F.MAX, with the corresponding condition holding for the fluorescence measurements after the stimulation.  F.MIN should be chosen to be greater than 1/rf, where rf is the dynamic range. By examining the equation for converting normalized fluorescence to free calcium i.e. ca = kd*(f - 1/rf)/(1-f), one can see that for f < 1/rf the calculated calcium concentrations are negative. Similarly for normalized fluorescence exactly equal to one, the implied calcium concentration is infinite, and negative for values greater than one. 
F.MIN <- .01
F.MAX <- .999

# max.iterations specifies the maximum number of stochastic gradient descent iterations to perform. In the current verstion, reaching max.iterations is the only method for terminating the fit, and exactly max.iterations will always be performed.
max.iterations <- 10000

# The full model contains nine total parameters:
# the dynamic range, rf
# the kd of the indicator dye, kd.dye
# the total concentration of the indicator dye, bt.dye
# the kd of endogenous buffers to be modeled, kd.1 and kd.2
# the total concentration of the corresponding buffers, bt.1 and bt.2
# a nonsaturable calcium buffer or a buffer with sufficiently low affinity that its buffering capacity does not change appreciably over the range of free calcium values investigated, kappa.nonsaturable
# the fraction of the intracellular volume that is accessible to calcium, accesible.volume
# all of the parameters are specified in a single vector, beta
# the user must specify starting parameters and upper and lower limits
# for a particular run, the user must specify whicy parameters are to be considered as fixed at there starting values, and which parameters the algorithm should modify to perform the fit. 

# parameter.is.fixed is a logical vector specifying which parameters are to be considered fixed. The convention is that values set to "TRUE" are treated as fixed parameters. Note that R is case sensitive, the logical values TRUE and FALSE must be in all capital letters. 
parameter.is.fixed <- c(rf = TRUE, kd.dye = FALSE, bt.dye = TRUE, kd.1 = FALSE, bt.1 = FALSE, kd.2 = TRUE, bt.2 = TRUE, kappa.nonsaturable = TRUE, accessible.volume = FALSE)

# The user must specify starting values for the parameters. The stochastic gradient descent algorithm is designed to escape local "false" optima in the parameter space that can be attributed to measurement noise. Therefore, it is not especially sensitive to starting estimates provided that a unique optimum exists. Extremely poor choices for starting parameters, however, may cause the algorithm to behave poorly. For parameters that are fixed, their values are fixed at the initial estimate provided in beta.
beta <- c(rf = 50, kd.dye = .380, bt.dye = 100, kd.1 = 0, bt.1 = 0, kd.2 = 0, bt.2 = 0, kappa.nonsaturable = 0, accessible.volume = 1)


#The algorithm requires lower and upper bounds on parameters. Poor performance can be guaranteed by improper specification of the upper and lower bounds on beta. Some simple physical considerations provide guidelines for choosing good boundaries. Clearly the kd for the dye and any endogenous buffers must be positive. If the kd of the endogenous buffer were too large then it would not even partially saturate over the range of observable free calcium values. Therefore, the upper bounds on the endogenous buffer kd's may be taken to be some multiple of the dye kd. Extremely low affinity buffers can be modeled as nonsaturable buffers. The nonsaturable buffering capacity, if present, must be strictly nonnegative and an upper limit must be estimated by considerations of the preparation. Strictly speaking, the fraction of accessible volume must be less than or equal to one, however, some allowance should be made for under estimation of the total volume.  The lower boundary for the fraction of accessible volume must be strictly greater than zero, and a suitable value may be estimated from knowledge of the preparation.
beta.lower <- c(rf = 20, kd.dye = .1, bt.dye = 50, kd.1 = 0, bt.1 = 0, kd.2 = 0, bt.2 = 0, kappa.nonsaturable = 0, accessible.volume = .05)

beta.upper <- c(rf = 200, kd.dye = 2, bt.dye = 200, kd.1 = 10, bt.1 = 1000, kd.2 = 10, bt.2 = 1000, kappa.nonsaturable = 200, accessible.volume = 1.1)

# num.steps: sets the initial learning rate eta = (beta.upper - beta.lower)/num.steps. Here eta is the learning rate. That is, the update rule for each iteration is beta = beta + eta * grad(Objective.Fxn)
num.steps <- 50
# decay.constant: controls the rate at which the learning rate eta "decays" for each iteration we have eta = eta * decay.constant
decay.constant <- .9995

# parameter.names specifies the names of the parameters, this is used only to write output files
parameter.names <- c("rf", "kd.dye", "bt.dye", "kd.endogenous.1", "bt.endogenous.1", "kd.endogenous.2", "bt.endogenous.2", "kappa.nonsaturable", "accessible.volume")

# The stochastic gradient descent algorithm does not intrinsicly generate reliable estimates for the standard errors of the parameters. By default no estimates of the standard error are provided. Estimates may be provided by setting do.bootstrap.estimate = TRUE. Note again that R is case sensitive, the logical values TRUE and FALSE must be all capital letters. If error estimates are desired, then the number of bootstrap replicates is specified by the value of bootstrap.replicates. A significant shortcoming of the total least squares estimation procedure is the computational time required, users are cautioned that calculating parameter estimates will require significant computational time. The fits to the bootstrapped data sets can be done concurrently on multicore processors. To use this facility, set n.threads to the desired number of threads. Due to the stochastic nature of the algorithm, it will occasionally get stuck at the boundaries of the parameter space, and may have difficulty moving back into the interior. The veracity of parameters under this case is suspect. By default, the program will discard estimates that are too close to the boundary. This behavior is controlled by the replace.on.boundary and boundary.margin variables. Boundary points are discarded if replace.on.boundary = TRUE, and retained otherwise. The criteria for discarding estimates is controlled by the bondary.margin variable. Estimates are discarded if the distance between the estimated parameter values and the parameter limits is less than boundary.margin.
do.bootstrap.estimate = FALSE
bootstrap.replicates = 8
n.threads = 4
replace.on.boundary = TRUE
boundary.margin = (beta.upper - beta.lower) * (.Machine$double.eps)^.5

#The data files may be specefied in one of two ways. If interactive.file.chooser = TRUE then, when this script is ran, the user will select the files interactively through a file chooser dialog box. If interactive.file.chooser = FALSE then the file names with full paths must be specefied below.

interactive.file.chooser = FALSE
#file path to the directory containing the data:
data.directory = "F:/2013/kd estimation/CalciumBufferingTLS"
#The fluorescence measurement data:
fluorescence.filename = "f.csv"
#Fluorescence measurement error estimates:  
fluorescence.se.filename = "f.sem.csv"
#Total calcium increment:
delta.ca.total.filename = "d.ca.csv"
#Total calcium increment error estimates:
delta.ca.total.se.filename = "d.ca.sem.csv"

#Advanced features
# parameter.boundary.margin: controls how close to the boundary the solution is allowed to be before restarting. If at any point we have abs(beta - beta.upper) < abs(beta.upper - beta.lower) * parameter.boundary.margin or abs(beta - beta.upper ) < abs(beta.upper - beta.lower) * parameter.boundary.margin, the beta is reset either randomly or the the best parameters so far
parameter.boundary.margin = .01
# p.restart: if the parameters are "too close to the boundary" as defined above, then the reset behavior described above is executed with probability p.restart
p.restart = 1
#p.restart.decay: if the parameters are reset, then p.restart = p.restart * p.restart.decay, this option is included to prevent the algorithm from restarting at every iteration, which may happen if the parameter limits are poorly chosen or on certain data sets
p.restart.decay = .5
# p.restart.grow: if the parameters are on the boundary and a restart action is not performed, then p.restart = p.restart * p.restart.grow
p.restart.grow = 1.1

# p.check.obj: evaluating the toal least squares objective function is very computationaly intensive, therefore we perform the evaluation with probability p.check.obj
p.check.obj = .05

# p.goto.best.par: with probably p.goto.best.par set beta to the best parameters observed so far
p.goto.best.par = .005

rescale.vector = c(1,.9,1,.9,.9,.9,.9,1,1.1111)
p.rescale = .005

#Done setting variables. Stop modifying here. When this script is executed, a sequence of file chooser dialog boxes will open. For each dialog box, the program will print to the R console a description of which file should be chosen. The program will request the location for the fluorescence data file, the fluorescence measurement standard errors, calcium increment, and calcium increment standard errors, in that order. After the data files are read, the fitting process will begin. A single run can be very time consuming, and calculating parameter standard errors from bootstraping will be proportionally longer by a factor of bootstrap.replicates/n.threads. 
#If do.bootstrap.estimates = FALSE, then the program will output status updates to the R console as described in the manual. If do.bootstrap.estimates = TRUE, then n.threads will be run in parallel. In this case, the parallel processing facility will supress output to the R console. Output is supressed even for n.threads=1.

#When the fitting process has finished, the results will be written to an output file. The output is qualitatively different depending on the do.bootstrap.estimate flag. For do.bootstrap.estimate = FALSE, the results are written to results.txt in human readable format. If do.bootstrap.estimate = TRUE, the output file is results.csv, with bootsrap.replicates rows. Each row corresponds to the solution for a single bootstrap data set.

#To run the script from inside RStudio: from the top menu, select code -> source.
source("Ca.Buffering.TLS.Sub.R")
