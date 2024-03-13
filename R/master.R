#' Sets up a Neuroimaging bootstrap-based simulation given a set of images and covariates.
#'
#' The images can be images from a real data set.
#' The bootstrap-based simulation conditions on the distribution of the sample by drawing subsets with replacement from the sample.
#' @param images String vector containing paths to the images used for the simulation.
#' @param data A dataframe with number of rows equal to the length of the images variable with corresponding covariates.
#' @param outdir A directory to save the output files that are used for the simulations.
#' @param nsim Number of simulations to setup for each sample size.
#' @param ns Number of subjects to evaluate for each simulation. Should be less than the number of images/subjects.
#' @param nMeas Approximate average number of measurements per subject.
#' @param by Specifications of the level (observation or individual) on which bootstrap is performed.
#' @param mask If performing simulations under an alternative, signal can be added to the images within the mask.
#' @param formres Formula used for the mixed effect model.
#' @param W weights assigned for each subject.
#' @param ncores number of cores for parallel commands.
#' @return Returns a data frame of directories where simulation setup files are stored. The rds file in each directory saves the data frame and the image locations in the variable "images".
#' @importFrom RNifti readNifti writeNifti
#' @importFrom pbmcapply pbmclapply
#' @export
# prepare the output directories
# All the randomization happens within this loop

simSetup = function(images, data, outdir, nsim=1000, ns=c(50, 100, 200, 400), nMeas=1, by='observations',
                    mask=NULL, formres, W, ncores=parallel::detectCores() ){

  if(by == 'observations'){
    # sample by observations with a fixed pre-specified sample size - fixed sample size
    sims = expand.grid(sim=1:nsim, n=ns, nMeas=nMeas, simdir=NA)
    data$images = images
    sims$simdir = do.call(c, pbmcapply::pbmclapply(1:nrow(sims), function(simind) {
      simdir = file.path(outdir, paste0('sim', sims[simind,'sim']), paste0('n', sims[simind,'n'], 'nMeas', sims[simind,'nMeas']) ) # create a path for each simulation setting (nMeas + n)
      n = sims[simind, 'n']
      dir.create(simdir, showWarnings=FALSE, recursive = TRUE)
      unlink(file.path(simdir, '*.nii.gz'), recursive=TRUE)
      # create random sample from demographics and roi data
      N4sim = round(n * nMeas) # total measurements
      tempdata = data[sample.int(nrow(data), N4sim, replace=TRUE),
                      c(all.vars(as.formula(formres)), W)] # random sample rows in data
      saveRDS(tempdata, file=file.path(simdir,'data.rds'))
      simdir
    }, mc.cores = ncores, mc.preschedule=FALSE) )
  }else if(by == 'individuals'){
    # sample by individuals with all their observations - fixed subject size
    id = unique(data$subjectid)

    sims = expand.grid(sim=1:nsim, n=ns, simdir=NA)
    data$images = images
    sims$simdir = do.call(c, pbmcapply::pbmclapply(1:nrow(sims), function(simind) {
      simdir = file.path(outdir, paste0('sim', sims[simind,'sim']), paste0('n', sims[simind,'n']) ) # create a path for each simulation setting (nMeas + n)
      n = sims[simind, 'n']
      dir.create(simdir, showWarnings=FALSE, recursive = TRUE)
      unlink(file.path(simdir, '*.nii.gz'), recursive=TRUE)
      # create random sample from demographics and roi data
      # N4sim = round(n * nMeas) # total measurements
      tempid = id[sample.int(length(id), n, replace=TRUE)]
      tempindex = unlist(sapply(tempid, function(x) which(data$subjectid == x)))
      tempdata = data[tempindex, c(all.vars(as.formula(formres)), W)] # random sample rows in data
      saveRDS(tempdata, file=file.path(simdir,'data.rds'))
      simdir
    }, mc.cores = ncores, mc.preschedule=FALSE) )
  }else{
    stop("Argument 'by' can be 'observations' or 'individuals'.")
  }
  sims
}


#' Creates parameter image in random locations for simulations
#'
#' Uses the mask variable to create length(rs) spheres of size rs with values equal to betas within the mask.
#' A random voxel is chosen within the mask and a sphere is placed at that location.
#' The spheres are masked by the mask image so that no parameter values exist outside the mask.
#' @param mask String or niftiImage object indication the area to select sphere locations from.
#' @param parameterImage output nifti for the parameter image.
#' @param rs vector of radii for signal spheres.
#' @param betas vector of parameters for signal spheres. Signal is constant throughout the sphere.
#' Model parameter values are set to params = betas * sd(y)/sd(x).
#' @return Returns the parameter image after writing it to file.path(outdir, 'signal.nii.gz').
#' @importFrom RNifti readNifti writeNifti
#' @export
parameterImage = function(mask, parameterImage, rs, betas){
  if(is.character(mask)) mask = readNifti(mask)
  outfile = mask
  inds = which(mask==1, arr.ind=TRUE)
  # random location in gray matter mask for center voxel
  centers = inds[sample.int(nrow(inds), length(rs)),, drop=FALSE]
  outfile[,,] = 0
  for(rind in 1:length(rs)){
    r = rs[rind]
    center = centers[rind,]
    inds = as.matrix(expand.grid(seq(-r,r), seq(-r,r), seq(-r,r)))
    inds = inds[ sqrt(rowSums(inds^2))<=r,]
    inds = sweep(inds, 2, center, '+') # sphere around center voxel
    # make sure all voxels are in the image
    inds = inds[apply(sapply(1:3, function(x) inds[,x]>0 & inds[,x]<=dim(outfile)[x]), 1, all), ]
    # It is possible that a small cluster can be captured inside a big cluster and then it won't exist in that simulation
    outfile[ inds ] = betas[rind]
    # signal only in the mask
  } # end for(rind)
  outfile = outfile * mask
  writeNifti(outfile, file=parameterImage)
  outfile
}


#' Runs simulations
#'
#' Uses the mask variable to create length(rs) spheres of size rs with values equal to betas within the mask.
#' A random voxel is chosen within the mask and a sphere is placed at that location.
#' The spheres are masked by the mask image so that no parameter values exist outside the mask.
#' @param simdirs Vector of simulation directories created by simSetup.
#' @param sims Vector of integers specifying simulation index. Does not have to be unique. Can be used for running a particular analysis in every kth simulation.
#' @param n Integer sample size.
#' @param simfunc Function to evaluate on the simulated data.
#' @param simfuncArgs List of arguments passed to simfunc.
#' @param mask mask image argument passed to genSimData. Only required for genSimData function if betaimg is not NULL or if method='synthetic'
#' @param method method argument to generate sample data in genSimData. Data can be generated by bootstrapping or synthetically.
#' @param syntheticSqrt If method is synthetic this is a square root of the covariance matrix used to generate the imaging data. Passed as an RDS file or a matrix.
#' @param ncores Number of cores to use.
#' @param ... Arguments passed to pbmclapply
#' @return Returns the parameter image after writing it to file.path(outdir, 'signal.nii.gz').
#' @importFrom RNifti readNifti writeNifti
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom pbj addSignal
#' @export

runSim = function(simdirs, n, sim, simfunc, simfuncArgs=NULL, mask=NULL, method=c('bootstrap', 'synthetic'), syntheticSqrt, ncores=parallel::detectCores(), ...){
  method = tolower(method)
  # simdir = simdirs$simdir[1]; sim = simdirs$sims[1]; nMeas=simdirs$nMeas[1]; n = simdirs$n[1]
  result = pbmcapply::pbmcmapply(function(simdir, n, simfunc, method, mask){
    # load data
    dat = readRDS(file.path(simdir, 'data.rds'))
    N4sim = dim(dat)[1] # number of observations in simulation
    # create names for temporary files
    dat$tmpfiles = tempfile(paste0(basename(dirname(simdir)), '/s', 1:nrow(dat)), fileext='.nii.gz')
    dir.create(dirname(dat$tmpfiles[1]), recursive=TRUE, showWarnings = FALSE)
    if(method=='bootstrap'){
      genSimData(files=dat$images, data = dat, outfiles=dat$tmpfiles, mask=mask, method=method)
    } else {
      dat$id = factor(genSimData(files=syntheticSqrt, data = dat, outfiles=dat$tmpfiles, mask=mask, form = simConfig$formres, method=method, n4sim=n, N4sim=N4sim)$id)
    }
    dat$images = dat$tmpfiles
    # simfuncArgs$sim = sim
    simfuncArgs$data = dat
    result = do.call(simfunc, args = simfuncArgs)
    unlink(dat$images)
    return(result)
  }, simdir=simdirs, n=n, sim = sim, MoreArgs=list(simfunc = simfunc, method = method, mask=mask), mc.cores = ncores)
  return(result)
}


#' List memory usage of all objects
#'
#' @param units character what units to use, passed to format for object.size.
#' @export
#' @importFrom utils object.size
memoryUse = function(units='MiB'){
  sort( sapply(ls(),function(x){format(object.size(get(x)), units=units)}))
}


#' Remove Random Effects from Mixed Model Formula
#'
#' This function takes a mixed model formula string and removes any random effect terms, leaving only the fixed effects.
#' It is designed to work with formula strings that follow the lme4 package notation for random effects.
#'
#' @param formula_string A character string representing a mixed model formula.
#' @return A character string of the formula with random effects removed.
#' @examples
#' mixed_model_formula_string <- "~ sex + group*time + ns(age, df=10) + (1|subjectid)"
#' fixed_effects_formula_string <- remove_random_effects(mixed_model_formula_string)
#' @export
remove_ranef <- function(formula_string) {
  # Use a regular expression to remove the random effects
  # This regex looks for the pattern (1|something) and removes it
  fixed_only_string <- gsub("\\+\\s*\\([^)]+\\|[^)]+\\)", "", formula_string, perl = TRUE)

  # Replace any double spaces that may have been created with a single space
  fixed_only_string <- gsub("\\s+", " ", fixed_only_string)
  # Trim leading and trailing whitespace
  fixed_only_string <- trimws(fixed_only_string)

  # Ensure the formula starts with the tilde (~)
  if (substr(fixed_only_string, 1, 1) != "~") {
    fixed_only_string <- paste0("~ ", fixed_only_string)
  }
  return(fixed_only_string)
}

