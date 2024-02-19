#' Generates data for simulation
#'
#' This function loads images specified in filelist and adds signal to the
#'  images. The amount and structure of signal is determined by betaimg and X.
#'  The output images are equal the images specified by filelist plus betaimg
#'  times the column of X that is not in Xred. This column is first
#'  residualized to Xred.
#' @param files a vector of .nii or .nii.gz images.
#' @param data a simulation dataset from simSetup function used for simulating y
#' @param n4sim total number of subjects that we're simulating
#' @param N4sim total number of measuremetns that we're simulating
#' @param nMeas Approximate average number of measurements per subject. Will not be used when by = 'individuals'.
#' @param formres formula used in linear mixed effect model
#' @param mask image mask where data exist.
#' @param method Method to generate data, either "synthetic" (i.e. multivariate normal) or "bootstrap."
#' @param lambda Simulated data are a convex coombination of normally distributed data with bootstrapped data. lambda=1 is fully bootstrapped data, lambda=0 is fully synthetic.
#' @param outfiles a vector of images to save the output.
#' @keywords  null power simulation
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom RNifti writeNifti
#' @export
# @examples

genSimData = function(files, data, outfiles=NULL, n4sim, nMeas=NULL, N4sim, formres, mask=NULL, method=c('bootstrap', 'synthetic'), lambda=0.5 ){
  if(tolower(method[1])=='synthetic'){
    # THIS CODE SIMULATES DATA FOR THE FULL SAMPLE, if n in simulation is larger than n in sample then it will create problems
    # contains residuals for entire study
    # files = readRDS(syntheticSqrt)
    # cov = estimates
    cov = if(length(files)==1) readRDS(files[1]) else files
    if(!is.list(cov)){ # list represents lmer
      y = (matrix(rnorm(nrow(cov)*ncol(cov)), nrow=ncol(cov), ncol=nrow(cov)) %*% cov)/sqrt(nrow(cov))
      rm(cov)
      temp = if(is.character(mask)) readNifti(mask) else mask
      trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
      RNifti::writeNifti(temp, outfiles[ind])
      })
      rm(y)
    } else {
      # This is a hybrid of a bootstrap and parametric simulation
      # generate y using sample residuals, coefficients, and random effects
      n = nrow(cov$resids)  # number of observations
      nRE = nrow(cov$ranefs) # number of subjects
      vare = apply(cov$resids,2,var)
      varb = apply(cov$ranefs,3,var)
      coef = cov$coefs # 15 * n
      sigma2 = cov$vary
      # N4sim = round(n4sim * nMeas)

      # simulation data is ordered by subjects
      id = sort(rep(1:n4sim, ceiling(N4sim/n4sim))[1:N4sim]) # ordered fake id (index)
      samp = sample(n, N4sim, replace=TRUE)
      reSamp_id = sample(nRE, n4sim, replace=TRUE)
      n_id = as.numeric(table(id))
      reSamp = rep(reSamp_id, n_id)

      # currently written only for random intercept
      # REs = matrix(rnorm(n4sim*nRE)/sqrt(nRE), nrow=n4sim) %*% cov$ranefs[,1,] # new random intercepts
      # epsilon = matrix(rnorm(N4sim * n)/sqrt(n), nrow=N4sim) %*% cov$resids # new residuals

      # convex combination of bootstrap resample and parametric resample
      # also assumes random intercept model only
      # REs = (1-lambda)* REs + lambda * cov$ranefs[reSamp,1,] # what is that combination
      # y = REs[id,] + sqrt(1-lambda) * epsilon + lambda * cov$resids[samp,]

      # bootstrap resample
      R_ranef = sweep(sweep(cov$ranefs[reSamp,1,], 2, sapply(varb, sqrt), FUN = "/"),
                      2, sapply(sigma2, sqrt), FUN = "*")
      R_resid = sweep(sweep(cov$resids[samp,], 2, sapply(vare, sqrt), FUN = "/"),
                      2, sapply(sigma2, sqrt), FUN = "*")

      # modify the formres to remove the random effect
      # formresfix = remove_ranef(form)
      formresfix = remove_ranef(simConfig$formres)
      # data = readRDS(file.path(simdir, 'data.rds'))
      y = model.matrix(as.formula(formresfix), data = data) %*% coef +
        sqrt(lambda) * R_ranef + sqrt(1-lambda) * R_resid

      # write out Nifti images
      temp = if(is.character(mask)) readNifti(mask) else mask
      trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
      RNifti::writeNifti(temp, outfiles[ind])
      })
      # return simulated data if requested
      list(outfiles=outfiles, simdata=y, id=id)
    }
  } else if(tolower(method[1])=='bootstrap'){
    # the bootstrapping has already been performed in the simulation Setup phase
    result = file.copy(files, outfiles)
  } else{
    stop('genSimData method is not correctly specified.')
  }
}
