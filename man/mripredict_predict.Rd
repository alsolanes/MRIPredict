\name{mripredict_predict}
\alias{mripredict_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict predict
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Apply a previously created mripredict model to new data
}
\usage{
mripredict_predict(mp, mri_data = NULL, mri_fu_paths_file = NULL, data_table_file = NULL, space, n_cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mp}{
%%     ~~Describe \code{mp} here~~
}
  \item{mri_data}{
%%     ~~Describe \code{mri_data} here~~
}
  \item{mri_fu_paths_file}{
%%     ~~Describe \code{mri_fu_paths_file} here~~
}
  \item{data_table_file}{
%%     ~~Describe \code{data_table_file} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{n_cores}{
%%     ~~Describe \code{n_cores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mp, mri_data = NULL, mri_fu_paths_file = NULL, data_table_file = NULL, 
    space, n_cores = 1) 
{
    .require("glmnet")
    .require("oro.nifti")
    .require("survival")
    .require("doParallel")
    .require("parallel")
    .require("logistf")
    SIGNIFICANCE_THRESHOLD = qnorm(0.975)
    N_ITERATIONS = mp$n_iterations
    N_M_IMPUTATIONS = mp$n_imputations
    if (n_cores == "auto") {
        n_cores = 2
    }
    .print_action("Setting new MRIPredict model")
    mp_test = mripredict(mri_data = mri_data, data_table_file = data_table_file, 
        response_var = mp$response_var, covariates = mp$covX_var, 
        predictor = mp$pred_var, response_family = mp$response_family, 
        modulation = mp$modulation)
    .print_action("Loading MRI data")
    a = Sys.time()
    mri = NULL
    if (!is.null(mp_test$mri_paths)) {
        mri = .load_mri(mp_test$mri_paths, mask = mp$mask$data, 
            space = space)
    }
    mri_fu = NULL
    if (!mp_test$modulation == "un" && !is.null(mp_test$mri_fu_paths)) {
        mri_fu = .load_mri(mp_test$mri_fu_paths, mp$mask_fu$data, 
            space = space)
    }
    mri_wm = NULL
    mri_wm_fu = NULL
    if (!is.null(mp_test$mri_wmmod_paths)) {
        mri_fu = .load_mri(mp_test$mri_fu_paths, mp$mask_fu$data, 
            space = space)
        mri_wm = .load_mri(mp_test$mri_wm_paths, mp$mask_wm$data, 
            space = space)
        mri_wm_fu = .load_mri(mp_test$mri_wmmod_paths, mp$mask_wmmod$data, 
            space = space)
    }
    cat("Time loading data:", difftime(Sys.time(), a, units = "mins"), 
        "mins.\n")
    .print_ok()
    .print_action("Creating response vector and covariate matrix")
    n_subjects = nrow(data_table_file)
    if (!is.null(mp$covX_transf)) 
        covX = .create_covX(data_table_file, mp$covX_transf)
    sites = NULL
    if (!is.null(mp$data_table_transf)) {
        data_informative_table = data.frame2glmnet.matrix(m = mp$data_table_transf, 
            x = data_table_file)
        if ("site" \%in\% colnames(data_table_file)) {
            sites = factor(x = data_table_file[, "site"])
        }
    }
    else {
        data_informative_table = NULL
        imp.data_informative_table_test = NULL
    }
    .print_ok()
    if (!is.null(mp$impute_obj)) {
        imp.data_informative_table_test = impute.glmnet.matrix(m = mp$impute_obj, 
            x = data_informative_table, nimp = N_M_IMPUTATIONS)
    }
    linPreds = c()
    i = 1
    for (iter in 1:N_ITERATIONS) {
        for (iter_imputation in 1:N_M_IMPUTATIONS) {
            if (!is.null(mp$impute_obj)) 
                data_informative_table = imp.data_informative_table_test[[iter_imputation]]
            model_list = mp$models[[i]]
            i = i + 1
            lasso_ijk = model_list$lasso_ijk
            signif_indx = model_list$signif_indx
            lasso_covB = model_list$lasso_covB
            lasso_covB_fu = model_list$lasso_fu_covB
            mask = model_list$mask
            scale_clinical = model_list$predictors_scale
            if (!is.null(data_informative_table)) {
                if (!is.null(mp$covX_var)) {
                  covX_test = cbind(1, matrix(data_informative_table[, 
                    mp$covX_var], nrow = nrow(data_informative_table)))
                }
                else {
                  covX_test = data_informative_table
                }
            }
            else {
                if (is.null(data_table_file)) {
                  covX_test = matrix(1, nrow = mri$n)
                }
                else {
                  covX_test = matrix(1, nrow = nrow(data_table_file))
                }
            }
            if (!is.null(mp$pred_var)) {
                predX_test = matrix(data_informative_table[, 
                  mp$pred_var], nrow = nrow(data_informative_table))
            }
            tipett_take_un = model_list$tipett_take_un
            preds = apply_model(mp = mp, mri_data = mri$data, 
                mri_fu_data = mri_fu$data, mri_wm_data = mri_wm$data, 
                mri_wm_fu_data = mri_wm_fu$data, covX_test = covX_test, 
                signif_indx = model_list$signif_indx, lasso_covB = lasso_covB, 
                lasso_covB_fu = lasso_covB_fu, mask = model_list$mask, 
                masks_3d = model_list$masks_3d, predX_test = predX_test, 
                scale_clinical = model_list$scale_predictors, 
                scale_mri = model_list$scale_mri, lasso = model_list$lasso, 
                lasso_predX_indx = model_list$lasso_predX_indx, 
                covX_site = sites, n_voxels_mask = model_list$n_voxels_mask, 
                combat = model_list$combat)
            linPreds = cbind(linPreds, preds)
        }
    }
    linPred = rowMeans(linPreds)
    pred = c()
    if (mp$response_family == "binomial") {
        prob = 1/(1 + exp(-linPred))
        pred = prob
        label = pred
        label[pred <= mp$bin_threshold] = mp$response_ref
        label[pred > mp$bin_threshold] = mp$response_event
        names(pred) = label
    }
    else if (mp$response_family == "gaussian") {
        pred = as.matrix(linPred)
    }
    else if (mp$response_family == "cox") {
        pred = linPred
    }
    pred
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
