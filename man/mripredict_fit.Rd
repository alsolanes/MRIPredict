\name{mripredict_fit}
\alias{mripredict_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mripredict_fit(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1, use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, ide_shiny = FALSE, standardize_images = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mp}{
%%     ~~Describe \code{mp} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{save_name}{
%%     ~~Describe \code{save_name} here~~
}
  \item{preloaded_covB_path}{
%%     ~~Describe \code{preloaded_covB_path} here~~
}
  \item{preloaded_covB_fu_path}{
%%     ~~Describe \code{preloaded_covB_fu_path} here~~
}
  \item{folds_file}{
%%     ~~Describe \code{folds_file} here~~
}
  \item{n_cores}{
%%     ~~Describe \code{n_cores} here~~
}
  \item{use_significant_voxels}{
%%     ~~Describe \code{use_significant_voxels} here~~
}
  \item{use_ensemble_learning}{
%%     ~~Describe \code{use_ensemble_learning} here~~
}
  \item{use_ensemble_voxels}{
%%     ~~Describe \code{use_ensemble_voxels} here~~
}
  \item{use_ensemble_subjects}{
%%     ~~Describe \code{use_ensemble_subjects} here~~
}
  \item{ide_shiny}{
%%     ~~Describe \code{ide_shiny} here~~
}
  \item{standardize_images}{
%%     ~~Describe \code{standardize_images} here~~
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
function (mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, 
    preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1, 
    use_significant_voxels = FALSE, use_ensemble_learning = FALSE, 
    use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, 
    ide_shiny = FALSE, standardize_images = FALSE) 
{
    .require("glmnet")
    .require("oro.nifti")
    .require("survival")
    .require("doParallel")
    .require("logistf")
    data_informative_table = NULL
    imp.data_informative_table_train = NULL
    imp.data_informative_table_test = NULL
    SIGNIFICANCE_THRESHOLD = qnorm(0.975)
    if (use_ensemble_learning) {
        N_ITERATIONS = 18
    }
    else {
        N_ITERATIONS = 1
    }
    n_multiple_imputations = 1
    N_M_IMPUTATIONS = 20
    if (n_cores == "auto") {
        n_cores = max(round(detectCores()/2), 1)
    }
    mp$mask$data <- mp$mask_fu$data <- mp$mask_wm$data <- mp$mask_wmmmod$data <- NULL
    if (!is.null(mp$mask_path)) 
        mp$mask = .load_mri(mp$mask_path, space = "NO_CHECK")
    if (!is.null(mp$mask_fu_path)) 
        mp$mask_fu = .load_mri(mp$mask_fu_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wm_path)) 
        mp$mask_wm = .load_mri(mp$mask_wm_path, space = "NO_CHECK")
    if (!is.null(mp$mask_wmmod_path)) 
        mp$mask_wmmod = .load_mri(mp$mask_wmmod_path, space = "NO_CHECK")
    .print_action("Loading MRI data")
    a = Sys.time()
    mri = NULL
    if (!is.null(mp$mri_paths)) {
        mri = .load_mri(mp$mri_paths, mask = mp$mask$data, space = space)
    }
    mri_fu = NULL
    if (!mp$modulation == "un" && !is.null(mp$mri_fu_paths)) {
        mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
            space = space)
    }
    mri_wm = NULL
    mri_wm_fu = NULL
    if (!is.null(mp$mri_wmmod_paths)) {
        mri_fu = .load_mri(mp$mri_fu_paths, mp$mask_fu$data, 
            space = space)
        mri_wm = .load_mri(mp$mri_wm_paths, mp$mask_wm$data, 
            space = space)
        mri_wm_fu = .load_mri(mp$mri_wmmod_paths, mp$mask_wmmod$data, 
            space = space)
    }
    cat("Time loading data:", difftime(Sys.time(), a, units = "mins"), 
        "mins.\n")
    .print_ok()
    .print_action("Creating response vector and covariate matrix")
    Y = c()
    if (mp$response_family == "cox") {
        Y = mp$data_table[, match(mp$response_var, colnames(mp$data_table))]
        Y = Surv(time = as.numeric(Y[, 1]), event = as.numeric(Y[, 
            2]))
    }
    else {
        Y = .create_Y(mp$data_table, mp$response_var, mp$response_event)
    }
    n_subjects = nrow(mp$data_table)
    if (!is.null(mp$covX_transf)) 
        covX = .create_covX(mp$data_table, mp$covX_transf)
    if (!is.null(mp$data_table_transf)) 
        data_informative_table = data.frame2glmnet.matrix(mp$data_table_transf, 
            mp$data_table)
    .print_ok()
    if ("site" \%in\% colnames(mp$data_table)) {
        sites = factor(x = mp$data_table[, "site"])
    }
    else {
        sites = NULL
    }
    n_subjects = min(length(Y), nrow(Y))
    subjects_used = matrix(nrow = n_subjects, ncol = N_ITERATIONS)
    subjects_used = list()
    z_cuts = c()
    train_thresholds_all = c()
    mp$models = list()
    model_counter = 1
    time_points = c(30, 60, 180, 360, 720)
    if (!is.null(data_informative_table)) {
        data_table_imputed_train = data_informative_table
        if (any(is.na(data_informative_table))) {
            n_multiple_imputations = N_M_IMPUTATIONS
            impute_obj = impute.glmnet.matrix_fit(x = data_table_imputed_train, 
                n_cores = 4)
            imp.data_informative_table_train = impute.glmnet.matrix(m = impute_obj, 
                x = data_table_imputed_train, nimp = N_M_IMPUTATIONS)
            mp$impute_obj = impute_obj
        }
        else {
            n_multiple_imputations = 1
            mp$impute_obj = NULL
        }
    }
    else {
        data_table_imputed_train = NULL
        imp.data_informative_table_train = NULL
    }
    start_fold.time <- Sys.time()
    mp$n_iterations = N_ITERATIONS
    mp$n_imputations = n_multiple_imputations
    for (iter in 1:N_ITERATIONS) {
        if (N_ITERATIONS > 1 && use_ensemble_subjects) {
            subjects_used[[iter]] = training[sample(1:length(training), 
                replace = TRUE)]
            training = subjects_used[[iter]]
        }
        if (mp$response_family == "cox") {
            trainY = Y
        }
        else {
            trainY = Y
            if (mp$response_family == "binomial") {
                mp$bin_threshold = sum(trainY == 1)/length(trainY)
            }
        }
        linPreds = c()
        if (mp$response_family == "cox") {
            train_thresholds = NULL
        }
        start.time <- Sys.time()
        if (N_ITERATIONS > 1) 
            cat("Iteration ", iter, "of", N_ITERATIONS, "\n", 
                save_name)
        sprintf("/nIteration: \%d", iter)
        preds = c()
        for (iter_imputation in 1:n_multiple_imputations) {
            if (n_multiple_imputations > 1 & !is.null(imp.data_informative_table_train)) {
                data_table_imputed_train = imp.data_informative_table_train[[iter_imputation]]
            }
            if (n_multiple_imputations > 1) 
                .print_action(paste("Imputation", iter_imputation, 
                  "of", n_multiple_imputations, "\n", save_name))
            internal_folds = c()
            internal_folds = assign.folds(y = trainY, family = mp$response_family, 
                nfolds = 10, site = sites)
            model_list = fit_model(mp = mp, data_informative_table = data_table_imputed_train, 
                Y = trainY, mri = mri$data, mri_fu = mri_fu$data, 
                mri_wm = mri_wm$data, mri_fu_wm = mri_wm_fu$data, 
                iter = ifelse(use_ensemble_voxels == FALSE, -1, 
                  iter), internal_folds = internal_folds, n_cores = n_cores, 
                use_significant_voxels = use_significant_voxels, 
                covX_site = sites, standardize_images = standardize_images)
            if (mp$modulation \%in\% c("cl", "clinical")) {
                n_voxels_mask = 0
            }
            tipett_take_un = model_list$tipett_take_un
            scale_clinical = model_list$scale_predictors
            scale_mri = model_list$scale_mri
            mp$models[[model_counter]] <- model_list
            model_counter <- model_counter + 1
            .print_ok()
            if (ide_shiny) 
                incProgress(1/(n_folds * N_ITERATIONS * n_multiple_imputations), 
                  detail = paste("Doing crossvalidation. Fold", 
                    fold, ". Iteration:", iter, ". Imputation:", 
                    iter_imputation))
        }
        .print_ok()
        end.time <- Sys.time()
        if (iter > 1) 
            message("Time per iteration:", difftime(end.time, 
                start.time, units = "mins"), " mins.")
    }
    cat("\n[Saving] - Saving MRIPredict object to:", sprintf("\%s_mp.rds", 
        save_name))
    saveRDS(mp, file = sprintf("\%s_mp.rds", save_name))
    mp
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
