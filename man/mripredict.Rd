\name{mripredict}
\alias{mripredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Create an object of type "mripredict". It contains the needed configuration for mripredict_cv, mripredict_fit and mripredict_predict.
}
\usage{
mripredict(mri_data = NULL, data_table_file, response_var, covariates = NULL, predictor = NULL, response_family, modulation, information_variables = NULL, mask_files_path = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mri_data}{
%%     ~~Describe \code{mri_data} here~~
}
  \item{data_table_file}{
%%     ~~Describe \code{data_table_file} here~~
}
  \item{response_var}{
%%     ~~Describe \code{response_var} here~~
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
}
  \item{predictor}{
%%     ~~Describe \code{predictor} here~~
}
  \item{response_family}{
%%     ~~Describe \code{response_family} here~~
}
  \item{modulation}{
%%     ~~Describe \code{modulation} here~~
}
  \item{information_variables}{
%%     ~~Describe \code{information_variables} here~~
}
  \item{mask_files_path}{
%%     ~~Describe \code{mask_files_path} here~~
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
function (mri_data = NULL, data_table_file, response_var, covariates = NULL, 
    predictor = NULL, response_family, modulation, information_variables = NULL, 
    mask_files_path = NULL) 
{
    if ((is.null(mri_data)) && (is.null(predictor))) {
        .check_stop("No MRI files or predictor variables found.")
    }
    .print_action("Setting new MRIPredict model")
    mri_paths = NULL
    mri_fu_paths = NULL
    mri_wm_paths = NULL
    mri_wmmod_paths = NULL
    mask_path = NULL
    mask_fu_path = NULL
    mask_wm_path = NULL
    mask_wmmod_path = NULL
    if (is.null(mri_data) & modulation != "clinical") {
        stop("No MRI data is selected, and the modulation is not set to 'clinical'. Please introduce MRI data or change modulation to 'clinical'.")
    }
    if (modulation != "clinical") {
        if (!is.null(mri_data)) 
            mri_data = as.data.frame(mri_data)
        if (length(mri_data) == 1 && !all(grepl(".nii.gz|.nii", 
            mri_data))) {
            if (all(mri_data != "")) {
                if (file.exists(mri_data)) {
                  info = file.info(mri_data)
                  if (info$size > 0) {
                    mri_paths = .read_1col_file(mri_data)
                    if (!is.null(mask_files_path)) 
                      mask_path = mask_files_path[1]
                    if (modulation != "un") {
                      mri_fu_paths = .read_2col_file(mri_data)
                      if (!is.null(mask_files_path)) 
                        mask_fu_path = mask_files_path[2]
                    }
                    if (modulation == "all") {
                      mri_wm_paths = .read_3col_file(mri_data)
                      mri_wmmod_paths <- .read_4col_file(mri_data)
                      if (!is.null(mask_files_path)) {
                        mask_wm_path = mask_files_path[3]
                        mask_wmmod_path = mask_files_path[4]
                      }
                    }
                  }
                  else {
                    mri_data = NULL
                  }
                }
                else {
                  stop("MRI paths file does not exist.")
                }
            }
        }
        else if (dim(mri_data)[2] <= 4) {
            mri_data <- as.data.frame(mri_data)
            if (all(mri_data != "")) {
                mri_paths = mri_data[, 1]
                if (ncol(mri_data) > 1 && modulation != "un") 
                  mri_fu_paths = mri_data[, 2]
                if (ncol(mri_data) > 2 && modulation == "all") {
                  mri_wm_paths = mri_data[, 3]
                  mri_wmmod_paths = mri_data[, 4]
                }
            }
            if (!is.null(mask_files_path)) 
                mask_path = mask_files_path[1]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                1) 
                mask_fu_path = mask_files_path[2]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                2) 
                mask_wm_path = mask_files_path[3]
            if (!is.null(mask_files_path) && length(mask_files_path) > 
                3) 
                mask_wm_fu_path = mask_files_path[4]
        }
        else if (is.numeric(mri_data) && dim(mri_data)[2] > 4) {
            mri_paths <- mri_data
        }
        else if (is.list(mri_data)) {
            mri_paths <- mri_data$gm_un
            mri_fu_paths <- mri_data$gm_fu
            mri_wm_paths <- mri_data$wm_un
            mri_wmmod_paths <- mri_data$wm_fu
        }
        else {
            stop("Wrong MRI data format.")
        }
    }
    if (!is.null(data_table_file) && any(data_table_file != "")) {
        if (length(data_table_file) == 1) {
            data_table = .read_data_table(data_table_file)
        }
        else {
            data_table = data_table_file
        }
    }
    else {
        data_table = NULL
    }
    if (!is.null(data_table_file)) 
        response = .check_response_var(data_table, response_var)
    else response = NULL
    response_levels = c()
    if (response_family == "cox") {
        response_levels = unique(response[, 2])
        if (length(response_levels) != 2) {
        }
    }
    else if ((response_family != "gaussian")) {
        response_levels = unique(response)
        if (length(response_levels) != 2) {
            mp = list(mri_paths = mri_paths, mri_fu_paths = mri_fu_paths, 
                mri_wm_paths = mri_wm_paths, mri_wmmod_paths = mri_wmmod_paths, 
                data_table = data_table, response_var = response_var, 
                response_ref = response_levels[1], response_event = response_levels[2], 
                modulation = modulation)
            attr(mp, "status") <- paste("Response variable", 
                response_var, "should have two levels")
            return(mp)
        }
    }
    mp = list(mri_paths = mri_paths, mri_fu_paths = mri_fu_paths, 
        mri_wm_paths = mri_wm_paths, mri_wmmod_paths = mri_wmmod_paths, 
        mask_path = mask_path, mask_fu_path = mask_fu_path, mask_wm_path = mask_wm_path, 
        mask_wmmod_path = mask_wmmod_path, data_table = data_table, 
        response_var = response_var, response_ref = response_levels[1], 
        response_event = response_levels[2], modulation = modulation)
    if (!is.null(covariates) & length(covariates) != 0) {
        data_cov = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            covariates)
        covX_info_fit = data.frame2glmnet.matrix_fit(data_cov)
        covX_info = data.frame2glmnet.matrix(covX_info_fit, data_cov)
        mp$covX_transf = covX_info_fit
        mp$covX_var = colnames(covX_info)
    }
    else {
        mp$covX_transf = NULL
        mp$covX_var = NULL
    }
    if (is.null(information_variables) & (!is.null(covariates) | 
        !is.null(predictor))) {
        information_variables = c(covariates, predictor)
    }
    if (!is.null(information_variables) & (length(information_variables) != 
        0)) {
        data_info = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            information_variables)
        data_info_fit = data.frame2glmnet.matrix_fit(data_info)
        data_info = data.frame2glmnet.matrix(data_info_fit, data_info)
        mp$data_table_transf = data_info_fit
        mp$data_table_var = information_variables
    }
    else {
        mp$data_table_transf = NULL
        mp$data_table_var = NULL
    }
    if (!is.null(predictor) & length(predictor) != 0) {
        data_pred = seleccio_covariables_data.frame2glmnet.matrix_fit(data_table, 
            predictor)
        pred_info_fit = data.frame2glmnet.matrix_fit(data_pred)
        pred_info = data.frame2glmnet.matrix(pred_info_fit, data_pred)
        mp$pred_transf = pred_info_fit
        mp$pred_var = colnames(pred_info)
    }
    else {
        mp$pred_transf = NULL
        mp$pred_var = NULL
    }
    mp$response_family = response_family
    if (!is.null(mri_paths)) {
        mp$mri_params = .load_mri(mri_paths = mri_paths[1])
        mp$mri_params$data <- NULL
    }
    class(mp) = "mripredict"
    .print_ok()
    attr(mp, "status") <- "OK"
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
