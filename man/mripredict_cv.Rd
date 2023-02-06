\name{mripredict_cv}
\alias{mripredict_cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict Crossvalidation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to perform a crossvalidation.
}
\usage{
mripredict_cv(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = NULL, n_cores = 1, use_significant_voxels = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, n_folds = 10, ide_shiny = FALSE, standardize_images = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mp}{
%%     ~~Describe \code{mp} here~~
An MRIPredict object, as created by the function mripredict.
}
  \item{space}{
%%     ~~Describe \code{space} here~~
Images should be MNI aligned. This field can be "MNI" or "NO_CHECK" if images are not in MNI but the user decides to not check space alignment.
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
  \item{use_ensemble_voxels}{
%%     ~~Describe \code{use_ensemble_voxels} here~~
}
  \item{use_ensemble_subjects}{
%%     ~~Describe \code{use_ensemble_subjects} here~~
}
  \item{n_folds}{
%%     ~~Describe \code{n_folds} here~~
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
Joaquim Radua and Aleix Solanes
Maintainer: Aleix Solanes <email: solanes@recerca.clinic.cat>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
x = list(mri$unmodulated, mri$modulated)

# create an mripredict parameters object. Predict the sex covarying by age centered and considering the score1 predictor variable, along with mri data.
mp <- mripredict(mri_data=x, data_table_file = clinical, 
                 response_var="female", covariates=c("cent_age"), predictor=c("score1"), 
                 response_family="binomial",modulation="un")
                 
# perform a 10-fold cross-validation
mp <- mripredict_cv(mp, space = "NO_CHECK", 
                    save_name = "testing", use_ensemble_voxels = TRUE,
                    folds_file = NULL, 
                    n_cores = 1, 
                    n_folds = 10)
## End(**Not run**)
}
}
\keyword{ models }
\keyword{ crossvalidation }
\keyword{ regression }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
