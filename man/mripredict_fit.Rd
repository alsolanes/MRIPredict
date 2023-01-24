\name{mripredict_fit}
\alias{mripredict_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Train a model to be applied later on new data.
}
\usage{
mripredict_fit(mp, space = "MNI", save_name = "results_cv", preloaded_covB_path = NULL, preloaded_covB_fu_path = NULL, folds_file = "", n_cores = 1, use_significant_voxels = FALSE, use_ensemble_learning = FALSE, use_ensemble_voxels = FALSE, use_ensemble_subjects = FALSE, ide_shiny = FALSE, standardize_images = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mp}{
An MRIPredict object, as created by the function mripredict.
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
person("Aleix", "Solanes", role = "aut",
            email = "solanes@recerca.clinic.cat")
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

# create an object of class mripredict
  mp <- mripredict(mri_data = mri_train, data_table_file = clinical_train,
                 response_var="female", covariates=c("cent_age"), predictor=c("score1"),
                 response_family="binomial", modulation="fu")
# train a model trying to predict sex, covarying per centered age, and score1 as a predictor variable
  mp <- mripredict_fit(mp, space = "MNI", save_name = "model_fit")

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
