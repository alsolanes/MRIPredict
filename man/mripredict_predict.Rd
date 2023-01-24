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
An MRIPredict object, as created by the function mripredict.
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
                 response_var="female", covariates=c("cent_age"), predictor=c("score1"),# "clinical_scale1"),
                 response_family="binomial", modulation="fu")
# fit a model
mp <- mripredict_fit(mp, space = "MNI", save_name = "model_fit")
# apply previous model to test data
preds <- mripredict_predict(mp, mri_data = mri_test, data_table_file = clinical_test)

## End(**Not run**)
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
