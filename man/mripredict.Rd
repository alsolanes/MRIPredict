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
