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
mripredict_predict(mp, mri_data = NULL, clinical_data = NULL, space, 
n_cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mp}{
An MRIPredict object, as created by the function mripredict.
}
  \item{mri_data}{
%%     ~~Describe \code{mri_data} here~~
A list of paths, sorted in the same order as the clinical_data dataframe.
}
  \item{clinical_data}{
%%     ~~Describe \code{clinical_data} here~~
A data.frame or csv file containing the relevant clinical information for each 
subject. Each row correspond to a subject, and they are in the same order as the 
list of paths of MRI images.
}
  \item{space}{
%%     ~~Describe \code{space} here~~
A string that can be 'MNI' or 'NO_CHECK'. Images are supposed to be aligned to 
MNI, this can be skipped by setting this to 'NO_CHECK'. 
'NO_CHECK' should be used only if the user is sure that all images are correctly 
aligned to the same space. Also take into account, 
that the software will output MNI coordinates for the MRI images, so if images 
are not MNI-aligned these values will not be reliable.
}
  \item{n_cores}{
%%     ~~Describe \code{n_cores} here~~
Number of cores to be used by the software. Take into account that more cores 
means more RAM, so take care.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An MRIPredict object that can be accessed like a common list. It contains the 
individual predictions along with information about the model.
}
\references{
%% ~put references to the literature/web site here ~
Solanes, A. (2022), et al. Combining MRI and clinical data to detect high 
relapse risk after the first episode of psychosis. Schizophrenia (Heidelb). 
2022 Nov 17;8(1):100. doi: 10.1038/s41537-022-00309-w.
}
\author{
Joaquim Radua and Aleix Solanes
Maintainer: Aleix Solanes <email: solanes@recerca.clinic.cat>
}

\examples{
\dontrun{

# create an object of class mripredict
mp <- mripredict(mri_data = mri_train, data_table_file = clinical_train,
                 response_var="female", covariates=c("cent_age"), 
                 predictor=c("score1"),# "clinical_scale1"),
                 response_family="binomial", modulation="fu")
# fit a model
mp <- mripredict_fit(mp, space = "MNI", save_name = "model_fit")
# apply previous model to test data
preds <- mripredict_predict(mp, mri_data = mri_test, 
data_table_file = clinical_test)

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
