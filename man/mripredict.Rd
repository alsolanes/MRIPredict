\name{mripredict}
\alias{mripredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Create an object of type "mripredict". It contains the needed configuration for 
mripredict_cv, mripredict_fit and mripredict_predict.
}
\usage{
mripredict(mri_data = NULL, clinical_data, response_var,covariates = NULL, 
predictor = NULL, response_family, modulation, information_variables = NULL, 
mask_files_path = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{response_var}{
%%     ~~Describe \code{response_var} here~~
A string with the name of the column that contains the response variable per 
each subject in the clinical_data table. 
If the response_family is 'cox' the response_var should contain a list of two 
variables, one for 'time' and one for 'status'.
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
A list of variables to be treated as confounding variables.
}
  \item{predictor}{
%%     ~~Describe \code{predictor} here~~
A list of variables that are supposed to apport information to the prediction, 
treated as predictor variables.
}
  \item{response_family}{
%%     ~~Describe \code{response_family} here~~
It can be 'binomial' for binary outcomes, 'gaussian' for continuous outcomes, 
or 'cox' for survival risk-estimation models.
}
  \item{modulation}{
%%     ~~Describe \code{modulation} here~~
It can be 'un' for unmodulated data, 'fu' for fully modulated data, or 
'clinical' if there is only clinical data and no MRI images.
}
  \item{information_variables}{
%%     ~~Describe \code{information_variables} here~~
Which variables should be used for multiple imputation. By default this will be 
set to the list of covariates and predictor variables.
}
  \item{mask_files_path}{
%%     ~~Describe \code{mask_files_path} here~~
Path to a niftii image to be used as a mask.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An MRIPredict object that can be accessed like a common list. It can be passed 
as a parameter in other mripredict_* functions.
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

# create an object of class mripredict to perform binary outcome prediction 
# using clinical and unmodulated structural MRI images.
mp <- mripredict(mri_data = mri_train, clinical_data = clinical_train,
                 response_var="female", covariates=c("cent_age"), 
                 predictor=c("score1"),
                 response_family="binomial", modulation="un")
# 
# perform a 10-fold cross-validation
mp <- mripredict_cv(mp, space = "NO_CHECK", 
                    save_name = "testing", 
                    folds_file = NULL, 
                    n_cores = 1, 
                    n_folds = 5)
# Train on training data, and apply it to new data                    
# train a model
mp <- mripredict_fit(mp, space = "MNI", save_name = "model_fit")
# apply a model on new data
predictions <- mripredict_predict(mp, mri_data = mri_test, 
clinical_data = clinical_test)


# perform a survival cox proportional hazards analysis
mp <- mripredict(mri_data = mri_train, clinical_data = clinical,
                 response_var=c("time","status"), covariates=c("cent_age"), 
                 predictor=c("score1"),
                 response_family="cox", modulation="un")

mp <- mripredict_cv(mp, space = "NO_CHECK", 
                    save_name = "testing", 
                    folds_file = NULL, 
                    n_cores = 1, 
                    n_folds = 5)

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
