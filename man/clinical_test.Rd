\name{clinical_test}
\alias{clinical_test}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
clinical testing demo data for MRIPredict
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Sample data only for demonstration purposes
}
\usage{data("clinical_test")}
\format{
  A data frame with 10 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{female}}{a character vector}
    \item{\code{cent_age}}{a numeric vector}
    \item{\code{score0}}{a numeric vector}
    \item{\code{score1}}{a numeric vector}
    \item{\code{score2}}{a numeric vector}
    \item{\code{site}}{a character vector}
    \item{\code{clinical_scale1}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{status}}{a numeric vector}
  }
}

\examples{
data(clinical_test)
## maybe str(clinical_test) ; plot(clinical_test) ...
}
\keyword{datasets}
