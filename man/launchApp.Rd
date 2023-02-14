\name{launchApp}
\alias{launchApp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
MRIPredict GUI shiny app
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
GUI launcher for MRIPredict package. An easy to use interface of the MRIPredict package.
}
\usage{
launchApp()
}
%- maybe also 'usage' for other objects documented here.

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Output files will be defined by the user in the interface.
}
\references{
%% ~put references to the literature/web site here ~
Solanes, A. (2022), et al. Combining MRI and clinical data to detect high relapse risk after the first episode of psychosis. Schizophrenia (Heidelb). 2022 Nov 17;8(1):100. doi: 10.1038/s41537-022-00309-w.
}
\author{
%%  ~~who you are~~
Joaquim Radua and Aleix Solanes
Maintainer: Aleix Solanes <email: solanes@recerca.clinic.cat>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# this will launch the interface. A shiny application.
\dontrun{
library(shiny)
launchApp()
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
}