\name{mri}
\alias{mri}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
MRI demo images for MRIPredict
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Sample data only for demonstration purposes
}
\usage{data("mri")}
\format{
  The format is:
List of 2
 $ unmodulated:List of 8
  ..$ dim       : int [1:3] 23 28 23
  ..$ sto.xyz   : num [1:4, 1:4] -8 0 0 0 0 8 0 0 0 0 ...
  ..$ scl.slope : num 1
  ..$ scl.inter : num 0
  ..$ sform.code: chr "NIFTI.XFORM.MNI.152"
  ..$ data      : num [1:23, 1:28, 1:23, 1:120] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ sto.ijk   : num [1:4, 1:4] -0.125 0 0 0 0 0.125 0 0 0 0 ...
  ..$ n         : int 120
  ..- attr(*, "class")= chr "mripredict_data"
 $ modulated  :List of 8
  ..$ dim       : int [1:3] 23 28 23
  ..$ sto.xyz   : num [1:4, 1:4] -8 0 0 0 0 8 0 0 0 0 ...
  ..$ scl.slope : num 1
  ..$ scl.inter : num 0
  ..$ sform.code: chr "NIFTI.XFORM.MNI.152"
  ..$ data      : num [1:23, 1:28, 1:23, 1:120] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ sto.ijk   : num [1:4, 1:4] -0.125 0 0 0 0 0.125 0 0 0 0 ...
  ..$ n         : int 120
  ..- attr(*, "class")= chr "mripredict_data"
}

\examples{
data(mri)
## maybe str(mri) ; plot(mri) ...
}
\keyword{datasets}
