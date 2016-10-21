\name{PPN.cores}
\alias{PPN.cores}
\title{PPN cores}
\description{Set mc.cores option from an environment variable.}
\usage{PPN.cores(variable = "PBS_NUM_PPN")}
\arguments{
  \item{variable}{The PBS_NUM_PPN variable is defined by PBS,
for example it will be 4 when qsub -l nodes=1:ppn=4}
}

\value{The new value of options(mc.cores).}

\author{Toby Dylan Hocking}





