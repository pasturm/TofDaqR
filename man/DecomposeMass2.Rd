% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DecomposeMass2}
\alias{DecomposeMass2}
\title{Calculates possible sum formulas from a list of fragment masses.}
\usage{
DecomposeMass2(
  targetMass,
  tolerance,
  atomMass,
  atomLabel,
  elementIndex1,
  elementIndex2,
  filterMinVal,
  filterMaxVal,
  maxHits,
  maxSearch
)
}
\arguments{
\item{targetMass}{Target mass.}

\item{tolerance}{Tolerance of target mass.}

\item{atomMass}{Numeric vector with fragment masses.}

\item{atomLabel}{String Vector of fragment labels.}

\item{elementIndex1}{Element index for count filters, first element index for ratio filters.}

\item{elementIndex2}{-1 for count filters, second element for ratio filters.}

\item{filterMinVal}{Counts or ratios that are smaller than this value are filtered out.}

\item{filterMaxVal}{Counts or ratios that are larger than this value are filtered out.}

\item{maxHits}{Maximum number of candidate hits to report.}

\item{maxSearch}{Maximum number of candidate formula to search.}
}
\value{
List with sum formula strings and the mass and mass errors of the sum
formulas.
}
\description{
\code{DecomposeMass2} calculates possible sum formulas from a list of fragment
masses.
}
\details{
Same as \code{\link{DecomposeMass}} but with additional parameters \code{maxHits}
and \code{maxSearch} to better fine-tune when to abort a compomer search.
Calculates possible sum formulas that amount to target mass +/- tolerance
given a target mass, mass tolerance and a list of fragment masses. Filters
specifying min/max values for absolute counts of a fragment or for ratios
between fragments can be specified in order to reduce the amount of results
and restrict hits to chemically reasonable sum formulae. Typically atomMass
and atomLabels are the masses/labels of elements but you are free to use
whatever you like (isotopes, amino acids, common fragments etc.).
}
\examples{
targetMass = 314
tolerance = 0.5
atomLabel = c("C", "H", "O")
n = length(atomLabel)
atomMass = rep(0, n)
for (i in 1:n) {
 atomMass[i] = GetMoleculeMass(atomLabel[i])
}
elementIndex1 = seq(along.with = atomLabel)-1
elementIndex2 = rep(-1, n)
filterMinVal = c(20, 20, 0)
filterMaxVal = c(22, 40, 5)
maxHits = 10
maxSearch = 1000
DecomposeMass2(targetMass, tolerance, atomMass, atomLabel, elementIndex1,
              elementIndex2, filterMinVal, filterMaxVal, maxHits, maxSearch)
}
\seealso{
Other Chemistry functions: 
\code{\link{DecomposeMass}()},
\code{\link{GetIsotopePattern2}()},
\code{\link{GetIsotopePattern}()},
\code{\link{GetMoleculeMass}()}
}
\concept{Chemistry functions}
