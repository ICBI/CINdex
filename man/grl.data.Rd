\docType{data}
\name{grl.data}
\alias{grl.data}
\title{Output of segmentation algorithm}
\format{A GRangesList}
\usage{
data(grl.data)
}
\value{An example output of segmentation algorithm}
\description{
To mathematically and quantitatively describe these alternations we first locate their genomic positions and measure their ranges. Such algorithms are referred to as segmentation algorithms. Bioconductor has several copy number segmentation algorithms. There are many copy number segmentation algorithms outside of Bioconductor as well, examples are Fused Margin Regression (FMR) and Circular Binary Segmentation (CBS). 

Segmentation results are typically have information about the start position and end position in the genome, and the segment value. The algorithms typically covers chromosomes 1 to 22 without any gaps, sometimes sex chromosomes are also included.

For more details refer tutorial in the accompanying vignette in the CINdex package
}
\keyword{datasets}
