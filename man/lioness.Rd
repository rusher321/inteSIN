% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lioness.R
\name{lioness}
\alias{lioness}
\title{Linear Interpolation to Obtain Network Estimates for Single Samples}
\usage{
lioness(matrix, n_method = .pearsonF, transF = F)
}
\arguments{
\item{matrix, }{row is gene/species, col is smaple/cell}

\item{n_method, }{function correlation method for feature-feature}

\item{transF, }{list output or matrix}
}
\value{
specific sample network
}
\description{
Linear Interpolation to Obtain Network Estimates for Single Samples
}
\examples{
num_samples <- 5
num_genes <- 6
gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
x <- as.data.frame(gene_matrix )
rownames(x) <- paste0("gene_", 1:num_genes)
colnames(x) <- paste0("sample_", 1:num_samples)
lioness <- lioness(t(x))
}
\references{
Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019 Apr 26;14:226-240. doi: 10.1016/j.isci.2019.03.021 Add to Citavi project by DOI. Epub 2019 Mar 28. PMID: 30981959 Add to Citavi project by Pubmed ID; PMCID: PMC6463816.
}
