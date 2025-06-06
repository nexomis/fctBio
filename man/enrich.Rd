% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.r
\name{enrich}
\alias{enrich}
\title{Simple Gene Set Enrichment}
\usage{
enrich(
  gene_input,
  ann,
  gene_universe = NULL,
  lim_pmin = 0.05,
  classic = FALSE,
  log_level = "WARN",
  hard_pmin_filter = TRUE
)
}
\arguments{
\item{gene_input}{Input gene list converted as an INTEGER vector.}

\item{ann}{Modified annotation table as returned from \link{load_ann_space}
but with}

\item{gene_universe}{Optional, Universe gene list as an INTEGER vector.}

\item{lim_pmin}{Minimum value for pmin (lowest possible p-value).
This limit is used for filtering before p-value correction.}

\item{classic}{Boolean, if TRUE performs analysis relative to the universe
(vs. parent-child method).}

\item{log_level}{Logging level (see logging package). Default is WARN.}

\item{hard_pmin_filter}{Boolean, if TRUE, terms with pmin > lim_pmin are
removed from outputs.}
}
\description{
This function performs enrichment analysis using topology-based annotations.
It applies the parent-child method for over-representation analysis.
Multi-testing correction is performed using BH and Bonferroni
methods for terms with a minimal p-value below 1e-4 by default. P-values
are calculated using the hypergeometric test. This function is intended
for internal use due to its reliance on integer IDs for genes and terms;
later a new implementation may allow direct use of this function.
}
