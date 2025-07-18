% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density_surv_analysis.r
\name{density_surv_analysis}
\alias{density_surv_analysis}
\title{Post-PIPA analysis - prognostic association analysis based on cell phenotype density}
\usage{
density_surv_analysis(
  surv_data = NULL,
  input_dir = NULL,
  output_dir = NULL,
  KM_xlab = "Survivals",
  density_fnm = NULL
)
}
\arguments{
\item{surv_data}{A data.frame: tumors in rows (tumor_ids in row.names), columns must include: "time", "cens";
if there are additional columns, they are treated as confounders}

\item{input_dir}{A file.path string of input directory storing the phenotype cell density data}

\item{output_dir}{A file.path string of output directory}

\item{KM_xlab}{A character string for labelling the x-axis in Kaplan-Meier curves}

\item{density_fnm}{A character string of input file name for cell density data}
}
\description{
Post-PIPA analysis - prognostic association analysis based on cell phenotype density
Patients/tumors are categorized into phenotype quartiles or zero+tertiles (if zero density cases >25% of total cases),
which are then used as (i) indicator variable in estimating Kaplan-Meier curves & log-rank test,
or (ii) as ordinal variable in Cox PH analysis
}
