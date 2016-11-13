#' Non-enriched dataset
#'
#' A dataset containing sequences, modifications and abundances of about 17000 peptides 
#' measured over 5 samples with 3 technical replicates each. 
#' 
#'
#' @format A data frame with 16982 rows and 17 variables, all samples are median normalized:
#' \describe{
#'    \item{Sequence}{The sequence of the peptide}
#'    \item{Modification}{The modification and its location}
#'    \item{gcNorm.ctrl2.1}{Sample: Control 2 Technical Replicate: 1}
#'    \item{gcNorm.ctrl2.2}{Sample: Control 2 Technical Replicate: 2}
#'    \item{gcNorm.ctrl2.3}{Sample: Control 2 Technical Replicate: 3}
#'    \item{gcNorm.ctrl1.1}{Sample: Control 1 Technical Replicate: 1}
#'    \item{gcNorm.ctrl1.2}{Sample: Control 1 Technical Replicate: 2}
#'    \item{gcNorm.ctrl1.3}{Sample: Control 1 Technical Replicate: 3}
#'    \item{gcNorm.CIP2A.1}{Sample: CIP2A Technical Replicate: 1}
#'    \item{gcNorm.CIP2A.2}{Sample: CIP2A Technical Replicate: 2}
#'    \item{gcNorm.CIP2A.3}{Sample: CIP2A Technical Replicate: 3}
#'    \item{gcNorm.RAS.1}{Sample: RAS Technical Replicate: 1}
#'    \item{gcNorm.RAS.2}{Sample: RAS Technical Replicate: 2}
#'    \item{gcNorm.RAS.3}{Sample: RAS Technical Replicate: 3}
#'    \item{gcNorm.OA.1}{Sample: OA Technical Replicate: 1}
#'    \item{gcNorm.OA.2}{Sample: OA Technical Replicate: 2}
#'    \item{gcNorm.OA.3}{Sample: OA Technical Replicate: 3}
#'    ...
#' }
#' @source \url{http://www.nature.com/articles/srep13099}
#' @return Example Non-enriched dataset
"non.enriched.rd"
