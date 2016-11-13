#' Pairwise Normalization of MS-based phosphoproteomic data
#'
#' @description  This function compensates for the bias introduced in global 
#' phosphorylation in the sample after using median normalization.
#' @param enriched The enriched dataframe with the sequence and modifications 
#' as the first and the second column in format of charachter, followed by the 
#' samples (and possibly their technical replicates) in numeric. 
#' Note that, exactly the same column order in the non-enriched dataset 
#' must be considered. The dataset is expected to be median normalized.
#' @param non.enriched The non-enriched dataframe with the sequence and modifications 
#' as the first and the second column in format of charachter, followed by the samples 
#' (and possibly their technical replicates) in numeric. Note that, exactly the same 
#' column order in the enriched dataset must be considered. The dataset is expected to 
#' be median normalized.
#' @param phospho a string that shows the term that represents phosphorylation in the 
#' modification column of the data. If it is not assigned, "Phospho" will be used as 
#' the default value
#' @param techRep a factor that holds information about columns order and the technical
#'  replicates of the samples
#' @details It is shown that global median normalization can introduce bias in the fold change of global phosphorylation between samples. 
#' It is suggested that by taking the non-enriched data into consideration, this bias could be compensated(Kauko et al. 2015).
#' @keywords Phosphoproteomics, Normalization, Mass-spectrometry
#' @import plyr
#' @importFrom stats median
#' @importFrom graphics boxplot
#' @return A dataframe with the normalize values
#' @export
#' @references \url{http://www.nature.com/articles/srep13099}
#' @examples
#' #Load the library
#' library(phosphonormalizer)
#' #The samples and their technical replicates
#' techRep <- factor(x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
#' #Call the function
#' norm <- normalizePhospho(enriched.rd, non.enriched.rd, techRep = techRep)
#' head(norm)
#' 
normalizePhospho <- function(enriched, non.enriched, phospho = NULL, techRep)
{
    stopifnot(!missing(x = enriched))
    stopifnot(!missing(x = non.enriched))
    stopifnot(!missing(x = techRep))
    
    mod = seq = NULL
    
    if(!inherits(x = enriched[,1], what = "character") & 
        !inherits(x = enriched[,2], what = "character") & 
        !inherits(x = non.enriched[,1], what = "character") &
        !inherits(x = non.enriched[,2], what = "character"))
        stop("The first two columns must be the sequence and modification!")
    
    if(ncol(enriched) != ncol(non.enriched))
        stop("The number of columns in enriched and non.enriched is different!")
    
    enriched.original.mat <- as.matrix(enriched[, -c(1,2)])
    seqMod <- enriched[, c(1,2)]
    
    if(!all((apply(X = enriched[, -c(1,2)], MARGIN = 2, 
                    function(x) inherits(x = x, what = "numeric")))) &
       !all((apply(X = non.enriched[, -c(1,2)], MARGIN = 2, 
                    function(x) inherits(x = x, what = "numeric")))))
        stop(paste0("Columns: 3 to ", ncol(enriched), "must be of type numeric!"))
    
    enriched <- enriched[apply(X = enriched[,-c(1,2)], MARGIN = 1, 
                                function(x) all(x != 0)),]
    non.enriched <- non.enriched[apply(X = non.enriched[,-c(1,2)], MARGIN = 1, 
                                        function(x) all(x != 0)),]
    
    colnames(enriched)[1:2] <- c("seq", "mod")
    colnames(non.enriched)[1:2] <- c("seq", "mod")
    
    if(is.null(phospho))
    {
        enriched <- enriched[grepl(pattern = "Phospho", x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = "Phospho", x = non.enriched$mod),]
    } else {
        enriched <- enriched[grepl(pattern = phospho, x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = phospho, x = non.enriched$mod),]
    }
    
    enriched <- plyr::ddply(enriched, plyr::.(seq, mod),
                            function(df) apply(df[, -c(1,2)], 2, sum))
    non.enriched <- plyr::ddply(non.enriched, plyr::.(seq, mod),
                                function(df) apply(df[, -c(1,2)], 2, sum))
    
    enriched[,"modSeq"] <- paste(enriched$seq, enriched$mod,sep = ", ")
    non.enriched[,"modSeq"] <- paste(non.enriched$seq, non.enriched$mod,sep = ", ")
    
    inter <- intersect(non.enriched$modSeq, enriched$modSeq)
    stopifnot(length(inter) > 0)
    enriched.olp.idx <- which(enriched$modSeq %in% inter)
    non.enriched.olp.idx <- which(non.enriched$modSeq %in% inter)
    
    enriched.mat <- enriched[enriched.olp.idx, ]
    non.enriched.mat <- non.enriched[non.enriched.olp.idx, ]
    
    enriched.mat <- enriched.mat[order(enriched.mat$modSeq),]
    non.enriched.mat <- non.enriched.mat[order(non.enriched.mat$modSeq),]
    
    enriched.mat <- as.matrix(enriched.mat[, -c(1,2,ncol(enriched))])
    non.enriched.mat <- as.matrix(non.enriched.mat[, -c(1,2,ncol(enriched))])
    
    ratios <- non.enriched.mat/enriched.mat
    colnames(ratios) <- as.numeric(techRep)
    ratios.avg <- matrix(nrow = nrow(ratios), ncol = length(levels(techRep)))
    
    for (tr in levels(techRep)) {
        ratios.avg[,as.numeric(tr)] <- rowMeans(ratios[,colnames(ratios) == levels(techRep)[as.numeric(tr)]])
    }
    
    max.fc <- log2(ratios.avg[,1]) - log2(ratios.avg[,2])
    
    for (i in 2:(ncol(ratios.avg)-1)) {
        for (j in (i+1):(ncol(ratios.avg))) {
            max.fc <- as.numeric(
                unlist(
                    apply(cbind(max.fc, log2(ratios.avg[,i]) - log2(ratios.avg[,j])),
                          1, max),
                    use.names = FALSE)
            )
        }
    }
    
    boxp <- boxplot(max.fc, plot = FALSE)
    ratios <- ratios[!(max.fc >= max(boxp$stats)),]
    
    ratios <- log10(ratios)
    
    col.sub <- rowMeans(ratios)
    
    ratios.norm <- as.matrix(apply(ratios, 2, 
                                   function(x) x - col.sub))
    
    factors <- 10^(apply(ratios.norm,2,median))
    
    data.frame(seqMod, t(t(enriched.original.mat) * factors)) 
}