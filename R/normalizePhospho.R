# Pairwise Normalization of MS-based phosphoproteomic data
# Sohrab Saraei
# This function compensates for the bias introduced in global 
# phosphorylation in the sample after using median normalization.

normalizePhospho <- function(enriched, non.enriched, phospho = NULL, techRep)
{
    
    stopifnot(!missing(x = enriched))
    stopifnot(!missing(x = non.enriched))
    stopifnot(!missing(x = techRep))
    stopifnot(class(enriched) == "data.frame")
    stopifnot(class(non.enriched) == "data.frame")
    stopifnot(class(phospho) == "character" | phospho == NULL)
    stopifnot(class(techRep) == "factor" & (length(techRep) == unique(ncol(enriched), ncol(non.enriched)) - 2))
    
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