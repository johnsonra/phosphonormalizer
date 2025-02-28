####################################################################
## Pairwise Normalization of MS-based phosphoproteomic data       ##
## Sohrab Saraei                                                  ##
## This function compensates for the bias introduced in global    ##
## phosphorylation in the sample after using median normalization.##
####################################################################

normalizePhospho <- function(enriched, non.enriched, phospho = NULL, samplesCols, modseqCols, techRep, plot.fc=NULL)
{
    #Check if all the necessary arguments are present
    if(missing(enriched)) stop("The function parameter (enriched) is missing!")
    if(missing(non.enriched)) stop("The function parameter (non.enriched) is missing!")
    if(missing(samplesCols)) stop("The function parameter (samplesCols) is missing!")
    if(missing(modseqCols)) stop("The function parameter (modseqCols) is missing!")
    if(missing(techRep)) stop("The function parameter (techRep) is missing!")
    
    #A function that checks the types of the columns of a data.frame
    
    #Check the type of the arguments
    class.dfCols <- function(df, types)
    {
        all(apply(df, MARGIN = 2, function(x) class(x)[1]) %in% types)
    }
    if(!is(enriched, "data.frame") & !is(enriched, "MSnSet"))
        stop("The enriched parameter must be of either types of data.frame or MSnSet")
    
    if(!is(non.enriched, "data.frame") & !is(non.enriched, "MSnSet"))
        stop("The non.enriched parameter must be of either types of data.frame or MSnSet")
    
    if(is(enriched, "MSnSet"))
    {
        enriched.eset <- enriched
        enriched <- cbind(MSnbase::fData(enriched)[, modseqCols$enriched], MSnbase::exprs(enriched)[, samplesCols$enriched])
        non.enriched <- cbind(MSnbase::fData(non.enriched)[, modseqCols$non.enriched], MSnbase::exprs(non.enriched)[, samplesCols$non.enriched])
        modseqCols <- data.frame(enriched = c(1,2), non.enriched = c(1,2))
        samplesCols <- data.frame(enriched = 3:ncol(enriched), non.enriched = 3:ncol(non.enriched))
    }
    
    if(!is.null(phospho) & !methods::is(phospho, "character"))
        stop("The phospho parameter must be of type of character")
    
    if(!methods::is(samplesCols, "data.frame") |
       ncol(samplesCols) != 2 | !class.dfCols(samplesCols, c("integer","numeric")) |
       all(!colnames(samplesCols) %in% c("enriched", "non.enriched")))
        stop("The samplesCols must be data.frame with two columns, with the column names enriched and non.enriched
             , of type numeric or integer, which must contain the column number of samples that hold the abundances")
    
    if(!methods::is(modseqCols, "data.frame") |
       ncol(modseqCols) != 2 | !class.dfCols(modseqCols, c("integer","numeric")) | 
       all(!colnames(modseqCols) %in% c("enriched", "non.enriched")))
        stop("The modseqCols must be data.frame with two columns, with the names enriched and non.enriched
             , of type numeric or integer, which must contain the column number of samples that hold the sequence 
             and modifications of the peptides")
    
    if(any(!class.dfCols(enriched[, modseqCols$enriched], "character") ,
           !class.dfCols(non.enriched[, modseqCols$non.enriched], "character")))
        stop("The sequence and modification columns that is specified are not of type charachter!")
    
    if(any(!class.dfCols(enriched[, samplesCols$enriched], c("numeric", "integer")) ,
           !class.dfCols(non.enriched[, samplesCols$non.enriched], c("numeric", "integer"))))
        stop("The samples specified are not of type of numeric")
    
    mod = seq = NULL
    enriched.original.mat <- as.matrix(enriched[, samplesCols$enriched])
    seqMod <- enriched[, modseqCols$enriched]
    #Removing peptides with non-quantified values across the runs
    enriched <- enriched[apply(X = enriched[ ,samplesCols$enriched], MARGIN = 1, 
                                function(x) all(x != 0)),]
    non.enriched <- non.enriched[apply(X = non.enriched[,samplesCols$non.enriched], MARGIN = 1, 
                                        function(x) all(x != 0)),]
    
    colnames(enriched)[modseqCols$enriched] <- c("seq", "mod")
    colnames(non.enriched)[modseqCols$non.enriched] <- c("seq", "mod")
    
    if(is.null(phospho))
    {
        enriched <- enriched[grepl(pattern = "Phospho", x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = "Phospho", x = non.enriched$mod),]
    } else {
        enriched <- enriched[grepl(pattern = phospho, x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = phospho, x = non.enriched$mod),]
    }
    
    enriched <- plyr::ddply(enriched, plyr::.(seq, mod),
                            function(df) colSums(df[, samplesCols$enriched]))
    non.enriched <- plyr::ddply(non.enriched, plyr::.(seq, mod),
                                function(df) colSums(df[, samplesCols$non.enriched]))
    
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
    if(length(inter) > 1) {
        ratios <- non.enriched.mat/enriched.mat

        colnames(ratios) <- as.numeric(techRep)
        ratios.avg <- matrix(nrow = nrow(ratios), ncol = length(levels(techRep)))
        for (tr in levels(techRep)) {
            tr_num <- techRep[techRep == tr] |> unique() |> as.numeric()
            if(nrow(ratios.avg) == 1) {
                ratios.avg[,tr_num] <- mean(ratios[,colnames(ratios) == tr_num])
            } else {
                ratios.avg[,tr_num] <- rowMeans(ratios[,colnames(ratios) == tr_num])
            }
        }


        max.fc <- log2(ratios.avg[,1]) - log2(ratios.avg[,2])

        for (i in 2:(ncol(ratios.avg)-1)) {
            for (j in (i+1):(ncol(ratios.avg))) {
                max.fc <- matrixStats::rowMaxs(cbind(max.fc, log2(ratios.avg[,i]) - log2(ratios.avg[,j])))
            }
        }

        boxp <- boxplot(max.fc, plot = FALSE)
        ratios <- ratios[!(max.fc > max(boxp$stats)),]

        ratios <- log10(ratios)
        if(methods::is(ratios, "matrix") | methods::is(ratios, "data.frame")) {
            col.sub <- rowMeans(ratios)
        } else {
            col.sub <- mean(ratios)
        }

        ratios.norm <- ratios - col.sub

        if(methods::is(ratios, "matrix") | methods::is(ratios, "data.frame")) {
            factors <- 10^(matrixStats::colMedians(ratios.norm))
        } else {
            factors <- ratios.norm
        }
    } else {
        factors <- as.numeric(non.enriched.mat/enriched.mat)
    }

    enriched.normalized.mat <- t(t(enriched.original.mat) * factors)
    if(!is.null(plot.fc)) {
        for(i in plot.fc$control) {
            tr_i <- as.numeric(techRep) == i
            for(j in plot.fc$samples) {
                tr_j <- as.numeric(techRep) == j
                a.original <- rowMeans(log2(enriched.original.mat[,tr_i]+1),na.rm=TRUE)
                b.original <- rowMeans(log2(enriched.original.mat[,tr_j]+1),na.rm=TRUE)
                fc.original <- a.original - b.original
                a.normnalized <- rowMeans(log2(enriched.normalized.mat[,tr_i]+1),na.rm=TRUE)
                b.normnalized <- rowMeans(log2(enriched.normalized.mat[,tr_j]+1),na.rm=TRUE)
                fc.normnalized <- a.normnalized - b.normnalized
                boxplot(cbind(fc.original, fc.normnalized), range=1.5, outline=FALSE, main=paste0("Peptide log fold changes", " (sample ", j, " vs sample ", i, ")"), names=c("Median normalized","Pairwise normalized"))
                abline(h=0, lty=2)}
        }

    }
    message(paste0("The number of peptides in the intersect is: ", length(inter)))
    message(paste0(length(plot.fc$control) * length(plot.fc$samples), " plots generated. Browse through them."))
    data.frame(seqMod, t(t(enriched.original.mat) * factors))
}
