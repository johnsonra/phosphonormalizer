context("phosphonormalizer function")

test_that("Basic tests", {
    samplesCols <- data.frame(enriched=3:17, non.enriched=3:17)
    modseqCols <- data.frame(enriched = 1:2, non.enriched = 1:2)
    techRep <- factor(x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
    phonor.out <- normalizePhospho(enriched = enriched.rd, non.enriched = non.enriched.rd, 
                             samplesCols = samplesCols, modseqCols = modseqCols, techRep = techRep)
    expect_that(phonor.out, is_a("data.frame"))
    expect_that(ncol(phonor.out), equals(ncol(enriched.rd)))
    expect_that(nrow(phonor.out), equals(nrow(enriched.rd)))
})

