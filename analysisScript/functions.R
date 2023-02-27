calculateQCMetric <- function(AllSampleSGLevel, negativeControls, panLethals ) {
  
  
  zPrimeScore <-
    data.frame(
      SampleName = c(),
      SampleID = c(),
      ZPrimeScore_treat = c(),
      ZPrimeScore_ctrl = c(),
      ttest_treat = c(),
      ttest_ctrl = c(),
      ZPrimeScore_treat_sg = c(),
      ZPrimeScore_ctrl_sg = c(),
      ttest_treat_sg = c(),
      ttest_ctrl_sg = c()
    )
  set = 0
  samples <- unique(AllSampleSGLevel$SAMPLEID)
  
  for (s in samples) {
    #s <- unique(SampleSGLevel$SAMPLEID)
    SampleSGLevel <- AllSampleSGLevel %>%
                      filter(SAMPLEID == s)
    
    if (!is.na(s)) {
      sname <- unique(SampleSGLevel[, "SAMPLE_NAME"])
      neg_data <-
        SampleSGLevel %>% 
        filter(Gene %in% negativeControls)
        #filter(!(Gene %in% unique(essential_genes_selected)))

      pos_data <- SampleSGLevel %>% 
        filter(Gene %in% panLethals)
        #filter((Gene %in% unique(essential_genes_selected)))

      neg_sum_treat <-
        neg_data %>%
        group_by(Gene) %>%
        dplyr::summarize(medianGene = median(X_PlasmidNorm.Treatment, na.rm=T)) 

      pos_sum_treat <-
        pos_data %>% 
        group_by(Gene) %>%
        dplyr::summarize(medianGene = median(X_PlasmidNorm.Treatment, na.rm=T))
      
      neg_sum_ctrl <-
        neg_data %>%
        group_by(Gene) %>%
        dplyr::summarize(medianGene = median(X_PlasmidNorm.Control, na.rm=T))
      
      pos_sum_ctrl <-
        pos_data %>% 
        group_by(Gene)%>%
        dplyr::summarize(medianGene = median(X_PlasmidNorm.Control, na.rm=T))

      neg_data_median_treat_sum <- median(neg_sum_treat$medianGene)
      neg_data_median_ctrl_sum <- median(neg_sum_ctrl$medianGene)
      pos_data_median_treat_sum <- median(pos_sum_treat$medianGene)
      pos_data_median_ctrl_sum <- median(pos_sum_ctrl$medianGene)
      
      
      neg_data_median_treat <-
        median(neg_data$X_PlasmidNorm.Treatment)
      neg_data_median_ctrl <- median(neg_data$X_PlasmidNorm.Control)
      pos_data_median_treat <-
        median(pos_data$X_PlasmidNorm.Treatment)
      pos_data_median_ctrl <- median(pos_data$X_PlasmidNorm.Control)
      
      #neg_data_sd_treat <- mad(neg_data$X_PlasmidNorm.Treatment, na.rm=T)
      #neg_data_sd_ctrl <- mad(neg_data$X_PlasmidNorm.Control, na.rm=T)
      #pos_data_sd_treat <- mad(pos_data$X_PlasmidNorm.Treatment, na.rm=T)
      #pos_data_sd_ctrl <- mad(pos_data$X_PlasmidNorm.Control, na.rm=T)
      
      #neg_data_sd_treat_sum <- mad(neg_sum_treat$medianGene, na.rm=T)
      #neg_data_sd_ctrl_sum <- mad(neg_sum_ctrl$medianGene, na.rm=T)
      #pos_data_sd_treat_sum <- mad(pos_sum_treat$medianGene, na.rm=T)
      #pos_data_sd_ctrl_sum <- mad(pos_sum_ctrl$medianGene, na.rm=T)
      
      neg_data_sd_treat <-
        sd(neg_data$X_PlasmidNorm.Treatment, na.rm = T)
      neg_data_sd_ctrl <-
        sd(neg_data$X_PlasmidNorm.Control, na.rm = T)
      pos_data_sd_treat <-
        sd(pos_data$X_PlasmidNorm.Treatment, na.rm = T)
      pos_data_sd_ctrl <-
        sd(pos_data$X_PlasmidNorm.Control, na.rm = T)
      
      neg_data_sd_treat_sum <- sd(neg_sum_treat$medianGene, na.rm = T)
      neg_data_sd_ctrl_sum <- sd(neg_sum_ctrl$medianGene, na.rm = T)
      pos_data_sd_treat_sum <- sd(pos_sum_treat$medianGene, na.rm = T)
      pos_data_sd_ctrl_sum <- sd(pos_sum_ctrl$medianGene, na.rm = T)
      
      zprime_treat <-
        1 - (((
          pos_data_sd_treat + neg_data_sd_treat
        )) / abs((
          pos_data_median_treat - neg_data_median_treat
        )))
      zprime_ctrl <-
        1 - (((
          pos_data_sd_ctrl + neg_data_sd_ctrl
        )) / abs((
          pos_data_median_ctrl - neg_data_median_ctrl
        )))
      
  
      ttest_treat <-
         t.test(
           neg_data$X_PlasmidNorm.Treatment,
           pos_data$X_PlasmidNorm.Treatment,
           alternative = "less"
         )
       ttest_ctrl <-
         t.test(
           neg_data$X_PlasmidNorm.Control,
           pos_data$X_PlasmidNorm.Control,
           alternative = "less"
         )

      ssmd_treat <-
        (pos_data_median_treat - neg_data_median_treat) / (1.4826 * (sqrt(
          pos_data_sd_treat ^ 2 + neg_data_sd_treat ^ 2
        )))
      ssmd_ctrl <-
        (pos_data_median_ctrl - neg_data_median_ctrl) / (1.4826 * (sqrt(
          pos_data_sd_ctrl ^ 2 + neg_data_sd_ctrl ^ 2
        )))


      zprime_treat_sum <-
        1 - (((
          pos_data_sd_treat_sum + neg_data_sd_treat_sum
        )) / abs((pos_data_median_treat_sum - neg_data_median_treat_sum)
        ))
      zprime_ctrl_sum <-
        1 - (((
          pos_data_sd_ctrl_sum + neg_data_sd_ctrl_sum
        )) / abs((
          pos_data_median_ctrl_sum - neg_data_median_ctrl_sum
        )))

       ttest_treat_sum <-
         t.test(neg_sum_treat$medianGene,
                pos_sum_treat$medianGene,
                alternative = "less")
       ttest_ctrl_sum <-
         t.test(neg_sum_ctrl$medianGene,
                pos_sum_ctrl$medianGene,
                alternative = "less")
       
      ssmd_treat_sum <-
        (pos_data_median_treat_sum - neg_data_median_treat_sum) / (1.4826 * (
          sqrt(pos_data_sd_treat_sum ^ 2 + neg_data_sd_treat_sum ^ 2)
        ))
      ssmd_ctrl_sum <-
        (pos_data_median_ctrl_sum - neg_data_median_ctrl_sum) / (1.4826 * (sqrt(
          pos_data_sd_ctrl_sum ^ 2 + neg_data_sd_ctrl_sum ^ 2
        )))

      if (set == 0) {
        zPrimeScore <-
          data.frame(
            SampleName = sname,
            SampleID = s,
            ZPrimeScore_treat_sg = zprime_treat,
            ZPrimeScore_ctrl_sg = zprime_ctrl,
            ttest_treat_sg = ttest_treat$statistic,
            ttest_ctrl_sg = ttest_ctrl$statistic,
            ssmd_treat_sg = ssmd_treat,
            ssmd_ctrl_sg = ssmd_ctrl,
            ZPrimeScore_treat = zprime_treat_sum,
            ZPrimeScore_ctrl = zprime_ctrl_sum,
            ttest_treat = ttest_treat_sum$statistic,
            ttest_ctrl = ttest_ctrl_sum$statistic,
            ssmd_treat = ssmd_treat_sum,
            ssmd_ctrl = ssmd_ctrl_sum
          )
        set = 1
      } else{
        zPrimeScore <-
          rbind(
            zPrimeScore,
            data.frame(
              SampleName = sname,
              SampleID = s,
              ZPrimeScore_treat_sg = zprime_treat,
              ZPrimeScore_ctrl_sg = zprime_ctrl,
              ttest_treat_sg = ttest_treat$statistic,
              ttest_ctrl_sg = ttest_ctrl$statistic,
              ssmd_treat_sg = ssmd_treat,
              ssmd_ctrl_sg = ssmd_ctrl,
              ZPrimeScore_treat = zprime_treat_sum,
              ZPrimeScore_ctrl = zprime_ctrl_sum,
              ttest_treat = ttest_treat_sum$statistic,
              ttest_ctrl = ttest_ctrl_sum$statistic,
              ssmd_treat = ssmd_treat_sum,
              ssmd_ctrl = ssmd_ctrl_sum
            )
          )
      }
    }
  }

  return(zPrimeScore)
}


qcdataclean <- function(qcdatain){
  QC <- pivot_longer(qcdatain, !starts_with("Sample"), names_to  = "QCMeasure", values_to="QCValue")
  
  
  qcdataret <- QC[which(QC$QCMeasure%in%c("ssmd_treat", "ssmd_ctrl")),]
  
  
  
  
  
  p <- ggplot(qcdataret,aes(x=QCValue))+
    geom_histogram(binwidth=0.1)+facet_grid(QCMeasure~.) +
    scale_x_continuous(limits = c(min(qcdataret$QCValue)-0.5, 0))
  dat <-  ggplot_build(p)$data[[1]]
  dat$SampleType <- ifelse(dat$PANEL==1, "Ctrl", "Treat")
  
  
  
  
  return(list(qcdataret = qcdataret, dat = dat))
}


##
#  function to calcuate the correlation of a large matrix
#
##

bigcor <- function(x, nblocks = 10, verbose = TRUE, ...)
{
  library(ff, quietly = TRUE)
  NCOL <- ncol(x)
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(MAT[, G1], MAT[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  gc()
  return(corMAT)
}