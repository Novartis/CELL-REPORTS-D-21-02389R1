
#' ---
#' title:  QC analysis of screen of BRAF mutant cell lines
#' author: Grainne Kerr
#' output:
#'    html_document:
#'        toc: true
#'        toc_depth: 2
#'    md_document:
#'        variant: gfm
#' always_allow_html: true
#' ---

# The above header specifies that both an html document and a markdown
# document should be rendered from this R script. The markdown file can be
# viewed natively within Bitbucket.
# Turn off messages and warning when generating html / markdown files,
# and save plots to a plots directory (located in the working directory)
knitr::opts_chunk$set(message = FALSE,
                      warning = FALSE,
                      fig.path = "plots/braf/",
                      echo=FALSE)

# The following line renders both the html and the md document.
# run in console outside script
# rmarkdown::render('main_brafmut.R', output_format = 'html_document', clean = TRUE)

# Set the working directory to the top-level directory of this repo
setwd(here::here())



library(tidyverse, quietly = TRUE)
library(plyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(hrbrthemes, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(UpSetR, quietly = TRUE)
library(GGally, quietly=TRUE)
library(cowplot, quietly = TRUE)

# Currently the following two lines are necessary to load the gdsx library
#dyn.load('/usr/prog/HDF5/1.10.4-foss-2018a/lib/libhdf5.so')
#dyn.load('/usr/prog/HDF5/1.10.4-foss-2018a/lib/libhdf5_hl.so.100')
#library(gdsx, quietly = TRUE)

source("functions.R")


#' 
#' # Summary of Analysis Method
#' 
#' For each sample (and corresponding plasmid library), barcodes targeting each sgRNA were counted and normalized using a TMM normalization method available using the edgeR Bioconductor package [1]. The log fold change of sgRNA abundance compared to plasmid library was calculated using the general linear model log likelihood ratio test (glmLRT) method in edgeR. This corresponds to the log fold change 10 days. The log fold change of sgRNA abundance after 10 days compound treatment compared to DMSO control was calculated, also using a glmLRT method. This log fold change corresponds to a measure of either enhancement of drug response or activation cell growth after compound treatment (the sgRNA score). A gene score for each library pool was taken as the median log fold change of the >=10 sgRNAs targeting it. The significance of the gene knockdown synergism or activation was also assessed using the redundant siRNA concept [2]. This is a statistical score that models the probability of a gene 'hit' based on the activities of multiple sgRNAs per gene. Briefly, all sgRNAs are ranked according to their log fold change signal (treatment vs. control). Then, for each gene, the rank distribution of sgRNAs targeting it is examined and a p-value is assigned, based on an iterative hypergeometric function. The –log10  of this p-value indicates the statistical significance of all sgRNAs targeting a single gene being unusually distributed toward the top ranking slots. A separate calculation is made for drug synergistic targets (unusually distributed to the left of the distribution) and activators (unusually distributed to the right of the distribution) [2]. 
#' 
#' ![Fig 1. Analysis workflow](../meta_data/workflow.png)
#'
#'
#' # Screen Quality Assessment
#'

# Some global variables to extract the data
OMIT_CELL_LINES <- c("OUMS-23")
INVESTIGATION_NUM <- 200

#' ## Pan lethal and negative control genes
#' A large number of sgRNAs targeting essential positive genes and negative control genes were included in each pool. These genes were selected from the screen across CCLE [3].
#' If the median RSA score was < -7 across all CCLE then it was selected as a pan lethal.
#' if the RSA score was > -1.5 and it showed a normal distribution then it was selected as a negative control.
#' This resulted in 378 pan lethal genes and 1221 negative controls.

# Read in the list of essential genes and control genes

panLethals <- read.table("~/Projects//cEPIC/data/lethals.txt", header=F)
negativeControls <- read.table("~/Projects//cEPIC/data/negativeControls.txt", header=F)




# Read the data for cEPIC investigation meta info from the webservice?
X <-
  read.delim("http://nrchbs-slp0083:8080/GDSx/rest/api/spotfire/cEPIC/getInvestigations")
X <-
  X %>% filter(INVESTIGATION_ID == INVESTIGATION_NUM) %>% unique()
X <- X %>% filter(!is.na(CONTROL_SAMPLE_ID))
X <- X[-grep(OMIT_CELL_LINES, X$SAMPLE_NAME), ]


# Need the sample and control sample ids to retrieve the gene knock out scores
samp.list <- X %>% select(SAMPLE_ID) %>% pull(SAMPLE_ID)

ctrl.samp.list <-
  X %>% select(CONTROL_SAMPLE_ID) %>% pull(CONTROL_SAMPLE_ID)

selectedSamplesAndControls <-
  paste(samp.list, ctrl.samp.list, sep = "_")

# Now retrieve the read-out data from the screens dropout values using the ids
# Retrieving in a for loop, easier to catch a file/sample pairing that does not exist
logFC = lapply(selectedSamplesAndControls, function(x) {
  out <- tryCatch({
    url = paste0(
      "http://nrchbs-slp0083:8080/GDSx/rest/api/spotfire/cEPIC/getInvestigationLogFC/",
      x
    )
    read.delim2(url)
  },
  error = function(cond) {
    message(paste("URL does not seem to exist:", url))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning = function(cond) {
    message(paste("URL caused a warning:", url))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)
  })
  out
})

# Change the list to a dataframe and change the type to a double
logFC = rbind.fill(logFC)
logFC %<>% dplyr::mutate(X_PlasmidNorm.Control = as.double(X_PlasmidNorm.Control))
logFC %<>% dplyr::mutate(X_PlasmidNorm.Treatment = as.double(X_PlasmidNorm.Treatment))
logFC %<>% dplyr::mutate(XTreatment_Control_LogFC = as.double(XTreatment_Control_LogFC))
logFC %<>% dplyr::mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1"))
logFC %<>% dplyr::mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1"))
logFC %<>% dplyr::mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*) EF1a.*", "\\1"))
logFC %<>% dplyr::mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)-CMV-.*", "\\1"))
logFC %<>% dplyr::mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)-cas9.*", "\\1"))
logFC %<>% dplyr::mutate(CLEANNAME = str_replace(tolower(BIOSPECIMEN), "-", ""))

# Extract names of biospecimens
biospecimen_names <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  pull(BIOSPECIMEN) %>%
  unique()
biospecimen_cleannames <- str_replace(tolower(biospecimen_names), "-", "")

# New retrieve the rsa information - this may be unnecessary.
rsaLogFC = lapply(selectedSamplesAndControls, function(x) {
  out <- tryCatch({
    url = paste0(
      "http://nrchbs-slp0083:8080/GDSx/rest/api/spotfire/cEPIC/getInvestigationRSALogFC/",
      x
    )
    read.delim2(url)
  },
  error = function(cond) {
    message(paste("URL does not seem to exist:", url))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning = function(cond) {
    message(paste("URL caused a warning:", url))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)
  })
  out
})
rsaLogFC = rbind.fill(rsaLogFC)


# The control genes and the dropout score in each of the samples and the corresponding non-essential gene dropout
summarizedLogFC_ctrl = logFC  %>%
  filter(Gene %in% c(panLethals[[1]], negativeControls[[1]])) %>%
  mutate(CONTROL_GENE = ifelse(
    Gene %in% panLethals[[1]],
    "PAN.LETHAL",
    "NEG.CONTROL")
  ) %>%
  group_by(SAMPLE_NAME, CONTROL_GENE) %>%
  dplyr::summarise(
    X_PlasmidNorm.Control.median = median(X_PlasmidNorm.Control, na.rm = TRUE),
    X_PlasmidNorm.Treatment.median = median(X_PlasmidNorm.Treatment, na.rm = TRUE),
    BIOSPECIMEN_NAME = unique(BIOSPECIMEN_NAME)
  ) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, pattern = "-Cas9-.*", replacement = "")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, pattern = "_EF1aL*", replacement = "")) %>%
  pivot_longer(
    c(
      X_PlasmidNorm.Control.median,
      X_PlasmidNorm.Treatment.median
    ),
    names_to = "CONDITION",
    values_to = "MED.GENE.SCORE"
  ) %>%
  mutate(CONDITION = str_replace(CONDITION, "X_PlasmidNorm.", "")) %>%
  mutate(CONDITION = str_replace(CONDITION, ".median", "")) %>%
  pivot_wider(names_from = CONTROL_GENE,
              values_from = c(MED.GENE.SCORE))

#' The summary table below shows the median score for control genes in each sample. As expected, the pan lethal genes have a large negative score and the negative control genes have a score ~0. This illustrates that the essential lethal genes are causing cell death after ten days growth in both the treated and the untreated arms.
summarizedLogFC_ctrl  %>% ungroup() %>%
  mutate(POOL = str_replace(SAMPLE_NAME, ".* EPIC(Pool[123]) .*", "\\1")) %>%
  select(!SAMPLE_NAME) %>%
  select(BIOSPECIMEN, POOL, CONDITION, NEG.CONTROL, PAN.LETHAL) %>%
  kableExtra::kable() %>%
  kableExtra::kable_classic()




#'
#' The figure below illustrates the "window" for identifying hits observed by looking at the control gene information of pan lethal (green) vs negative control genes(red/brown). Note that TMM normalization was used in the analysis pipeline.
#' With the exception of MDST8, the dropout in the treated and control arms are similar range.
#'

#+ fig.height = 6, fig.width = 10, label = 'controlGEneDropout2'
summarizedLogFC_ctrl %>%
  ggplot() +
  geom_segment(aes(
    x = BIOSPECIMEN,
    xend = BIOSPECIMEN,
    y = PAN.LETHAL,
    yend = NEG.CONTROL
  ),
  color = "grey") +
  geom_point(aes(
    x = BIOSPECIMEN,
    y = PAN.LETHAL,
    color = rgb(0.2, 0.7, 0.1, 0.5)
  ), size = 3) +
  geom_point(aes(
    x = BIOSPECIMEN,
    y = NEG.CONTROL,
    color = rgb(0.7, 0.2, 0.1, 0.5)
  ), size = 3) +
  facet_wrap( ~ CONDITION, ncol = 5) +
  scale_color_manual(
    name = "",
    values = c(rgb(0.2, 0.7, 0.1, 0.5), rgb(0.7, 0.2, 0.1, 0.5)) ,
    labels = c("PanLethal", "NegativeControl")
  ) +
  theme_ipsum() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust =  0
    ),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Drop out of control genes") +
  xlab("") +
  ylab("Gene Knockdown Score")


#' ## SSMD score
#' To quantify the effect size and hence potential to identify synergistic effects, the strictly standardized mean difference measure was used. This measure evaluates the difference between the positive and negative controls in the data, while taking into consideration the variablity of the data. Here, we accepted a score of -1 as inferring an adequate effect size (a more negative score infers a better separation). The table below shows the QC score for each sample using the positive and negative control genes. 

zPrimeScore_panLethalControls <- calculateQCMetric(logFC, negativeControls[[1]], panLethals[[1]])

plotDataQC_panLethalControls  <- qcdataclean(zPrimeScore_panLethalControls )


plotDataQC_panLethalControls[["qcdataret"]] %>%
  mutate(Pool = str_replace(SampleName, ".* EPIC(Pool[123]) .*", "\\1")) %>%
  mutate(SampleName = str_replace(SampleName, "(.*)_EF1aL.*", "\\1")) %>%
  mutate(SampleName = str_replace(SampleName, "(.*)-Cas9.*", "\\1"))%>%
  select(-SampleID) %>% select(SampleName, Pool, QCMeasure, QCValue)%>%
  arrange(QCMeasure, SampleName, Pool) %>%
  kableExtra::kable() %>%
  kableExtra::kable_classic()

#' This cumulative plot shows that all of the control arms have a SSMD score of -1 or better. In the treatment arm 13/15 have a score of -1 or better.
#' 
#+ fig.height = 6, fig.width = 10, label = 'cumulativePlotZPrimeScorePanLethals'


ggplot(data.table::setDT(plotDataQC_panLethalControls[["dat"]])[,y2:=cumsum(y),"PANEL"],aes(x=x)) +
  geom_bar(aes(y=y2, fill=SampleType ), stat="identity", position="dodge")+
  labs(title="Strictly standardized mean difference score", 
       x="SSMD", 
       y="Cumulative Count")+
  theme_ipsum() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust =  0
        ),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5))






#' ## Technical reproduction
#' We considered the guides' fold-change profiles (with respect to initial representation) along the two branches as pseudo replicates of the same experiment thus estimating pipeline robustness and experiment reproducibility. We focused the analysis to guides where we expected to see an effect irrespective of the treatment condition, i.e. the pan-lethal/essential genes. Including all gene guides, the majority of which show 0 log-fold change, provides no additional information to assess the reproducibility.
#' 
#' 
#' #### Correlation between control and treated arms
#' 
#' We examined the knockdown of the essential genes and observe a high correlation between pan-lethal genes in both the treated and control arms of the experiment (the effect of lethal genes should be similar in both arms), the lowest correlation observed in the MDST8 cell line.
#' 
#' 

#+ fig.height = 6, fig.width = 10, label = 'panLethalGeneCorrelation'
logFC  %>%
  filter(Gene %in% c(panLethals[[1]])) %>%
  mutate(CONTROL_GENE = ifelse(
    Gene %in% panLethals[[1]],
    "PAN.LETHAL",
    "NEG.CONTROL"))%>%
  mutate(POOL_NAME = gsub("EPIC", "", POOL_NAME)) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  group_by(Gene, BIOSPECIMEN, POOL_NAME) %>% 
  dplyr::summarize(treatment.plasmidnorm = median(X_PlasmidNorm.Treatment), control.plasmidnorm = median(X_PlasmidNorm.Control)) %>%
  pivot_wider(id_cols=c(Gene, BIOSPECIMEN, POOL_NAME), values_from=c(treatment.plasmidnorm, control.plasmidnorm)) %>%
  ungroup() %>% 
  select(-Gene) %>% 
  ggscatter(x="treatment.plasmidnorm_", y="control.plasmidnorm_", add="reg.line", conf.int=TRUE) +
  stat_cor(method = "pearson", size=4)+
  facet_grid(rows=vars(BIOSPECIMEN), cols=vars(POOL_NAME))+
  theme_ipsum() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust =  0
        ),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5))




#' #### Correlation between cell lines
#'
#' As shown below we see high correlation of knockdown of pan lethal genes in all experiment arms between cell lines. ('T.' = Treatment, 'C.' = Control)
#'
#'
#'

#+ fig.height = 6, fig.width = 10, label = 'corAcrossCellLines'
cor_pre_mat <- logFC %>%
  filter(Gene %in% c(panLethals[[1]])) %>%
  pivot_wider(id_cols=c(Gene, Gene_Pool, CRISPR_NAME), names_from=c(BIOSPECIMEN, POOL_NAME), values_from = c(X_PlasmidNorm.Treatment, X_PlasmidNorm.Control), values_fill = NA)
colnames(cor_pre_mat) <- str_replace(colnames(cor_pre_mat), "X_PlasmidNorm.", "")

cormat_p3 <- cor(cor_pre_mat %>% select(-Gene, -Gene_Pool, -CRISPR_NAME) %>% select(ends_with("Pool3")), use="complete.obs")
colnames(cormat_p3) <- str_replace(colnames(cormat_p3), "_EPICPool3", "")
rownames(cormat_p3) <- str_replace(rownames(cormat_p3), "_EPICPool3", "")
colnames(cormat_p3) <- str_replace(colnames(cormat_p3), "Control_", "C.")
rownames(cormat_p3) <- str_replace(rownames(cormat_p3), "Control_", "C.")
colnames(cormat_p3) <- str_replace(colnames(cormat_p3), "Treatment_", "T.")
rownames(cormat_p3) <- str_replace(rownames(cormat_p3), "Treatment_", "T.")
cormat_p3 <- cormat_p3[order(rownames(cormat_p3)),order(colnames(cormat_p3))]

cormat_p1 <- cor(cor_pre_mat %>% select(-Gene, -Gene_Pool, -CRISPR_NAME) %>% select(ends_with("Pool1")), use="complete.obs")
colnames(cormat_p1) <- str_replace(colnames(cormat_p1), "_EPICPool1", "")
rownames(cormat_p1) <- str_replace(rownames(cormat_p1), "_EPICPool1", "")
colnames(cormat_p1) <- str_replace(colnames(cormat_p1), "Control_", "C.")
rownames(cormat_p1) <- str_replace(rownames(cormat_p1), "Control_", "C.")
colnames(cormat_p1) <- str_replace(colnames(cormat_p1), "Treatment_", "T.")
rownames(cormat_p1) <- str_replace(rownames(cormat_p1), "Treatment_", "T.")
cormat_p1 <- cormat_p1[order(rownames(cormat_p1)),order(colnames(cormat_p1))]


cormat_p2 <- cor(cor_pre_mat %>% select(-Gene, -Gene_Pool, -CRISPR_NAME) %>% select(ends_with("Pool2")), use="complete.obs", method="pearson")
colnames(cormat_p2) <- str_replace(colnames(cormat_p2), "_EPICPool2", "")
rownames(cormat_p2) <- str_replace(rownames(cormat_p2), "_EPICPool2", "")
colnames(cormat_p2) <- str_replace(colnames(cormat_p2), "Control_", "C.")
rownames(cormat_p2) <- str_replace(rownames(cormat_p2), "Control_", "C.")
colnames(cormat_p2) <- str_replace(colnames(cormat_p2), "Treatment_", "T.")
rownames(cormat_p2) <- str_replace(rownames(cormat_p2), "Treatment_", "T.")
cormat_p2 <- cormat_p2[order(rownames(cormat_p2)),order(colnames(cormat_p2))]


p1 <- ggcorr(data=NULL, 
             cor_matrix = cormat_p1, geom = "blank", label = TRUE, label_round = 2) +  
  geom_point(aes(size = coefficient, color = coefficient > 0, alpha = ifelse(abs(coefficient) > 0.8, "HIGH", ifelse(abs(coefficient) > 0.5, "MED","LOW")))) +
  scale_size_continuous(range = c(5, 12)) + 
  scale_alpha_manual(values = c("HIGH" = 0.5, "MED"=0.25, "LOW" = 0)) +  
  guides(color = "none", alpha = "none")+ 
  theme(legend.position = "none")


p2 <- ggcorr(data=NULL, 
             cor_matrix = cormat_p2, geom = "blank", label = TRUE, label_round = 2, ) +  
  geom_point(aes(size = coefficient, color = coefficient > 0, alpha = ifelse(abs(coefficient) > 0.8, "HIGH", ifelse(abs(coefficient) > 0.5, "MED","LOW")))) +
  scale_size_continuous(range = c(5, 12)) + 
  scale_alpha_manual(values = c("HIGH" = 0.5, "MED"=0.25, "LOW" = 0)) +  
  guides(color = "none", alpha = "none")+ 
  theme(legend.position = "none")
 
p3 <- ggcorr(data=NULL, 
             cor_matrix = cormat_p3, geom = "blank", label = TRUE, label_round = 2) +  
  geom_point(aes(size = coefficient, color = coefficient > 0, alpha = ifelse(abs(coefficient) > 0.8, "HIGH", ifelse(abs(coefficient) > 0.5, "MED","LOW")))) +
  scale_size_continuous(range = c(5, 12)) + 
  scale_alpha_manual(values = c("HIGH" = 0.5, "MED"=0.25, "LOW" = 0)) +  
  guides(color = "none", alpha = "none")+ 
  theme(legend.position = "none") 

#+ fig.height = 10, fig.width = 16, label = 'cellLineCorrelation'
plot_grid(p1,p2, p3, ncol = 3, labels=c("A) Correlation Panlethal genes Pool1", "B) Correlation Panlethal genes Pool2", "C) Correlation Panlethal genes Pool3"))


cormat_p1 %>%
  kableExtra::kable(caption = "Correlation essential/pan lethal genes pool 1") %>%
  kableExtra::kable_classic()

cormat_p2 %>%
  kableExtra::kable(caption = "Correlation essential/pan lethal genes pool 2") %>%
  kableExtra::kable_classic()

cormat_p3 %>%
  kableExtra::kable(caption = "Correlation essential/pan lethal genes pool 3") %>%
  kableExtra::kable_classic()


#' # Synergistic Genes

#'
#'
#' The gene score distribution plots below show the distribution of the logFC essential/panlethal genes of the treated/Untreated. This corresponds the the treatment effect, and we would expect both to be ~0 (no difference in effect of pan lethal genes in the treated and untreated arm.) The distributions all center around 0, with the exception of the MDST8 cell line. In this cell line the pan-lethal genes are skewed to the right, indicating that the knock out of lethal genes are more effective in the untreated arm compared to the treated arm (“pseudo-rescuer”).  The dashed line indicates threshold used to select synergistic genes. 
#'

#+ fig.height = 6, fig.width = 10, label = 'logFCDistribution'
logFC%>%
  filter(Gene %in% c(panLethals[[1]], negativeControls[[1]])) %>%
  mutate(CONTROL_GENE = ifelse(
    Gene %in% panLethals[[1]],
    "PAN.LETHAL",
    "NEG.CONTROL"))%>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  group_by(Gene_Pool, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T), CONTROL_GENE = unique(CONTROL_GENE)) %>%
  ungroup() %>% 
  ggplot()+
  geom_density(aes(x=logFC, fill=CONTROL_GENE), alpha = 0.4) +
  geom_vline(aes(xintercept = -0.5), linetype = "dashed", colour="grey") +
  geom_text(aes(x=-0.5, label="-0.5", y=5), angle=90) +
  facet_grid(cols=vars(BIOSPECIMEN))+
  theme_ipsum() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust =  0
        ),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
 
# Create a list for the input to upsetR
listInput_rko <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  filter(BIOSPECIMEN == "RKO") %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  pull(Gene)

listInput_ht29 <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  filter(BIOSPECIMEN == "HT-29") %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  pull(Gene)

listInput_ls411n <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  filter(BIOSPECIMEN == "LS411N") %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  pull(Gene)

listInput_mdst8 <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  filter(BIOSPECIMEN == "MDST8") %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  pull(Gene)

listInput_snuc5 <- logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  filter(BIOSPECIMEN == "SNU-C5") %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  pull(Gene)


listInput <- list(RKO=unique(listInput_rko), HT29 = listInput_ht29, LS411N = listInput_ls411n, MDST8 = listInput_mdst8, SNUC5 = listInput_snuc5)

#'
#'  The upset plot shows the number synergistic genes in each cell line (set) (logFC less than or equal to - 0.5), and the number in each intersection thereof.
#' 
upset(fromList(listInput), order.by = "freq")




#' 
#' The genes which were synergistic in 2 or more cell lines.
#' 
  logFC  %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN_NAME, "(.*)-Cas9.*", "\\1")) %>%
  mutate(BIOSPECIMEN = str_replace(BIOSPECIMEN, "(.*)_EF1aL.*", "\\1")) %>%
  group_by(Gene, BIOSPECIMEN) %>% 
  dplyr::summarize(logFC = median(XTreatment_Control_LogFC, na.rm=T)) %>%
  filter(logFC < -0.5) %>%
  dplyr::count(Gene, sort=TRUE) %>%
    filter(n > 1) %>%
    kableExtra::kable() %>%
    kableExtra::kable_classic()
  
  res <- list(GENE=c())
  for(i in seq_along(listInput)){
    names(listInput[[i]]) <- listInput[[i]]
    res <- bind_rows(res, listInput[[i]])
  }
  rownames(res) <- names(listInput)
  kableExtra::kable(t(as.data.frame(res)), caption = "Gene hits per cell line", longtable = TRUE, row.names = FALSE) %>%
    kableExtra::kable_styling(font_size = 10, full_width = FALSE)
  
  
 


#' # References
#' 1. https://bioconductor.org/packages/release/bioc/html/edgeR.html McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), 4288-4297. doi: 10.1093/nar/gks042.
#' 2. https://doi.org/10.1038/nmeth1089 König, R., Chiang, Cy., Tu, B. et al. A probability-based approach for the analysis of large-scale RNAi screens. Nat Methods 4, 847–849 (2007). 
#' 3. https://doi.org/10.1016/j.cell.2017.07.005 "McDonald ER 3rd, de Weck A, Schlabach MR, Billy E, Mavrakis KJ, Hoffman GR, Belur D, Castelletti D, Frias E, Gampa K, Golji J, Kao I, Li L, Megel P, Perkins TA, Ramadan N, Ruddy DA, Silver SJ, Sovath S, Stump M, Weber O, Widmer R, Yu J, Yu K, Yue Y, Abramowski D, Ackley E, Barrett R, Berger J, Bernard JL, Billig R, Brachmann SM, Buxton F, Caothien R, Caushi JX, Chung FS, Cortés-Cros M, deBeaumont RS, Delaunay C, Desplat A, Duong W, Dwoske DA, Eldridge RS, Farsidjani A, Feng F, Feng J, Flemming D, Forrester W, Galli GG, Gao Z, Gauter F, Gibaja V, Haas K, Hattenberger M, Hood T, Hurov KE, Jagani Z, Jenal M, Johnson JA, Jones MD, Kapoor A, Korn J, Liu J, Liu Q, Liu S, Liu Y, Loo AT, Macchi KJ, Martin T, McAllister G, Meyer A, Mollé S, Pagliarini RA, Phadke T, Repko B, Schouwey T, Shanahan F, Shen Q, Stamm C, Stephan C, Stucke VM, Tiedt R, Varadarajan M, Venkatesan K, Vitari AC, Wallroth M, Weiler J, Zhang J, Mickanin C, Myer VE, Porter JA, Lai A, Bitter H, Lees E, Keen N, Kauffmann A, Stegmeier F, Hofmann F, Schmelzle T, Sellers WR. Project DRIVE A Compendium of Cancer Dependencies and Synthetic Lethal Relationships Uncovered by Large-Scale, Deep RNAi Screening. Cell. 2017 Jul 27;170(3):577-592.e10. doi: 10.1016/j.cell.2017.07.005. PMID: 28753431."
#' 4. https://doi.org/10.1038/nature23270. "Manguso, R., Pope, H., Zimmer, M. et al. In vivo CRISPR screening identifies Ptpn2 as a cancer immunotherapy target. Nature 547, 413–418 (2017). https://doi.org/10.1038/nature23270"
 
  
 

# Print the R and package versions used for this analysis
#' # R version and attached packages for analysis
sessionInfo()
