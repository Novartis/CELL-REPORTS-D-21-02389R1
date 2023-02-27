
 
 avana.data <- load_data("CRISPR_depMAP", "1.1.1")
 
 avana.colset <- avana.data$colAnnotation %>% 
   filter(CLEANNAME %in% logFC$CLEANNAME)
 
 avana.rowset <- avana.data$rowAnnotation %>% 
   filter(gene_symbol %in% unique(logFC$Gene))
 
 avana.dataset <- avana.data$data[rownames(avana.rowset), rownames(avana.colset)]
 colnames(avana.dataset) <- avana.colset[colnames(avana.dataset),]$CCLE_Name
 
 avana.dataset.essential <- avana.dataset %>%
   as.data.frame() %>%
   filter(rownames(avana.dataset) %in% panLethals[[1]])
 
 tmp <- avana.colset %>% select(CCLE_Name, CLEANNAME )
 avana.dataset.essential.long <- avana.dataset.essential %>% 
   rownames_to_column(var="GENE") %>% 
   pivot_longer(!GENE) %>% 
   inner_join(tmp, by = c("name" = "CCLE_Name"))
 
 corAvanaData <- logFC %>%
   filter(Gene %in% c(panLethals[[1]])) %>%
   mutate(CONTROL_GENE = ifelse(
     Gene %in% panLethals[[1]],
     "PAN.LETHAL",
     "NEG.CONTROL")) %>%
   group_by(Gene, CLEANNAME) %>% 
   dplyr::summarize( control.plasmidnorm = median(X_PlasmidNorm.Control)) %>% 
   inner_join(avana.dataset.essential.long, by=c("CLEANNAME" = "CLEANNAME", "Gene" = "GENE"))
 