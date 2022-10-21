
############## Function for meta-analysis plot
mres_f = function(mts, range, nstuds, protcomp) {
  if (protcomp == "All") {
    nstuds = as.numeric(as.character(nstuds))
    sel1 = which(mts$I2<range/100)
    sel2 = which(mts$k.study>=nstuds)
    sel = intersect(sel1, sel2)
    meta_sel = mts[sel,]
  } else {
    meta_sel = mts[which(mts$gene_name %in% protcomp),]
  }
  
  meta_sel = meta_sel[order(meta_sel$TE.random, decreasing = T),]
  meta_sel$order = 1:nrow(meta_sel)
  
  error_y = list(
    type = "data",
    symmetric = FALSE,
    array = c(meta_sel$upper.random-meta_sel$TE.random),
    arrayminus = c(meta_sel$TE.random-meta_sel$lower.random),
    thickness=0.5,
    width=2)
  
  plot_ly(data=meta_sel, x=meta_sel$order, y=meta_sel$TE.random, 
               text = meta_sel$gene_name, color=meta_sel$I2*100, error_y=error_y,
               type="scatter", mode="markers", size=1, 
               colors=c("#56B4E9", "#0072B2", "#F0E442", "#E69F00", "#D55E00")) %>%
    layout(legend=list(title=list(text='Heterogeneity - I2 %')), 
           xaxis = list(title = list(text ='MA: Statistically-significant Proteins')),
           yaxis = list(title = list(text ='MA: Standardised mean difference + 95% CI'))) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = paste("Protein_summary", nstuds, "studies", "below_", range, "heterogeneity", sep="_"),
      width = 900,
      height = 600, dpi=300
    ))
  
}

####### Function for heterogeneity plot
f_phet = function(mts, range, nstuds) {
  nstuds = as.numeric(as.character(nstuds))
  
  sel1 = which(mts$I2<range/100)
  sel2 = which(mts$k.study>=nstuds)
  
  sel = intersect(sel1, sel2)
  
  plot_ly(data = mts[sel,], x = mts[sel,]$TE.random, y = mts[sel,]$I2*100, 
                  text = mts[sel,]$gene_name, mode = "markers",
                  color=as.factor(mts[sel,]$k.study)) %>% 
    layout(legend=list(title=list(text='N of studies')), 
           xaxis = list(title = list(text ='MA: Standardised mean difference')),
           yaxis = list(title = list(text ='I2 % heterogeneity'))) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = "Heterogeneity plot",
      width = 800,
      height = 600, dpi=300
    ))
  
}