
##### Function for DOR plot
plot_dor = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = pprot[[z]][,-c(1,6)]
  names = pprot[[z]]$study
  names = paste(str_to_title(gsub("_.+", "", names)), "et al.")
  sel = grep("-", names)
  names[sel] = paste(gsub("-[A-Z][a-z]+", "", names[sel]), 
        str_extract(names[sel], "-[A-Z][a-z]+"), sep="")
  
  if (any(study!="None")) {
    if(any(study %in% c("Plasma", "Plasma and serum", "Serum"))) {
      nam1 = names
      nam2 = t2$Author[which(t2$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    } else {
      nam1 = names
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    }
  }
  
  if (ms!="None") {
    nam1 = names
    nam2 = t2$Author[which(t2$Type.of.MS.acquisition==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
    names = names[-rem]
  }
  
  #### Forest plots
  mf = madauni(mtdatr)
  smf = summary(mf)
  
  ##### DOR - univariate
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = names
  fmf = fmf[order(fmf$DOR, decreasing = F),]
  addit = c(smf$CIcoef[2,], NA, "Summary")
  fmf = as.data.frame(rbind(addit, fmf))
  fmf$names = ordered(fmf$names, levels=fmf$names)
  fmf[,1:4] = apply(fmf[,1:4], 2, function(x) round(as.numeric(x), 3))
  n2 = nrow(fmf)-1
  
  error_fmf_x = list(type = "data",
                     symmetric = FALSE,
                     array = c(fmf$`97.5%`-fmf$DOR),
                     arrayminus = c(fmf$DOR-fmf$`2.5%`),
                     thickness=1,
                     width=2)
  plot_ly(fmf, x = fmf$DOR, y=fmf$names, 
          error_x=error_fmf_x,
          type = 'scatter', mode="markers",
          marker=list(size=c(15, rep(12,n2)), 
                      symbol=c(2, rep(1, n2)),
                      color=c("#0072B2", rep("black", n2)))) %>%
    layout(title=paste("Forest plot for protein:", z),
           xaxis = list(
             range=c(round((min(fmf$`2.5%`)-0.5), 1), 
                     round((max(fmf$`97.5%`)+0.5), 1)),
             title='Diagnostic Odds Ratio (DOR), with 95% CI (log)',
             font=list(size=15),
             zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=1),
           yaxis = list(
             font=list(size=15),
             zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1)) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = paste("DOR_protein", protein, sep="_"),
      width = 600,
      height = 400, dpi=300
    ))
  
}

###### Function for sensitivity
plot_glmm_sens = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = pprot[[z]][,-c(1,6)]
  names = pprot[[z]]$study
  names = paste(str_to_title(gsub("_.+", "", names)), "et al.")
  sel = grep("-", names)
  names[sel] = paste(gsub("-[A-Z][a-z]+", "", names[sel]), 
                     str_extract(names[sel], "-[A-Z][a-z]+"), sep="")
  
  if (any(study!="None")) {
    if(any(study %in% c("Plasma", "Plasma and serum", "Serum"))) {
      nam1 = names
      nam2 = t2$Author[which(t2$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    } else {
      nam1 = names
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    }
  }
  
  if (ms!="None") {
    nam1 = names
    nam2 = t2$Author[which(t2$Type.of.MS.acquisition==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
    names = names[-rem]
  }
  
  
  #### Forest plots
  mf = madauni(mtdatr)
  
  ##### DOR - univariate (To keep the same order in the plot)
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = names
  fmf = fmf[order(fmf$DOR, decreasing = F),]
  n2 = nrow(fmf)
  
  ##### Sensitivity plot
  mf1 = as.data.frame(do.call(cbind, mf$descr$sens))
  smf = as.numeric(tabs[[z]][1,1:3])
  mf1$names = names
  
  addit = c(smf, "Mean Sensitivity")
  mf1 = as.data.frame(rbind(addit, mf1))
  mf1$names = as.factor(mf1$names)
  mf1$names = ordered(mf1$names, levels=c("Mean Sensitivity", as.character(fmf$names)))
  mf1[,1:3] = apply(mf1[,1:3], 2, function(x) round(as.numeric(x), 3))
  
  error_mf1_x = list(type = "data",
                     symmetric = FALSE,
                     array = c(mf1$`97.5%`-mf1$sens),
                     arrayminus = c(mf1$sens-mf1$`2.5%`),
                     thickness=1,
                     width=2)
  
  plot_ly(mf1, x = mf1$sens, y=mf1$names, 
          error_x=error_mf1_x,
          type = 'scatter',
          marker=list(size=c(15, rep(12,n2)), 
                      symbol=c(2, rep(1, n2)),
                      color=c("#0072B2", rep("black", n2)))
          ) %>%
    layout(title=paste("Forest plot for protein:", z),
           xaxis = list(
             title='Sensitivity, with 95% CI',
             font=list(size=15),
             range=c(0,1),
             #zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1),
           yaxis = list(
             font=list(size=15),
             #zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1)) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = paste("Sensitivity_protein", protein, sep="_"),
      width = 600,
      height = 400, dpi=300
    ))
}

##### Function for specificity
plot_glmm_spec = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = pprot[[z]][,-c(1,6)]
  names = pprot[[z]]$study
  names = paste(str_to_title(gsub("_.+", "", names)), "et al.")
  sel = grep("-", names)
  names[sel] = paste(gsub("-[A-Z][a-z]+", "", names[sel]), 
                     str_extract(names[sel], "-[A-Z][a-z]+"), sep="")
  
  if (any(study!="None")) {
    if(any(study %in% c("Plasma", "Plasma and serum", "Serum"))) {
      nam1 = names
      nam2 = t2$Author[which(t2$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    } else {
      nam1 = names
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    }
  }
  
  if (ms!="None") {
    nam1 = names
    nam2 = t2$Author[which(t2$Type.of.MS.acquisition==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
    names = names[-rem]
  }
  
  
  #### Forest plots
  mf = madauni(mtdatr)
  
  ##### DOR - univariate (To keep the same order in the plot)
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = names
  fmf = fmf[order(fmf$DOR, decreasing = F),]
  n2 = nrow(fmf)
  
  ##### Specificity plot
  mf2 = as.data.frame(do.call(cbind, mf$descr$spec))
  smf = as.numeric(tabs[[z]][2,1:3])
  mf2$names = names
  
  addit = c(smf, "Mean Specificity")
  mf2 = as.data.frame(rbind(addit, mf2))
  mf2$names = as.factor(mf2$names)
  mf2$names = ordered(mf2$names, levels=c("Mean Specificity", as.character(fmf$names)))
  mf2[,1:3] = apply(mf2[,1:3], 2, function(x) round(as.numeric(x), 3))
  
  error_mf2_x = list(type = "data",
                     symmetric = FALSE,
                     array = c(mf2$`97.5%`-mf2$spec),
                     arrayminus = c(mf2$spec-mf2$`2.5%`),
                     thickness=1,
                     width=2)
  
  plot_ly(mf2, x = mf2$spec, y=mf2$names, 
          error_x=error_mf2_x,
          type = 'scatter',
          marker=list(size=c(15, rep(12,n2)), 
                      symbol=c(2, rep(1, n2)),
                      color=c("#0072B2", rep("black", n2)))
          ) %>%
    layout(title=paste("Forest plot for protein:", z),
           xaxis = list(
             title='Specificity, with 95% CI',
             range=c(0,1),
             font=list(size=15),
             #zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1),
           yaxis = list(
             font=list(size=15),
             #zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1)) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = paste("Specificity_protein", protein, sep="_"),
      width = 600,
      height = 400, dpi=300
    ))
  
}

source("functions_for_glmm_srocs.R")

##### Function for SROC
plot_glmm_sroc = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = pprot[[z]][,-c(1)]
  names = pprot[[z]]$study
  names = paste(str_to_title(gsub("_.+", "", names)), "et al.")
  sel = grep("-", names)
  names[sel] = paste(gsub("-[A-Z][a-z]+", "", names[sel]), 
                     str_extract(names[sel], "-[A-Z][a-z]+"), sep="")
  
  if (any(study!="None")) {
    if(any(study %in% c("Plasma", "Plasma and serum", "Serum"))) {
      nam1 = names
      nam2 = t2$Author[which(t2$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    } else {
      nam1 = names
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
    }
  }
  
  if (ms!="None") {
    nam1 = names
    nam2 = t2$Author[which(t2$Type.of.MS.acquisition==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
    names = names[-rem]
  }
  
  long_data = reshape.data(mtdatr)
  model <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|study), 
                  data=long_data, 
                  family = binomial(link='logit'),           
                  nAGQ=0)
  table = tab.glmer(model)
  fitted = fit.glmer(model, mtdatr)
  
  ##### Make a plotly SROC plot #####
  #### Extract data for the SROCs
  sr = sroc.glmer(fitted) ## To get the SROC curve data
  re = ROCellipse.glmer(fitted, level=0.95)
  rdf = as.data.frame(re$ROCellipse)
  
  sr = as.data.frame(sr)
  s = table[,1:3]
  s[2,] = 1-s[2,]
  sspr = as.data.frame(s)
  tsspr = as.data.frame(t(sspr))
  colnames(tsspr) = c("y", "x")
  a = AUC_calc(sr$sens)
  
  # To figure out how to calculate DOR and LR
  sf = summary(SummaryPts.glmer(fitted))
  sf[4, c(1,3,4)]
  
  ##### Convert to plotly
  plot_ly(x = sr$fpr, y = sr$sens, type = 'scatter', mode = 'lines', name="SROC") %>%
    add_lines(y = c(0, 1), x = c(0, 1), name="Chance", line=list(dash="dash")) %>%
    add_trace(x = s[2,1], y = s[1,1], 
              name = 'MA Estimate', mode = 'markers', 
              marker=list(color="#0072B2", size=15, symbol=2)) %>%
    add_lines(y=c(sspr$est[1], sspr$est[1]), 
              x=c(sspr$`2.5 %`[2], sspr$`97.5 %`[2]),
              name="95% CI FPR")  %>%
    add_lines(x=c(sspr$est[2], sspr$est[2]), 
              y=c(sspr$`2.5 %`[1], sspr$`97.5 %`[1]),
              name="95% CI TPR") %>%
    add_trace(x=rdf$V1, y=rdf$V2, name="95% CI \nellipse around \nmean values") %>%
    add_annotations(x=1, y=0.1,
                    text = paste("AUC: ", round(a, 3)),
                    xanchor = 'right', font=list(size=16),
                    showarrow = F
    ) %>%
    layout(title=paste("SROC curve for protein:", z),
           xaxis = list(
             title='FPR (1-Specificity)',
             zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1),
           yaxis = list(
             title='Sensitivity',
             zerolinecolor = 'gray',
             zerolinewidth = 1,
             #gridcolor = 'ffff',
             dtick=0.1)) %>%
    add_annotations(x=1, y=0.3,
                    text = paste("Mean DOR: ", round(sf[4,1], 3)),
                    xanchor = 'right', font=list(size=16),
                    showarrow = F
    )  %>%
    add_annotations(x=1, y=0.25,
                    text = paste("95% CI:", round(sf[4,3], 3), "-", round(sf[4,4], 3)),
                    xanchor = 'right', font=list(size=16),
                    showarrow = F
    ) %>%
    config(toImageButtonOptions = list(
      format = "svg",
      filename = paste("SROC_protein", protein, sep="_"),
      width = 600,
      height = 400, dpi=300
    ))
  
}

tab_glmm_sroc = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = pprot[[z]][,-c(1,6)]
  names1 = pprot[[z]]$study
  names = paste(str_to_title(gsub("_.+", "", names1)), "et al.")
  sel = grep("-", names)
  names[sel] = paste(gsub("-[A-Z][a-z]+", "", names[sel]), 
                     str_extract(names[sel], "-[A-Z][a-z]+"), sep="")
  
  if (any(study!="None")) {
    if(any(study %in% c("Plasma", "Plasma and serum", "Serum"))) {
      nam1 = names
      nam2 = t2$Author[which(t2$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
      names1 = names1[-rem]
    } else {
      nam1 = names
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
      names = names[-rem]
      names1 = names1[-rem]
    }
  }
  
  if (ms!="None") {
    nam1 = names
    nam2 = t2$Author[which(t2$Type.of.MS.acquisition==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
    names = names[-rem]
    names1 = names1[-rem]
  }

  res = madauni(mtdatr)
  ress = cbind(cbind(res$descr$sens$sens, res$descr$sens$sens.ci), 
    cbind(res$descr$spec$spec, res$descr$spec$spec.ci))
  ress = apply(ress, 2, as.numeric)
  
  mtdatr$study = names
  mtdatr = mtdatr[,c(5,1:4)]
  sum_tab = cbind(mtdatr, ress)
  colnames(sum_tab)[grep("V", colnames(sum_tab))] = c("Sensitivity", "Specificity")
  
  long_data = reshape.data(mtdatr)
  model <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|study), 
                 data=long_data, 
                 family = binomial(link='logit'),           
                 nAGQ=0)
  table = tab.glmer(model)
  fitted = fit.glmer(model, mtdatr)
  sr = sroc.glmer(fitted) ## To get the SROC curve data
  sr = as.data.frame(sr)
  s = table[,1:3]
  sspr = as.data.frame(s)
  a = round(AUC_calc(sr$sens), 3)
  sf = summary(SummaryPts.glmer(fitted))
  
  alpha = selres$alphamin[which(selres$protein==z)]
  names(alpha) = selres$study[which(selres$protein==z)]
  C1 = selres$c1[which(selres$protein==z)]
  names(C1) = selres$study[which(selres$protein==z)]
  
  sum_tab$Alpha = alpha[match(names1, names(alpha))]
  sum_tab$C1 = C1[match(names1, names(C1))]

  val = list()
  for(i in 1:8) {
    val[[i]] = rep("", 4)
  }
  
  s = apply(s, 2, function(x) round(x, 3))
  add1 = cbind(rownames(s), s[,1], rep("",2), s[,2:3], rep("", 2), 
               c("AUC:", a), do.call(cbind, val)[1:2,-c(1:2)])
  
  add2 = cbind(rownames(sf), sf, do.call(cbind, val))
  colnames(add1) = colnames(sum_tab)
  colnames(add2) = colnames(sum_tab)
  sum_tab[,6:13] = apply(sum_tab[,6:13], 2, function(x) round(x, 3))
  
  sroctab = rbind(sum_tab, rep("", 13), 
                  c("SUMMARY RESULTS", rep("", 12)),
                  c("Estimate", "Mean", "", colnames(s)[2:3], rep("", 8)), add1, 
                  c("Estimate", colnames(sf), rep("", 8)), add2)
  sroctab 
}
