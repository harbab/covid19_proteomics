
##### Function for DOR plot
plot_dor = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = mtdats[[z]]
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
    } else {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = studies2[as.numeric(rownames(mtdatr))]
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
  }
  
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  mf = madauni(mtdatr)
  smf = summary(mf)
  
  ##### DOR - univariate
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = nam
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
plot_sens = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = mtdats[[z]]
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
    } else {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = studies2[as.numeric(rownames(mtdatr))]
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
  }
  
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  mf = madauni(mtdatr)
  smf = summary(mf)
  
  ##### DOR - univariate
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = nam
  fmf = fmf[order(fmf$DOR, decreasing = F),]
  addit = c(smf$CIcoef[2,], NA, "Summary")
  fmf = as.data.frame(rbind(addit, fmf))
  fmf$names = ordered(fmf$names, levels=fmf$names)
  n2 = nrow(fmf)-1
  
  ##### Sensitivity plot
  mf1 = as.data.frame(do.call(cbind, mf$descr$sens))
  msf = summary(fit.reitsmas[[protein]])
  smf = msf$coefficients
  mf1$names = nam
  
  addit = c(smf[3,c(1,5,6)], "Mean Sensitivity")
  mf1 = as.data.frame(rbind(addit, mf1))
  mf1$names = as.factor(mf1$names)
  mf1$names = ordered(mf1$names, levels=c("Mean Sensitivity", as.character(fmf$names[-1])))
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
plot_spec = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = mtdats[[z]]
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
    } else {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = studies2[as.numeric(rownames(mtdatr))]
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
  }
  
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  mf = madauni(mtdatr)
  smf = summary(mf)
  
  ##### DOR - univariate
  fmf = as.data.frame(log(do.call(cbind, mf$descr$DOR)))
  fmf$names = nam
  fmf = fmf[order(fmf$DOR, decreasing = F),]
  addit = c(smf$CIcoef[2,], NA, "Summary")
  fmf = as.data.frame(rbind(addit, fmf))
  fmf$names = ordered(fmf$names, levels=fmf$names)
  n2 = nrow(fmf)-1
  
  ##### Specificity plot
  mf2 = as.data.frame(do.call(cbind, mf$descr$spec))
  msf = summary(fit.reitsmas[[z]])
  smf = msf$coefficients
  mf2$names = nam
  
  addit = c(1-smf[4,c(1,6,5)], "Mean Specificity")
  mf2 = as.data.frame(rbind(addit, mf2))
  mf2$names = as.factor(mf2$names)
  mf2$names = ordered(mf2$names, levels=c("Mean Specificity", as.character(fmf$names[-1])))
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

##### Function for SROC
plot_sroc = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = mtdats[[z]]
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
    } else {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = studies2[as.numeric(rownames(mtdatr))]
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
  }
  
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  fit.reitsma = reitsma(mtdatr)
  
  #### Extract data for the SROCs
  s = summary(fit.reitsma)
  sr = sroc(fit.reitsma) ## To get the SROC curve data
  msr = mcsroc(fit.reitsma)
  re = ROCellipse(fit.reitsma, level=0.95)
  rdf = as.data.frame(re$ROCellipse)
  
  sr = as.data.frame(sr)
  sspr = as.data.frame(s$coefficients[3:4,c(1,5:6)])
  tsspr = as.data.frame(t(sspr))
  colnames(tsspr) = c("y", "x")
  
  a = AUC(fit.reitsma)
  a = a$AUC
  
  sfit = SummaryPts(fit.reitsma)
  sf = summary(sfit)
  
  ##### Convert to plotly
  plot_ly(sr, x = sr$fpr, y = sr$V2, type = 'scatter', mode = 'lines', name="SROC") %>%
    add_lines(y = c(0, 1), x = c(0, 1), name="Chance", line=list(dash="dash")) %>%
    add_trace(x = s$coefficients[4], y = s$coefficients[3], 
              name = 'MA Estimate', mode = 'markers', 
              marker=list(color="#0072B2", size=15, symbol=2)) %>%
    add_lines(y=c(sspr$Estimate[1], sspr$Estimate[1]), 
              x=c(sspr$`95%ci.lb`[2], sspr$`95%ci.ub`[2]),
              name="95% CI FPR")  %>%
    add_lines(x=c(sspr$Estimate[2], sspr$Estimate[2]), 
              y=c(sspr$`95%ci.lb`[1], sspr$`95%ci.ub`[1]),
              name="95% CI TPR") %>%
    add_trace(x=rdf$V1, y=rdf$V2, name="95% CI \nellipse around \nmean values") %>%
    add_annotations(x=1, y=0.1,
                    text = paste("AUC: ", round(a, 3)),
                    xanchor = 'right', font=list(size=16),
                    showarrow = F
    ) %>%
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
    config(toImageButtonOptions = list(
        format = "svg",
        filename = paste("SROC_protein", protein, sep="_"),
        width = 600,
        height = 400, dpi=300
      ))
  
}

tab_sroc = function(protein, study, ms) {
  #### Choose a protein of interest
  z = protein
  mtdatr = mtdats[[z]]
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      mtdatr = mtdatr[-rem,]
    } else {
      nam1 = studies2[as.numeric(rownames(mtdatr))]
      rem = which(nam1 %in% study)
      mtdatr = mtdatr[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = studies2[as.numeric(rownames(mtdatr))]
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    mtdatr = mtdatr[-rem,]
  }
  
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  #### Forest plots
  nam = studies2[as.numeric(rownames(mtdatr))]
  fit.reitsma = reitsma(mtdatr)
  
  mtdatr$Author = nam
  mtdatr = mtdatr[,c(5,1:4)]
  mtdatr
  
  #### Extract data for the SROCs
  s = summary(fit.reitsma)
  sspr = as.data.frame(s$coefficients[3:4,c(1,5:6)])
  a = as.numeric(as.character(sspr[1,]))
  b = 1-as.numeric(as.character(sspr[2,]))
  b = b[c(1,3,2)]
  sroctab = rbind(mtdatr, c("Parameter", "Estimate", "95%CI_li", "95%CI_ui", ""),
        c("Sensitivity", round(a, 3), ""),
        c("Specificity", round(b, 3), ""),
        c("AUC", round(s$AUC$AUC, 3), rep("", 3)))
  
  sroctab
  
}

