#### Functions for summary of SROC curves
plot_glmm_srocsum = function(bivglm_sign, nstuds2) {
  nstuds = as.numeric(as.character(nstuds2))
  
  sel = which(bivglm_sign$N_study>=nstuds)
  
  plot_ly(data=bivglm_sign[sel,], x=1-bivglm_sign[sel,]$Specificity, y=bivglm_sign[sel,]$Sensitivity,
          text = bivglm_sign[sel,]$protein, mode = "markers", type="scatter", 
          size=1,
          #error_y=error_y, error_x=error_x,
          color=as.factor(bivglm_sign[sel,]$N_study)) %>%
    add_lines(y = c(0, 1), x = c(0, 1), name="Chance", line=list(dash="dash"), inherit=F) %>%
    layout(title=paste("SROC protein summary"),
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
      filename = "SROC_summary",
      width = 600,
      height = 400, dpi=300
    ))
  
}


plot_glmm_sroc_3d = function(bivglm_sign, nstuds2) {
  nstuds = as.numeric(as.character(nstuds2))

  sel = which(bivglm_sign$N_study>=nstuds)
  
  plot_ly(bivglm_sign[sel,], x=bivglm_sign[sel,]$AUC, y=bivglm_sign[sel,]$Sensitivity, z=bivglm_sign[sel,]$Specificity,
        text = bivglm_sign[sel,]$protein,
        color=as.factor(bivglm_sign[sel,]$N_study)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'AUC'),
                      yaxis = list(title = 'Sensitivity', range=c(0.5,1)),
                      zaxis = list(title = 'Specificity', range=c(0.5,1))))
  
}


