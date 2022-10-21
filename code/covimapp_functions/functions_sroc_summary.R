#### Functions for summary of SROC curves
plot_srocsum = function(dorf, nstuds2) {
  nstuds = as.numeric(as.character(nstuds2))
  
  sel = which(dorf$k.study>=nstuds)
  
  # error_y = list(
  #   type = "data",
  #   symmetric = FALSE,
  #   array = c(dorf$upper.Sensitivity-dorf$Sensitivity),
  #   arrayminus = c(dorf$Sensitivity-dorf$lower.Sensitivity),
  #   thickness=0.3,
  #   width=2)
  # 
  # error_x = list(
  #   type = "data",
  #   symmetric = FALSE,
  #   array = c(dorf$upper.Specificity-dorf$Specificity),
  #   arrayminus = c(dorf$Specificity-dorf$lower.Specificity),
  #   thickness=0.3,
  #   width=2)
  
  plot_ly(data=dorf[sel,], x=1-dorf[sel,]$Specificity, y=dorf[sel,]$Sensitivity,
          text = dorf[sel,]$Gene_name, mode = "markers", type="scatter", 
          size=1,
          #error_y=error_y, error_x=error_x,
          color=as.factor(dorf[sel,]$k.study)) %>%
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


plot_sroc_3d = function(dorf, nstuds2) {
  nstuds = as.numeric(as.character(nstuds2))

  sel = which(dorf$k.study>=nstuds)
  
  plot_ly(dorf[sel,], x=dorf[sel,]$Mean.DOR, y=dorf[sel,]$Sensitivity, z=dorf[sel,]$Specificity,
        text = dorf[sel,]$Gene_name,
        color=as.factor(dorf[sel,]$k.study)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DOR'),
                      yaxis = list(title = 'Sensitivity', range=c(0.5,1)),
                      zaxis = list(title = 'Specificity', range=c(0.5,1))))
  
}


