
############## Function for forest plots
forest_f = function(protein, study, ms) {
  param = list()
  eff_size = list()
  g_eff = list()
  
  for (z in 1:ns) {
    x = df[[z]]
    param[[z]] = x[which(x$gene_name==protein),1:6]
    names(param)[[z]] = studies[[z]]
  }
  
  sume = as.data.frame(do.call(rbind, param))
  
  sume$author = gsub("_.+", "", rownames(sume))
  sume$publication_year = str_extract(rownames(sume), "[0-9]+")
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = sume$author
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      sume = sume[-rem,]
    } else {
      rem = which(sume$author %in% study)
      sume = sume[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = sume$author
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    sume = sume[-rem,]
  }
  
  sume_calc <- esc_mean_sd(grp1m = sume$mean_covid,
                           grp1sd = sume$sd_covid,
                           grp1n = sume$number_covid,
                           grp2m = sume$mean_healthy,
                           grp2sd = sume$sd_healthy,
                           grp2n = sume$number_healthy,
                           study = sume$author,
                           es.type = "g") %>% 
    as.data.frame()
  
  m.gen <- metagen(TE = es,
                   seTE = se,
                   studlab = study,
                   data = sume_calc,
                   sm = "SMD",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   #hakn = TRUE,
                   title = paste("Change in protein:", protein))
  
  m.gen <- update.meta(m.gen, prediction = TRUE)
  
  forest.meta(m.gen, 
              sortvar = TE,
              prediction=TRUE,
              print.tau2 = FALSE,
              leftlabs = c("Author", "n", "g", "SE"))
  
}


############## Function for forest plots
mgen_f = function(protein) {
  param = list()
  eff_size = list()
  g_eff = list()
  
  for (z in 1:ns) {
    x = df[[z]]
    param[[z]] = x[which(x$gene_name==protein),1:6]
    names(param)[[z]] = studies[[z]]
  }
  
  sume = as.data.frame(do.call(rbind, param))
  
  sume$author = gsub("_.+", "", rownames(sume))
  sume$publication_year = str_extract(rownames(sume), "[0-9]+")
  
  sume_calc <- esc_mean_sd(grp1m = sume$mean_covid,
                           grp1sd = sume$sd_covid,
                           grp1n = sume$number_covid,
                           grp2m = sume$mean_healthy,
                           grp2sd = sume$sd_healthy,
                           grp2n = sume$number_healthy,
                           study = sume$author,
                           es.type = "g") %>% 
    as.data.frame()
  
  m.gen <- metagen(TE = es,
                   seTE = se,
                   studlab = study,
                   data = sume_calc,
                   sm = "SMD",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   #hakn = TRUE,
                   title = paste("Change in protein:", protein))
  
  tbl = do.call(cbind, m.gen[c(1:8,10,11)])
  
  f = c("Summary estimate - fixed effect", unlist(m.gen[c(13:19)]), "", "")
  r = c("Summary estimate - random effect", unlist(m.gen[c(20:26)]), "", "")
  tbl2 = as.data.frame(rbind(tbl, f, r))
  
  tbl2[,-1] = apply(tbl2[,-1], 2, as.numeric)
  
  tbl2
  
}

############## Function for funnel plots
funnel_f = function(protein, study, ms) {
  param = list()
  eff_size = list()
  g_eff = list()
  
  for (z in 1:ns) {
    x = df[[z]]
    param[[z]] = x[which(x$gene_name==protein),1:6]
    names(param)[[z]] = studies[[z]]
  }
  
  sume = as.data.frame(do.call(rbind, param))
  
  sume$author = gsub("_.+", "", rownames(sume))
  sume$publication_year = str_extract(rownames(sume), "[0-9]+")
  
  if (study!="None") {
    if(study %in% c("Plasma", "Serum")) {
      nam1 = sume$author
      nam2 = t1$Author[which(t1$Samples==study)]
      rem = which(nam1 %in% nam2)
      sume = sume[-rem,]
    } else {
      rem = which(sume$author %in% study)
      sume = sume[-rem,]
    }
  }
  
  if (ms!="None") {
    nam1 = sume$author
    nam2 = t2$Author[which(t2$`Type of MS acquisition`==ms)]
    rem = which(nam1 %in% nam2)
    sume = sume[-rem,]
  }
  
  sume_calc <- esc_mean_sd(grp1m = sume$mean_covid,
                           grp1sd = sume$sd_covid,
                           grp1n = sume$number_covid,
                           grp2m = sume$mean_healthy,
                           grp2sd = sume$sd_healthy,
                           grp2n = sume$number_healthy,
                           study = sume$author,
                           es.type = "g") %>% 
    as.data.frame()
  
  m.gen <- metagen(TE = es,
                   seTE = se,
                   studlab = study,
                   data = sume_calc,
                   sm = "SMD",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "REML",
                   #hakn = TRUE,
                   title = paste("Change in protein:", protein))
  
  m.gen <- update.meta(m.gen, prediction = TRUE)
  
  col.contour = c("#165a72", "#7cc0d8", "#d3eaf2")
  
  funnel.meta(m.gen,
              xlim = c(min(m.gen$TE)-0.5, max(m.gen$TE)+0.5),
              studlab = TRUE,
              contour = c(0.9, 0.95, 0.99),   ### For contour plots
              col.contour = col.contour
  )
  legend(x = max(m.gen$TE)-0.1, y = 0.01, 
         legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
         fill = col.contour)
  
}

############## Function for summary tables
sume_f = function(protein) {
  param = list()
  eff_size = list()
  g_eff = list()
  
  for (z in 1:ns) {
    x = df[[z]]
    param[[z]] = x[which(x$gene_name==protein),1:6]
    names(param)[[z]] = studies[[z]]
  }
  
  sume = as.data.frame(do.call(rbind, param))
  
  sume$author = rownames(sume)
  sume$publication_year = t1$Year[match(sume$author, t1$Author)]
  
  tbl = sume[,c(7,8,1:6)]
  tbl$nsum = tbl$number_covid+tbl$number_healthy
  rownames(tbl) = NULL
  colnames(tbl) = c("Author", "Publication Year", "Mean COVID19", "Mean Controls", 
                    "SD COVID19", "SD Controls", "Number of cases", "Number of Controls", "Sum")
  tbl
  
  tbl[,3:6] = apply(tbl[,3:6], 2, function(x) round(x, 3))
  sr = c("All studies", "2020-2022", rep("", 4), sum(tbl$`Number of cases`), 
         sum(tbl$`Number of Controls`), sum(tbl$Sum))
  tbl = rbind(tbl, sr)
  
  tbl
}