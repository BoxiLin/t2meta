library(ggplot2)
library(dplyr)


### Read in t1e #####
rm(list = ls())

e <- c("norm")
m <- c("additive")
r <- c("linear")


ee <- "norm"
mm <- "additive"
rr <- "linear"

for (ee in e) {
  for (mm in m) {
    
    glb_err_dist <- ee
    glb_model <- mm
    
    setwd(paste("path_to_output/",glb_err_dist,"/",sep = ""))
    file_list<- list.files()[grep("MC",list.files(), fixed = TRUE)] 
    
    for (file in file_list){
      print(file)
      # if the merged data set doesn't exist, create it
      if (!exists("dataset_raw")){
        dataset_raw <- readRDS(file_list[1])
      }
      # if the merged data set does exist, append to it
      else {
        temp_dataset <-readRDS(file)
        dataset_raw<-rbind(dataset_raw, temp_dataset)
        rm(temp_dataset)
      }
    }
    
    for (rr in r) {
      print(paste("e:",ee,"m:",mm,"r:",rr))
      glb_response_type <- rr
      dataset <- dataset_raw %>% filter(response==glb_response_type)
      dataset$setting <-paste0(dataset$maf,dataset$r, dataset$err_dist, dataset$model)
      
      colnames(dataset)[colnames(dataset)=="p.TSG.f"]<- "TSG.F"
      colnames(dataset)[colnames(dataset)=="p.TSG.min"]<- "TSG.min"
      colnames(dataset)[colnames(dataset)=="p.TSG.1"]<- "TSG.L"
      colnames(dataset)[colnames(dataset)=="p.TSG.q"]<- "TSG.Q"
      colnames(dataset)[colnames(dataset)=="p.TI.joint"]<- "TI.joint"
      colnames(dataset)[colnames(dataset)=="p.T.main"]<- "T.main"
      
      
      test_list = c("TSG.F","TSG.min", "TSG.L", "T.main" , "TSG.Q","TI.joint")
      id_columns <- colnames(dataset)[!(colnames(dataset)%in% test_list | 
                                          colnames(dataset)=="batch")]
      data_avg = aggregate(dataset[, test_list],
                           list(dataset$setting), mean) %>% 
        left_join(unique(dataset[,id_columns]), 
                  by = c("Group.1"="setting")) %>%
        select(-Group.1)
      
      
      df <- reshape2::melt(data_avg, variable.name = 'Method', 
                           id.vars = id_columns[id_columns!="setting"]) 
        
      df$lower =  df$value- 1.96*sqrt(df$value*(1-df$value)/(10^7))
      df$upper =  df$value+ 1.96*sqrt(df$value*(1-df$value)/(10^7))
      
      err_dist_labs = c("Standard Normal", "Student-t, df=4", "Chi-square, df=4")
      names(err_dist_labs) <-  c("norm","t","chisq")
      model = paste("Simulation model =", unique(df$model))
      names(model) <-  unique(df$model)
      
      maf = paste("maf =", unique(dataset$maf))
      names(maf) <-  unique(dataset$maf)
      
      r = paste("r =", unique(dataset$r))
      names(r) <- unique(dataset$r)
      
      df1 <- filter(df, (err_dist == glb_err_dist) & (model == glb_model))
      # Main T1E figures for correctly-specified models
      ggplot2::ggplot(df1, aes_string(y="value",x="r",
                                      group = "Method", color = "Method", 
                                      shape = "Method",size="Method"))+
        geom_errorbar(aes(ymin=lower, ymax=upper, width=0.1))+
        geom_point(size = 3) + 
        geom_line()+
        geom_hline(yintercept=1e-5, linetype="dashed", color = "black")+
        facet_grid(Method~maf, labeller = labeller(maf = maf))+
        scale_color_manual(name = "Method", 
                           values=c("TSG.F" = "black", "TSG.min" = "darkgoldenrod3",
                                    "TSG.L" = "blue","T.main" = "purple",
                                    "TSG.Q" = "hotpink","TI.joint" = "red"))+
        scale_size_manual(name = "Method",
                          values = c("TSG.F"=0.8,"TSG.min"=0.8,
                                     "TSG.L"=0.8,"T.main"=0.8,
                                     "TSG.Q"=0.8,"TI.joint"=1.5))+
        scale_shape_manual(name = "Method", 
                           values = c("TSG.F"=1,"TSG.min"=1,
                                      "TSG.L"=4,"T.main"=1,
                                      "TSG.Q"=4,"TI.joint"=1))+
        labs(shape = "Type")+
        
        xlab("Female-to-female residual error ratio (r)") + ylab("Empirical type I error rate")+
        ylim(c(0.000001,0.00002))+
        theme_bw()+
        ggtitle("Scenario Null: No genetic effect in female and male")
      ggsave(paste("../../../plot/scen0_a1e-5_",glb_model,"-", glb_err_dist,"-",glb_response_type,"heteroskedasticity.jpeg", sep = ""), width = 9, height = 11, dpi = 300)
    }
  }
}


### Power #####

rm(list = ls())

e <- c("norm","t","chisq")
m <- c("additive","dominant")
r <- c("linear","binary")

e <- c("norm")
m <- c("additive")
r <- c("linear")

xvals = c("beta.m", "beta.f","beta.m")
xlabs = c("Genetic effect size","Female genetic effect size","Male genetic effect size")
mains = c("Scenario A1: Homogeneous genetic effects between female and male",
          "Scenario A2: Female only genetic effects",
          "Scenario A3: Heterogeneous genetic effects")

for (ee in e) {
  print(ee)
  glb_err_dist <- ee
  
  setwd(paste("/Volumes/struglis/boxi/frontier3/out/power/",glb_err_dist,"/",sep = ""))
  file_list<- list.files()[grep("20221120",list.files(), fixed = TRUE)] 
  for (file in file_list){
    #    print(file)
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset_raw")){
      dataset_raw <- readRDS(file_list[1])
    }
    # if the merged dataset does exist, append to it
    else {
      temp_dataset <-readRDS(file)
      dataset_raw<-rbind(dataset_raw, temp_dataset)
      rm(temp_dataset)
    }
  }
}

dataset_raw$scenario[dataset_raw$scenario==4] <- 3
dataset<- dplyr::filter(dataset_raw, ((beta.m <= 0.3) &(beta.m >= -0.3)))

colnames(dataset)[colnames(dataset)=="p.TSG.f"]<- "TSG.F"
colnames(dataset)[colnames(dataset)=="p.TSG.min"]<- "TSG.min"
colnames(dataset)[colnames(dataset)=="p.TSG.1"]<- "TSG.L"
colnames(dataset)[colnames(dataset)=="p.TSG.q"]<- "TSG.Q"
colnames(dataset)[colnames(dataset)=="p.TI.joint"]<- "TI.joint"
colnames(dataset)[colnames(dataset)=="p.T.main"]<- "T.main"

test_list = c("TSG.F","TSG.min", "TSG.L", "T.main" , "TSG.Q","TI.joint")

maf = paste("maf =", unique(dataset$maf))
names(maf) <-  unique(dataset$maf)

k = paste("k =", unique(dataset$k))
names(k) <-  unique(dataset$k)

r = paste("r =", unique(dataset$r))
names(r) <-  unique(dataset$r)

err_dist = paste("Residual dist. =", unique(dataset$err_dist))
names(err_dist) <-  unique(dataset$err_dist)

### Main Figure: Normal responses #####
r = "linear"
e = "norm"
m = "additive"
hrr = c(0.5,1,2,5)

for (hr in hrr) {
for (ee in e) {
  for (mm in m) {
    for (rr in r) {
      for (ss in 1:3) {
        glb_model <- mm
        glb_response_type <- rr
        print(paste("error dist:",ee,"model:",mm,"response type:",rr))
        
        dataset_sub <- dataset %>% filter(response==glb_response_type, 
                                          model == glb_model,
                                          err_dist == ee,
                                          scenario == ss, r==hr)
        
        df <- reshape2::melt(dataset_sub, variable.name = 'Method', 
                             id = colnames(dataset_sub)[!colnames(dataset_sub)%in% test_list])
        
        
        # df$lower =  df$value- 1.96*sqrt(df$value*(1-df$value)/(10^7))
        # df$upper =  df$value+ 1.96*sqrt(df$value*(1-df$value)/(10^7))
        
        xval = xvals[ss]
        xlab = xlabs[ss]
        main = mains[ss]
        
        # Main T1E figures for correctly-specified models
        ggplot2::ggplot(df, aes_string(y="value",x=xval,
                                       group = "Method", color = "Method", 
                                       shape = "Method",size="Method"))+
          #      geom_errorbar(aes(ymin=lower, ymax=upper, width=0.1))+
          geom_point(size = 3) + 
          geom_line()+
          facet_grid(.~maf, labeller = labeller(maf = maf))+
          scale_color_manual(name = "Method", 
                             values=c("TSG.F" = "black", "TSG.min" = "darkgoldenrod3",
                                      "TSG.L" = "blue","T.main" = "purple",
                                      "TSG.Q" = "hotpink","TI.joint" = "red"))+
          scale_size_manual(name = "Method",
                            values = c("TSG.F"=0.8,"TSG.min"=0.8,
                                       "TSG.L"=0.8,"T.main"=0.8,
                                       "TSG.Q"=0.8,"TI.joint"=1.5))+
          scale_shape_manual(name = "Method", 
                             values = c("TSG.F"=1,"TSG.min"=1,
                                        "TSG.L"=4,"T.main"=1,
                                        "TSG.Q"=4,"TI.joint"=1))+
          labs(shape = "Type")+
          
          xlab(xlab) + ylab("Empirical power")+
          #      ylim(c(0.000001,0.00003))+
          # scale_x_continuous(breaks=seq(-0.3,0.3,0.1), 
          #                    limits = c(-0.3, 0.3),
          #                    labels=c("-0.3", "-0.2", "-0.1","0","0.1","0.2","0.3"))+
          theme_bw()+
          ggtitle(main)
        ggsave(paste("../../../plot/scene",ss,"_a5e-8_",glb_model,"-",glb_response_type,"-r-",hr,".jpeg", sep = ""), width = 10, height = 3, dpi = 300)
      }
    }
  }
}
}

### Supp Figure: Normal responses: maf*k #####
r = "linear"
e = "norm"
m = "additive"

for (ee in e) {
  for (mm in m) {
    for (rr in r) {
      for (ss in 1:3) {
        glb_model <- mm
        glb_response_type <- rr
        print(paste("error dist:",ee,"model:",mm,"response type:",rr))
        
        dataset_sub <- dataset %>% filter(response==glb_response_type, 
                                          model == glb_model,
                                          err_dist == ee,
                                          scenario == ss)
        
        df <- reshape2::melt(dataset_sub, variable.name = 'Method', 
                             id = colnames(dataset_sub)[!colnames(dataset_sub)%in% test_list])
        
        
        xval = xvals[ss]
        xlab = xlabs[ss]
        main = mains[ss]
        
        # Main T1E figures for correctly-specified models
        ggplot2::ggplot(df, aes_string(y="value",x=xval,
                                       group = "Method", color = "Method", 
                                       shape = "Method",size="Method"))+
          #      geom_errorbar(aes(ymin=lower, ymax=upper, width=0.1))+
          geom_point(size = 3) + 
          geom_line()+
          facet_grid(k~maf, labeller = labeller(maf = maf,k=k))+
          scale_color_manual(name = "Method", 
                             values=c("TSG.F" = "black", "TSG.min" = "darkgoldenrod3",
                                      "TSG.L" = "blue","T.main" = "purple",
                                      "TSG.Q" = "hotpink","TI.joint" = "red"))+
          scale_size_manual(name = "Method",
                            values = c("TSG.F"=0.8,"TSG.min"=0.8,
                                       "TSG.L"=0.8,"T.main"=0.8,
                                       "TSG.Q"=0.8,"TI.joint"=1.5))+
          scale_shape_manual(name = "Method", 
                             values = c("TSG.F"=1,"TSG.min"=1,
                                        "TSG.L"=4,"T.main"=1,
                                        "TSG.Q"=4,"TI.joint"=1))+
          labs(shape = "Type")+
          
          xlab(xlab) + ylab("Empirical power")+
          theme_bw()+
          ggtitle(main)
        ggsave(paste("../../../image/power/supp/scen",ss,"_a5e-8_",glb_model,"-",glb_response_type,"-normal_mafxk.jpeg", sep = ""), 
               width = 10, height = 8, dpi = 300)
      }
    }
  }
}


### Supplementary Figures:quantitative traits ##########
e <- c("norm","t","chisq")
m <- c("additive","dominant")
r <- c("linear","binary")
for (ee in e) {
  for (mm in m) {
    for (rr in r) {
      for (ss in 1:3) {
        glb_model <- mm
        glb_response_type <- rr
        print(paste("error dist:",ee,"model:",mm,"response type:",rr))
        
        dataset_sub <- dataset %>% filter(response==glb_response_type, 
                                          model == glb_model,
                                          scenario == ss, k==1,maf==0.25)
        
        df <- reshape2::melt(dataset_sub, variable.name = 'Method', 
                             id = colnames(dataset_sub)[!colnames(dataset_sub)%in% test_list])
        
        
        # df$lower =  df$value- 1.96*sqrt(df$value*(1-df$value)/(10^7))
        # df$upper =  df$value+ 1.96*sqrt(df$value*(1-df$value)/(10^7))
        
        xval = xvals[ss]
        xlab = xlabs[ss]
        main = mains[ss]
        
        # Main T1E figures for correctly-specified models
        ggplot2::ggplot(df, aes_string(y="value",x=xval,
                                       group = "Method", color = "Method", 
                                       shape = "Method",size="Method"))+
          #      geom_errorbar(aes(ymin=lower, ymax=upper, width=0.1))+
          geom_point(size = 3) + 
          geom_line()+
          facet_grid(.~err_dist, labeller = labeller(maf = maf,k=k))+
          scale_color_manual(name = "Method", 
                             values=c("TSG.F" = "black", "TSG.min" = "darkgoldenrod3",
                                      "TSG.L" = "blue","T.main" = "purple",
                                      "TSG.Q" = "hotpink","TI.joint" = "red"))+
          scale_size_manual(name = "Method",
                            values = c("TSG.F"=0.8,"TSG.min"=0.8,
                                       "TSG.L"=0.8,"T.main"=0.8,
                                       "TSG.Q"=0.8,"TI.joint"=1.5))+
          scale_shape_manual(name = "Method", 
                             values = c("TSG.F"=1,"TSG.min"=1,
                                        "TSG.L"=4,"T.main"=1,
                                        "TSG.Q"=4,"TI.joint"=1))+
          labs(shape = "Type")+
          
          xlab(xlab) + ylab("Empirical power")+
          theme_bw()+
          ggtitle(main)
        ggsave(paste("../../../image/power/supp/scen",ss,"_a5e-8_",glb_model,"-",glb_response_type,".jpeg", sep = ""), width = 10, height = 3, dpi = 300)
      }
    }
  }
}

r = "binary"
e = "norm"
for (ee in e) {
  for (mm in m) {
    for (rr in r) {
      for (ss in 1:3) {
        glb_model <- mm
        glb_response_type <- rr
        print(paste("error dist:",ee,"model:",mm,"response type:",rr))
        
        dataset_sub <- dataset %>% filter(response==glb_response_type, 
                                          model == glb_model,
                                          scenario == ss,maf==0.25)
        
        
        df <- reshape2::melt(dataset_sub, variable.name = 'Method', 
                             id = colnames(dataset_sub)[!colnames(dataset_sub)%in% test_list])
        
        
        xval = xvals[ss]
        xlab = xlabs[ss]
        main = mains[ss]
        
        # Main T1E figures for correctly-specified models
        ggplot2::ggplot(df, aes_string(y="value",x=xval,
                                       group = "Method", color = "Method", 
                                       shape = "Method",size="Method"))+
          #      geom_errorbar(aes(ymin=lower, ymax=upper, width=0.1))+
          geom_point(size = 3) + 
          geom_line()+
          facet_grid(.~k, labeller = labeller(k = k))+
          scale_color_manual(name = "Method", 
                             values=c("TSG.F" = "black", "TSG.min" = "darkgoldenrod3",
                                      "TSG.L" = "blue","T.main" = "purple",
                                      "TSG.Q" = "hotpink","TI.joint" = "red"))+
          scale_size_manual(name = "Method",
                            values = c("TSG.F"=0.8,"TSG.min"=0.8,
                                       "TSG.L"=0.8,"T.main"=0.8,
                                       "TSG.Q"=0.8,"TI.joint"=1.5))+
          scale_shape_manual(name = "Method", 
                             values = c("TSG.F"=1,"TSG.min"=1,
                                        "TSG.L"=4,"T.main"=1,
                                        "TSG.Q"=4,"TI.joint"=1))+
          labs(shape = "Type")+
          
          xlab(xlab) + ylab("Empirical power")+
          theme_bw()+
          ggtitle(main)
        ggsave(paste("../../../image/power/supp/scen",ss,"_a5e-8_",glb_model,"-",glb_response_type,".jpeg", sep = ""), width = 10, height = 3, dpi = 300)
      }
    }
  }
}

