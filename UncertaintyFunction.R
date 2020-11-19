# Bayesian Method to estimate EY.
library(MCMCpack)
library(dplyr)
library(ggplot2)
library(poLCA)
library(grid)

Data2EY.bayes = function(data, Z.sims = 50, N.sims = 400){
  # input data, and size, 
  # output: a dataframe 
  
  data = as.data.frame(data)
  # lcm to get posteiror distribution
  lcm = poLCA(data = data+1, cbind(int.ha, int.hai, int.na) ~ 1, 
        nclass = 2, na.rm = F)
  
  df.exp.prob = data.frame()
  
  for (sim in 1:Z.sims){
    cat("sample Z from posteiror distribution at ", sim, "time. \n")
    # simulate Z
    data$deZ = apply(lcm$posterior, 1, function(x) sample(c(0,1), size = 1, prob = x))
    
    # fit bayesian logistic model
    m1 = glm(h1 ~ int.ha + int.hai + int.na + deZ , family = binomial, data = data)
    glm.bayes= MCMClogit(h1 ~ int.ha + int.hai + int.na + deZ, 
                         data = data, burnin = 5000, mcmc = 5000, b0=0, B0=0.01) # 
    
    # sample beta from posteiror, 400, 40000, every 1000 
    design.X = model.matrix(m1)
    sample.index = floor(seq(2000,5000, length.out = N.sims))
    ##  EY = design.mat  %*% t(glm.bayes[sample.index,])
    ## plot(glm.bayes)
    
    
    
    HL.ha = max(data$int.ha)
    HL.hai = max(data$int.hai)
    HL.na = max(data$int.na)
    ## For single
    ##------
    # # EY.HAI.
    # for (i in 0:HL.hai){ 
    #   new.X = design.X 
    #   new.X[, (1 + 1):(HL.hai + 1)] = 0
    #   if( i > 0 ) {new.X[,i+1] = 1}
    #   EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
    #   EY = colMeans(exp(EY.expit) / (1+exp(EY.expit))) 
    #   df.exp.prob = rbind(df.exp.prob, 
    #                       matrix(c(rep(i, length(EY)), EY, 
    #                                rep("All", length(EY)), rep("HAI", length(EY))), 
    #                              byrow = F,ncol = 4 ) )
    # }
    # # EY.HA
    # for (i in 0:HL.ha){ 
    #   new.X = design.X 
    #   new.X[, (1 + HL.hai + 1):(1 + HL.hai + HL.ha)] = 0
    #   if( i > 0 ) {new.X[,1 + HL.hai + i] = 1}
    #   EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
    #   EY = colMeans(exp(EY.expit) / (1+exp(EY.expit))) 
    #   df.exp.prob = rbind(df.exp.prob, 
    #                       matrix(c(rep(i, length(EY)), EY, 
    #                                rep("All", length(EY)), rep("HA", length(EY))), 
    #                              byrow = F,ncol = 4 ) )
    # }
    # # EY.NA
    # for (i in 0:HL.na){ 
    #   new.X = design.X 
    #   new.X[, (1 + HL.hai + HL.ha + 1):(1 + HL.hai + HL.ha + HL.na)] = 0
    #   if( i > 0 ) {new.X[,1 + HL.hai + HL.ha + i] = 1}
    #   EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
    #   EY = colMeans(exp(EY.expit) / (1+exp(EY.expit))) 
    #   df.exp.prob = rbind(df.exp.prob, 
    #                       matrix(c(rep(i, length(EY)), EY, 
    #                                rep("All", length(EY)), rep("NA", length(EY))), 
    #                              byrow = F,ncol = 4 ) )
    # }
    ##-----
    
    ## For multiple cause
    for (i in 0:HL.ha){
      new.X = design.X 
      new.X[, (1 + 1):(1 + HL.hai + HL.ha + HL.na)] = 0
      if( i > 0 ) {new.X[,i+1] = 1}
      
      for (j in 0:HL.hai){
        if( j > 0 ) {new.X[,1 + HL.ha + j] = 1}
        
        for(k in 0:HL.na){
          if( k > 0 ) {new.X[,1 + HL.ha + HL.hai + k] = 1}
          
          # calculate corresponding EY for i,j,k
          EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
          EY = colMeans(exp(EY.expit) / (1+exp(EY.expit))) 
          df.exp.prob = rbind(df.exp.prob, 
                              matrix(c(rep(4*i + 2*j + k + 1, length(EY)),
                                       rep(i, length(EY)), 
                                       rep(j, length(EY)),
                                       rep(k, length(EY)),
                                       EY),
                                     byrow = F,ncol = 5 )
                              )
                              
          
        }
      }
    }
    
    
  }
  
  # Add name to the dataframe
  colnames(df.exp.prob) = c("Group", "HA", "HAI", "int.NA", "EY")
  df.exp.prob$EY = as.numeric(as.character(df.exp.prob$EY))
  
  return(df.exp.prob)
  
}

Data2EY.bayes.marginal = function(data, Z.sims=50, N.sims=400){
  # input data, and size, 
  # output: a dataframe 
  
  data = as.data.frame(data)
  # lcm to get posteiror distribution
  lcm = poLCA(data = data+1, cbind(int.ha, int.hai, int.na) ~ 1, 
              nclass = 2, na.rm = F)
  
  df.exp.prob = data.frame()
  
  for (sim in 1:Z.sims){
    cat("sample Z from posteiror distribution at ", sim, "time. \n")
    # simulate Z
    data$deZ = apply(lcm$posterior, 1, function(x) sample(c(0,1), size = 1, prob = x))
    
    # fit bayesian logistic model
    m1 = glm(h1 ~ int.ha + int.hai + int.na + deZ , family = binomial, data = data)
    glm.bayes= MCMClogit(h1 ~ int.ha + int.hai + int.na + deZ, 
                         data = data, burnin = 5000, mcmc = 5000, b0=0, B0=0.01) # 
    
    # sample beta from posteiror, 400, 40000, every 1000 
    design.X = model.matrix(m1)
    sample.index = floor(seq(2000,5000, length.out = N.sims))
    ##  EY = design.mat  %*% t(glm.bayes[sample.index,])
    ## plot(glm.bayes)
    
    
    
    HL.ha = max(data$int.ha)
    HL.hai = max(data$int.hai)
    HL.na = max(data$int.na)
    ## For single
    ##------
    # EY.HA
    for (i in 0:HL.ha){
      new.X = design.X
      new.X[, (1 + 1):(1+HL.ha)] = 0
      if( i > 0 ) {new.X[,i+1] = 1}
      EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
      EY = colMeans(exp(EY.expit) / (1+exp(EY.expit)))
      df.exp.prob = rbind(df.exp.prob,
                          matrix(c(rep(i, length(EY)), EY,
                                   rep("All", length(EY)), rep("HA", length(EY))),
                                 byrow = F,ncol = 4 ) )
    }
    
    # EY.HAI
    for (i in 0:HL.hai){
      new.X = design.X
      new.X[, (1 + HL.ha + 1):(1 + HL.ha + HL.hai)] = 0
      if( i > 0 ) {new.X[,1 + HL.ha + i] = 1}
      EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
      EY = colMeans(exp(EY.expit) / (1+exp(EY.expit)))
      df.exp.prob = rbind(df.exp.prob,
                          matrix(c(rep(i, length(EY)), EY,
                                   rep("All", length(EY)), rep("HAI", length(EY))),
                                 byrow = F,ncol = 4 ) )
    }
    # EY.NA
    for (i in 0:HL.na){
      new.X = design.X
      new.X[, (1 + HL.hai + HL.ha + 1):(1 + HL.hai + HL.ha + HL.na)] = 0
      if( i > 0 ) {new.X[,1 + HL.hai + HL.ha + i] = 1}
      EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
      EY = colMeans(exp(EY.expit) / (1+exp(EY.expit)))
      df.exp.prob = rbind(df.exp.prob,
                          matrix(c(rep(i, length(EY)), EY,
                                   rep("All", length(EY)), rep("NA", length(EY))),
                                 byrow = F,ncol = 4 ) )
    }
    ##-----
    
    ## For multiple cause
    ##--------
    # for (i in 0:HL.ha){
    #   new.X = design.X 
    #   new.X[, (1 + 1):(1 + HL.hai + HL.ha + HL.na)] = 0
    #   if( i > 0 ) {new.X[,i+1] = 1}
    #   
    #   for (j in 0:HL.hai){
    #     if( j > 0 ) {new.X[,1 + HL.ha + j] = 1}
    #     
    #     for(k in 0:HL.na){
    #       if( k > 0 ) {new.X[,1 + HL.ha + HL.hai + k] = 1}
    #       
    #       # calculate corresponding EY for i,j,k
    #       EY.expit = new.X  %*% t(glm.bayes[sample.index,]) #X^T %*% beta
    #       EY = colMeans(exp(EY.expit) / (1+exp(EY.expit))) 
    #       df.exp.prob = rbind(df.exp.prob, 
    #                           matrix(c(rep(4*i + 2*j + k + 1, length(EY)),
    #                                    rep(i, length(EY)), 
    #                                    rep(j, length(EY)),
    #                                    rep(k, length(EY)),
    #                                    EY),
    #                                  byrow = F,ncol = 5 )
    #       )
    #       
    #       
    #     }
    #   }
    # }
    ##-------
    
    
  }
  
  # Add name to the dataframe
  colnames(df.exp.prob) = c("Level", "EY", "Age.Grp", "Biomarker")
  df.exp.prob$EY = as.numeric(as.character(df.exp.prob$EY))
  
  return(df.exp.prob)
  
}

Data2EY.bootstrap.marginal = function(data, bootstime = 500){
  
  df.exp.prob = data.frame()
  n.data = nrow(data)
  data = as.data.frame(data)
  for (boots in 1:bootstime){
    cat("Boots Time:", boots, "\n")
    data.boots = data[sample(1:n.data, n.data, replace = T),]
    
    lcm.boots = poLCA(data = as.data.frame(data.boots+1), cbind(int.ha, int.hai, int.na) ~ 1, 
                      nclass = 2, na.rm = F)
    
    data.boots$deZ = lcm.boots$posterior[,1]
    glm0 = glm(data = data.boots , 
               h1 ~ + int.ha + int.hai + int.na + deZ, 
               family = binomial)
    
    ## Aftering resampling, some level may not occur, so we need ti adjust the level.
    hai.levels = max(data.boots$int.hai)
    #if("NA's" %in% names(hai.levels) ) hai.levels = hai.levels[-length(hai.levels)]
    #zero_ind = which(hai.levels == 0)
    #if( !is_empty(zero_ind) ) hai.levels = hai.levels[-zero_ind]
    
    
    # All ages
    ## HAI
    HL.ha = max(data.boots$int.ha)
    HL.hai = max(data.boots$int.hai)
    HL.na = max(data.boots$int.na)
    EY.hai.3 = rep(0,HL.hai+1)
    newdata2 = data.boots
    for (hai in 0:HL.hai){
      newdata2$int.hai = rep(hai,n.data)
      EY.hai.3[hai+1] = mean(predict(glm0, newdata= newdata2, type="response"),na.rm=T)
    }
    df.exp.prob = rbind(df.exp.prob, 
                        matrix(c(0:HL.hai, 
                                 EY.hai.3, 
                                 rep("All", length(EY.hai.3)), rep("HAI", length(EY.hai.3))), 
                               byrow = F,ncol = 4 ) )
  
    
    EY.ha2.3 = rep(0,HL.ha+1)
    newdata2 = data.boots
    for (ha2 in 0:HL.ha){
      newdata2$int.ha = rep(ha2,n.data)
      EY.ha2.3[ha2+1] = mean(predict(glm0, newdata= newdata2, type="response"),na.rm=T)
    }
    df.exp.prob = rbind(df.exp.prob, 
                        matrix(c(0:HL.ha, 
                                 EY.ha2.3, 
                                 rep("All", length(EY.ha2.3)), 
                                 rep("HA", length(EY.ha2.3))), 
                               byrow = F,ncol = 4 ) )
    
    EY.na.3 = rep(0,HL.na+1)
    newdata2 = data.boots
    for (na in 0:HL.na){
      newdata2$int.na = rep(na,n.data)
      EY.na.3[na+1] =  mean(predict(glm0, newdata= newdata2, type="response"),na.rm=T)
    }
    df.exp.prob = rbind(df.exp.prob, 
                        matrix(c(0:HL.na, 
                                 EY.na.3, 
                                 rep("All", length(EY.na.3)), 
                                 rep("NA", length(EY.na.3))), 
                               byrow = F,ncol = 4 ) )
    
    
    
  }
  # Add name to the data.frame
  colnames(df.exp.prob) = c("Level", "EY", "Age.Grp", "Biomarker")
  df.exp.prob$EY = as.numeric(as.character(df.exp.prob$EY))
  
  return(df.exp.prob)

  
}
  
Data2EY.bootstrap = function(data, bootstime = 500){
  
  df.exp.prob = data.frame()
  n.data = nrow(data)
  
  data = as.data.frame(data)
  for (boots in 1:bootstime){
    cat("Boots Time:", boots, "\n")
    data.boots = data[sample(1:n.data, n.data, replace = T),]
    
    lcm.boots = poLCA(data = as.data.frame(data.boots+1), cbind(int.ha, int.hai, int.na) ~ 1, 
                           nclass = 2, na.rm = F)
    
    data.boots$deZ = lcm.boots$posterior[,1]
    glm0 = glm(data = data.boots , 
               h1 ~ + int.ha + int.hai + int.na + deZ, 
               family = binomial)
    
    ## Aftering resampling, some level may not occur, so we need ti adjust the level.
    ## For multiple cause
    HL.ha = max(data.boots$int.ha)
    HL.hai = max(data.boots$int.hai)
    HL.na = max(data.boots$int.na)
    for (i in 0:HL.ha){
      for (j in 0:HL.hai){
        for(k in 0:HL.na){
          newdata2 = data.boots
          newdata2$int.ha = rep(i, n.data)
          newdata2$int.hai = rep(j, n.data)
          newdata2$int.na = rep(k, n.data)
          EY = mean(predict(glm0, newdata= newdata2, type="response"), na.rm=T)
          df.exp.prob = rbind(df.exp.prob, 
                              c(4*i + 2*j + k + 1,
                                i, j, k, EY))
        }
      }
    }
  }
  
  # Add name to the data.frame
  colnames(df.exp.prob) = c("Group", "HA", "HAI", "int.NA", "EY")
  df.exp.prob$EY = as.numeric(as.character(df.exp.prob$EY))
  
  return(df.exp.prob)
}

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
plot.marginal = function(data, df.result, Z.sims, N.sims, method = "Bayes"){
  # Input: data is the simulated dataset
  # df.result 
  # n is the dataset size; 
  # Z.sims is the sampled z times from posterior distribution
  # N.sims is the sampled beta times from posterior distribution
  # caculate marginal 
  data = as.data.frame(data)
  n = nrow(data)
  marginal.HA = data %>% group_by(int.ha) %>% summarise(EY = mean(h1))
  marginal.HAI = data %>% group_by(int.hai) %>% summarise(EY = mean(h1))
  marginal.NA = data %>% group_by(int.na) %>% summarise(EY = mean(h1))
  
  # Insert marginal EY into dataframe
  df.Bayes.marginal = df.result
  df.Bayes.marginal$TrueEY = 0
  for(i in 1:nrow(marginal.HA)){
    df.Bayes.marginal[which(df.Bayes.marginal$Biomarker == "HA" & 
                              df.Bayes.marginal$Level == marginal.HA$int.ha[i]), 
                      "TrueEY"] = marginal.HA$EY[i]
  }
  for(i in 1:nrow(marginal.NA)){
    df.Bayes.marginal[which(df.Bayes.marginal$Biomarker == "NA" & 
                              df.Bayes.marginal$Level == marginal.NA$int.na[i]), 
                      "TrueEY"] = marginal.NA$EY[i]
  }
  for(i in 1:nrow(marginal.HAI)){
    df.Bayes.marginal[which(df.Bayes.marginal$Biomarker == "HAI" & 
                              df.Bayes.marginal$Level == marginal.HAI$int.hai[i]), 
                      "TrueEY"] = marginal.HAI$EY[i]
  }
  
  # Plot
  p1 = df.Bayes.marginal %>% filter(Age.Grp == "All", Biomarker == "HAI")  %>%
    ggplot(aes(x = Level, y = EY)) + geom_violin() + 
    geom_point(aes(x = Level, y = TrueEY), shape =4, color = "#56B4E9") + 
    ylab("") + xlab("Log-transformed Level") +
    ggtitle(paste(method,",N=", n, ", #Z = ", Z.sims, "\n#beta = ", N.sims, ", HAI")) +
    stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red", size = 0.2) + ylim(0,1) 
  
  p2 = df.Bayes.marginal %>% filter(Age.Grp == "All", Biomarker == "HA")  %>%
    ggplot(aes(x = Level, y = EY)) + geom_violin() + 
    geom_point(aes(x = Level, y = TrueEY), shape =4, color = "#56B4E9") + 
    ylab("") + xlab("Log-transformed Level") +
    ggtitle("HA Stalk") +
    stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red", size = 0.2) + ylim(0,1) 
  
  p3 = df.Bayes.marginal %>% filter(Age.Grp == "All", Biomarker == "NA")  %>%
    ggplot(aes(x = Level, y = EY)) + geom_violin() + 
    geom_point(aes(x = Level, y = TrueEY), shape =4, color = "#56B4E9") + 
    ylab("") + xlab("Log-transformed Level") +
    ggtitle("NA") +
    stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red", size = 0.2) + ylim(0,1) 
  
  png(paste("sim1.",method,".n",n,"_Z",Z.sims,"_b",N.sims,".png", sep = ""),
      res = 100, width=788,height=534)
  grid.newpage()  ##新建页面
  pushViewport(viewport(layout = grid.layout(1,3))) ####将页面分成2*2矩阵
  print(p1, vp = vplayout(1,1))
  print(p2, vp = vplayout(1,2))
  print(p3, vp = vplayout(1,3))
  dev.off()
  
}
plot.combined = function(data, df.result, Z.sims, N.sims, method = "Bayes"){
  # Input: data is the simulated dataset
  # df.result 
  # n is the dataset size; 
  # Z.sims is the sampled z times from posterior distribution
  # N.sims is the sampled beta times from posterior distribution
  
  # caculate EY at combined level
  data = as.data.frame(data)
  n = nrow(data)
  df.result[which(df.result$Group %in% c(1,2)), "TrueEY"] = 0.9
  df.result[which(df.result$Group %in% c(3,4)), "TrueEY"] = 0.7
  df.result[which(df.result$Group %in% c(5,6)), "TrueEY"] = 0.5
  df.result[which(df.result$Group %in% c(7,8)), "TrueEY"] = 0.2
  
  #
  p1.sim = df.result%>% 
    ggplot(aes(x = as.factor(Group), y = EY)) + 
    geom_violin(aes(x = as.factor(Group), y = EY)) + 
    geom_point(aes(x = as.factor(Group), y = TrueEY), shape = 4, color = "#56B4E9") + 
    xlab("Group") +
    ggtitle(paste(method, ".N: ", n, ", Z.sims: ", Z.sims, ", Beta sims: ", N.sims, sep =""))+ 
    ylab("Estimated Probability of PCR-confirmed Postive Result") + 
    stat_summary(fun.data=mean_sdl, 
                 geom="pointrange", color="red", size = 0.2)
  
  png(paste("sim1.",method,".combined.n",n,"_Z",Z.sims,"_b",N.sims,".png", sep= ""),
      res = 100, width=788,height=534)
  grid.newpage()
  print(p1.sim)
  dev.off()
  
}


# n200 - z50 * beta400; boostrap
# n400 - z50 * beta400; boostrap 500
# 1500 - 50
# 10000 - 400
Z.sims = 50
N.sims = 400
df.Bayes.marginal = Data2EY.bayes.marginal(sim.data, Z.sims = Z.sims, N.sims = N.sims)
df.Bayes = Data2EY.bayes(sim.data, Z.sim = Z.sims, N.sims = N.sims)
df.Bootstrap.marginal = Data2EY.bootstrap.marginal(sim.data)
df.Boostrap = Data2EY.bootstrap(sim.data)


# save data.
#df.marginal.n10000.z500.b300 = df.Bayes.marginal
#df.n10000.z500.b300 = df.Bayes
plot.marginal(sim.data, df.Bayes.marginal, Z.sims, N.sims, method = "Bayes")
plot.combined(sim.data, df.Bayes, Z.sims, N.sims, method = "Bayes")
plot.marginal(sim.data, df.Bootstrap.marginal, Z.sims, N.sims, method = "Boots")
plot.combined(sim.data, df.Boostrap, Z.sims, N.sims, method = "Boots")

# Save Data
##------
df.n1500.bayes = df.Bayes
df.n1500.bayes.marginal = df.Bayes.marginal
df.n1500.Boostrap = df.Boostrap
df.n1500.bootstrap.marginal = df.bootstrap.marginal
save.image("../SaveData/n1500.sim1.RData")


## Compare 
df.Bayes = Data2EY.bayes(sim.data.2, Z.sim = 50)
df.Boostrap = Data2EY.bootstrap(sim.data.2, bootstime = 500)

##---------

