
#helper functions 


getROC_prob <- function(biomark, groupID, criteria) {
  ROC = data.frame('crit' = criteria, 
                   'hits' = criteria,
                   'misses'=criteria,
                   'CRs' = criteria, 
                   'FAs' = criteria,
                   'TP' = criteria, 
                   'TN' = criteria,
                   'FP' = criteria, 
                   'FN' = criteria, 
                   'accRaw' = criteria)
  
  AUC = 0
  for(ci in 1:length(criteria)){
    guess = rep(1,length(biomark))
    guess[biomark<criteria[ci]] = 2
    ROC$hits[ci] = sum(guess==2 & groupID==2)
    ROC$misses[ci] = sum(guess==1 & groupID==2)
    ROC$CRs[ci] = sum(guess==1 & groupID==1)
    ROC$FAs[ci] = sum(guess==2 & groupID==1)
    ROC$TP[ci] = sum(guess==2 & groupID==2) / sum(groupID==2)
    ROC$TN[ci] = sum(guess==1 & groupID==1) / sum(groupID==1)
    ROC$FP[ci] = sum(guess==2 & groupID==1) / sum(groupID==1)
    ROC$FN[ci] = sum(guess==1 & groupID==2) / sum(groupID==2)
    ROC$accRaw[ci] = (sum(guess==2 & groupID==2) + sum(guess==1 & groupID==1)) / length(groupID)
    if(ci>1){
      baseDist = ROC$FP[ci] - ROC$FP[ci-1]
      #top triangle area:
      tria = ((ROC$TP[ci] - ROC$TP[ci-1]) * baseDist) / 2
      #base rectangle area:
      reca = ROC$TP[ci-1] * baseDist
      AUC = AUC + tria + reca
    } else {
      baseDist = ROC$FP[ci]
      #top triangle area:
      tria = ((ROC$TP[ci]) * baseDist) / 2
      AUC = AUC + tria
    }
    
  }
  baseDist = 1 - ROC$FP[ci]
  #top triangle area: 
  tria = ((1 - ROC$TP[ci]) * baseDist) / 2
  #base rectangle area: 
  reca = ROC$TP[ci] * baseDist
  AUC = AUC + tria + reca
  ROC$AUC = AUC
  
  
  #sensitivity, true positive rate
  ROC$TPR = ROC$hits/(ROC$hits + ROC$misses)
  #fallout, false positive rate
  ROC$FPR = ROC$FAs/(ROC$FAs + ROC$CRs)
  #specificity, true negative rate
  ROC$TNR = ROC$CRs/(ROC$CRs + ROC$FAs)
  #accuracy
  ROC$acc = (ROC$TP + ROC$TN) / (ROC$TP + ROC$TN + ROC$FP + ROC$FN)
  #what's the optimal criterion? 
  ci = which(ROC$acc == max(ROC$acc))
  if(length(ci)>1){
    ci = ci[1]
  }
  
  baseRate = sum(groupID==2) / length(groupID)
  #bayesian posterior: 
  post = ROC$TP[ci]*baseRate / (ROC$TP[ci]*baseRate + ROC$FP[ci]*(1-baseRate))
  ROC$base = baseRate
  ROC$post = post
  
  #false positive index
  guess = rep(1,length(biomark))
  guess[biomark>criteria[ci]] = 2
  ROC$FPi = sum(guess==2 & groupID==1) / sum(guess==2 & groupID==2)
  #precision, positive predictive value
  ROC$PPV = ROC$hits[ci] / (ROC$hits[ci] + ROC$FAs[ci])
  #negative predictive value
  ROC$NPV = ROC$CRs[ci] / (ROC$CRs[ci] + ROC$misses[ci])
  
  
  
  return(ROC)
}

getResults <- function(dat, varNames, eyeCons, diagGroups){
  #results for this bootstrap sample only 
  results = data.frame('dependentVariable' = rep('A', 4*length(varNames)),
                       'age' = rep(0, 4*length(varNames)),
                       'sex' = rep(0, 4*length(varNames)),
                       'IQ' = rep(0, 4*length(varNames)),
                       'para' = rep(0, 4*length(varNames)),
                       'paranoia_gvif'= rep(0, 4*length(varNames)),
                       'para_p' = rep(0, 4*length(varNames)),
                       'age_b' = rep(0, 4*length(varNames)),
                       'sex_b' = rep(0, 4*length(varNames)),
                       'IQ_b' = rep(0, 4*length(varNames)),
                       'para_b' = rep(0, 4*length(varNames)),
                       'n' = rep(0, 4*length(varNames)),
                       'out' = rep(0, 4*length(varNames)),
                       'outlierFirst' = rep(0, 4*length(varNames)),
                       'eyeCon' = rep('a', 4*length(varNames)),
                       'diagGroup' = rep('a', 4*length(varNames)),
                       'dvMean' = rep(0, 4*length(varNames)),
                       'dvSD' = rep(0, 4*length(varNames)),
                       'nonLin' = rep(0, 4*length(varNames)))
  
  
  
  for(ii in 1:length(eyeCons)){
    for(jj in 1:length(diagGroups)){
      
      
      #down select for eye condition and diag group 
      curDat = dat[dat$group == diagGroups[jj] & dat$eyes==eyeCons[ii], ]
      
      if(diagGroups[jj] == 'CON'){
        IV = curDat$par_total
      } else { 
        IV = curDat$pos_p6
      }
      
      
      for(tt in 1:length(varNames)){
        if(!grepl('PACmi', varNames[tt] ) & #skipping the non z-scored raw PAC values
           !grepl('phase', varNames[tt])) { #skipping phase preferences of PAC
          
          # print(paste(eyeCons[ii], diagGroups[jj], tt))
          #### outlier removal ####
          tmpOut <- getModDat(dat, IV, varNames, tt,  eyeCons, ii, 
                              diagGroups, jj)
          modDat <- tmpOut[[1]]
          dat <- tmpOut[[2]]
          
          
          
          #optional plot for this variable
          # modelPlot(modDat, varNames, tt,  eyeCons, ii,
          #           diagGroups, jj)
          
          
          
          
          modDat$sex <- factor(modDat$sex, levels = c( "F","M"))  
          
          
          
          #fit multiple linear regression
          curLM = lm(dv ~ age + sex + iv  +IQ , data = modDat)
          
          
          
          
          #is a non-linear fit better for age? 
          curLM_nonLin = lm(dv~ age + I(age^2), data = modDat)
          curLM_ageAlone = lm(dv~age, data = modDat)
          
          nonLinTest = lrtest(curLM_ageAlone, curLM_nonLin)
          
          if( nonLinTest$`Pr(>Chisq)`[2] <.05) { #use the non linear if it's better
            curLM = lm(dv ~ age + sex + iv  +IQ + I(age^2), data = modDat)
          }
          
          aovTab = Anova(curLM, type = 3)
          
          #adjust factor for storing results: 
          ai = ((ii-1)*2+(jj-1))*length(varNames)
          
          
          results$nonLin[tt+ai] = nonLinTest$`Pr(>Chisq)`[2] #non lin better? 
          
          sink(file = "nul")
          
          # Call gvif and store its output
          vifTab <- gvif(curLM)  # Ensure to replace 'curLM' with your actual model variable
          
          # Restore console output
          sink()
          
          ri = which(grepl('iv', rownames(vifTab)))
          results$paranoia_gvif[tt+ai] = vifTab[ri,3]
          
          
          results[tt+ai,1] = varNames[tt]
          ei = which(row.names(aovTab) == 'Residuals')
          ri = which(row.names(aovTab) == 'iv')
          results$para[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          results$para_p[tt+ai] = aovTab$`Pr(>F)`[ri] 
          
          ri = which(row.names(aovTab) == 'sex')
          results$sex[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          
          ri = which(row.names(aovTab) == 'IQ')
          results$IQ[tt+ai] = aovTab$`Sum Sq`[ri] / (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          ri = which(row.names(aovTab) == 'age')
          results$age[tt+ai] = aovTab$`Sum Sq`[ri] / 
            (aovTab$`Sum Sq`[ri] + aovTab$`Sum Sq`[ei])
          
          
          ri = which(names(curLM$coefficients) == 'age')
          results$age_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'sexM')
          results$sex_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'IQ')
          results$IQ_b[tt+ai] = curLM$coefficients[[ri]]
          ri = which(names(curLM$coefficients) == 'iv')
          results$para_b[tt+ai] = curLM$coefficients[[ri]]
          
          
          
          results$n[tt+ai] = length(modDat[,1])
          results$out[tt+ai] = length(curDat[,1])-length(modDat$dv)
          results$dvMean[tt+ai] = modDat$dv_mean[1]
          results$eyeCon[tt+ai] = eyeCons[ii]
          results$diagGroup[tt+ai] = diagGroups[jj]
          results$dvSD[tt+ai] = modDat$dv_sd[1]
          
          
        }
        
        
      }
      
    }
  }
  results <- results %>% filter(dependentVariable != 'A')
  
  results$type = 'A'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('pow', x, ignore.case = T))] = 'power'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('slope', x, ignore.case = T))] = 'slope'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('relalpha', x, ignore.case = T) | 
                       grepl('logalpha', x, ignore.case = T))] = 'alpha'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('sampEnt', x, ignore.case = T))] = 'sampEnt'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('ispc', x, ignore.case = T))] = 'ispc'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('PAC', x, ignore.case = T))] = 'PAC'
  results$type[apply(as.matrix(results$dependentVariable), 
                     1, function(x) grepl('phase', x, ignore.case = T))] = 'phase'
  
  
  
  return(list(results, dat))
}



makeDemoTable <- function(dat){
  dataSumTable = data.frame('group' = rep(c('AD', 'ASD', 'CON')), 
                            'data set' = as.vector(apply(as.matrix(unique(dat$dataSet)), 1, function(x) rep(x, 3))),
                            'IQ' = rep('', 15),
                            'IQ metric' = c(rep('DAS GCA', 3), rep('MSEL', 3), rep('WTAR', 3), rep('DAS GCA', 3), rep('WTAR', 3)),
                            'age' = rep('', 15),
                            'n female' = rep('', 15),
                            'n total' = rep('', 15),
                            'orig channels' = rep('', 15),
                            'final channels' = rep('', 15),
                            'orig epochs' = rep('', 15),
                            'final epochs' = rep('', 15))
  
  
  
  for(ii in 1:15) { 
    cur = filter(dat, group == dataSumTable$group[ii], dataSet == dataSumTable$data.set[ii])
    if(length(cur$age) > 0){
      #replace na IQ values
      dat$IQ[dat$group == dataSumTable$group[ii] & dat$dataSet == dataSumTable$data.set[ii] & is.na(dat$IQ)] = mean(cur$IQ, na.rm = T)
      dataSumTable$IQ[ii] = paste(round(mean(cur$IQ, na.rm = T)), ' (', round(sd(cur$IQ, na.rm = T)),'; ', min(cur$IQ), '-', max(cur$IQ), ')', sep='')
      dataSumTable$age[ii] = paste(round(mean(cur$age)), ' (', round(sd(cur$age)),'; ', min(cur$age), '-', max(cur$age), ')', sep='')
      dataSumTable$n.female[ii] = sum(cur$sex == 'F')
      dataSumTable$n.total[ii] = length(cur$sex)
      dataSumTable$orig.channels[ii] = paste(round(mean(cur$nbChanOrig)), ' (', round(sd(cur$nbChanOrig)), ')', sep='')
      dataSumTable$final.channels[ii] = paste(round(mean(cur$nbChanFinal)), ' (', round(sd(cur$nbChanFinal)), ')', sep='')
      dataSumTable$orig.epochs[ii] = paste(round(mean(cur$nbTrialOrig[cur$nbTrialOrig<999])), ' (', round(sd(cur$nbTrialOrig[cur$nbTrialOrig<999])),'; ', min(cur$nbTrialOrig[cur$nbTrialOrig<999]), '-', max(cur$nbTrialOrig[cur$nbTrialOrig<999]), ')', sep='')
      dataSumTable$final.epochs[ii] = paste(round(mean(cur$nbTrialFinal)), ' (', round(sd(cur$nbTrialFinal)),'; ', min(cur$nbTrialFinal), '-', max(cur$nbTrialFinal), ')', sep='')
      
      
    }
  }
  
  write.csv(dataSumTable, "G:\\My Drive\\Milne\\clusterDifferences\\table1.csv")
  
  dataSumTable %>% 
    kbl(align = 'c') %>% 
    kable_classic(full_width = F, 
                  font_size = 20) %>%
    row_spec(1, align = 'c')%>%
    footnote(general = "DAS GCA = Differential Ability Scales General Conceptual Ability
                      MSEL = Mullen Scales of Early Learning
                      WTAR = Weschler Test of Adult Reading
                      biomarkCon = The Autism Biomarkers Consortium for Clinical Trials
                      biomarkDev = Biomarkers of Developmental Trajectories and Treatment in ASD
                      bpSZ = Bipolar & Schizophrenia Consortium for Parsing Intermediate Phenotypes 
                      femaleASD = Multimodal Developmental Neurogenetics of Females with ASD
                      socBrain = The Social Brain in Schizophrenia and Autism Spectrum Disorders
                      Numeric values indicate mean and (standard deviation).",
             general_title = "Table 1: ",
             footnote_as_chunk = T, title_format = c("italic", "underline")
    )  
}


myNMI <- function(set1, set2){
  
  # Get unique cluster IDs and remove -1 (noise or outliers)
  clustIDs <- unique(set1[set1 != 0])
  clustIDs2 <- unique(set2[set2 != 0])
  
  # Calculate the normalized mutual information
  nominSum <- 0
  denomh <- 0
  denoml <- 0
  n <- length(set1)
  
  # Iterate over each cluster in set1
  for (c in clustIDs) {
    vec1 <- as.numeric(set1 == c)
    nh <- sum(vec1)
    
    # Iterate over each cluster in set2
    for (cc in clustIDs2) {
      vec2 <- as.numeric(set2 == cc)
      nl <- sum(vec2)
      vec3 <- as.numeric(vec1 == 1 & vec2 == 1)
      nhl <- sum(vec3)
      
      # Correction for zero overlap
      if (nhl == 0) {
        nhl <- .000000001
      }
      
      nominSum <- nominSum + nhl * log10((n * nhl) / (nh * nl))
    }
    
    denomh <- denomh + nh * log10(nh / n)
  }
  
  # Calculate the denominator for clusters in set2
  for (cc in clustIDs2) {
    vec2 <- as.numeric(set2 == cc)
    nl <- sum(vec2)
    denoml <- denoml + nl * log10(nl / n)
  }
  
  # Calculate NMI
  nmi <- nominSum / sqrt(denomh * denoml)
  
  # Check for values outside [0,1] range
  if (nmi < 0) {
    # warning("Check inputs for cluster scheme of all -1")
    nmi <- 0
  } else if (nmi > 1) {
    # warning("Check inputs for cluster scheme of all -1")
    nmi <- 1
  }
  
  # Print the NMI result
  return(nmi)
}

getModDat <- function(dat, IV, varNames, tt,  eyeCons, ii, 
                      diagGroups, jj){

    groupIDX = (dat$group == diagGroups[jj] & dat$eyes == eyeCons[ii])
    groupIDX = which(groupIDX)
    deleteIDX = !is.na(dat[groupIDX,varNames[tt]] ) & !is.na(IV)
    groupIDX = groupIDX[deleteIDX]
    IV = IV[deleteIDX]
   
    temp = dat[groupIDX,]
    dv = temp[[varNames[tt]]]
    if(is.character(dv)){
      dv <- as.numeric(dv)
    }
    n = length(dv)
    lims = quantile(dv, c(.1, .90))
    dv_mean = mean(dv[dv>lims[1] & dv<lims[2]])
    dv_sd = sd(dv[dv>lims[1] & dv<lims[2]])
    dv_z = (dv - dv_mean) / dv_sd
    nonOut = which(abs(dv_z)<5) #find the values of the dv less than 5 sds in trimmed data
    dv_mean = mean(dv[nonOut])
    dv_sd = sd(dv[nonOut])
    dv_z = (dv - dv_mean) / dv_sd
    
    
    
    #as run for paper: 
    dv = dv[abs(dv_z)<5]
    age = temp$age[abs(dv_z)<5]
    sex = temp$sex[abs(dv_z)<5]
    IQ = temp$IQ[abs(dv_z)<5]
    IV = IV[abs(dv_z)<5]
    #standardizing the dependent variable prior to model calculation
    #it was discovered after data analysis began that many values, particularly
    #PAC values could not easily be fit because values were so small that floating
    #point precision had a problem. To get around this, all variables were z-scored
    #prior to model fitting. Mean and standard deviation are stored so that beta
    #weights can be converted back later if necessary
    dv_mean = mean(dv)
    dv_sd = sd(dv)
    dv = (dv - dv_mean) / dv_sd
    # if(length(keep) != length(dv_z)){
    # print(paste(length(keep), length(dv_z), tt))
    # }
    keepz = abs(dv_z)<5
    
 
  
  dat[groupIDX[!keepz], varNames[tt]] = NA #replace outlier values
  
  
  
  modDat = data.frame('dv' = dv, 'age' = age, 'sex' = factor(sex), 'IQ' = IQ, 'iv' = IV)
  modDat$dv_sd = dv_sd
  modDat$dv_mean = dv_mean
  return(list(modDat, dat))
  
}

modelPlot <- function(modDat, varNames, tt,  eyeCons, ii, 
                      diagGroups, jj){
  ggplot(modDat, aes(x = iv, y = dv*dv_sd + dv_mean)) +
    geom_jitter() +
    theme_classic() +
    theme(axis.line = element_line(color = 'black', size = 3),
          axis.ticks = element_line(colour = "black", size = 2),
          axis.ticks.length=unit(-.25, "cm"),
          text = element_text(size = 20)) +
    ylab(varNames[tt]) +
    xlab('paranoia measure') + 
    ggtitle(paste(varNames[tt], '; diag group: ', 
                  diagGroups[jj], '; eyes: ', eyeCons[ii], sep = ''))
}


randomSplit <- function(dat, prop = .5){
 
  oneHalfIDX = c()
  
  #key variables
  S = dat$sex
  E = dat$eyes
  D = dat$group

  for(si in unique(S)){
    for(ei in unique(E)){
      for(di in unique(D)){
        curIDX = which(S == si & ei == E & di == D)
        oneHalfIDX = c(oneHalfIDX, sample(curIDX, size = round(prop*length(curIDX))))
      }
    }
  }
  allIDX = seq(1,length(dat$name))
  
  heldOut = dat[allIDX[-oneHalfIDX], ]
  dat = dat[oneHalfIDX,]
  return(list(heldOut, dat))
  
}






# 
# dat <- dat %>% arrange(age)
# allIdx = c()
# for(ii in 2:length(ageGroups)){
#   temp = which(dat$age>ageGroups[ii-1] & dat$age<=ageGroups[ii] )
#   tempS = dat$sex[temp]
#   tempG = dat$group[temp]
# 
#   keepers = c()
#   #do removal separately for male and female
#   temp_f = temp[tempS=='F']
#   temp_m = temp[tempS=='M']
#   #bring diagnosis info along
#   temp_fG = tempG[tempS=='F']
#   temp_mG = tempG[tempS=='M']
#   #do removal for each diagnosis
#   for(di in c('AD', 'ASD', 'CON')){
#     curM = temp_m[temp_mG == di]
#     curF = temp_f[temp_fG == di]
#     keepers = c(keepers, sample(curM, round(length(curM)*.5)))
#     keepers = c(keepers, sample(curF, round(length(curF)*.5)))
#   }
#   allIdx = c(allIdx, is.element(temp, keepers))
# 
# }
# heldOut = dat[allIdx == F, ]
# dat = dat[allIdx,]
# write.csv(dat, "C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\trainSet.csv")
# write.csv(heldOut, "C:\\Users\\pc1aod\\Documents\\GitHub\\SheffieldAutismBiomarkers\\testSet.csv")
  


getBalanced <- function(modDat){
  modDat$ID = seq(1,length(modDat$group), 1)
  #go through the model Data and create a new dataFrame that has matched samples
  outDat <- modDat[0,]
  
  ADgroup = modDat %>% filter(group == 'AD')
  ASDgroup = modDat %>% filter(group == 'ASD')
  CONgroup = modDat %>% filter(group == 'CON')
  oi = 1
  for(jj in 1:length(ADgroup$group)){
   
    curAD = ADgroup[jj,]
    sampTmp = CONgroup %>% filter(sex == curAD$sex, 
                                  age>curAD$age-5,
                                  age<curAD$age+5, 
                                  IQ>curAD$IQ-10,
                                  IQ<curAD$IQ+10)
    #check that we have at least one match option: 
    if(length(sampTmp$group)>0){
      outDat[oi,] = sample_n(sampTmp, 1)
      CONgroup <- CONgroup %>% filter(ID != outDat$ID[oi])
      oi = oi+1
      outDat[oi,] = curAD
      oi = oi+1
    }
    
    
  }
  
  for(jj in 1:length(ASDgroup$group)){
    curAD = ASDgroup[jj,]
    sampTmp = CONgroup %>% filter(sex == curAD$sex, 
                                  age>curAD$age-5,
                                  age<curAD$age+5, 
                                  IQ>curAD$IQ-10,
                                  IQ<curAD$IQ+10)
    #check that we have at least one match option: 
    if(length(sampTmp$group)>0){
      outDat[oi,] = sample_n(sampTmp, 1)
      CONgroup <- CONgroup %>% filter(ID != outDat$ID[oi])
      oi = oi+1
      outDat[oi,] = curAD
      oi = oi+1
    }
    
    
  }
  outDat<- outDat[ , !(names(outDat) %in% c('ID'))]
  return(outDat)
  
  
}






getROC <- function(biomark, groupID, criteria, flipLogic = F) {
  ROC = data.frame('crit' = criteria,
                   'hits' = criteria,
                   'misses'=criteria,
                   'CRs' = criteria,
                   'FAs' = criteria,
                   'TP' = criteria,
                   'TN' = criteria,
                   'FP' = criteria,
                   'FN' = criteria,
                   'accRaw' = criteria)
  
  if(flipLogic){
    criteria = rev(criteria)
    ROC$crit = criteria
  }

  AUC = 0
  for(ci in 1:length(criteria)){
    guess = rep(1,length(biomark))
    if(flipLogic){
      guess[biomark<criteria[ci]] = 2
    } else {
      guess[biomark>criteria[ci]] = 2
    }
    ROC$hits[ci] = sum(guess==2 & groupID==2)
    ROC$misses[ci] = sum(guess==1 & groupID==2)
    ROC$CRs[ci] = sum(guess==1 & groupID==1)
    ROC$FAs[ci] = sum(guess==2 & groupID==1)
    ROC$TP[ci] = sum(guess==2 & groupID==2) / sum(groupID==2)
    ROC$TN[ci] = sum(guess==1 & groupID==1) / sum(groupID==1)
    ROC$FP[ci] = sum(guess==2 & groupID==1) / sum(groupID==1)
    ROC$FN[ci] = sum(guess==1 & groupID==2) / sum(groupID==2)
    ROC$accRaw[ci] = (sum(guess==2 & groupID==2) + sum(guess==1 & groupID==1)) / length(groupID)
    if(ci>1){
      baseDist = ROC$FP[ci] - ROC$FP[ci-1]
      #top triangle area:
      tria = ((ROC$TP[ci] - ROC$TP[ci-1]) * baseDist) / 2
      #base rectangle area:
      reca = ROC$TP[ci-1] * baseDist
      AUC = AUC + tria + reca
    } else {
      baseDist = ROC$FP[ci]
      #top triangle area:
      tria = ((ROC$TP[ci]) * baseDist) / 2
      AUC = AUC + tria
    }

  }
  baseDist = 1 - ROC$FP[ci]
  #top triangle area:
  tria = ((1 - ROC$TP[ci]) * baseDist) / 2
  #base rectangle area:
  reca = ROC$TP[ci] * baseDist
  AUC = AUC + tria + reca
  ROC$AUC = AUC


  #sensitivity, true positive rate
  ROC$TPR = ROC$hits/(ROC$hits + ROC$misses)
  #fallout, false positive rate
  ROC$FPR = ROC$FAs/(ROC$FAs + ROC$CRs)
  #specificity, true negative rate
  ROC$TNR = ROC$CRs/(ROC$CRs + ROC$FAs)
  #accuracy
  ROC$acc = (ROC$TP + ROC$TN) / (ROC$TP + ROC$TN + ROC$FP + ROC$FN)
  #what's the optimal criterion?
  ci = which(ROC$acc == max(ROC$acc))
  ci = ci[1]

  baseRate = sum(groupID==2) / length(groupID)
  #bayesian posterior:
  post = ROC$TP[ci]*baseRate / (ROC$TP[ci]*baseRate + ROC$FP[ci]*(1-baseRate))
  ROC$base = baseRate
  ROC$post = post

  #false positive index
  guess = rep(1,length(biomark))
  guess[biomark>criteria[ci]] = 2
  ROC$FPi = sum(guess==2 & groupID==1) / sum(guess==2 & groupID==2)
  #precision, positive predictive value
  ROC$PPV = ROC$hits[ci] / (ROC$hits[ci] + ROC$FAs[ci])
  #negative predictive value
  ROC$NPV = ROC$CRs[ci] / (ROC$CRs[ci] + ROC$misses[ci])



  return(ROC)
}