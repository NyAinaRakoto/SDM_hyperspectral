#########Library 
#Readapted Based on spectrait function and 
#https://github.com/Antguz/PLSR_leaf-traits/blob/main/data_processing/optimal_number.R
# Also readapted from the guide of Burnett, A. C., Anderson, J., Davidson, K. J., Ely, K. S.,
#Lamour, J., Li, Q., ... & Serbin, S. P. (2021). 
#A best-practice guide to predicting plant traits from leaf-level hyperspectral 
#data using partial least squares regression. Journal of Experimental Botany, 72(18), 6175-6189.


library(spectrolab)
library(parallel)
library(readxl)
#Read all the data 
#Read spectra
dir <- "directory"
airref_51 <- read.csv(paste(dir, 
                            "/Mean_Reflectance.csv", sep=""), check.names = FALSE)

airref_spec = as_spectra(airref_51[,-1])/10000
dim(airref_spec)

plot(airref_spec, lwd = 0.75, lty = 1, col = "grey25", main = "All Spectra")

# Make a matrix from a `spectra` object
spec_as_mat = as.matrix(airref_spec, fix_names = "none")

# set the desired spectra wavelength range to include
wv <- colnames(airref_51)[-1]
colnames(spec_as_mat) <- c(paste0("Wave_",colnames(airref_51)[-1]))
rownames(spec_as_mat) <- airref_51$X

###read the traits ####
dir2 <- "D:/Sericea_Data/Tallgrass2022/Impact_sericea/Traits_Inventory"

CWMtraits <- data.frame(read.csv(paste(dir2,"/NPK_2022_100plots.csv", sep="")))#IN CONCENTRATION
summary(CWMtraits)

rownames(CWMtraits) <- CWMtraits$Quadrat
CWMtraits <- CWMtraits[,-1]

#Join the dataset traits and reflectance
CWM_spec <- merge(CWMtraits,spec_as_mat, by = "row.names")

# Make a matrix out of the reflectance alone, not with the traits
CWM_spec_only <- as.matrix(CWM_spec[,-c(1:4)]) 
rownames(CWM_spec_only) <- CWM_spec$Row.names

"%notin%" <- Negate(`%in%`)

# Create a data frame of the data to be used for the PLSR analysis
plsr.data <- data.frame(CWM_spec[, which(names(CWM_spec) %notin% paste0("Wave_",wv))],
                        Spectra=I(CWM_spec_only))

###############

##### TOOOOOO RUNNNNN The optimal function check 

####################################### function to find best number of component
optimal_number <- function(model, ref_train, int_ncomp =10, iterations = 30, length.seg = 10, threads = 12) {
  
  #Frame of performance
  frame <- data.frame(Spectra = NA, 
                      Iteration = NA, 
                      Abs_min = NA, 
                      Optimal = NA)
  
  frame <- frame[0,]
  final <- frame
  
  #RMSE and RMSEP
  fRMSEP <- matrix(NA, ncol = (int_ncomp + 3), nrow = 1)
  names <- c("Spectra", "Iteration", "Intercept", paste0("Comp_", 1:int_ncomp))
  fRMSEP <- as.data.frame(fRMSEP)
  colnames(fRMSEP) <- names
  fRMSEP <- fRMSEP[0,]
  final_fRMSEP <- fRMSEP
  
  fPRESS <- matrix(NA, ncol = (int_ncomp + 2), nrow = 1)
  names <- c("Spectra", "Iteration", paste0("Comp_", 1:int_ncomp))
  fPRESS <- as.data.frame(fPRESS)
  colnames(fPRESS) <- names
  fPRESS <- fPRESS[0,]
  final_fPRESS <- fPRESS
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = iterations, style = 3) 
  
  #Parallel 
  cl <- makeCluster(threads)
  pls.options(parallel = cl)
  
  #Loop
  for(i in 1:iterations) {
    
    setTxtProgressBar(pb, i)
    
    ###Frames copy
    frame_iteration <- frame
    fRMSEP_iteration <- fRMSEP
    fPRESS_iteration <- fPRESS
    
    ###Cross validation
    #CVseg <- CVSeg_lifeforms(train_ID$Life_form,  length.seg =  length.seg)
    
    ###---------------------------------Reflectance
    #Creation of the model
    pls_components_ref <- plsr(model, 
                               data= ref_train, 
                               # scale = TRUE,
                               # center = TRUE,
                               ncomp= int_ncomp,
                               validation = "CV", 
                               # segments = CVseg,
                               jackknife= TRUE,
                               trace= FALSE, 
                               method = "oscorespls")
    
    #Spectra
    frame_iteration[1, 1] <- "Reflectance"
    
    #Iteration
    frame_iteration[1, 2] <- i
    
    ###Selection of n components
    #Min absolute
    frame_iteration[1, 3] <- which.min(as.vector(pls_components_ref$validation$PRESS))
    
    #Onesigma
    frame_iteration[1, 4] <- selectNcomp(pls_components_ref, "onesigma", plot=TRUE)
    
    ###RMSEP
    vec <- as.vector(RMSEP(pls_components_ref)$val[1,,])
    
    fRMSEP_iteration[1, 1] <- "Reflectance"
    fRMSEP_iteration[1, 2] <- i
    fRMSEP_iteration[1, 3:(int_ncomp + 3)] <- vec
    
    ###PRESS
    vec <- as.vector(pls_components_ref$validation$PRESS)
    
    fPRESS_iteration[1, 1] <- "Reflectance"
    fPRESS_iteration[1, 2] <- i
    fPRESS_iteration[1, 3:(int_ncomp + 2)] <- vec
    
    ###---------------------------------Save results
    
    final <- rbind(final, frame_iteration)
    final_fRMSEP <- rbind(final_fRMSEP,fRMSEP_iteration)
    final_fPRESS <- rbind(final_fPRESS, fPRESS_iteration)
    
  }
  
  stopCluster(cl)
  
  return(list(frame = final, RMSEP = final_fRMSEP, PRESS = final_fPRESS))
}

###########Run opt number
names(plsr.data)

# Name the traits of interest 
inVar <- "Height" # example Height 

library(pls)

# Find the number of components
a <- optimal_number(as.formula(paste(inVar,"~","Spectra")),
                    ref_train=plsr.data
                    , int_ncomp = 15, iterations = 30,
                    length.seg = 10, threads = 2)

# Check PRESS statistics
boxplot(a$PRESS[,-c(1,2)], xlab = "Number of components",
        ylab = "PRESS",
        col = "burlywood1",
        border = "brown",
        horizontal = F,
        notch = F, cex.lab=1.5)
par(cex.axis= 1)
par(cex.lab= 1)

# write down the chosen number of component here
nComps <- 5 # change according to results from PRESS

##############PLS permutation  function
###########################
pls_permutation <- function(model, dataset, ncomp, iterations, prop,
                            verbose, outdir) {
  coefs <- array(0,dim=c(239,iterations,ncomp)) # This number (239) needs to be changed depending on the number of bands, 
  # In this example, 239 is the band numbers (238) + 1 extra column
  press.out <- array(data=NA, dim=c(iterations,ncomp))
  print("*** Running permutation test.  Please hang tight, this can take awhile ***")
  print("Options:")
  print(paste("Max Components:",ncomp, "Iterations:", iterations, 
              "Data Proportion (percent):", prop*100, sep=" "))
  expr <- data.frame(matrix(nrow=iterations, ncol= 2))
  val_output <- data.frame(matrix(ncol=iterations, nrow= length(plsr.data[,inVar])))
  exprdata <- data.frame(matrix(nrow=iterations, ncol= 2))
  valdata_output <- data.frame(matrix(ncol=iterations, nrow= 40))#40))
  
  if (verbose) {
    j <- 1 # <--- Numeric counter for progress bar
    pb <- utils::txtProgressBar(min = 0, max = iterations, 
                                char="*",width=70,style = 3)
  }
  
  for (i in seq_along(1:iterations)) {
    rows <- sample(1:nrow(dataset),floor(prop*nrow(dataset)))
    sub.data <- dataset[rows,]
    val.sub.data <- dataset[-rows,]
    plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),
                     # scale=TRUE,
                     # center= TRUE, 
                     ncomp=nComps,
                     trace=FALSE,
                     data= sub.data)
    saveRDS(plsr.out, paste(outdir, "saveRDS/", inVar, "_model_", i, sep=""))
    pred_val <- predict(plsr.out,newdata=val.sub.data)
    valdata.plsr.output <- data.frame(rownames(pred_val), val.sub.data[,inVar],
                                      PLSR_Predicted=as.vector(predict(plsr.out,
                                                                       newdata = val.sub.data,
                                                                       ncomp=nComps, type="response")[,,1]))
    valdata_output[,i] <- valdata.plsr.output$PLSR_Predicted
    
    write.csv(valdata.plsr.output,paste(outdir, "valdata/", inVar,"_", i, "_estimated_alldata.csv", sep=""))
    
    valdata.R2 <- RSQUARE(valdata.plsr.output$val.sub.data...inVar., valdata.plsr.output$PLSR_Predicted)
    valdata.RMSEP <- round(pls::RMSEP(plsr.out, newdata=val.sub.data , intercept=F)[[1]][nComps],2)
    
    exprdata[i,1] <- valdata.R2
    exprdata[i,2] <- valdata.RMSEP
    
    # create a new list with PRESS and permuted coefficients x wavelength x component number
    # print("*** Providing PRESS and coefficient array output ***")
    # output <- list(PRESS=press.out, coef_array=coefs)
    
    coefs[,i,] <- coef(plsr.out,  ncomp = nComps,intercept = TRUE)
    
    #Final prediction for the whole dataset 
    val.plsr.output <- data.frame(plsr.data[,inVar],
                                  PLSR_Predicted=as.vector(predict(plsr.out,
                                                                   newdata = plsr.data,
                                                                   ncomp=nComps, type="response")[,,1]))
    val_output[,i] <- val.plsr.output$PLSR_Predicted
    
    
    val.R2 <- RSQUARE(val.plsr.output$plsr.data...inVar.,val.plsr.output$PLSR_Predicted)
    val.RMSEP <- round(pls::RMSEP(plsr.out, newdata=dataset, intercept=F)[[1]][nComps],2)
    
    expr[i,1] <- val.R2
    expr[i,2] <- val.RMSEP
    
    ### Display progress to console
    if (verbose) {
      setTxtProgressBar(pb, j)    # show progress bar
      j <- j+1                    # <--- increase counter by 1
      flush.console()             #<--- show output in real-time
    }
  }
  if (verbose) {
    close(pb)
  }
  colnames(expr) <- c("R2", "RMSEP")
  write.csv(expr,paste(outdir, inVar, "_exp_metric_all.csv", sep=""))
  val_output$field <- val.plsr.output[,1]
  write.csv(val_output, paste(outdir, inVar, "_val_plsr_output", i, ".csv", sep=""))
  
  colnames(exprdata) <- c("R2", "RMSEP")
  write.csv(exprdata,paste(outdir, inVar, "_valdata_exp_metric_all.csv", sep=""))
  valdata_output$field <- valdata.plsr.output[,1]
  write.csv(valdata_output, paste(outdir, inVar, "_valdata_val_plsr_output", i, ".csv", sep=""))
  
  output <- list(expr, coef_array=coefs)
  return(output)
  
}

###############################
#########Fit the model now 
RSQUARE = function(y_actual,y_predict){
  cor(y_actual,y_predict)^2
}

nComps <-5
names(plsr.data)
inVar <- "Height"

# update the output directory 
dir3 <- "output_directory"

# create some needed output folders
dir3 <- "C:/Documents"
setwd("dir3")
dir.create("saveRDS", recursive = TRUE)
dir.create("valdata", recursive = TRUE)

# Run pls permutation now
plsr_perm <- pls_permutation(as.formula(paste(inVar,"~","Spectra")), 
                             plsr.data, ncomp=nComps, 50, 0.6,
                             FALSE,
                             # "C:/Users/nrakoto/Desktop/Ny/New_run_cal_val_2021/")
                             dir3)

# get bootstrap intercept and coefficient values
bootstrap_intercept <- plsr_perm$coef_array[1,,nComps]
bootstrap_coef <- plsr_perm$coef_array[2:length(plsr_perm$coef_array[,1,nComps]),
                                       ,nComps]
###########Output data
# Bootstrap Coefficients
out.jk.coefs <- data.frame(Iteration=seq(1,length(bootstrap_intercept),1),
                           Intercept=bootstrap_intercept,t(bootstrap_coef))
head(out.jk.coefs)[1:6]

outdir <- dir3
write.csv(out.jk.coefs,file=file.path(outdir,paste0(inVar,
                                                    '_Bootstrap_PLSR_Coefficients.csv')),
          row.names=FALSE)

#Bootstrap plot
##############
#function 
f.plot.coef <- function(
    Z,                  ## Coefficient matrix with each row corresponding to the coefficients and wavelength in columns
    wv,                 ## vector of wavelengths 
    xlim=NULL,          ## vector to change the default xlim of the plots (ex xlim = c(500, 2400))
    position,## Position of the legend (see base function legend for help)
    type='Coefficient', ## Name of the y axis and of the legend
    plot_label=NULL     ## optional label for plot
){
  
  if(is.null(xlim)){xlim=c(min(wv),max(wv))}
  mean_spec <- colMeans(Z)
  spectra_quantiles <- apply(Z,2,quantile,na.rm=T,probs=c(0,0.01,0.025,0.05,0.5,0.95,0.975,0.99,1))
  
  plot(x=NULL,y=NULL,xlim=xlim,ylim=c(min(Z),max(Z)),xlab="Wavelength (nm)",
       ylab=type,main=plot_label)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[9,], rev(spectra_quantiles[1,])),
          col="grey60",border=NA)
  polygon(c(wv ,rev(wv)),c(spectra_quantiles[6,], rev(spectra_quantiles[4,])),
          col="#99CC99",border=NA)
  lines(wv,mean_spec,lwd=2, lty=1, col="black")
  lines(wv,spectra_quantiles[1,], lty=3, col="grey60")
  lines(wv,spectra_quantiles[9,], lty=3, col="grey60")
  legend(position,legend=c(paste("Mean",type),"Min/Max (range)", "95% CI"),lty=c(1,1,1),
         lwd=c(2,10,10),col=c("black","grey50","#99CC99"),bty="n")
  box(lwd=2.2)
}

#################

f.plot.coef(Z = t(bootstrap_coef), wv = as.numeric(wv),
            plot_label="Bootstrap regression coefficients",position = 'topright')
abline(h=0,lty=2,col="grey50")
box(lwd=2.2)

setwd(outdir)

##########################Prepare for plot of predicted vs. observed traits
#########Read all files
names(plsr.data)

inVar <- "Height"
Var <- "Height"
units <- "(cm)"

##############Read metrics for the validation 
metrics <- read.csv(paste(inVar, "_valdata_exp_metric_all.csv", sep=""))
metrics_mean <- colMeans(metrics)
metrics_sd <- apply(metrics,2,sd)

######### plot now 
#Read all plot from the valdata
dr <- paste(dir3, "valdata/", sep="")
temp_csv <- list.files( dr, pattern=inVar) 
names_csv= lapply(paste(dr, "/", temp_csv, sep=""), read.csv)

df_all <- Reduce(function(...) merge(..., all = TRUE, by="rownames.pred_val."), names_csv)
library(dplyr) 

# Get the data, it has 50 iterations
seq2 <- seq(2, 151, 3)
df_all_n1 <- df_all[, -c(seq2)]
names(df_all_n1)
seq3 <- seq(2, 101, 2) 
df_all_n2 <- df_all_n1[, -c(seq3)]
names(df_all_n2)
pred_val <- rowMeans(df_all_n2[,-1], na.rm=TRUE)
pred_val_sd <- apply(df_all_n2[,-1], 1, sd, na.rm=TRUE)
pred_val_only <- data.frame(df_all_n2$rownames.pred_val., pred_val, pred_val_sd)
colnames(pred_val_only) <- c("rownames", "pred", "sd")

orig_data <- data.frame(plsr.data[,1], plsr.data[,inVar])
colnames(orig_data) <- c('rownames', "field")
pred_field_val <- merge(pred_val_only, orig_data, by="rownames")

# Save the data from iterations
write.csv(pred_field_val, paste(inVar, "_predicted_validation_toplot.csv", sep=""))


pred_field_val <- read.csv(paste(inVar, "_predicted_validation_toplot.csv", sep=""))

#######################PLOT FOR JUST THE VALIDATION 
r2 <- paste0("R\u00B2 = ",round(metrics_mean[2],2), " (", round(metrics_sd[2],2), ")")
rmse <- paste0("RMSEP = ",round(metrics_mean[3],2), " (", round(metrics_sd[3],2), ")")
library(ggplot2)
ggplot(pred_field_val , aes(y=pred, x=field)) +theme_classic() +
  geom_point(colour = "black", size = 3) + 
  geom_errorbar(aes(ymin=pred-sd, ymax=pred+sd), width=.2,
                position=position_dodge(0.05), 
                color="dark grey")+
  geom_abline(intercept = 0, slope = 1, color="dark green",
              linetype="dashed", size=1.5) +
  labs(y=paste0("Predicted ", paste(Var), " ", units),
       x=paste0("Observed ", paste(Var), " ", units))+
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"),
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))+
  ggtitle("Observed vs. predicted validation dataset")+
  #change the following annotated depending on the x-axis and y-axis value range
  annotate("text", x = 50, y=c(120, 110), label = c(r2,rmse), size=5) 