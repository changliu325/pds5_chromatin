



options(stringsAsFactors=F, scipen=999)
library(data.table)
library(shape)
library("splus2R")

max_pos <- c(30427671,19698289,23459830,18585056,26975502)
cen_start <- c(11500001, 1100001, 10300001,  1500001,  9000001)
cen_end <-   c(17700000, 7200000, 17300000,  6300000, 16000000)

window <- 2000

chromosome_arms <- matrix(0, nrow=10, ncol=3)
chromosome_arms[,1] <- rep(1:5, each=2)
chromosome_arms[c(1,3,5,7,9),2] <- 1; chromosome_arms[c(1,3,5,7,9),3] <- ceiling(cen_start/window)
chromosome_arms[c(2,4,6,8,10),3] <- ceiling(max_pos/window); chromosome_arms[c(2,4,6,8,10),2] <- ceiling(cen_end/window)

#select arms:
chromosome_arms <- chromosome_arms[c(2,4,5,6,8,10),]

insulation_mat <- c()

for(j in 1:nrow(chromosome_arms)){
  chr <- chromosome_arms[j,1]; start_bin <- chromosome_arms[j,2]; end_bin <- chromosome_arms[j,3]
  file_name_WT <- paste("ICE_col0_ctr_rep1rep2_chr_Chr",chr,"_binsize_2kb", sep="")
  file_name_mu <- paste("ICE_pds5_rep1rep2_chr_Chr",chr,"_binsize_2kb", sep="")
  file_name_mu_wapl <- paste("ICE_wapl_double_rep1rep2_chr_Chr",chr,"_binsize_2kb", sep="")
  
  WT <- data.matrix(fread(file_name_WT, header=F, sep='auto', data.table=F))[start_bin:end_bin,start_bin:end_bin]
  mu <- data.matrix(fread(file_name_mu, header=F, sep='auto', data.table=F))[start_bin:end_bin,start_bin:end_bin]
  mu_wapl <- data.matrix(fread(file_name_mu_wapl, header=F, sep='auto', data.table=F))[start_bin:end_bin,start_bin:end_bin]
  
  blocked_bins <- unique(c(which(apply(WT,1,sum)==0), which(apply(mu,1,sum)==0)))
  WT[blocked_bins,] <- NA; WT[,blocked_bins] <- NA
  mu[blocked_bins,] <- NA; mu[,blocked_bins] <- NA
  mu_wapl[blocked_bins,] <- NA; mu_wapl[,blocked_bins] <- NA
  
  look_up <- 25 #number of bins upstream/downstream to check
  HiC_WT_exp <- matrix(NA, nrow=nrow(WT), ncol=ncol(WT)) #expected interaction values
  for(i in 1:(2*look_up-1)){
    band = (row(HiC_WT_exp)==col(HiC_WT_exp)-i)
    HiC_WT_exp[band] = mean(WT[band], na.rm=T)
  }
  
  #calculate insulation scores:
  insulation_scores_WT <- rep(NA, nrow(WT))
  for(i in 1:nrow(WT)){
    if(i<look_up|(i+look_up)>nrow(WT)){next}
    sub_mat <- WT[(i-look_up):(i-1),(i+1):(i+look_up)]/HiC_WT_exp[(i-look_up):(i-1),(i+1):(i+look_up)]
    sub_mat_5 <- WT[(i-look_up):(i-1),(i-look_up):(i-1)]/HiC_WT_exp[(i-look_up):(i-1),(i-look_up):(i-1)]
    sub_mat_3 <- WT[(i+1):(i+look_up),(i+1):(i+look_up)]/HiC_WT_exp[(i+1):(i+look_up),(i+1):(i+look_up)]
    #  insulation_scores_WT[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)
    #if consider 5' and 3' regions:
    insulation_scores_WT[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)/mean(c(sub_mat_5[upper.tri(sub_mat_5)], sub_mat_3[upper.tri(sub_mat_3)]),na.rm=T)
    
  }
  
  HiC_mu_exp <- matrix(NA, nrow=nrow(mu), ncol=ncol(mu)) #expected interaction values
  for(i in 1:(2*look_up-1)){
    band = (row(HiC_mu_exp)==col(HiC_mu_exp)-i)
    HiC_mu_exp[band] = mean(mu[band], na.rm=T)
  }
  
  HiC_wapl_exp <- matrix(NA, nrow=nrow(mu_wapl), ncol=ncol(mu_wapl)) #expected interaction values
  for(i in 1:(2*look_up-1)){
    band = (row(HiC_wapl_exp)==col(HiC_wapl_exp)-i)
    HiC_wapl_exp[band] = mean(mu_wapl[band], na.rm=T)
  }
  
  #calculate insulation scores:
  insulation_scores_mu <- rep(NA, nrow(mu))
  for(i in 1:nrow(mu)){
    if(i<look_up|(i+look_up)>nrow(mu)){next}
    sub_mat <- mu[(i-look_up):(i-1),(i+1):(i+look_up)]/HiC_mu_exp[(i-look_up):(i-1),(i+1):(i+look_up)]
    sub_mat_5 <- mu[(i-look_up):(i-1),(i-look_up):(i-1)]/HiC_mu_exp[(i-look_up):(i-1),(i-look_up):(i-1)]
    sub_mat_3 <- mu[(i+1):(i+look_up),(i+1):(i+look_up)]/HiC_mu_exp[(i+1):(i+look_up),(i+1):(i+look_up)]
    #  insulation_scores_mu[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)
    #if consider 5' and 3' regions:
    insulation_scores_mu[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)/mean(c(sub_mat_5[upper.tri(sub_mat_5)], sub_mat_3[upper.tri(sub_mat_3)]),na.rm=T)
  }
  
  #calculate insulation scores:
  insulation_scores_wapl <- rep(NA, nrow(mu_wapl))
  for(i in 1:nrow(mu_wapl)){
    if(i<look_up|(i+look_up)>nrow(mu_wapl)){next}
    sub_mat <- mu_wapl[(i-look_up):(i-1),(i+1):(i+look_up)]/HiC_wapl_exp[(i-look_up):(i-1),(i+1):(i+look_up)]
    sub_mat_5 <- mu_wapl[(i-look_up):(i-1),(i-look_up):(i-1)]/HiC_wapl_exp[(i-look_up):(i-1),(i-look_up):(i-1)]
    sub_mat_3 <- mu_wapl[(i+1):(i+look_up),(i+1):(i+look_up)]/HiC_wapl_exp[(i+1):(i+look_up),(i+1):(i+look_up)]
    #  insulation_scores_mu[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)
    #if consider 5' and 3' regions:
    insulation_scores_wapl[i] <- mean(data.matrix(sub_mat[!upper.tri(sub_mat)]), na.rm=T)/mean(c(sub_mat_5[upper.tri(sub_mat_5)], sub_mat_3[upper.tri(sub_mat_3)]),na.rm=T)
  }
  
  local_min_WT <- peaks(-log(insulation_scores_WT,2), span=10, strict=TRUE, endbehavior=0)
  local_min_mu <- peaks(-log(insulation_scores_mu,2), span=10, strict=TRUE, endbehavior=0)
  local_min_wapl <- peaks(-log(insulation_scores_wapl,2), span=10, strict=TRUE, endbehavior=0)
  
  cut_off <- -0.25 # only consider regions with strong insulation
  local_min_mu[log(insulation_scores_mu,2)>cut_off] <- F
  local_min_WT[log(insulation_scores_WT,2)>cut_off] <- F
  local_min_wapl[log(insulation_scores_wapl,2)>cut_off] <- F
  
  insu_mat <- data.frame(chr=chr,bin_5k=start_bin:end_bin, log2WT_score=log(insulation_scores_WT,2),
                         log2mu_score=log(insulation_scores_mu,2), log2wapl_score=log(insulation_scores_wapl,2), WT_local_min=as.numeric(local_min_WT), 
                         mu_local_min=as.numeric(local_min_mu), wapl_local_min=as.numeric(local_min_wapl))
  insulation_mat <- rbind(insulation_mat, insu_mat)
}  

write.table(insulation_mat, file = "local_insulation_2kb_lookup_50kb.txt", row.names=F, quote=F, sep="\t")



