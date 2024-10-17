

##Load normalized matrix:
options(stringsAsFactors=F)
options(scipen=999)
library(data.table)
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
registerDoParallel(cores=9)

max_pos <- c(30427671,19698289,23459830,18585056,26975502)
cen_start <- c(11500001, 1100001, 10300001,  1500001,  9000001)
cen_end <-   c(17700000, 7200000, 17300000,  6300000, 16000000)
window <- 2000
chromosome_arms <- matrix(0, nrow=10, ncol=3)
chromosome_arms[,1] <- rep(1:5, each=2)
chromosome_arms[c(1,3,5,7,9),2] <- 1; chromosome_arms[c(1,3,5,7,9),3] <- ceiling(cen_start/window)
chromosome_arms[c(2,4,6,8,10),3] <- ceiling(max_pos/window); chromosome_arms[c(2,4,6,8,10),2] <- ceiling(cen_end/window)

chromosome_arms <- chromosome_arms[c(2,4,5,6,8,10),]


TAD_border_genome <- foreach(line_id=1:nrow(chromosome_arms))%dopar%{

chr <- chromosome_arms[line_id, 1]
fileName <- paste("ICE_pds5_rep1rep2_chr_Chr",chr,"_binsize_2kb",sep="")
normalized_16 <- data.matrix(fread(fileName, header=F, sep='auto', data.table=F))
ROI_left <- chromosome_arms[line_id, 2]; ROI_right <- chromosome_arms[line_id, 3]

normalized_16 <- normalized_16[ROI_left:ROI_right,ROI_left:ROI_right]

#large matrix takes long time to calculate, but normally TADs are small (<500kb), so we break the matrix into serveral overlapping
# smaller matrices, with each having size around 5MB and 2.5 MB for bin size of 10 kb and 5 kb, respectively
#define: query region:
#get bin ID:
sum_sub_mat <- ceiling(dim(normalized_16)[1]/250)
sub_mat_all <- matrix(0, nrow=sum_sub_mat-1, ncol=2)
sub_mat_all[,1] <- seq(1, dim(normalized_16)[1],250)[1:(sum_sub_mat-1)]
sub_mat_all[,2] <- seq(1, dim(normalized_16)[1]+250,250)[3:(sum_sub_mat+1)]
sub_mat_all[sum_sub_mat-1,2] <- dim(normalized_16)[1]

TAD_border_all <- numeric()

for(sub_mat_index in 1:nrow(sub_mat_all)){
  
view_start_bin <- sub_mat_all[sub_mat_index,1]
view_end_bin <- sub_mat_all[sub_mat_index,2]

z_combined=normalized_16[view_start_bin:view_end_bin,view_start_bin:view_end_bin]

# Arrowhead matrix transformation:
Arrowhead <- z_combined
bad_id <- which(apply(Arrowhead,1, "sum")==0) #these bins has been masked during normalization
Arrowhead[,] <- NA
z_combined_2 <- z_combined+2
for(i in 1:(nrow(z_combined)-1)){
  for(d in 1:(nrow(z_combined)-i)){
    if(d>=i){next}else
      if((i+d)%in%bad_id | (i-d)%in%bad_id ){next}else
        Arrowhead[i+d, i] <- Arrowhead[i, i+d] <- (z_combined[i, i-d]-z_combined[i, i+d])/(z_combined[i, i-d]+z_combined[i, i+d])
       # Arrowhead[i+d, i] <- Arrowhead[i, i+d] <- (log(z_combined_2[i, i-d],2)-log(z_combined_2[i, i+d],2))/(log(z_combined_2[i, i-d],2)+log(z_combined_2[i, i+d],2))
      
  }
}

# S matrix to store "domain corner score"
# the score is composed of three parts: difference in sum, difference in sign, total variance

S_mat <- vector(mode="list",length=3)
S_mat[[1]] <- S_mat[[2]] <- S_mat[[3]] <- S_mat[[4]] <- matrix(0, nrow(Arrowhead), nrow(Arrowhead))
# matrix1 stores difference in sum, matrix2 stores difference in sign, matrix3 stores total variance
# matrix4 stores the final score
# To check of entry[a,b] is the domain conner:

pick_up_tri <- function(mat, a, b, part){
  # pick up the lower-left, or the upper-right "triangle" from an arrowhead matrix with given coordinates
  # Note this function only applies to arrow head matrix calculation
  # mat, matrix, with a, b given, which are row_id and col_id
  # part, specify "lower" or "upper"
  # this funcion returns a vector of entries selected.
  
  result <- numeric()
  
  if(part=="upper"){
    row_no <- (1:(b-a)+1:(b-a)%%2)/2
    for(i in 1:(b-a)){
      result <- c(result, mat[a:(a+row_no[i]-1),a+i])
    }
    return(result[!is.na(result)])
  }
  
  if(part=="lower"){
    row_no <- rev((1:(b-a)+1:(b-a)%%2)/2)
    for(i in 1:(b-a)){
      result <- c(result, mat[b:(b-row_no[i]+1),b+i])}
    return(result[!is.na(result)])
  }
}

for(a in 1:(nrow(Arrowhead)-1)){
  for(b in (a+1):(nrow(Arrowhead))){
    if(is.na(Arrowhead[a,b])){next} # close to matrix left border or bad points, skip
    if((2*b-a)>nrow(Arrowhead)){next} # close to matrix right border
    
    #extract pixel from upper triangle "U":
    U <- pick_up_tri(Arrowhead, a,b, "upper")
    #and lower triangle "L":
    L <- pick_up_tri(Arrowhead,a,b, "lower")
    
    S_mat[[1]][a,b] <- (sum(U<0) + sum(L>0))/length(c(U,L)) # "average sign"
    S_mat[[2]][a,b] <- (sum(-U) + sum(L))/length(c(U,L)) # "average"
    S_mat[[3]][a,b] <- var(U) + var(L)
  }
  S_mat[[3]][is.na(S_mat[[3]])] <- 0
}

#S_mat[[4]] <- S_mat[[1]] #seems work best for bin=10kb
S_mat[[4]] <- S_mat[[1]]/max(S_mat[[1]])+S_mat[[2]]/max(S_mat[[2]])-S_mat[[3]]/max(S_mat[[3]]) #Rao et al., 2014

image_filtered <- S_mat[[4]]
# to remove score of points close to the diagonal (we don't consider very small-sized TAD)
# So here we don't consider TAD if is has size less than 40k, which means excluding values 
# from lines within three steps away from the diagnal
cutoff_TAD_size <- 4
for(i in 1:(cutoff_TAD_size-1)){
  image_filtered[(row(image_filtered)+i)==col(image_filtered)] <- 0
}

image_filtered[image_filtered<0.95] <- 0 # an arbiterily set cutoff (works well for rice at bin=5kb)

# identify clusters with 25-connectivity (two layers around). Assume that the real sigal come from the two TAD borders
# is a cluster consisting of at least 6 pixels, then for each retained cluster, the one pixel with 
# the highest score defines TAD borders
image_cluster_index <- image_cluster <- matrix(0, nrow(image_filtered)+4, nrow(image_filtered)+4) #to avoid out-of-boundary errors
image_cluster[3:(nrow(image_filtered)+2),3:(nrow(image_filtered)+2)] <- image_filtered
cluster_id <- 0
for(i in 3:(nrow(image_filtered)+2)){
  for(j in i:(nrow(image_filtered)+2)){
    if(image_cluster[i,j]==0){next}else{
      # new cluster:
      if(all(c(image_cluster_index[(i-2):(i+2),(j-2):(j+2)])==0)){cluster_id <- cluster_id+1; image_cluster_index[i,j] <- cluster_id; next}else{
        
        # find out the smallest index from neighbor, and sign this value to all connected clusters/pixel
        indices <- c(image_cluster_index[(i-2):(i+2),(j-2):(j+2)])
        none_0_index <- indices[indices!=0]
        image_cluster_index[i,j] <- min(none_0_index)
        for(k in none_0_index){
          image_cluster_index[image_cluster_index==k] <- min(none_0_index)
        }
        
      } 
      
    }
    
  }
}

image_cluster_index <- image_cluster_index[3:(nrow(image_filtered)+2),3:(nrow(image_filtered)+2)]
image_cluster_index[image_cluster_index==0] <- NA
image_filtered[image_filtered==0] <- NA
TAD_border <- numeric()
  #if no clusters found, move on to next
  if(all(is.na(image_cluster_index))){next}else{
  for(i in 1:max(image_cluster_index, na.rm=T)){
    if(sum(image_cluster_index==i, na.rm=T)>=6){
      max_score <- max(image_filtered[which(image_cluster_index==i,arr.ind=T)])
      border <- which(image_cluster_index==i&image_filtered==max_score,arr.ind=T)[1,]
     #second filter:
      if(mean(c(image_filtered[border[1],(border[2]-1):border[2]], 
                image_filtered[border[1]:(border[1]+1),border[2]]), na.rm=T)>1.05){TAD_border <- rbind(TAD_border,border)}
    }
  }
  if(length(TAD_border)>0){TAD_border_all <- rbind(TAD_border_all, TAD_border+view_start_bin-1)}
  }
}

TAD_border <- TAD_border_all
TAD_border <- TAD_border+ROI_left-1 #so that we obtain the correct bin position of the chromosome
TAD_border <- cbind(rep(chr, nrow(TAD_border)), TAD_border, rep(1, nrow(TAD_border)))

print(paste("finished calling TADs from chromosome: -- ", chr))
return(TAD_border)
}

TAD_final <- do.call("rbind", TAD_border_genome)
colnames(TAD_final) <- c("chr","start_2k_bin_id", "end_2k_bin_id","flag")
TAD_final <- unique(TAD_final)
write.table(TAD_final, file="TAD_pds5_bin_2kb_0.95_1.05", row.names=F, quote=F, sep="\t")

