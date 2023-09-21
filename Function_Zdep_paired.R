## Function to compute Z_dep statistic and compare its values between two groups 
## of sites paires (a single distance matrix with a factor of group apprtenance)

## -> Input a distance matrix of community composition, a distance matrix of 
## community appartenance and a geographical distance matrix

Z_test_paired <- function(distcom, envidist, group, host, nperm, byhost){

  options(scipen = 999)

  ####
  distcom   = as.dist(distcom)
  envidist  = as.dist(envidist)
  group     = as.dist(group)
  host     = as.dist(host)
  
  ## GLM
  # distcom_vec <- as.vector(distcom)
  # dist_vec    <- as.vector(envidist)
  # group_vec   <- as.numeric(as.vector(group))
  # host_vec   <- as.numeric(as.vector(host))
  # 
  # data <- data.frame(distcom = distcom_vec, geodist = dist_vec, samegroup = group_vec, same_host = host_vec)
  
  data <- data.frame(distcom = as.vector(distcom), 
                     #com_mat_bin = as.vector(com_mat_bin), 
                     geodist = as.vector(envidist), 
                     samegroup = as.factor(as.vector(group)),
                     same_host = as.factor(host))
  
  ## measure distance between plots intra and inter
  same_ER_dist <- data |> 
    filter(samegroup == 1) |> 
    summarize(mean = mean(geodist), max = max(geodist), min = min(geodist))
  diff_ER_dist <- data |> 
    filter(samegroup == 0) |> 
    summarize(mean = mean(geodist), max = max(geodist), min = min(geodist))
  
  ## filter by host if specified
  if(byhost == TRUE) {
    data = data |>
      filter(same_host == 1)#2116 observations
  }
  # keeping only the ones within the same geographic distance boundaries
  data <-  data |>     
  filter(geodist >= diff_ER_dist$min & geodist  <= same_ER_dist$max)
  
  
  # Community 1
  mod1  <- glm(distcom ~ geodist*samegroup, 
               data = data,
               family = quasibinomial(link = "log"))

  int_1   <- exp(mod1$coefficients[1]) # transform intercept with exp but not slope
  int_2   <- exp(mod1$coefficients[1] + mod1$coefficients[3]) # transform intercept with exp but not slope
  slp_1   <- mod1$coefficients[2]
  slp_2   <- mod1$coefficients[2] + mod1$coefficients[4] 
  
#################################################  
  #### Site - Block permutation ####
  nsites = nrow(as.matrix(envidist))
  
  # Convert triangular matrix into square matrix
  xdis.matrix   <- as.matrix(envidist) 
  ydis.matrix   <- as.matrix(distcom)  
  group.matrix  <- as.matrix(group)  
  host.matrix   <- as.matrix(host)
  
  blocks.xdis   <- matrix(rep(0, nsites*(nsites - 1)), ncol = nsites)    # Empty matrix, to store site-block spatial distances 
  blocks.ydis   <- matrix(rep(0, nsites*(nsites - 1)), ncol = nsites)    # Empty matrix, to store site-block community similarities
  blocks.group  <- matrix(rep(0, nsites*(nsites - 1)), ncol = nsites)
  blocks.host   <- matrix(rep(0, nsites*(nsites - 1)), ncol = nsites)
  
  for (i in 1:nsites){
    
    blocks.xdis[,i]  <- xdis.matrix[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.ydis[,i]  <- ydis.matrix[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.group[,i] <- group.matrix[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.host[,i] <- host.matrix[i,][-i]
    }
  
  # RESAMPLING
  
  int_v1 <- vector()
  int_v2 <- vector()
  slp_v1 <- vector()
  slp_v2 <- vector()
  
perm_data = data.frame()

  for (i in 1:nperm){
    
    # Block-resample with replacement: The matrices with site-blocks are resampled ensuring that the same cases are selected in both matrices
    b.i.resample<-sample(x = 1:nsites, size=nsites, replace=T)                           # Index vector for resampling with replacement (DS1)
    
    b.xdis.resample  <- blocks.xdis[, b.i.resample]                                     # Resampling (with replacement) of blocks of spatial distance values (DS1)
    b.ydis.resample  <- blocks.ydis[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS1)
    b.group.resample <- blocks.group[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS2)
    b.host.resample  <- blocks.host[, b.i.resample]  
    
    # Downsizing of resampled site-blocks: To preserve the size of the observed distance/similarity matrices, only the number of cases in the original distance/similarity matrix is randomly selected
    boot.index <- sample(1:length(b.xdis.resample), size = nsites*(nsites - 1)/2)
    
    boot.b.xdis   <- b.xdis.resample[boot.index] # Index vector for downsizing of resampled site-blocks (DS1)
    boot.b.ydis   <- b.ydis.resample[boot.index]                                         # Downsize of site-blocks of community similarity (DS1)
    boot.b.group  <- b.group.resample[boot.index]                                         # Downsize of site-blocks of community similarity (DS2)
    boot.b.host   <- b.host.resample[boot.index]
      
    data_boot <- data.frame(boot.b.xdis, boot.b.ydis, boot.b.group, boot.b.host, perm = rep(paste(i)))
    ## If choose to keep only by host pairs of site, this should be included here
    ## so only these pairwise distances are kept in the multi-dataset to draw the graph
    if(byhost == TRUE) {
      
      data_boot = data_boot|>
        filter(boot.b.host == 1) 
    }
    
    ## keep data within the shared range of geographical distances
    data_boot <- data_boot |>
      filter(boot.b.xdis >= diff_ER_dist$min & boot.b.xdis  <= same_ER_dist$max)
    
    
    ## total gataset to plot the graph
    perm_data <- rbind(perm_data, data_boot)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2 <- glm(boot.b.ydis ~ boot.b.xdis * boot.b.group, 
                data = data_boot,
                family = quasibinomial(link = "log"))
    
    int_v1[i] <- exp(mod2$coefficients[1]) # transform intercept with exp but not slope
    int_v2[i] <- exp(mod2$coefficients[1] + mod2$coefficients[3])
    slp_v1[i] <- mod2$coefficients[2]
    slp_v2[i] <- mod2$coefficients[2] + mod2$coefficients[4]
  }
  param_permut <- data.frame(int_1 = int_v1, int_2 = int_v2,
                             slp_1 = slp_v1, slp_1 = slp_v1)
  ### Z_dep
  test_int <- (int_1 - int_2)/ sqrt((var(int_v1) + var(int_v2)) - 2 * cov(int_v1, int_v2))
  test_slp <- (slp_1 - slp_2)/ sqrt((var(slp_v1) + var(slp_v2)) - 2 * cov(slp_v1, slp_v2))
  
  ## pvalues calculation
  # intercept
  if (test_int<=0) {
    pval_int <-2*pnorm(test_int, mean=0, sd=1, lower.tail =T)
  } else if (test_int>0){
    pval_int<-2*pnorm(test_int, mean=0, sd=1, lower.tail =F)
  }
  
  #slope
  if (test_slp<=0) {
    pval_slp <-2*pnorm(test_slp, mean=0, sd=1, lower.tail =T)
  } else if (test_slp>0){
    pval_slp<-2*pnorm(test_slp, mean=0, sd=1, lower.tail =F)
  }
  
  res <- tibble(var = c("Intercept", "Slope"),
                Est_1 = c(int_1, slp_1),
                Est_2 = c(int_2, slp_2),
                Var_1 = round(c(var(int_v1), var(slp_v1)), 5),
                Var_2 = round(c(var(int_v2), var(slp_v2)), 5),
                Zdep = c(test_int, test_slp), 
                p_value = (c(if(pval_int > 0.001){round(pval_int, 4)} else {"<0.001"}, 
                             if(pval_slp > 0.001){round(pval_slp, 4)} else {"<0.001"})))
  res
  return(list(res, param_permut, perm_data))
}

