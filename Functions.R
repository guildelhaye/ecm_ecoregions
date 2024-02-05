## Function to compute distance decay model and model fit.
## -> Input a distance matrix of community composition, a distance matrix of 
## community appartenance and a geographical distance matrix
dd_mod_perm <- function(distcom, envidist, host, nperm, byhost){
  set.seed(666)
  options(scipen = 999)
  
  ####
  distcom   = as.dist(distcom)
  envidist  = as.dist(envidist)
  host     = as.dist(host)
  
  data <- data.frame(distcom = as.vector(distcom), 
                     geodist = as.vector(envidist), 
                     same_host = as.factor(host))
  
  ## filter by host if specified
  if(byhost == TRUE) {
    data = data |>
      filter(same_host == 1)
  }
  
  # null model 
  null.x <- rep(1, nrow(data))
  null.nlm <- nlsLM(distcom ~ a * (null.x), 
                    start = list(a = 1), 
                    data = data)
  
  null.deviance <- deviance(null.nlm)
  
  # Model of distance decay
  mod1.1 <- nlsLM(distcom ~ a * exp(b * geodist), 
                  data = data, 
                  start = c(a = 1, b = 0))
  
  int_1 <- summary(mod1.1)$parameters[1, 1]
  slp_1 <- summary(mod1.1)$parameters[2, 1]
  
  pseudo.r.squared <- 1 - deviance(mod1.1)/null.deviance
  
  #################################################  
  #### Site - Block permutation ####
  # RESAMPLING
  
  int_v1 <- vector()
  slp_v1 <- vector()
  
  perm_data = data.frame()
  perm.dev.R2.distr <- numeric(nperm)
  data_boot <- data
  
  for (i in 1:nperm){
    
    data_boot$distcom <- sample(data_boot$distcom)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2.1 <- nlsLM(distcom ~ a * exp(b * geodist), 
                    data = data_boot, 
                    start = c(a = 1, b = 0))
    
    int_v1[i] <- summary(mod2.1)$parameters[1, 1]
    slp_v1[i] <- summary(mod2.1)$parameters[2, 1]
    
    # R2 distribution
    perm.dev.R2.distr[i] <- 1 - deviance(mod2.1)/null.deviance
    
    if ((i / 50)  %in% 1:50) {print(paste0("Permutation = ", i))}
  }
  param_permut <- data.frame(int_1 = int_v1, slp_1 = slp_v1)
  
  p.value <- mean(perm.dev.R2.distr > pseudo.r.squared)
  
  parameters <- list(first.parameter = int_1, 
                     second.parameter = slp_1,
                     pseudo.r.squared = pseudo.r.squared, 
                     p.value = ifelse(p.value == 0, 1/nperm, p.value))
  parameters
  return(list(parameters, param_permut, perm_data, data))
}

# Function to compute Z_dep statistic and compare its values between two groups 
## of sites paires (a single distance matrix with a factor of group apprtenance)

## -> Input a distance matrix of community composition, a distance matrix of 
## community appartenance and a geographical distance matrix

Z_test_paired <- function(distcom, envidist, group, host, nperm, byhost){
  set.seed(666)
  options(scipen = 999)
  
  ####
  distcom   = as.dist(distcom)
  envidist  = as.dist(envidist)
  group     = as.dist(group)
  host     = as.dist(host)
  
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
  
  # Model of distance decay
  mod1.1 <- nlsLM(distcom ~ a * exp(b * geodist), 
                  data = data |> filter(samegroup == 0), 
                  start = c(a = 1, b = 0))
  
  int_1 <- summary(mod1.1)$parameters[1, 1]
  slp_1 <- summary(mod1.1)$parameters[2, 1]
  
  mod1.2 <- nlsLM(distcom ~ a * exp(b * geodist), 
                  data = data |> filter(samegroup == 1), 
                  start = c(a = 1, b = 0))
  int_2 <- summary(mod1.2)$parameters[1, 1]
  slp_2 <- summary(mod1.2)$parameters[2, 1]
  
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
    mod2.1 <- nlsLM(boot.b.ydis ~ a * exp(b * boot.b.xdis), 
                    data = data_boot |> filter(boot.b.group == 0), 
                    start = c(a = 1, b = 0))
    
    mod2.2 <- nlsLM(boot.b.ydis ~ a * exp(b * boot.b.xdis), 
                    data = data_boot |> filter(boot.b.group == 1), 
                    start = c(a = 1, b = 0))
    
    int_v1[i] <- summary(mod2.1)$parameters[1, 1]
    int_v2[i] <- summary(mod2.2)$parameters[1, 1]
    slp_v1[i] <- summary(mod2.1)$parameters[2, 1]
    slp_v2[i] <- summary(mod2.2)$parameters[2, 1]
    
    if ((i / 50)  %in% 1:50) {print(paste0("Permutation = ", i))}
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
  
  res <- data.frame(var = c("Intercept", "Slope"),
                    Est_1 = round(c(int_1, slp_1), 2),
                    Est_2 = round(c(int_2, slp_2), 2),
                    sd_1 = round(c(sqrt(var(int_v1)), sqrt(var(slp_v1))), 4),
                    sd_2 = round(c(sqrt(var(int_v2)), sqrt(var(slp_v2))), 4),
                    Zdep = round(c(test_int, test_slp), 2), 
                    p_value = c(ifelse(pval_int == 0, 1/nperm, pval_int), 
                                ifelse(pval_slp == 0, 1/nperm, pval_slp))
  )
  res
  return(list(res, param_permut, perm_data))
}

#### function to compare the parameters of two distance decay curve
#### following a negative-exponential model (based on Martin-Devasa et al 2022)

#distcom1 = distance matrix for community 1
#distcom2 = distance matrix for community 2
#envidist = distance matrix between plots
#nperm = number of permutation to calculate variance in parameters

Z_test <- function(distcom1, distcom2, 
                   envidist1, envidist2, 
                   nperm){
  set.seed(666)
  options(scipen=999)
  distcom1 = as.dist(distcom1)
  distcom2 = as.dist(distcom2)
  envidist1 = as.dist(envidist1)
  envidist2 = as.dist(envidist2)
  
  ## GLM
  dist1_vec <- as.vector(envidist1)
  dist2_vec <- as.vector(envidist2)
  distcom1_vec <- as.vector(distcom1)
  distcom2_vec <- as.vector(distcom2)
  
  # ## Using method Millar 2011 -> abandonned 
  # # Community 1
  # mod1.1 <- glm(distcom1_vec ~ dist1_vec, 
  #               family = quasibinomial(link = "log"))
  # int_1 <- exp(mod1.1$coefficients[1]) # transform intercept with exp but not slope
  # slp_1 <- mod1.1$coefficients[2]
  # 
  # # Community 2
  # mod1.2 <- glm(distcom2_vec ~ dist2_vec, 
  #               family = quasibinomial(link = "log"))
  # int_2 <- exp(mod1.2$coefficients[1]) # transform intercept with exp but not slope
  # slp_2 <- mod1.2$coefficients[2]
  
  # Using nlsLM (cf betapart package)
  # Community 1
  mod1.1 <- nlsLM(distcom1_vec ~ a * exp(b*dist1_vec), 
                  start = c(a = 1, b = 0))
  
  int_1 <- summary(mod1.1)$parameters[1, 1]
  slp_1 <- summary(mod1.1)$parameters[2, 1]
  
  # Community 2
  mod1.2 <- nlsLM(distcom2_vec ~ a * exp(b*dist2_vec), 
                  start = c(a = 1, b = 0))
  
  int_2 <- summary(mod1.2)$parameters[1, 1]
  slp_2 <- summary(mod1.2)$parameters[2, 1]
  
  ## Site - Block permutation
  # forcommunity 1
  nsites1 = nrow(as.matrix(envidist1))
  xdis.matrix1 <- as.matrix(envidist1)# Convert triangular matrix into square matrix
  ydis.matrix1 <- as.matrix(distcom1)  
  
  # Convert triangular matrix into square matrix
  
  blocks.xdis1  <- matrix(rep(0, nsites1*(nsites1 - 1)), ncol = nsites1)    # Empty matrix, to store site-block spatial distances 
  blocks.ydis1  <- matrix(rep(0, nsites1*(nsites1 - 1)), ncol = nsites1)    # Empty matrix, to store site-block community similarities
  
  for (i in 1:nsites1){
    
    blocks.xdis1[,i]  <- xdis.matrix1[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.ydis1[,i]  <- ydis.matrix1[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
  }
  
  # RESAMPLING
  
  int_v1 <- vector()
  slp_v1 <- vector()
  
  for (i in 1:nperm){
    
    # Block-resample with replacement: The matrices with site-blocks are resampled ensuring that the same cases are selected in both matrices
    b.i.resample <- sample(x = 1:nsites1, size=nsites1, replace=T)                           # Index vector for resampling with replacement (DS1)
    
    b.xdis.resample1 <- blocks.xdis1[, b.i.resample]                                     # Resampling (with replacement) of blocks of spatial distance values (DS1)
    b.ydis.resample1 <- blocks.ydis1[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS1)
    
    # Downsizing of resampled site-blocks: To preserve the size of the observed distance/similarity matrices, only the number of cases in the original distance/similarity matrix is randomly selected
    boot.index <- sample(1:length(b.xdis.resample1), 
                         size = nsites1*(nsites1 - 1)/2)
    
    boot.b.xdis1 <- b.xdis.resample1[boot.index] # Index vector for downsizing of resampled site-blocks (DS1)
    boot.b.ydis1 <- b.ydis.resample1[boot.index]                                         # Downsize of site-blocks of community similarity (DS1)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2.1 <- nlsLM(boot.b.ydis1 ~ a * exp(b*boot.b.xdis1), 
                    start = c(a = 1, b = 0))
    
    int_v1[i] <- summary(mod2.1)$parameters[1, 1]
    slp_v1[i] <- summary(mod2.1)$parameters[2, 1]
    
  }
  
  # forcommunity 2
  nsites2      <- nrow(as.matrix(envidist2))
  xdis.matrix2 <- as.matrix(envidist2)# Convert triangular matrix into square matrix
  ydis.matrix2 <- as.matrix(distcom2)  
  
  # Convert triangular matrix into square matrix
  
  blocks.xdis2  <- matrix(rep(0, nsites2*(nsites2 - 1)), ncol = nsites2)    # Empty matrix, to store site-block spatial distances 
  blocks.ydis2  <- matrix(rep(0, nsites2*(nsites2 - 1)), ncol = nsites2)    # Empty matrix, to store site-block community similarities
  
  for (i in 1:nsites2){
    
    blocks.xdis2[,i]  <- xdis.matrix2[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.ydis2[,i]  <- ydis.matrix2[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
  }
  
  # RESAMPLING
  
  int_v2 <- vector()
  slp_v2 <- vector()
  
  for (i in 1:nperm){
    
    # Block-resample with replacement: The matrices with site-blocks are resampled ensuring that the same cases are selected in both matrices
    b.i.resample<-sample(x = 1:nsites2, 
                         size=nsites2,
                         replace=T)  # Index vector for resampling with replacement (DS1)
    
    b.xdis.resample2 <- blocks.xdis2[, b.i.resample]                                     # Resampling (with replacement) of blocks of spatial distance values (DS1)
    b.ydis.resample2 <- blocks.ydis2[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS1)
    
    # Downsizing of resampled site-blocks: To preserve the size of the observed distance/similarity matrices, only the number of cases in the original distance/similarity matrix is randomly selected
    boot.index <- sample(1:length(b.xdis.resample2), 
                         size = nsites2*(nsites2 - 1)/2)
    
    boot.b.xdis2 <- b.xdis.resample2[boot.index] # Index vector for downsizing of resampled site-blocks (DS1)
    boot.b.ydis2 <- b.ydis.resample2[boot.index]                                         # Downsize of site-blocks of community similarity (DS1)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2.2 <- nlsLM(boot.b.ydis2 ~ a * exp(b*boot.b.xdis2), 
                    start = c(a = 1, b = 0))
    
    int_v2[i] <- summary(mod2.2)$parameters[1, 1]
    slp_v2[i] <- summary(mod2.2)$parameters[2, 1]
    
  }
  
  ### Z_dep
  test_int <- (int_1 - int_2)/ sqrt((var(int_v1) + var(int_v2)) - 2 * cov(int_v1, int_v2))
  test_slp <- (slp_1 - slp_2)/ sqrt((var(slp_v1) + var(slp_v2)) - 2 * cov(slp_v1, slp_v2))
  
  ## pvalues calculation
  # intercept
  if (test_int<=0) {
    pval_int <- 2*pnorm(test_int, mean=0, sd=1, lower.tail=T)
  } else if (test_int>0){
    pval_int <- 2*pnorm(test_int, mean=0, sd=1, lower.tail=F)
  }
  
  #slope
  if (test_slp<=0) {
    pval_slp <- 2*pnorm(test_slp, mean=0, sd=1, lower.tail=T)
  } else if (test_slp>0){
    pval_slp <- 2*pnorm(test_slp, mean=0, sd=1, lower.tail=F)
  }
  
  res <- data.frame(var = c("Intercept", "Slope"),
                Est_1 = c(int_1, slp_1),
                Est_2 = c(int_2, slp_2),
                sd_1 = round(c(sqrt(var(int_v1)), sqrt(var(slp_v1))), 5),
                sd_2 = round(c(sqrt(var(int_v2)), sqrt(var(slp_v2))), 5),
                Zdep = c(test_int, test_slp), 
                p_value = c(ifelse(pval_int == 0, 1/nperm, pval_int), 
                            ifelse(pval_slp == 0, 1/nperm, pval_slp)))
  return(res)
}

## Compares the distance decay parameters between host trees 
Z_test_tree <- function(tree1, tree2, nperm){
  set.seed(666)
  options(scipen=999)
  ## select data 
  # Tree 1
  site_tree1 <- site_val |> 
    rownames_to_column("code") |> 
    filter(tree_species == tree1)
  
  ## tree community
  com_tree1 <-  site_tree1 |>
    select(code) |>
    left_join(OTU |> rownames_to_column("code")) |>
    column_to_rownames("code") %>%
    select(which(!colSums(.) %in% 0))
  
  # Community similarity 
  distcom1 <- as.dist(1 - vegdist(com_tree1, method = "horn", binary = FALSE))
  
  # Geographic distance by tree 
  distmat_tree1 <- as.dist(round(geodist(x = site_tree1[,c("longitude", "latitude")], 
                                         measure = "geodesic")/1000, 1)) ## transform in km and in distance matrix
  envidist1 <- distmat_tree1/1000 # unit = 1 000km
  
  # Tree 2
  site_tree2 <- site_val |> 
    rownames_to_column("code") |> 
    filter(tree_species == tree2)
  
  ## tree community
  com_tree2 <-  site_tree2 |>
    select(code) |>
    left_join(OTU |> rownames_to_column("code")) |>
    column_to_rownames("code") %>%
    select(which(!colSums(.) %in% 0))
  
  # Community similarity 
  distcom2 <- as.dist(1 - vegdist(com_tree2, method = "horn", binary = FALSE))
  
  # Geographic distance by tree 
  distmat_tree2 <- as.dist(round(geodist(x = site_tree2[,c("longitude", "latitude")], 
                                         measure = "geodesic")/1000, 1)) ## transform in km and in distance matrix
  envidist2 <- distmat_tree2/1000 # unit = 1 000km
  
  
  distcom1 = as.dist(distcom1)
  distcom2 = as.dist(distcom2)
  envidist1 = as.dist(envidist1)
  envidist2 = as.dist(envidist2)
  
  ## GLM
  dist1_vec <- as.vector(envidist1)
  dist2_vec <- as.vector(envidist2)
  distcom1_vec <- as.vector(distcom1)
  distcom2_vec <- as.vector(distcom2)
  
  # ## Using method Millar 2011 -> abandonned 
  # # Community 1
  # mod1.1 <- glm(distcom1_vec ~ dist1_vec, 
  #               family = quasibinomial(link = "log"))
  # int_1 <- exp(mod1.1$coefficients[1]) # transform intercept with exp but not slope
  # slp_1 <- mod1.1$coefficients[2]
  # 
  # # Community 2
  # mod1.2 <- glm(distcom2_vec ~ dist2_vec, 
  #               family = quasibinomial(link = "log"))
  # int_2 <- exp(mod1.2$coefficients[1]) # transform intercept with exp but not slope
  # slp_2 <- mod1.2$coefficients[2]
  
  # Using nlsLM (cf betapart package)
  # Community 1
  mod1.1 <- nlsLM(distcom1_vec ~ a * exp(b*dist1_vec), 
                  start = c(a = 1, b = 0))
  
  int_1 <- summary(mod1.1)$parameters[1, 1]
  slp_1 <- summary(mod1.1)$parameters[2, 1]
  
  # Community 2
  mod1.2 <- nlsLM(distcom2_vec ~ a * exp(b*dist2_vec), 
                  start = c(a = 1, b = 0))
  
  int_2 <- summary(mod1.2)$parameters[1, 1]
  slp_2 <- summary(mod1.2)$parameters[2, 1]
  
  ## Site - Block permutation
  # forcommunity 1
  nsites1 = nrow(as.matrix(envidist1))
  xdis.matrix1 <- as.matrix(envidist1)# Convert triangular matrix into square matrix
  ydis.matrix1 <- as.matrix(distcom1)  
  
  # Convert triangular matrix into square matrix
  
  blocks.xdis1  <- matrix(rep(0, nsites1*(nsites1 - 1)), ncol = nsites1)    # Empty matrix, to store site-block spatial distances 
  blocks.ydis1  <- matrix(rep(0, nsites1*(nsites1 - 1)), ncol = nsites1)    # Empty matrix, to store site-block community similarities
  
  for (i in 1:nsites1){
    
    blocks.xdis1[,i]  <- xdis.matrix1[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.ydis1[,i]  <- ydis.matrix1[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
  }
  
  # RESAMPLING
  
  int_v1 <- vector()
  slp_v1 <- vector()
  
  for (i in 1:nperm){
    
    # Block-resample with replacement: The matrices with site-blocks are resampled ensuring that the same cases are selected in both matrices
    b.i.resample <- sample(x = 1:nsites1, size=nsites1, replace=T)                           # Index vector for resampling with replacement (DS1)
    
    b.xdis.resample1 <- blocks.xdis1[, b.i.resample]                                     # Resampling (with replacement) of blocks of spatial distance values (DS1)
    b.ydis.resample1 <- blocks.ydis1[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS1)
    
    # Downsizing of resampled site-blocks: To preserve the size of the observed distance/similarity matrices, only the number of cases in the original distance/similarity matrix is randomly selected
    boot.index <- sample(1:length(b.xdis.resample1), 
                         size = nsites1*(nsites1 - 1)/2)
    
    boot.b.xdis1 <- b.xdis.resample1[boot.index] # Index vector for downsizing of resampled site-blocks (DS1)
    boot.b.ydis1 <- b.ydis.resample1[boot.index]                                         # Downsize of site-blocks of community similarity (DS1)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2.1 <- nlsLM(boot.b.ydis1 ~ a * exp(b*boot.b.xdis1), 
                    start = c(a = 1, b = 0))
    
    int_v1[i] <- summary(mod2.1)$parameters[1, 1]
    slp_v1[i] <- summary(mod2.1)$parameters[2, 1]
    
  }
  
  # forcommunity 2
  nsites2      <- nrow(as.matrix(envidist2))
  xdis.matrix2 <- as.matrix(envidist2)# Convert triangular matrix into square matrix
  ydis.matrix2 <- as.matrix(distcom2)  
  
  # Convert triangular matrix into square matrix
  
  blocks.xdis2  <- matrix(rep(0, nsites2*(nsites2 - 1)), ncol = nsites2)    # Empty matrix, to store site-block spatial distances 
  blocks.ydis2  <- matrix(rep(0, nsites2*(nsites2 - 1)), ncol = nsites2)    # Empty matrix, to store site-block community similarities
  
  for (i in 1:nsites2){
    
    blocks.xdis2[,i]  <- xdis.matrix2[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
    blocks.ydis2[,i]  <- ydis.matrix2[i,][-i]    # Store block corresponding to i site while excluding the diagonal value
  }
  
  # RESAMPLING
  
  int_v2 <- vector()
  slp_v2 <- vector()
  
  for (i in 1:nperm){
    
    # Block-resample with replacement: The matrices with site-blocks are resampled ensuring that the same cases are selected in both matrices
    b.i.resample<-sample(x = 1:nsites2, 
                         size=nsites2,
                         replace=T)  # Index vector for resampling with replacement (DS1)
    
    b.xdis.resample2 <- blocks.xdis2[, b.i.resample]                                     # Resampling (with replacement) of blocks of spatial distance values (DS1)
    b.ydis.resample2 <- blocks.ydis2[, b.i.resample]                                     # Resampling (with replacement) of blocks of community similarity values (DS1)
    
    # Downsizing of resampled site-blocks: To preserve the size of the observed distance/similarity matrices, only the number of cases in the original distance/similarity matrix is randomly selected
    boot.index <- sample(1:length(b.xdis.resample2), 
                         size = nsites2*(nsites2 - 1)/2)
    
    boot.b.xdis2 <- b.xdis.resample2[boot.index] # Index vector for downsizing of resampled site-blocks (DS1)
    boot.b.ydis2 <- b.ydis.resample2[boot.index]                                         # Downsize of site-blocks of community similarity (DS1)
    
    ### C.2 DISTANCE-DECAY MODEL AND PARAMETER ESTIMATION FOR EACH RESAMPLE 
    mod2.2 <- nlsLM(boot.b.ydis2 ~ a * exp(b*boot.b.xdis2), 
                    start = c(a = 1, b = 0))
    
    int_v2[i] <- summary(mod2.2)$parameters[1, 1]
    slp_v2[i] <- summary(mod2.2)$parameters[2, 1]
    
  }
  
  ### Z_dep
  test_int <- (int_1 - int_2)/ sqrt((var(int_v1) + var(int_v2)) - 2 * cov(int_v1, int_v2))
  test_slp <- (slp_1 - slp_2)/ sqrt((var(slp_v1) + var(slp_v2)) - 2 * cov(slp_v1, slp_v2))
  
  ## pvalues calculation
  # intercept
  if (test_int<=0) {
    pval_int <- 2*pnorm(test_int, mean=0, sd=1, lower.tail=T)
  } else if (test_int>0){
    pval_int <- 2*pnorm(test_int, mean=0, sd=1, lower.tail=F)
  }
  
  #slope
  if (test_slp<=0) {
    pval_slp <- 2*pnorm(test_slp, mean=0, sd=1, lower.tail=T)
  } else if (test_slp>0){
    pval_slp <- 2*pnorm(test_slp, mean=0, sd=1, lower.tail=F)
  }
  
  res <- data.frame(var = c("Intercept", "Slope"),
                Est_1 = c(int_1, slp_1),
                Est_2 = c(int_2, slp_2),
                sd_1 = round(c(sqrt(var(int_v1)), sqrt(var(slp_v1))), 5),
                sd_2 = round(c(sqrt(var(int_v2)), sqrt(var(slp_v2))), 5),
                Zdep = c(test_int, test_slp), 
                p_value = c(ifelse(pval_int == 0, 1/nperm, pval_int), 
                            ifelse(pval_slp == 0, 1/nperm, pval_slp)))
  return(res)
}

##Function to plot distance decay model 
plot_dd <- function(model){
  
  r2 <- round(model[[1]]$pseudo.r.squared, 2)
  p_val <- round(model[[1]]$p.value, 4)
  
  model[[4]] |>
    ggplot(aes(x = geodist, y = distcom)) + 
    geom_point(alpha = 0.05) +
    geom_smooth(method = "glm", 
                formula = y ~ exp(-x), 
                se = F, lwd = 3) +
    theme_light() +
    #xlim(0.07, 1.25) +
    ylim(0, 1) +
    labs(x = "Geographic distance (x 1000 km)", 
         y = "Community similarity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 3.25)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    geom_text(x=1.5, y=0.9, label=paste0("Intercept = ", round(model[[1]]$first.parameter, 2), 
                                         "   Slope = ", round(model[[1]]$second.parameter, 2))) +
    geom_text(x=1.5, y=0.8, label=paste0("pseudo R2 = ", r2, 
                                         "   p-value = ", p_val))
}

# helper function to plot distance decay betweentwo groups with spagheti plots 
plot_dist_dec <- function(data, nperm, ncurves){
  # number of curves on graph
  perm_n <- sample(1:nperm, size = ncurves, replace = F)
  
  g <-  data.frame(data) |>
    filter(perm %in% perm_n) |>
    ggplot(aes(
      y = boot.b.ydis, x = boot.b.xdis)) + 
    geom_line(stat = "smooth", method = "glm",
              formula = y ~ exp(-x), se = F, fullrange = T,
              linewidth = 0.5, alpha = 0.1,
              aes(linetype = as.factor(boot.b.group),
                  colour = as.factor(boot.b.group), fill = perm),
              show.legend = F) +
    geom_smooth(method = "glm",
                formula = y ~ exp(-x), se = F, fullrange = T,
                linewidth = 2, alpha = 1,
                aes(linetype = as.factor(boot.b.group),
                    colour = as.factor(boot.b.group)),
                show.legend = F) +
    coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.7)) +
    labs(x = "Geographic distance (km)", 
         y = "Community similarity",
         colour = "") +
    theme_classic() +
    theme(plot.margin = unit(c(3, 10, 3, 3), "mm")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA), 
                       labels = scales::label_number(scale = 1000)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_color_brewer(palette = "Dark2")
  return(g)
}

# calculate distance decay for subgroups of fungi
dd_group <- function(type_var, var, nperm){
  
  if(type_var == "phylum") {
    traits_com <- trait %>% 
      filter(phylum == var)} else {
        if(type_var == "order"){
          traits_com <- trait %>% 
            filter(order == var)} else {
              if(type_var == "fruitbody"){
                traits_com <- trait %>% 
                  filter(fruitbody == var)} else {
                    if(type_var == "type"){
                      traits_com <- trait %>% 
                        filter(type == var)}}}}
  
  otu <- OTU %>% 
    select(any_of(traits_com$SH))
  
  ## check communities where the subgroup is present
  pres <- otu %>% 
    transmute(sum = as.numeric(rowSums(across(where(is.numeric)))>=1))
  
  ## select the communities
  otu <- otu[pres$sum==1,]
  
  ## subset the matrices except for community, recalculate
  com <- 1 - vegdist(otu, method = "horn", binary = FALSE)
  dist <- dist_subset(distmat, which(pres$sum == 1))
  sem <- dist_subset(same_ecoreg_matrix, which(pres$sum == 1))
  tre <- dist_subset(same_host_matrix, which(pres$sum == 1))
  
  # general model
  mod_group <- dd_mod_perm(
    distcom = com,
    envidist = dist,
    host = tre,
    nperm = nperm,
    byhost = TRUE)
  mod_group <- mod_group[[1]]
  
  ## intercept comparison
  zres <- Z_test_paired(distcom = com, 
                        envidist = dist, 
                        byhost = TRUE,
                        group = sem, 
                        host = tre,
                        nperm = nperm)
  
  res <- list(Intercept = mod_group$first.parameter, 
              Slope = mod_group$second.parameter,
              pseudo_R2 = mod_group$pseudo.r.squared,
              p_value = mod_group$p.value, 
              ZDep_tes = zres[[1]])
  return(res)
}

# calculate distance decay for different host trees
dd_tree <- function(tree, nperm){
  
  site_tree <- site_val |> 
    rownames_to_column("code") |> 
    filter(tree_species == tree)
  
  ## tree community
  com_tree <-  site_tree |>
    select(code) |>
    left_join(OTU |> rownames_to_column("code")) |>
    column_to_rownames("code") %>%
    select(which(!colSums(.) %in% 0))
  
  # Community similarity 
  com_mat_tree <- as.dist(1 - vegdist(com_tree, method = "horn", binary = FALSE))
  
  # Geographic distance by tree 
  distmat_tree <- as.dist(round(geodist(x = site_tree[,c("longitude", "latitude")], 
                                        measure = "geodesic")/1000, 1)) ## transform in km and in distance matrix
  distmat_tree <- distmat_tree/1000 # unit = 1 000km
  #distmat_tree <- dist_subset(distmat, which(site_val$tree_species %in% tree))
  
  # Same ecoregion
  sem_tree <- matrix(NA, nrow = nrow(site_tree), ncol = nrow(site_tree))
  for(i in 1:nrow(site_tree)) {
    for(j in 1:nrow(site_tree)) {
      sem_tree[i,j]=as.numeric(site_tree$ecoregion[i]==site_tree$ecoregion[j])
    }
  }
  
  # matrix of "same host" for next function
  dime <- dim(sem_tree)
  mat_host = as.dist(matrix(1, nrow = dime[1], ncol = dime[2]))
  
  #General model
  mod_tree <- dd_mod_perm(
    distcom = com_mat_tree,
    envidist = distmat_tree,
    host = mat_host,
    nperm = nperm,
    byhost = TRUE)
  mod_tree <- mod_tree[[1]]
  
  ### intercept comparison
  Z_tree <-Z_test_paired(distcom = com_mat_tree, 
                         envidist = distmat_tree, 
                         group = sem_tree,
                         host = mat_host, 
                         byhost= TRUE,
                         nperm = nperm)
  
  res <- list(Intercept = mod_tree$first.parameter, 
              Slope = mod_tree$second.parameter,
              pseudo_R2 = mod_tree$pseudo.r.squared,
              p_value = mod_tree$p.value, 
              ZDep_tes = Z_tree[[1]])
  return(res)
}

## Test indicators species for each host tree separately 
## -> easy to observe in the field and avoid the non-generalisability 
## because some regions have just conifers and others have just broadleaf 
## and not balanced proportions of hosts
indval_group <- function(tree){
  
  ## select subset of sites
  site_tree <- site_val |> 
    rownames_to_column("code") |> 
    filter(tree_species == tree)
  
  eco_sup2 <- site_tree |> 
    count(ecoregion) |>
    filter(n >= 2)
  
  site_tree <- site_tree |>
    filter(ecoregion %in% eco_sup2$ecoregion)
  
  ## community
  com_tree <-  site_tree |>
    select(code) |>
    left_join(OTU |> rownames_to_column("code")) |>
    column_to_rownames("code") %>%
    select(which(!colSums(.) %in% 0))
  
  com_tree <- com_tree |> 
    select_if(colSums(com_tree != 0) >= 2 & colSums(com_tree) >= 10)
  
  rownames(com_tree) == site_tree$code
  
  ## calculate Indval
  indval_tree = multipatt(com_tree, 
                          site_tree$ecoregion, 
                          func = "IndVal.g", 
                          duleg = TRUE, # to only test single ER (no combination)
                          control = how(nperm=50000))
  summary(indval_tree)
  
  ## Correct p-values (FDR) and bind the results to original output
  fdr.p.value <- p.adjust(indval_tree$sign[,"p.value"], method="fdr")
  fdr.p.value
  
  indval_tree_p <- cbind(indval_tree$sign, fdr.p.value)
  
  ## Omit NA values and print only those with p <= 0.05
  attach(indval_tree_p) ## allows sorting by header
  indval.tree.fdr <- indval_tree_p[order(fdr.p.value, p.value, na.last=NA),] ## sort output by fdr, then p-value, omit "NA"
  detach(indval_tree_p)

  ## the A component is the specificity. 
  ## Species with a high specificity are found mostly 
  ## or only in that group 
  indval_a <- data.frame(indval_tree$A) 
  
  specificity <- data.frame(ecoregion=NA, specificity=NA)
  
  for(i in 1:ncol(indval_tree$A)){
    specificity[i, 1] <- colnames(indval_a)[i]
    specificity[i, 2] <- sum(indval_a[,i] > 0.8)
  }
  res <- list(indval.tree.fdr, specificity)
  
  saveRDS(res, file = paste0("output_indval/indval_", tree, ".Rdata"))
  return(res)
}