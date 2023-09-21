library(tidyverse)
library(vegan)
library(geodist)
library(usedist)
library(RColorBrewer)
library(pairwiseAdonis)
library(ggpubr)
library(indicspecies)
library(iNEXT)
source("Function_Zdep_paired.R")

#### Data preparation -------------------------------------------------------
### Load the data 
site_val <- read.csv("data/sites.csv") |>
  arrange(code) |>
  column_to_rownames("code") |>
  #mutate_if(is.character, as.factor) |>
  mutate(MAP_BC = as.numeric(MAP_BC))
  
## summary of the ecoregions
site_val |> count(ecoregion)

site_val %>% 
  group_by(ecoregion, tree_species)%>% 
  count() %>% 
  pivot_wider(names_from="tree_species", values_from = "n")

### Load the OTU table
OTU <- read.csv("data/OTU.csv") |> 
  arrange(code) |>
  column_to_rownames("code") |> 
  mutate_if(is.integer, as.numeric)

rownames(site_val) == rownames(OTU)

### Load the trait table for trait specific analyses
trait <- read_csv("data/trait_phylo.csv")

### Calculate the disimilarity matrices between communities and distance
## Distance matrix using morista-horn distance
dist_com     <- vegdist(OTU, method = "horn", binary = FALSE)
dist_com_bin <- vegdist(OTU, method = "horn", binary = TRUE)

## Community similarity matrix 
sim_com     <- 1 - dist_com
sim_com_bin <- 1- dist_com_bin

#### Models of distance decay - Data preparation ----
# Matrix of appartenance to the same ecoregion 
same_ecoreg_matrix <- matrix(NA, nrow = nrow(site_val), ncol = nrow(site_val))  
for(i in 1:nrow(site_val)) {
  for(j in 1:nrow(site_val)) {
    same_ecoreg_matrix[i,j]=as.numeric(site_val$ecoregion[i]==site_val$ecoregion[j])
  }
}
same_ecoreg_matrix <- as.dist(same_ecoreg_matrix)

## Matrix of appartenance to the same host
same_host_matrix <- matrix(NA, nrow = nrow(site_val), ncol = nrow(site_val))  
for(i in 1:nrow(site_val)) {
  for(j in 1:nrow(site_val)) {
    same_host_matrix[i,j]=as.numeric(site_val$tree_species[i]==site_val$tree_species[j])
  }
}
same_host_matrix <- as.dist(same_host_matrix)

## Matrix of geographic distances 
distmat <- as.dist(round(geodist(x = site_val[,c("longitude", "latitude")], 
                                 measure = "geodesic")/1000, 1)) ## transform in km and in distance matrix
distmat <- distmat/1000 # unit = 1 000km

# # Create data frame
data <- data.frame(com_mat = as.vector(sim_com),
                   com_mat_bin = as.vector(sim_com_bin),
                   distmat = as.vector(distmat),
                   same_ecoreg_matrix = as.factor(as.vector(same_ecoreg_matrix)),
                   same_host_matrix = as.factor(same_host_matrix))

# data filtered by distance max between intra sites
#Mean distance int
same_ER_dist <- data |> filter(same_ecoreg_matrix == 1) |>
  summarize(mean = mean(distmat), max = max(distmat), min = min(distmat))
diff_ER_dist <- data |> filter(same_ecoreg_matrix == 0) |>
  summarize(mean = mean(distmat), max = max(distmat), min = min(distmat))

data <- data |>
  filter(distmat >= diff_ER_dist$min & distmat  <= same_ER_dist$max)

## dataset keeping only the pairs of similar hosts remove host difference effect
data_same_host <- data |>
  filter(same_host_matrix == 1) |>#2116 observations
  # while keeping only the ones within the same geographic distance boundaries
  filter(distmat >= diff_ER_dist$min & distmat  <= same_ER_dist$max) # 1573 observations

data_same_host |>
  group_by(same_ecoreg_matrix) |>
  count() # 1224 in different ecoregions and 349 in the same

data # for raw data
data_same_host # same host species comparison and distance range

## General model -----
## Helper function that calculate the model, plot the graph and 
## adds on the graph models parameters and deviance relative to a model
## with simple intercept only
plot_dd <- function(data, y, x){

    summary(mod <- glm(y ~ x, 
                       data = data, 
                       family = quasibinomial(link = "log"))) 
  
  coef <- c(exp(mod$coefficients[1]), mod$coefficients[2])
  
  # null model
  null = glm(y ~ 1, 
             data = data, 
             family = quasibinomial(link = "log"))
  test <- anova(mod, null, test="Chisq")
  
  data |>
  ggplot(aes(x, y)) + 
    geom_point(alpha = 0.05) +
    geom_smooth(method = "glm", 
                formula = y ~ exp(-x), 
                se = F, lwd = 3) +
    theme_light() +
    xlim(0.07, 1.25) +
    ylim(0, 1) +
    labs(x = "Geographic distance (x 1000 km)", 
         y = "Community similarity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.065, 1.25)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    geom_text(x=0.6, y=0.9, label=paste0("Intercept = ", round(coef[1], 2), 
                                         "   Slope = ", round(coef[2], 2))) +
    geom_text(x=0.6, y=0.8, label=paste0("Deviance = ", round(test$Deviance[2], 2), 
                                         "   p-value = ", round(test$`Pr(>Chi)`[2], 2)))
}

## abundance all plots
dd_all <- plot_dd(data = data, 
                  y = data$com_mat, 
                  x = data$distmat)

## abundance - same tree species
dd_same_tree <- plot_dd(data = data_same_host, 
                        y = data_same_host$com_mat, 
                        x = data_same_host$distmat)

## occurence all plots
dd_all_bin <- plot_dd(data = data, 
                  y = data$com_mat_bin, 
                  x = data$distmat)

## occurence - same tree species
dd_same_tree_bin <- plot_dd(data = data_same_host, 
                        y = data_same_host$com_mat_bin, 
                        x = data_same_host$distmat)
options(repr.plot.width = 4, repr.plot.height =3)
ggarrange(dd_all, dd_same_tree, dd_all_bin, dd_same_tree_bin, 
          nrow = 2, ncol = 2, labels = "auto")

# ggsave("figures/distance_decay_general.png")

## Test the difference in distance decay parameters (with graphs) ####
## Based on the method of Martin-Devasa et al 2022
# helper function to plot distance decay with spagheti plots 
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


nperm <- 1000 # number of permutations
ncurves = 300 # number of curves in spaghetti plot
bytree <- TRUE # only keep pairs of plots with the same host

# with abundance
Z_t <- Z_test_paired(distcom = sim_com, 
                     envidist = distmat, 
                     group = same_ecoreg_matrix, # Change eco or bioreg
                     host = same_host_matrix,
                     byhost=bytree, 
                     nperm = nperm)

# with occurrence 
Z_t_bin <- Z_test_paired(distcom = sim_com_bin, 
                         envidist = distmat, 
                         group = same_ecoreg_matrix, # Change eco or bioreg
                         host = same_host_matrix,
                         byhost=bytree,
                         nperm = nperm)
# Results
Z_t[[1]]
Z_t_bin[[1]]

## plot figure with incertitude lines
DD <- plot_dist_dec(data = Z_t[3], 
                    nperm = nperm, 
                    ncurves = ncurves)

DD.bin <- plot_dist_dec(data = Z_t_bin[3], 
                        nperm = nperm, 
                        ncurves = ncurves)

ggarrange(DD + rremove("xlab") + rremove( "x.text"), DD.bin, 
          nrow = 2, labels = "auto")
# ggsave("figures/ecoregion_DD_horn_per_tree.jpg")

## Distance decay parameters for subgroups of fungi ----
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
traits_com

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

data <- data.frame(com = as.vector(com), 
                   dist = as.vector(dist), 
                   sem = as.vector(as.factor(sem)),
                   tre = as.vector(tre)) |>
  filter(tre == 1)

# data filtered by distance max between intra sites
#Mean distance int
same_ER_dist <- data |> filter(sem == 1) |>
  summarize(mean = mean(dist), max = max(dist), min = min(dist))
diff_ER_dist <- data |> filter(sem == 0) |>
  summarize(mean = mean(dist), max = max(dist), min = min(dist))

data <- data |>
  filter(dist >= diff_ER_dist$min & dist  <= same_ER_dist$max)

## Model global 
summary(mod_group <- glm(com ~ dist, 
                   data = data, 
                   family = quasibinomial(link = "log"))) 

coef <- c(exp(mod_group$coefficients[1]), mod_group$coefficients[2])

# null model
null = glm(com ~ 1, 
           data = data, 
           family = quasibinomial(link = "log"))
test <- anova(mod_group, null, test="Chisq")

## intercept comparison
zres <- Z_test_paired(distcom = com, 
                      envidist = dist, 
                      byhost = TRUE,
                      group = sem, 
                      host = tre,
                      nperm = nperm)

res_mod <- paste0("Intercept = ", round(coef[1], 2), 
       "   Slope = ", round(coef[2], 2),
       "   Deviance = ", round(test$Deviance[2], 2), 
       "   p-value = ", round(test$`Pr(>Chi)`[2], 2))
res <- list(zres[[1]], res_mod)

return(res)
}

# Phylum
dd_group(type_var = "phylum", 
         var = "Basidiomycota", 
         nperm = nperm)

dd_group(type_var = "phylum", 
         var = "Ascomycota", 
         nperm = nperm)

# order
dd_group(type_var = "order", 
         var = "Agaricales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Russulales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Atheliales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Boletales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Cantharellales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Pezizales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Thelephorales", 
         nperm = nperm)

dd_group(type_var = "order", 
         var = "Sebacinales", 
         nperm = nperm)

# fruitbody
dd_group(type_var = "fruitbody", 
         var = "epigeous", 
         nperm = nperm)

dd_group(type_var = "fruitbody", 
         var = "hypogeous", 
         nperm = nperm)

# type
dd_group(type_var = "type", 
         var = "mushroom", 
         nperm = nperm)

dd_group(type_var = "type", 
         var = "sclerotium", 
         nperm = nperm)

dd_group(type_var = "type", 
         var = "truffle", 
         nperm = nperm)

dd_group(type_var = "type", 
         var = "crust", 
         nperm = nperm)

## Distance decay per host tree ------------------------------------------
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
  sem_tree <- as.dist(sem_tree)
  ###  
  data <- data.frame(com = as.vector(com_mat_tree), 
                     dist = as.vector(distmat_tree), 
                     sem = as.vector(sem_tree))
  
  same_ER_dist <- data |> filter(sem == 1) |>
    summarize(mean = mean(dist), max = max(dist), min = min(dist))
  diff_ER_dist <- data |> filter(sem == 0) |>
    summarize(mean = mean(dist), max = max(dist), min = min(dist))
  
  data <- data |> filter(dist >= diff_ER_dist$min & dist  <= same_ER_dist$max)
  ###
  summary(mod_tree <- glm(com ~ dist, 
                          data = data, 
                          family = quasibinomial(link = "log")))
  
  
  coef <- c(exp(mod_tree$coefficients[1]), mod_tree$coefficients[2])
  
  # null model
  null = glm(com ~ 1, 
             data = data, 
             family = quasibinomial(link = "log"))
  
  
  test <- anova(mod_tree, null, test="Chisq")
  
  ###
  Z_tree <-Z_test_paired(distcom = com_mat_tree, 
                         envidist = distmat_tree, 
                         group = sem_tree,
                         host = sem_tree, # does not matter here because not used
                         byhost= FALSE,
                         nperm = nperm)
  
  res_mod <- paste0("Intercept = ", round(coef[1], 2), 
                    "   Slope = ", round(coef[2], 2),
                    "   Deviance = ", round(test$Deviance[2], 2), 
                    "   p-value = ", round(test$`Pr(>Chi)`[2], 2))
  res <- list(Z_tree[[1]], res_mod)
  return(res)
}

####
dd_tree(tree = "BEECH", nperm = nperm)
dd_tree(tree = "PINE", nperm = nperm)
dd_tree(tree = "SPRUCE", nperm = nperm)
dd_tree(tree = "OAK", nperm = nperm)



#### Permanova testing differences in composition overall ----
(a1 = adonis2(dist_com ~ ecoregion,
        permutations = 9999,
        strata = site_val$tree_species,
        data = site_val))

## pairwise post-hoc
(pairwise <- pairwise.adonis2(dist_com ~ as.factor(ecoregion), 
              data = site_val,
              strata = "tree_species",
              nperm = 5000))

site_val|> 
  select(ecoregion) |> 
  count(ecoregion) |> 
  mutate(er_code = as.numeric(as.factor(ecoregion)))

df_pair <- data.frame(pairwise) |>
  t() |> 
  data.frame() |> 
  select(ecoregion = as.factor.ecoregion.) |> 
  rownames_to_column("variable") |>
  filter(variable != "parent_call") |> 
  separate(variable, into = c("site1", "vs", "site2", "param")) |>
  unite(sites, c("site1", "vs", "site2"), sep = "_") |>
  mutate(ecoregion = as.numeric(ecoregion)) |>
  pivot_wider(names_from = param, values_from = ecoregion) |>
  select(-Df)#|>
  # write.csv("pairwise_permanova_output.csv")

#### Variance partition -----
VP <- vegan::varpart(dist_com, 
                     X = ~ecoregion, ~MAT_BC + MAP_BC + ATS_BC + APS_BC, ~tree_species,
                     data = site_val)

par(mar=c(1, 1, 1, 1))
plot(VP,Xnames = c("Ecoregion", "Climate", "Host"))

# Test of host total
anova(dbrda(dist_com ~ tree_species, 
            data = site_val), 
      permutations = how(nperm = 999))

# Test of ecoregion total
anova(dbrda(dist_com ~ecoregion, 
            data = site_val), 
      permutations = how(nperm = 999))

# Test of climate total
anova(dbrda(dist_com ~MAT_BC + MAP_BC + ATS_BC + APS_BC, data = site_val), 
      permutations = how(nperm = 999))

#### NMDS representation ----
nmds <- metaMDS(dist_com, k = 6)
plot(nmds)

d1 <- data.frame(nmds$points[,]) %>%
  rownames_to_column("code") %>%
  left_join(site_val |>  rownames_to_column("code")) %>%
  select(code, MDS1, MDS2, MDS3, 
         tree_species, ecoregion) 

d <- d1 %>%
  group_by(ecoregion) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()

d2 <- d %>%
  count(ecoregion) %>%
  select(-n) %>%
  mutate(ecoregion_abbr = c(
    "Alps", "Baltic", "Cantabrian", "Carpathian", "Celtic",
    "Central", "English", "Atlantic", "Pannonian", "Sarmatic",
     "Taiga", "Western"))

d <- d %>% left_join(d2)
str(d)

# Figure
ggplot(data = d, aes(y = MDS2, x = MDS1))  +
  facet_wrap(~ ecoregion_abbr) +
  #geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  scale_color_brewer(palette="Paired") +
  geom_polygon(colour = "red",alpha = 0, lwd = 1) +
  geom_point(data = d1, 
             aes(y = MDS2, x = MDS1), alpha = 0.3)
#ggsave("figure/nmds.png")

#### Cluster of ecoregions ----
site_val <- site_val |> rownames_to_column("code")
otu_clust <- OTU %>%
  rownames_to_column("code") %>%
  left_join(site_val[,c("code","ecoregion")]) %>%
  dplyr::select(-code) %>%
  group_by(ecoregion) %>%
  summarise_all(sum) %>%
  column_to_rownames("ecoregion")

# Number of species in each ecoregion 
otu_ER_occ <- otu_clust
otu_ER_occ[otu_ER_occ > 0] <- 1
otu_ER_occ |> rowSums()

# Ecoregion completeness using inext ----
sac_er <- iNEXT(t(otu_clust), datatype = "abundance")
sac_er$iNextEst$coverage_based |> filter(Method == "Observed")

## distance and tree building
dist_clust <- vegdist(otu_clust, method = "horn", binary = FALSE)
dend <- as.dendrogram(hclust(dist_clust))

par(mar=c(5,3,1,15))
# plot(dend, horiz = T,
#      xlab = "Moristia-Horn Distance")%>%
# set("branches_k_color", k = 4) %>%
#   rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2)

library(dendextend)
dend %>% set("branches_k_color", k = 4) %>% plot(horiz = TRUE, 
                                                 xlab = "Morisita-Horn Distance")
dend %>% 
  rect.dendrogram(k = 4, horiz = TRUE, border = 8, lty = 5, lwd = 2,
                  xlab = "Moristia-Horn Distance")

#### Indicator species of each ecoregion ----
## Calculation of IndVal.g 
 # -> the g is a correction factor for the fact that some species 
 # are present in more groups than others : Tichý and Chytrý 2006

# Keep species that are present in more than 2 sites -> 695 species
OTU_sup2 <- OTU |> 
  select_if(colSums(OTU != 0) >= 2) 

corresp_code <- tibble(code_ecoregion = c(1:12)) |>
  mutate(ecoregion_abbr = c(
  "Alps", "Baltic", "Cantabrian", "Carpathian", "Celtic",
  "Central", "English", "Atlantic", "Pannonian", "Sarmatic",
  "Taiga", "Western"))

site_val <- site_val |> 
  rownames_to_column("code") |>
  mutate(code_ecoregion = as.numeric(as.factor(ecoregion)))|>
  left_join(corresp_code)
  
#OTU_gen <- OTU_sup2[,colnames(OTU_sup2) %in% generalists$nms]

indval = multipatt(OTU_sup2, 
                   site_val$ecoregion_abbr, 
                   control = how(nperm=999))
#indval_gen = multipatt(OTU_gen, site_val$ecoregion_code, control = how(nperm=999))

summary(indval, indvalcomp = T)

saveRDS(indval, file = "output_indval/indVal.Rdata")

## Test indicators species for each host tree separately 
## -> easy to observe in the field and avoid the non-generalisability 
## because some regions have just conifers and others have just broadleaf 
## and not balanced proportions of hosts
indval_group <- function(tree){
  
  ## select subset of sites
    site_tree <- site_val |> 
   # rownames_to_column("code") |> 
    filter(tree_species == tree)
  
  ## community
  com_tree <-  site_tree |>
    select(code) |>
    left_join(OTU_sup2 |> rownames_to_column("code")) |>
    column_to_rownames("code") %>%
    select(which(!colSums(.) %in% 0))

rownames(com_tree) == site_tree$code

  ## calculate Indval
indval_tree = multipatt(com_tree, 
                         site_tree$ecoregion_abbr, 
                         control = how(nperm=999))
summary(indval_tree, indvalcomp = T)
saveRDS(indval_tree, file = paste0("output_indval/indval_", tree, ".Rdata"))
return(indval_tree)
}

iv_pine <- indval_group("PINE")
iv_spruce <- indval_group("SPRUCE")
iv_oak <- indval_group("OAK")
iv_beech <- indval_group("BEECH")

## the A component is the specificity. 
## Species with a high specificity are found mostly 
## or only in that group 
indval_a <- data.frame(indval$A[,1:12])
specificity <- data.frame(ecoregion=NA, specificity=NA)

for(i in 1:12){
specificity[i, 1] <- colnames(indval_a)[i]
specificity[i, 2] <- sum(indval_a[,i] == 1)
}

######## END #############