library(tidyverse)
library(vegan)
library(geodist)
library(usedist)
library(minpack.lm)
library(RColorBrewer)
library(pairwiseAdonis)
library(ggpubr)
library(indicspecies)
library(iNEXT)
source("Functions.R")

#### Data preparation -------------------------------------------------------
### Load the data 
data <- readRDS("data/data.rds")

site_val <- data[[2]] |>
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
OTU <- data[[1]] |> 
  arrange(code) |>
  column_to_rownames("code") |> 
  mutate_if(is.integer, as.numeric)

rownames(site_val) == rownames(OTU)
n_tips <- sum(OTU) # number of individual root tips
n_otu <- ncol(OTU) # number of OTUs

### Load the trait table for trait specific analyses
trait <- data[[3]]

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

## general parameters
nperm = 1000

## General model -----
## abundance all plots
dd_all <- dd_mod_perm(
  distcom = sim_com,
  envidist = distmat,
  host = same_host_matrix,
  nperm = nperm,
  byhost = FALSE)

## abundance - same tree species
dd_same_tree <- dd_mod_perm(
  distcom = sim_com,
  envidist = distmat,
  host = same_host_matrix,
  nperm = nperm,
  byhost = TRUE)

## occurence all plots
dd_all_bin <- dd_mod_perm(
  distcom = sim_com_bin,
  envidist = distmat,
  host = same_host_matrix,
  nperm = nperm,
  byhost = FALSE)

## occurence - same tree species
dd_same_tree_bin <- dd_mod_perm(
  distcom = sim_com_bin,
  envidist = distmat,
  host = same_host_matrix,
  nperm = nperm,
  byhost = TRUE)

## Test the difference in distance decay parameters (with graphs) ####
## Based on the method of Martin-Devasa et al 2022
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

## Distance decay parameters for subgroups of fungi ----
# Phylum
dd_basidio <- dd_group(type_var = "phylum", 
         var = "Basidiomycota", 
         nperm = nperm)

dd_asco <- dd_group(type_var = "phylum", 
         var = "Ascomycota", 
         nperm = nperm)

# order
dd_agaricales <- dd_group(type_var = "order", 
         var = "Agaricales", 
         nperm = nperm)

dd_russulales <- dd_group(type_var = "order", 
         var = "Russulales", 
         nperm = nperm)

dd_atheliales <- dd_group(type_var = "order", 
         var = "Atheliales", 
         nperm = nperm)

dd_boletales <- dd_group(type_var = "order", 
         var = "Boletales", 
         nperm = nperm)

dd_cantharellales <- dd_group(type_var = "order", 
         var = "Cantharellales", 
         nperm = nperm)

dd_pezizales <- dd_group(type_var = "order", 
         var = "Pezizales", 
         nperm = nperm)

dd_telephorales <- dd_group(type_var = "order", 
         var = "Thelephorales", 
         nperm = nperm)

# fruitbody
dd_epigeous <- dd_group(type_var = "fruitbody", 
         var = "epigeous", 
         nperm = nperm)

dd_hypogeous <- dd_group(type_var = "fruitbody", 
         var = "hypogeous", 
         nperm = nperm)

# type
dd_mushroom <- dd_group(type_var = "type", 
         var = "mushroom", 
         nperm = nperm)

dd_sclerotium <- dd_group(type_var = "type", 
         var = "sclerotium", 
         nperm = nperm)

dd_truffle <- dd_group(type_var = "type", 
         var = "truffle", 
         nperm = nperm)

dd_crust <- dd_group(type_var = "type", 
         var = "crust", 
         nperm = nperm)


# proportion of total species and root tips 
# for the different groups
russulales <- trait |> filter(order == "Russulales") #285 OTU
atheliales <- trait |> filter(order == "Atheliales") #114 OTU
otu_russ_ath <- OTU |> 
  select(any_of(
    c(russulales$SH, atheliales$SH)
    ))
sum(otu_russ_ath)/n_tips

# proportions of basidios vs ascos
basidio <- trait |> filter(phylum == "Basidiomycota")
basidio_otu <- OTU |> select(all_of(basidio$SH))
n_tips_basidio <- sum(basidio_otu) # 19509 root tips 1224 OTU
asco <- trait |> filter(phylum == "Ascomycota")
asco_otu <- OTU |> select(all_of(asco$SH))
n_tips_asco <- sum(asco_otu) # 4423 root tips 126 OTU
n_tips_basidio/(n_tips_basidio + n_tips_asco)

## Distance decay per host tree ----
dd_beech   <- dd_tree(tree = "BEECH", nperm = nperm)
dd_pine    <- dd_tree(tree = "PINE", nperm = nperm)
dd_spruce  <- dd_tree(tree = "SPRUCE", nperm = nperm)
dd_oak     <- dd_tree(tree = "OAK", nperm = nperm)
adj_p = p.adjust(c(dd_beech[[4]], dd_pine[[4]], dd_spruce[[4]], dd_oak[[4]]))

# test for differences in DD parameters between host trees
oak_beech    <-  Z_test_tree(tree1 = "OAK", tree2 = "BEECH", nperm = nperm)
oak_pine     <-  Z_test_tree(tree1 = "OAK", tree2 = "PINE", nperm = nperm)
oak_spruce   <-  Z_test_tree(tree1 = "OAK", tree2 = "SPRUCE", nperm = nperm)
beech_pine   <-  Z_test_tree(tree1 = "BEECH", tree2 = "PINE", nperm = nperm)
beech_spruce <-  Z_test_tree(tree1 = "BEECH", tree2 = "SPRUCE", nperm = nperm)
pine_spruce  <-  Z_test_tree(tree1 = "PINE", tree2 = "SPRUCE", nperm = nperm)



##### Table Results-----
results_table <- data.frame(
  Group = c("Abundance", "Occurence", "Epigeous", "Hypogeous", 
            "Mushroom", "Sclerotium", "Crust", "Truffle", 
            "Basidiomycota", "Ascomycota", 
            "Russulales", "Atheliales", "Telephorales", "Agaricales", 
            "Boletales", "Cantharellales", "Pezizales",
            "P. sylverstris", "P. abies", "F. sylvatica", "Quercus spp."), 
  Intercept = round(c(dd_same_tree[[1]]$first.parameter, 
                      dd_same_tree_bin[[1]]$first.parameter,
                      dd_epigeous[[1]], dd_hypogeous[[1]], dd_mushroom[[1]],
                      dd_sclerotium[[1]], dd_crust[[1]], dd_truffle[[1]],
                      dd_basidio[[1]], dd_asco[[1]], dd_russulales[[1]],
                      dd_atheliales[[1]], dd_telephorales[[1]], dd_agaricales[[1]],
                      dd_boletales[[1]], dd_cantharellales[[1]], dd_pezizales[[1]], 
                      dd_pine[[1]], dd_spruce[[1]], dd_beech[[1]], dd_oak[[1]]), 2), 
  Slope = round(c(dd_same_tree[[1]]$second.parameter, 
                  dd_same_tree_bin[[1]]$second.parameter,
                  dd_epigeous[[2]], dd_hypogeous[[2]], dd_mushroom[[2]],
                  dd_sclerotium[[2]], dd_crust[[2]], dd_truffle[[2]],
                  dd_basidio[[2]], dd_asco[[2]], dd_russulales[[2]],
                  dd_atheliales[[2]], dd_telephorales[[2]], dd_agaricales[[2]],
                  dd_boletales[[2]], dd_cantharellales[[2]], dd_pezizales[[2]],
                  dd_pine[[2]], dd_spruce[[2]], dd_beech[[2]], dd_oak[[2]]), 2),
  Intra = c(Z_t[[1]][1, 2], Z_t_bin[[1]][1, 2], 
            dd_epigeous[[5]][1, 2], dd_hypogeous[[5]][1, 2], dd_mushroom[[5]][1, 2],
            dd_sclerotium[[5]][1, 2], dd_crust[[5]][1, 2], dd_truffle[[5]][1, 2],
            dd_basidio[[5]][1, 2], dd_asco[[5]][1, 2], dd_russulales[[5]][1, 2],
            dd_atheliales[[5]][1, 2], dd_telephorales[[5]][1, 2], dd_agaricales[[5]][1, 2],
            dd_boletales[[5]][1, 2], dd_cantharellales[[5]][1, 2], dd_pezizales[[5]][1, 2],
            dd_pine[[5]][1, 2], dd_spruce[[5]][1, 2], dd_beech[[5]][1, 2],dd_oak[[5]][1, 2]),
  Inter = c(Z_t[[1]][1, 3], Z_t_bin[[1]][1, 3], 
            dd_epigeous[[5]][1, 3], dd_hypogeous[[5]][1, 3], dd_mushroom[[5]][1, 3],
            dd_sclerotium[[5]][1, 3], dd_crust[[5]][1, 3], dd_truffle[[5]][1, 3],
            dd_basidio[[5]][1, 3], dd_asco[[5]][1, 3], dd_russulales[[5]][1, 3],
            dd_atheliales[[5]][1, 3], dd_telephorales[[5]][1, 3], dd_agaricales[[5]][1, 3],
            dd_boletales[[5]][1, 3], dd_cantharellales[[5]][1, 3], dd_pezizales[[5]][1, 3],
            dd_pine[[5]][1, 3], dd_spruce[[5]][1, 3], dd_beech[[5]][1, 3],dd_oak[[5]][1, 3]),
  Zdep = c(Z_t[[1]][1, 6], Z_t_bin[[1]][1, 6], 
           dd_epigeous[[5]][1, 6], dd_hypogeous[[5]][1, 6], dd_mushroom[[5]][1, 6],
           dd_sclerotium[[5]][1, 6], dd_crust[[5]][1, 6], dd_truffle[[5]][1, 6],
           dd_basidio[[5]][1, 6], dd_asco[[5]][1, 6], dd_russulales[[5]][1, 6],
           dd_atheliales[[5]][1, 6], dd_telephorales[[5]][1, 6], dd_agaricales[[5]][1, 6],
           dd_boletales[[5]][1, 6], dd_cantharellales[[5]][1, 6], dd_pezizales[[5]][1, 6],
           dd_pine[[5]][1, 6], dd_spruce[[5]][1, 6], dd_beech[[5]][1, 6],dd_oak[[5]][1, 6]),
  p.value = c(Z_t[[1]][1, 7], Z_t_bin[[1]][1, 7], 
              dd_epigeous[[5]][1, 7], dd_hypogeous[[5]][1, 7], dd_mushroom[[5]][1, 7],
              dd_sclerotium[[5]][1, 7], dd_crust[[5]][1, 7], dd_truffle[[5]][1, 7],
              dd_basidio[[5]][1, 7], dd_asco[[5]][1, 7], dd_russulales[[5]][1, 7],
              dd_atheliales[[5]][1, 7], dd_telephorales[[5]][1, 7], dd_agaricales[[5]][1, 7],
              dd_boletales[[5]][1, 7], dd_cantharellales[[5]][1, 7], dd_pezizales[[5]][1, 7],
              dd_pine[[5]][1, 7], dd_spruce[[5]][1, 7], dd_beech[[5]][1, 7],dd_oak[[5]][1, 7])
)
## adjust the p value for multiple comparisons
results_table$p_adj <- p.adjust(results_table$p.value, method = "BH")
write.csv(results_table, "tables/dd_intra_inter.csv") 

## Table zdep test trees 
table_tree <- data.frame(
  trees = c("oak_beech", "oak_pine", "oak_spruce", 
            "beech_pine", "beech_spruce", "pine_spruce"), 
  intercept1 = c(oak_beech[1,2], oak_pine[1,2], oak_spruce[1,2], 
                 beech_pine[1,2], beech_spruce[1,2], pine_spruce[1,2]),
  intercept2 = c(oak_beech[1,3], oak_pine[1,3], oak_spruce[1,3], 
                 beech_pine[1,3], beech_spruce[1,3], pine_spruce[1,3]),
  zdep_intercept = c(oak_beech[1,6], oak_pine[1,6], oak_spruce[1,6], 
                     beech_pine[1,6], beech_spruce[1,6], pine_spruce[1,6]),
  p_adj_intercept = p.adjust(c(oak_beech[1,7], oak_pine[1,7], oak_spruce[1,7], 
                     beech_pine[1,7], beech_spruce[1,7], pine_spruce[1,7]), method = "BH"),
  slope1 = c(oak_beech[2,2], oak_pine[2,2], oak_spruce[2,2], 
             beech_pine[2,2], beech_spruce[2,2], pine_spruce[2,2]),
  slope2 = c(oak_beech[2,3], oak_pine[2,3], oak_spruce[2,3], 
             beech_pine[2,3], beech_spruce[2,3], pine_spruce[2,3]),
  zdep_slope = c(oak_beech[2,6], oak_pine[2,6], oak_spruce[2,6], 
                     beech_pine[2,6], beech_spruce[2,6], pine_spruce[2,6]),
  p_adj_slope = p.adjust(c(oak_beech[2,7], oak_pine[2,7], oak_spruce[2,7], 
                  beech_pine[2,7], beech_spruce[2,7], pine_spruce[2,7]), method = "BH")) 
write.csv(table_tree, "tables/dd_tree_comparison.csv")  

#### Figures ----
## General DD models 
dd_all_pl <- plot_dd(dd_all)
dd_same_tree_pl <- plot_dd(dd_same_tree)
dd_all_bin_pl <- plot_dd(dd_all_bin)
dd_same_tree_bin_pl <- plot_dd(dd_same_tree_bin)

# options(repr.plot.width = 5, repr.plot.height =3)
# ggarrange(dd_all_pl, dd_same_tree_pl, 
#           dd_all_bin_pl, dd_same_tree_bin_pl, 
#           nrow = 2, ncol = 2, labels = "auto")
# ggsave("figures/distance_decay_general.png")


## plot figure with incertitude lines
ncurves = 300 # number of curves in spaghetti plot
DD <- plot_dist_dec(data = Z_t[3],
                    nperm = nperm,
                    ncurves = ncurves)
DD.bin <- plot_dist_dec(data = Z_t_bin[3],
                       nperm = nperm,
                       ncurves = ncurves)
ggarrange(DD + rremove("xlab") + rremove( "x.text"), DD.bin,
         nrow = 2, labels = "auto")
# ggsave("figures/ecoregion_DD_horn_per_tree.jpg")

##############################################################
#### Permanova testing differences in composition overall ----
(a1 = adonis2(dist_com ~ ecoregion,
        permutations = 9999,
        strata = site_val$tree_species,
        data = site_val))

## pairwise post-hoc
site_val_code <- site_val |> mutate(ecoregion = as.numeric(as.factor(ecoregion)))
(pairwise <- pairwise.adonis2(dist_com ~ as.factor(ecoregion), 
              data = site_val_code,
              strata = "tree_species",
              nperm = 5000))

site_val|>
  select(ecoregion) |>
  count(ecoregion) |>
  mutate(er_code = as.numeric(as.factor(ecoregion)))

corresp_code <- tibble(code_ecoregion = c(1:12)) |>
  mutate(ecoregion_abbr = c(
    "Alps", "Baltic", "Cantabrian", "Carpathian", "Celtic",
    "Central", "English", "Atlantic", "Pannonian", "Sarmatic",
    "Taiga", "Western"))

df_pair <- data.frame(pairwise) |>
  t() |> 
  data.frame() |> 
  select(ecoregion = as.factor('ecoregion')) |> 
  rownames_to_column("variable") |>
  filter(variable != "parent_call") |> 
  separate(variable, into = c("site1", "vs", "site2", "param")) |> 
 # unite(sites, c("site1", "vs", "site2"), sep = "_") |>
  mutate(ecoregion = as.numeric(ecoregion)) |>
  pivot_wider(names_from = param, values_from = ecoregion) |>
  select(-Df) |>
  mutate(adj_p = p.adjust(Pr, method = "BH"))|>
  #select(site1, site2, R2, adj_p) |>
 mutate(site1 = as.integer(str_replace_all(site1, "X", "")), 
        site2 = as.integer(site2)) |>
  left_join(corresp_code, by = join_by(site1 == code_ecoregion)) |>
  left_join(corresp_code, by = join_by(site2 == code_ecoregion)) |>
  select(site_1 = ecoregion_abbr.x, site_2 = ecoregion_abbr.y, 
         SumOfSqs, R2, 'F', adj_p)
# write.csv(df_pair, "tables/pairwise_permanova_output.csv")

# Test for difference in the dispersion of each group 
# to test patterns observed on the NMDS 
multi_disper <- betadisper(
  d = dist_com, 
  group = site_val$ecoregion, 
  type = "centroid")

anova(multi_disper)
tukey <- TukeyHSD(multi_disper)
colnames(tukey$group) <- c("difference", "lower", "upper", "adj_p")
tukey <- data.frame(tukey$group) |> mutate(adj_p = round(adj_p, 9))
write.csv(tukey, "permdisp.csv")

#### Variance partition -----
VP <- vegan::varpart(dist_com, 
                     X = ~ecoregion, ~MAT_BC + MAP_BC + ATS_BC + APS_BC, ~tree_species,
                     data = site_val)

par(mar=c(1, 1, 1, 1))
plot(VP,
     Xnames = c("Ecoregion", "Climate", "Host"), 
     digits = 2)
savePlot("gg.jpeg")
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
#plot(nmds)

d1 <- data.frame(nmds$points[,]) %>%
  rownames_to_column("code") %>%
  left_join(site_val |>  rownames_to_column("code")) %>%
  select(code, NMDS1 = MDS1, NMDS2 = MDS2, NMDS3 = MDS3, 
         tree_species, ecoregion) 

d <- d1 %>%
  group_by(ecoregion) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()

d2 <- d %>%
  count(ecoregion) %>%
  select(-n) %>%
  mutate(ecoregion_abbr = c(
    "Alps", "Baltic", "Cantabrian", "Carpathian", "Celtic",
    "Central", "English", "Atlantic", "Pannonian", "Sarmatic",
     "Scandinavian", "Western"))

d <- d %>% left_join(d2)
view(d)

# Figure
ggplot(data = d, aes(y = NMDS2, x = NMDS1))  +
  facet_wrap(~ ecoregion, 
             labeller = label_wrap_gen(20)) +
  #geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = "none") +
  scale_color_brewer(palette="Paired") +
  geom_polygon(colour = "red",alpha = 0, lwd = 1) +
  geom_point(data = d1, 
             aes(y = NMDS2, x = NMDS1), alpha = 0.3)
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

corresp_code <- tibble(code_ecoregion = c(1:12)) |>
  mutate(ecoregion_abbr = c(
    "Alps", "Baltic", "Cantabrian", "Carpathian", "Celtic",
    "Central", "English", "Atlantic", "Pannonian", "Sarmatic",
    "Taiga", "Western"))

site_val_indval <- site_val |>
  rownames_to_column("code") |>
  mutate(code_ecoregion = as.numeric(as.factor(ecoregion)))|>
  left_join(corresp_code)

# Keep species that are present in three sites or more -> 481 species
OTU_sup2 <- OTU |> 
  select_if(colSums(OTU != 0) >= 2) 

# Keeps only OTUs with more than 5 tips in the whole dataset
OTU_sup_10_tips <- OTU_sup2 |> colSums() |> 
  data.frame() |> 
  rownames_to_column("SH") |>
  select(SH, n = colSums.OTU_sup2.) |>
  filter(n > 10)

OTU_indval <- OTU_sup2 |> 
  select(all_of(OTU_sup_10_tips$SH))
dim(OTU_indval)

indval = multipatt(OTU_indval, 
                   site_val_indval$code_ecoregion, 
                   func = "IndVal.g", 
                   # the .g is to correct the difference in site numbers between categories
                   duleg = TRUE, # to only test single ER (no combination)
                   control = how(nperm=50000), 
                   )

## Correct p-values (FDR) and bind the results to original output
fdr.p.value <- p.adjust(indval$sign[,"p.value"], method="fdr")
indval.all.fdr <- cbind(indval$sign, fdr.p.value)

## Omit NA values and print only those with p <= 0.05
attach(indval.all.fdr) ## allows sorting by header
indval.all.fdr.nona.sort <- indval.all.fdr[order(fdr.p.value, p.value, na.last=NA),] ## sort output by fdr, then p-value, omit "NA"
detach(indval.all.fdr)
view(indval.all.fdr.nona.sort |> 
       filter(stat > 0.5 & p.value <= 0.05)
     )

## the A component is the specificity. 
## Species with a high specificity (>0.8) are found mostly 
## or only in that group 
indval_a <- data.frame(indval$A[,1:12])
specificity <- data.frame(ecoregion=NA, specificity=NA)

for(i in 1:12){
  specificity[i, 1] <- colnames(indval_a)[i]
  specificity[i, 2] <- sum(indval_a[,i] > 0.8)
}

res_indval <- list(indval.all.fdr.nona.sort, specificity)
saveRDS(indval, file = "output_indval/indVal_OTUsup2.Rdata")

## Test for each tree 
iv_pine <- indval_group("PINE")
iv_spruce <- indval_group("SPRUCE")
iv_oak <- indval_group("OAK")
iv_beech <- indval_group("BEECH")
######## END #############