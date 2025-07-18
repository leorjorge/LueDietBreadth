library(tidyverse)
library(reshape2)
library(vegan)
library(picante)
library(ape)
library(hillR)
library(bbmle)
library(emmeans)
library(colorspace)
library(phytools)
library(dizzy)
library(cmdstanr)
library(phytools)
library(brms)
library(tidybayes)
library(patchwork)
library(ggtree)
library(beepr)
library(ggdist)
options(mc.cores = parallel::detectCores())
theme_set(theme_classic())


# Data entry ####
d.long <- read.csv("allIN_Vials_Clean.csv", row.names = 1)[,-c(8:10)]
d <- aggregate(HPD~Vial, data = d.long, FUN = function(x) sum (x == "P"))
colnames(d)[2] <- "P"
d$H <- aggregate(HPD~Vial, data = d.long, FUN = function(x) sum (x == "H"))$HPD
d$D <- aggregate(HPD~Vial, data = d.long, FUN = function(x) sum (x == "D"))$HPD
d <- cbind(d, d.long[match(d$Vial, d.long$Vial),c(1:3,5)])
C <- d[d$Parasitoid == "C",]
Hc <- aggregate(H ~ Host+Temp, data = C, FUN = mean)
d$Hc <- round(Hc$H[match(paste(d$Host,d$Temp), paste(Hc$Host, Hc$Temp))])
d <- d[d$Parasitoid != "C",]
d$Parasitoid[d$Parasitoid == "D"] <- "T"
d$Temp <- factor(d$Temp, levels = c("24", "20", "28"))

d$HcM <- pmax(d$H, d$P, d$Hc)
d$SP <- d$P/d$HcM
d$DI <- (d$HcM-d$H)/d$HcM

# Multi-parasitoid experiment
d.multi <- read.csv("Mix_ptoid_treatment.csv") 
d.multi.A <- aggregate(Ptoid~Rep+Host+Parasitoid+Temp, 
                       data = d.multi, FUN = function(x) sum (x == "MA"))
colnames(d.multi.A)[5] <- "P"
d.multi.A$Hc <- round(Hc$H[match(paste(d.multi.A$Host,d.multi.A$Temp), paste(Hc$Host, Hc$Temp))])
d.multi.A <- rbind(d.multi.A[,-1], d[d$Parasitoid == "A", colnames(d.multi.A[-1])])
d.multi.A$HcM <- pmax(d.multi.A$P, d.multi.A$Hc)

d.multi.G <- aggregate(Ptoid~Rep+Host+Parasitoid+Temp, 
                       data = d.multi, FUN = function(x) sum (x == "MG"))
colnames(d.multi.G)[5] <- "P"
d.multi.G$Hc <- round(Hc$H[match(paste(d.multi.G$Host,d.multi.G$Temp), paste(Hc$Host, Hc$Temp))])
d.multi.G <- rbind(d.multi.G[,-1], d[d$Parasitoid == "G", colnames(d.multi.G[-1])])
d.multi.G$HcM <- pmax(d.multi.G$P, d.multi.G$Hc)

d.multi.L <- aggregate(Ptoid~Rep+Host+Parasitoid+Temp, 
                       data = d.multi, FUN = function(x) sum (x == "ML"))
colnames(d.multi.L)[5] <- "P"
d.multi.L$Hc <- round(Hc$H[match(paste(d.multi.L$Host,d.multi.L$Temp), paste(Hc$Host, Hc$Temp))])
d.multi.L <- rbind(d.multi.L[,-1], d[d$Parasitoid == "L", colnames(d.multi.L[-1])])
d.multi.L$HcM <- pmax(d.multi.L$P, d.multi.L$Hc)

# Host phylogeny
Flies <- c("D.pallidifrons", "D.sulfurigastersulfurigaster", "D.rubida", "D.pseudoananassae", 
           "D.bipectinata", "D.pseudotakahashii", "D.birchii") # The phylogeny has multiple sulfurigaster "subspecies"
names(Flies) <- Flies
Phy <- read.tree(file = "Fig1_S1_data.treefile")
Phylo <- prune.missing(Flies, phylo = Phy)$tree
Phylo$tip.label <- c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB")
Phy.vcv <- vcv.phylo(Phylo)

# Control analysis ####
dC <- aggregate(HPD~Vial, data = d.long, FUN = function(x) sum (x == "H"))
dC[,3:6] <- d.long[match(dC$Vial, d.long$Vial),c(1:3,5)]
dC <- dC[dC$Parasitoid == "C",]
dC$Temp <- factor(dC$Temp, levels = c("24", "20", "28"))
dC$Host.Ph <- dC$Host
dC.Avail <- dcast(data = dC, formula = Temp ~ Host, value.var = "HPD", fun.aggregate = mean)
rownames(dC.Avail) <- dC.Avail$Temp
dC.Avail <- dC.Avail[,-1]/50

##Phylogenetic Model ####
m.ControlPhylo <- brm(HPD| trials(50) ~ Temp + (1 + Temp|Host) + 
                         (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                      family = binomial(),data = dC, data2 = list(Phy.vcv = Phy.vcv), 
                      chains = 8,iter = 2250, warmup = 1000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 20),
                      save_pars = save_pars(all = T))

#Contrasts
ControlPhylo_EMM <- emmeans(m.ControlPhylo, ~Temp|Host, type = "response", 
                            re_formula = NULL, ref = 1)
contrast(ControlPhylo_EMM, method = 'trt.vs.ctrl')

ControlPhylo_EMM_Temp <- emmeans(m.ControlPhylo, ~Temp, type = "response", 
                                 re_formula = NULL, ref = 1)
contrast(ControlPhylo_EMM_Temp, method = 'trt.vs.ctrl')

##Non-phylogenetic model ####
m.Control <- brm(HPD| trials(50) ~ Temp + (1 + Temp|Host),
                 family = binomial(),data = dC, 
                 chains = 8,iter = 2250, warmup = 1000, 
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 save_pars = save_pars(all = T))

#Contrasts
Control_EMM <- emmeans(m.Control, ~Temp|Host, type = "response", 
                       re_formula = NULL, ref = 1)
contrast(Control_EMM, method = 'trt.vs.ctrl')

Control_EMM_Temp <- emmeans(m.Control, ~Temp, type = "response", 
                            re_formula = NULL, ref = 1)
contrast(Control_EMM_Temp, method = 'trt.vs.ctrl')

# models for Parasitism Rate and Degree of Infestation ###################

d$Host.Ph <- d$Host
d.larval <- d |> 
   dplyr::filter(Parasitoid != "T")
d.pupal <- d |> 
   dplyr::filter(Parasitoid == "T")

## Parasitism Rate ####
### Phylotenetic models ####
#larval parasitoids
m.PRPhylo <- brm(P| trials(HcM) ~ Temp + (1 + Temp|Host) + (1 + Temp|Parasitoid) + 
                    (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                 family = binomial(),data = d.larval, data2 = list(Phy.vcv = Phy.vcv), 
                 chains = 8,iter = 2250, warmup = 1000, 
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 save_pars = save_pars(all = T))
m.PRPhylo <- add_criterion(m.PRPhylo, criterion = "loo", moment_match = TRUE)

PRPhylo_EMM <- emmeans(m.PRPhylo, ~Temp|Host+Parasitoid, type = "response", 
                       re_formula = NULL, ref = 1)
PR_contrast.larval <- summary(contrast(PRPhylo_EMM, method = 'trt.vs.ctrl'))
write.table(PR_contrast.larval, file = "Individual Contrasts PR larval.txt")

PR_EMM_Temp.larval <- emmeans(m.PRPhylo, trt.vs.ctrl~Temp, 
                              type = "response", re_formula = NULL, ref = 1)
# pupal parasitoids
m.PRPhylo.pupal <- brm(P| trials(HcM) ~ Temp + (1 + Temp|Host) + 
                          (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                       family = binomial(),data = d.pupal, data2 = list(Phy.vcv = Phy.vcv), 
                       chains = 8,iter = 2250, warmup = 1000, 
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       save_pars = save_pars(all = T))
m.PRPhylo.pupal <- add_criterion(m.PRPhylo.pupal, criterion = "loo", moment_match = TRUE)

PRPhylo_EMM.pupal <- emmeans(m.PRPhylo.pupal, ~Temp|Host, type = "response", 
                             re_formula = NULL, ref = 1)
PR_contrast.pupal <- summary(contrast(PRPhylo_EMM.pupal, method = 'trt.vs.ctrl'))
write.table(PR_contrast.pupal, file = "Individual Contrasts PR pupal.txt")
PR_EMM_Temp.pupal <- emmeans(m.PRPhylo.pupal, trt.vs.ctrl~Temp, type = "response", re_formula = NULL, ref = 1)
PR_contrasts <- PR_contrast.pupal |> 
   mutate(Parasitoid = "T") |> 
   bind_rows(PR_contrast.larval) |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

#Variance partitioning
summary(aov(lm(d$SP[d$Parasitoid != "T"]~d$Host[d$Parasitoid != "T"]+d$Parasitoid[d$Parasitoid != "T"] + d$Temp[d$Parasitoid != "T"])))
summary(aov(lm(d$SP[d$Parasitoid == "T"]~d$Host[d$Parasitoid == "T"]+ d$Temp[d$Parasitoid == "T"])))

###Non-Phylogenetic models ####
#larval parasitoids
m.PRnoPhylo <- brm(P| trials(HcM) ~ Temp + (1 + Temp|Host) + (1 + Temp|Parasitoid),
                   family = binomial(),data = d.larval, chains = 8,iter = 2250, warmup = 1000, 
                   control = list(adapt_delta = 0.99, max_treedepth = 20),
                   save_pars = save_pars(all = T))
m.PRnoPhylo <- add_criterion(m.PRnoPhylo, criterion = "loo", moment_match = TRUE)
loo_compare(m.PRPhylo, m.PRnoPhylo, criterion = "loo")

PRnoPhylo_EMM <- emmeans(m.PRnoPhylo, ~Temp|Host+Parasitoid, type = "response", 
                         re_formula = NULL, ref = 1)

PR_contrast.larvalnophylo <- summary(contrast(PRnoPhylo_EMM, method = 'trt.vs.ctrl'))
write.table(PR_contrast.larvalnophylo, file = "Individual Contrasts PR larval nophylo.txt")

# pupal parasitoid
m.PRnoPhylo.pupal <- brm(P| trials(HcM) ~ Temp + (1 + Temp|Host),
                         family = binomial(),data = d.pupal, chains = 8,iter = 2250, warmup = 1000, 
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         save_pars = save_pars(all = T))

m.PRnoPhylo.pupal <- add_criterion(m.PRnoPhylo.pupal, criterion = "loo", moment_match = TRUE)
loo_compare(m.PRPhylo.pupal, m.PRnoPhylo.pupal, criterion = "loo")

PRnoPhylo_EMM.pupal <- emmeans(m.PRnoPhylo.pupal, ~Temp|Host, type = "response", 
                               re_formula = NULL, ref = 1)
PR_contrast.pupal.noPhylo <- summary(contrast(PRnoPhylo_EMM.pupal, method = 'trt.vs.ctrl'))
write.table(PR_contrast.pupal.noPhylo, file = "Individual Contrasts PR pupal nophylo.txt")

PR_contrasts_nophylo <- PR_contrast.pupal.noPhylo |> 
   mutate(Parasitoid = "T") |> 
   bind_rows(PR_contrast.larvalnophylo) |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

##Degree of infestation####
###Phylogenetic models ####
##larval parasitoids
m.DIPhylo <- brm(HcM-H| trials(HcM) ~ Temp + (1 + Temp|Host) + (1 + Temp|Parasitoid) + 
                    (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                 family = binomial(),data = d.larval, data2 = list(Phy.vcv = Phy.vcv), 
                 chains = 8,iter = 2250, warmup = 1000, 
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 save_pars = save_pars(all = T), backend = "cmdstanr")
m.DIPhylo <- add_criterion(m.DIPhylo, criterion = "loo", moment_match = TRUE)

DIPhylo_EMM <- emmeans(m.DIPhylo, ~Temp|Host+Parasitoid, type = "response", 
                       re_formula = NULL, ref = 1)
DI_contrast.larval <- summary(contrast(DIPhylo_EMM, method = 'trt.vs.ctrl'))
write.table(DI_contrast.larval, file = "Individual Contrasts DI larval.txt")
DI_EMM_Temp.larval <- emmeans(m.DIPhylo, trt.vs.ctrl~Temp, 
                              type = "response", re_formula = NULL, ref = 1)

## pupal parasitoid
m.DIPhylo.pupal <- brm(HcM-H| trials(HcM) ~ Temp + (1 + Temp|Host) + 
                          (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                       family = binomial(),data = d.pupal, data2 = list(Phy.vcv = Phy.vcv), 
                       chains = 8,iter = 2250, warmup = 1000, 
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       save_pars = save_pars(all = T))
m.DIPhylo.pupal <- add_criterion(m.DIPhylo.pupal, criterion = "loo", moment_match = TRUE)

DIPhylo_EMM.pupal <- emmeans(m.DIPhylo.pupal, ~Temp|Host, type = "response", 
                             re_formula = NULL, ref = 1)
DI_contrast.pupal <- summary(contrast(DIPhylo_EMM.pupal, method = 'trt.vs.ctrl'))
write.table(DI_contrast.pupal, file = "Individual Contrasts DI pupal.txt")
DI_EMM_Temp.pupal <- emmeans(m.DIPhylo.pupal, trt.vs.ctrl~Temp, type = "response", re_formula = NULL, ref = 1)
DI_contrasts <- DI_contrast.pupal |> 
   mutate(Parasitoid = "T") |> 
   bind_rows(DI_contrast.larval)  |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

#Variance partitioning
summary(aov(lm(d$DI[d$Parasitoid != "T"]~d$Host[d$Parasitoid != "T"]+d$Parasitoid[d$Parasitoid != "T"] + d$Temp[d$Parasitoid != "T"])))
summary(aov(lm(d$DI[d$Parasitoid == "T"]~d$Host[d$Parasitoid == "T"]+ d$Temp[d$Parasitoid == "T"])))

###Non-Phylogenetic models ####
##larval parasitoids
m.DInoPhylo <- brm(HcM-H| trials(HcM) ~ Temp + (1 + Temp|Host) + (1 + Temp|Parasitoid),
                   family = binomial(),data = d.larval, chains = 8,iter = 2250, warmup = 1000, 
                   control = list(adapt_delta = 0.99, max_treedepth = 20),
                   save_pars = save_pars(all = T))
m.DInoPhylo <- add_criterion(m.DInoPhylo, criterion = "loo", moment_match = TRUE)
loo_compare(m.DIPhylo, m.DInoPhylo, criterion = "loo")

DInoPhylo_EMM <- emmeans(m.DInoPhylo, ~Temp|Host+Parasitoid, type = "response", 
                         re_formula = NULL, ref = 1)
DI_contrast.larval.noPhylo <- summary(contrast(DInoPhylo_EMM, method = 'trt.vs.ctrl'))
write.table(DI_contrast.larval.noPhylo, file = "Individual Contrasts DI larval nophylo.txt")

## pupal
m.DInoPhylo.pupal <- brm(HcM-H| trials(HcM) ~ Temp + (1 + Temp|Host),
                         family = binomial(),data = d.pupal, chains = 8,iter = 2250, warmup = 1000, 
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         save_pars = save_pars(all = T))
m.DInoPhylo.pupal <- add_criterion(m.DInoPhylo.pupal, criterion = "loo", moment_match = TRUE)
loo_compare(m.DIPhylo.pupal, m.DInoPhylo.pupal, criterion = "loo")

DInoPhylo_EMM.pupal <- emmeans(m.DInoPhylo.pupal, ~Temp|Host, type = "response", 
                               re_formula = NULL, ref = 1)
DI_contrast.pupal.nophylo <- summary(contrast(DInoPhylo_EMM.pupal, method = 'trt.vs.ctrl'))
write.table(DI_contrast.pupal.nophylo, file = "Individual Contrasts DI pupal nophylo.txt")

DI_contrasts.nophylo <- DI_contrast.pupal.nophylo |> 
   mutate(Parasitoid = "T") |> 
   bind_rows(DI_contrast.larval.noPhylo)  |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

## Multi-parasitoid models ####
### Phylogenetic models ####
d.multi.A$Host.Ph <- d.multi.A$Host
d.multi.G$Host.Ph <- d.multi.G$Host
d.multi.L$Host.Ph <- d.multi.L$Host

m.PRPhylo.multi.A <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host) + 
                            (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                         family = binomial(),data = d.multi.A, data2 = list(Phy.vcv = Phy.vcv), 
                         chains = 8,iter = 2250, warmup = 1000, 
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         save_pars = save_pars(all = T))
PRPhylo.multi.A_EMM <- emmeans(m.PRPhylo.multi.A, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                               re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

m.PRPhylo.multi.G <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host) + 
                            (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                         family = binomial(),data = d.multi.G, data2 = list(Phy.vcv = Phy.vcv), 
                         chains = 8,iter = 2250, warmup = 1000, 
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         save_pars = save_pars(all = T))
PRPhylo.multi.G_EMM <- emmeans(m.PRPhylo.multi.G, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                               re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

m.PRPhylo.multi.L <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host) + 
                            (1 + Temp|gr(Host.Ph, cov = Phy.vcv)),
                         family = binomial(),data = d.multi.L, data2 = list(Phy.vcv = Phy.vcv), 
                         chains = 8,iter = 2250, warmup = 1000, 
                         control = list(adapt_delta = 0.99, max_treedepth = 20),
                         save_pars = save_pars(all = T))
PRPhylo.multi.L_EMM <- emmeans(m.PRPhylo.multi.L, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                               re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

### Non-phylogenetic models ####
m.PR.multi.A <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host),
                    family = binomial(),data = d.multi.A, 
                    chains = 8,iter = 2250, warmup = 1000, 
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    save_pars = save_pars(all = T))
PR.multi.A_EMM <- emmeans(m.PR.multi.A, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                          re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

m.PR.multi.G <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host),
                    family = binomial(),data = d.multi.G, 
                    chains = 8,iter = 2250, warmup = 1000, 
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    save_pars = save_pars(all = T))
PR.multi.G_EMM <- emmeans(m.PR.multi.G, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                          re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

m.PR.multi.L <- brm(P| trials(HcM) ~ Temp * Parasitoid + (1 + Temp + Parasitoid|Host),
                    family = binomial(),data = d.multi.L, 
                    chains = 8,iter = 2250, warmup = 1000, 
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    save_pars = save_pars(all = T))
PR.multi.L_EMM <- emmeans(m.PR.multi.L, trt.vs.ctrl~Parasitoid|Host + Temp, type = "response", 
                          re_formula = NULL, ref = 1)$contrasts |> 
   summary() |> 
   mutate(Host = fct_relevel(Host, "RUB", "SUL", "PAL", "BIP", "PSA", "PST", "BIR"))

# Specialization ###########

d.Network <- dcast(data = d, formula = Temp + Parasitoid ~ Host, value.var = "P", fun.aggregate = sum)
Phylo.star <- compute.brlen(starTree(species = levels(as.factor(d$Host))))

Spec <- d |>  
   group_by(Host,Temp,Parasitoid) %>% 
   mutate(n = row_number()) %>%
   dcast(formula = Parasitoid + Temp + Host ~ n, value.var = "P") %>%
   group_by(Parasitoid, Temp) |> 
   nest() |> 
   mutate(data = lapply(data, function(df) {
      mat <- t(df[,c(2:8)])
      colnames(mat) <- df$Host
      return(mat)
   })) |> 
   mutate(Shannon = lapply(data, hill_phylo, tree = Phylo.star, q = 1),
          Phylo = lapply(data, hill_phylo, tree = Phylo, q = 1)) |> 
   select(-data) |> 
   unnest_longer(col = c(Shannon, Phylo), indices_include = F)

## Specialization models ####
### Phylogenetic diversity ####
# all species
PhyloSpec_Mod.full <- lm(Phylo ~ Parasitoid*Temp, data = Spec)
PhyloSpec_Mod.add <- lm(Phylo ~ Parasitoid+Temp, data = Spec)
PhyloSpec_Mod.Temp <- lm(Phylo ~ Temp, data = Spec)
PhyloSpec_Mod.Para <- lm(Phylo ~ Parasitoid, data = Spec)
PhyloSpec_Mod.null <- lm(Phylo ~ 1, data = Spec)
AICctab(PhyloSpec_Mod.full, PhyloSpec_Mod.add, PhyloSpec_Mod.Temp, PhyloSpec_Mod.Para, PhyloSpec_Mod.null)
emmeans(PhyloSpec_Mod.full, trt.vs.ctrl ~ Temp|Parasitoid)

#larval parasitoids
PhyloSpec_Mod.larval.full <- lm(Phylo ~ Parasitoid*Temp, data = Spec[Spec$Parasitoid != "T",])
emmeans(PhyloSpec_Mod.larval.full, trt.vs.ctrl ~ Temp|Parasitoid)

#pupal parasitoid
PhyloSpec_Mod.pupal.full <- lm(Phylo ~ Temp, data = Spec[Spec$Parasitoid == "T",])
emmeans(PhyloSpec_Mod.pupal.full, trt.vs.ctrl ~ Temp)

### Shannon diversity ####
# all species
DivSpec_Mod.full <- lm(Shannon ~ Parasitoid*Temp, data = Spec)
DivSpec_Mod.add <- lm(Shannon ~ Parasitoid+Temp, data = Spec)
DivSpec_Mod.Temp <- lm(Shannon ~ Temp, data = Spec)
DivSpec_Mod.Para <- lm(Shannon ~ Parasitoid, data = Spec)
DivSpec_Mod.null <- lm(Shannon ~ 1, data = Spec)
AICctab(DivSpec_Mod.full, DivSpec_Mod.add, DivSpec_Mod.Temp, DivSpec_Mod.Para, DivSpec_Mod.null)
emmeans(DivSpec_Mod.full, trt.vs.ctrl ~ Temp|Parasitoid)

#larval parasitoids
DivSpec_Mod.larval.full <- lm(Shannon ~ Parasitoid*Temp, data = Spec[Spec$Parasitoid != "T",])
emmeans(DivSpec_Mod.larval.full, trt.vs.ctrl ~ Temp|Parasitoid)

#pupal parasitoid
DivSpec_Mod.pupal.full <- lm(Shannon ~ Temp, data = Spec[Spec$Parasitoid == "T",])
emmeans(DivSpec_Mod.pupal.full, trt.vs.ctrl ~ Temp)

## sensitivity analysis of specialization replicates ####
Spec.rep <- vector(mode = "list", length = 1000)
for(i in 1:length(Spec.rep)){
   Spec.rep[[i]] <- d |>  
      group_by(Host,Temp,Parasitoid) %>% 
      mutate(n = sample(row_number())) %>%
      dcast(formula = Parasitoid + Temp + Host ~ n, value.var = "P") %>%
      group_by(Parasitoid, Temp) |> 
      nest() |> 
      mutate(data = lapply(data, function(df) {
         mat <- t(df[,c(2:8)])
         colnames(mat) <- df$Host
         return(mat)
      })) |> 
      mutate(Shannon = lapply(data, hill_phylo, tree = Phylo.star, q = 1),
             Phylo = lapply(data, hill_phylo, tree = Phylo, q = 1)) |> 
      select(-data) |> 
      unnest_longer(col = c(Shannon, Phylo), indices_include = F)
}

Spec.rep <- list_rbind(Spec.rep, names_to = "rep")

### Phylogenetic diversity ####
Spec.error.phylo <- Spec.rep |> 
   select(-Shannon) |>  
   group_by(Parasitoid, Temp) |> 
   summarise(mean.Phylo = median(Phylo),
             lower = quantile(Phylo, probs = 0.25),
             upper = quantile(Phylo, probs = 0.75)) |> 
   mutate(Temp = fct_relevel(Temp, c("20", "24", "28")))

### Shannon diversity ####
Spec.error.shannon <- Spec.rep |> 
   select(-Phylo) |>  
   group_by(Parasitoid, Temp) |> 
   summarise(mean.Shannon = median(Shannon),
             lower = quantile(Shannon, probs = 0.25),
             upper = quantile(Shannon, probs = 0.75)) |> 
   mutate(Temp = fct_relevel(Temp, c("20", "24", "28")))

## Available diversity in control ####
dC.W <- dC |> 
   arrange(Temp, Host) |> 
   group_by(Host,Temp) |> 
   mutate(n = row_number()) |> 
   select(HPD, Host, Temp, n) |> 
   dplyr::filter(n <= 7) |> 
   pivot_wider(names_from = c(Temp, Host), values_from = HPD, names_glue = "{Temp}_{Host}") %>%
   select(-1)
### Phylogenetic diversity ####
dC.phylo <- data.frame(Temp = levels(d$Temp))
dC.phylo[,2:8] <- NA
colnames(dC.phylo)[2:8] <- 1:7
for(i in 1:3){
   df <- dC.W[,(1+(i-1)*7):(i*7)]
   colnames(df) <- levels(as.factor(dC$Host))
   dC.phylo[i,2:8] <- hill_phylo(df, tree = Phylo, q = 1)
}

dC.Spec <- pivot_longer(dC.phylo, cols = -1, names_to = NULL, values_to = "Phylo") |> 
   mutate(Temp = fct_relevel(Temp, "20", "24", "28"))

dC.Spec.test <- emmeans(lm(Phylo~Temp, data = dC.Spec), ~Temp)
contrast(dC.Spec.test, method = "trt.vs.ctrl", ref = 2)

###Shannon diversity ####
dC.nophylo <- data.frame(Temp = levels(d$Temp))
dC.nophylo[,2:8] <- NA
colnames(dC.nophylo)[2:8] <- 1:7
for(i in 1:3){
   df <- dC.W[,(1+(i-1)*7):(i*7)]
   colnames(df) <- levels(as.factor(dC$Host))
   dC.nophylo[i,2:8] <- hill_phylo(df, tree = Phylo.star, q = 1)
}

dC.Spec.nophylo <- pivot_longer(dC.nophylo, cols = -1, names_to = NULL, values_to = "Spec") |> 
   mutate(Temp = fct_relevel(Temp, "20", "24", "28"))

dC.Spec.test.nophylo <- emmeans(lm(Spec~Temp, data = dC.Spec.nophylo), ~Temp)
contrast(dC.Spec.test.nophylo, method = "trt.vs.ctrl", ref = 2)

##DSI ####
### Phylogenetic ####
DSI <- d.Network[,1:2]
DSI$DSI.Phyl <- NA
DSI$DSI.Phyl[1:4] <- dsi(Int = d.Network[1:4,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo), diag = T), Abund = t(dC.Avail[1,]+1))$DSIstar
DSI$DSI.Phyl[5:8] <- dsi(Int = d.Network[5:8,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo), diag = T), Abund = t(dC.Avail[2,]+1))$DSIstar
DSI$DSI.Phyl[9:12] <- dsi(Int = d.Network[9:12,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo), diag = T), Abund = t(dC.Avail[3,]+1))$DSIstar

###Non-Phylogenetic ####
DSI$DSI.noPhyl <- NA
DSI$DSI.noPhyl[1:4] <- dsi(Int = d.Network[1:4,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo.star), diag = T), Abund = t(dC.Avail[1,]+1))$DSIstar
DSI$DSI.noPhyl[5:8] <- dsi(Int = d.Network[5:8,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo.star), diag = T), Abund = t(dC.Avail[2,]+1))$DSIstar
DSI$DSI.noPhyl[9:12] <- dsi(Int = d.Network[9:12,-c(1:2)], Dist = as.dist(cophenetic.phylo(Phylo.star), diag = T), Abund = t(dC.Avail[3,]+1))$DSIstar


# Encapsulation #################
encap <- read.csv("encapsulation.csv") |> 
   drop_na(flies_checked_n) |> 
   mutate(vial = 1:n())

encap.pupal <- read.csv("encapsulation.csv") |> 
   drop_na(flies_checked_n) |> 
   dplyr::filter(Parasitoid %in% c("C", "D")) |> 
   mutate(vial = 1:n())

mod.encap <- brm(melanized_n| trials(flies_checked_n) ~ Temp + (1 + Temp|Parasitoid),
                 family = binomial(), data = encap, 
                 chains = 8,iter = 2250, warmup = 1000, 
                 control = list(adapt_delta = 0.99, max_treedepth = 20),
                 save_pars = save_pars(all = T), backend = "cmdstanr")

mod.encap.pupal <- brm(melanized_n| trials(flies_checked_n) ~ Temp + (1 + Temp|Parasitoid),
                       family = binomial(), data = encap.pupal, 
                       chains = 8,iter = 2250, warmup = 1000, 
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       save_pars = save_pars(all = T), backend = "cmdstanr")

emm.encap <- emmeans(mod.encap, ~ Parasitoid+Temp, type = "response", re_formula = NULL, ref = 2)
encap.contrast.temp <- as.data.frame(contrast(contrast(
   emm.encap, method = "trt.vs.ctrl", ref = 2, by = "Temp", exclude = 3), 
   by = "contrast", method = "pairwise"))
print(encap.contrast.temp[,], digits = 4)

emm.encap.notemp <- emmeans(mod.encap, ~ Parasitoid, type = "response", re_formula = NULL, ref = 2)
encap.contrast.parasitoids <- as.data.frame(contrast(
   contrast(emm.encap.notemp, method = "trt.vs.ctrl", ref = 2, exclude = 3), method = "pairwise"))
print(encap.contrast.parasitoids[,], digits = 4)

emm.encap.pupal <- emmeans(mod.encap.pupal, ~ Parasitoid+Temp, type = "response", re_formula = NULL)
contrast(contrast(emm.encap.pupal, method = "trt.vs.ctrl", by = "Temp", exclude = 3), by = "contrast", method = "pairwise")

emm.encap.notemp.pupal <- emmeans(mod.encap.pupal, ~ Parasitoid, type = "response", re_formula = NULL)
contrast(emm.encap.notemp.pupal, method = "trt.vs.ctrl")

