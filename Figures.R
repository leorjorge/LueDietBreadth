## Phylogeny used in most figures
PhyloPlot <- ggtree(Phylo, ladderize = F) + 
   geom_tiplab(align = T)
ggsave("Phylo.svg", width = 7, height = 10)

# Fig. 1 Raw data heatmaps ####
#PR
PR_estimate <- d |> 
   select(P, Host, Parasitoid, Temp, HcM) |> 
   mutate(Temp = factor(Temp, levels = c(20,24,28))) |> 
   mutate(Host = factor(Host, levels = c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))) |> 
   group_by(Host, Parasitoid, Temp) |> 
   summarise(prob = mean(P/HcM, na.rm = T)) |> 
   mutate(prob = replace_na(prob,0)) |> 
   ungroup()

Heatmap_PR <- ggplot(PR_estimate, mapping = aes(x = Parasitoid, y = Host, fill = prob))+
   scale_x_discrete("parasitoids", expand = c(0, 0)) + 
   scale_y_discrete(NULL, expand = c(0, 0)) + 
   geom_tile(color = "black", linewidth = 0.5) +
   coord_fixed() + 
   facet_wrap(~Temp, labeller = as_labeller(c("20" = "cooling (20°)",
                                              "24" = "ambient (24°)",
                                              "28" = "warming (28°)"))) +
   scale_fill_continuous_sequential(limits = c(0,1), palette = "Reds 2",
                                    breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
   labs(fill = "PR") +
   theme(legend.position = "right",
         strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text.x = element_text(size = 14),
         #         axis.text.y = element_blank(),
         axis.title = element_text(size = 16))
ggsave("Fig1a_PRHeatmap.svg", plot = Heatmap_PR, height = 5.5, width = 8.5)

DI_estimate <- d |> 
   select(H, Host, Parasitoid, Temp, HcM) |> 
   mutate(Temp = factor(Temp, levels = c(20,24,28))) |> 
   mutate(Host = factor(Host, levels = c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))) |> 
   group_by(Host, Parasitoid, Temp) |> 
   summarise(prob = mean((HcM-H)/HcM, na.rm = T)) |> 
   mutate(prob = replace_na(prob,0)) |> 
   ungroup()

Heatmap_DI <- Heatmap_PR %+% DI_estimate + 
   labs(fill = "DI") +
   scale_fill_continuous_sequential(limits = c(0,1), palette = "Greens 2",
                                    breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"))

ggsave("Fig1b_DIHeatmap.svg", plot = Heatmap_DI, height = 5.5, width = 8.5)

# Fig. 2 Control survival ####
Ctrl_Heatmap <- dC |> 
   mutate(Host = factor(Host, levels = c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))) |> 
   group_by(Host, Temp) |> 
   summarise(surv = mean(HPD)/50) |> 
   mutate(Temp = factor(Temp, levels = c(20,24,28))) |> 
   ungroup() |> 
   ggplot(mapping = aes(x = Temp, y = Host, fill = surv))+
   scale_x_discrete(expand = c(0, 0)) + 
   scale_y_discrete(NULL, expand = c(0, 0)) + 
   geom_tile(color = "black", linewidth = 0.5) +
   coord_fixed() + 
   scale_fill_continuous_sequential(limits = c(0,1), palette = "Purples 2",
                                    breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
   labs(x = "temperature", fill = "") +
   theme(legend.position = "right",
         axis.text.x = element_text(size = 14),
         axis.text.y = element_blank(),
         axis.title = element_text(size = 16))
ggsave("Fig2_ControlHeatmap.svg", plot = Ctrl_Heatmap)

# Fig. 3 Specialization ####

Spec.plot <- Spec
Spec.plot$Temp <- factor(Spec.plot$Temp, levels = c("20", "24", "28"))

pPhylo.spec <- ggplot(Spec.plot, aes(x = Parasitoid, y = Phylo, fill = Temp)) + 
   geom_dots(side = "both", position  = position_dodge(width = 0.8), 
             linewidth = 0, binwidth = 0.04, alpha = 0.5) + 
   geom_boxplot(outliers = F, coef = 0, alpha = 0.5, position  = position_dodge(width = 0.8), 
                aes(colour = Temp)) +
   scale_fill_discrete_diverging(palette = "Berlin",) + 
   scale_color_discrete_diverging(palette = "Berlin",) + 
   labs(x = "parasitoids", 
        y = "phylogenetic diversity",
        fill = "temp",
        colour = "temp",
        title = "A: host diversity") +
   coord_cartesian(ylim = c(0.5, 3.5))

pAvail <- ggplot(dC.Spec, aes(x = Temp, y = Phylo, fill = Temp)) + 
   geom_dots(side = "both", linewidth = 0, alpha = 0.5) + 
   geom_boxplot(outliers = F, coef = 0, alpha = 0.5, aes(colour = Temp)) +
   scale_fill_discrete_diverging(palette = "Berlin",) + 
   scale_color_discrete_diverging(palette = "Berlin",) + 
   labs(x = "temperature", 
        y = "phylogenetic diversity",
        title = "B: available diversity") +
   theme(legend.position="none") +
   coord_cartesian(ylim = c(0.5, 3.5))
DSI.plot <- DSI |> 
   mutate(Temp = fct_relevel(Temp, "20", "24", "28")) |>  
   mutate(Parasitoid = fct_relevel(Parasitoid, "A", "G", "L", "T"))

pDSI <- ggplot(DSI.plot, mapping = aes(x = Parasitoid, y = DSI.Phyl, colour = Temp)) +
   geom_point(position = position_dodge(width = 0.8)) +
   scale_color_discrete_diverging(palette = "Berlin",) +
   labs(y = "DSI*",
        title = "C: DSI*") + 
   coord_cartesian(ylim = c(-1,1)) +
   theme(legend.position="none")

Fig3 <- pPhylo.spec+pAvail+pDSI + 
   plot_layout(guides = "collect", widths = c(2,1,1))

ggsave("Fig3.svg", plot = Fig3, height = 4, width = 8)

# Fig. S1 Contrasts ####
S1_PR <- ggplot(data = PR_contrasts, mapping = aes(x = Parasitoid, y = odds.ratio, colour = Host)) + 
   geom_point(position = position_dodge(width = 0.9))+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), position = position_dodge(width = 0.9)) +
   geom_hline(yintercept = 1, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   labs(y = "odds ratio") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16))
ggsave("S1_PR.svg", S1_PR, width = 9, height = 4.5)

S1_DI <- S1_PR %+% DI_contrasts
ggsave("S1_DI.svg", S1_DI, width = 9, height = 4.5)

#Fig. S2 Control contrasts####
Contrast_Control <- summary(contrast(ControlPhylo_EMM, method = 'trt.vs.ctrl')) |> 
   mutate(Host = factor(Host, levels = rev(c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))))

S2_Ctrl <- ggplot(data = Contrast_Control, mapping = aes(x = Host, y = odds.ratio)) + 
   geom_point()+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD)) +
   geom_hline(yintercept = 1, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   scale_x_discrete(labels = c(BIR = "*D. birchii*", 
                               PST = "*D. pseudotakahashii*" , 
                               PSA = "*D. pseudoananassae*", 
                               BIP = "*D. bipectinata*", 
                               PAL = "*D. pallidifrons*", 
                               SUL = "*D. sulfurigaster*", 
                               RUB = "*D. rubida*")) +
   labs(y = "odds ratio") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         axis.text.x =  ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1))
ggsave("S2_Control.svg", S2_Ctrl, width = 9, height = 4.5)

# Fig. S3 Sensitivity ####
Spec.error.phylo.plot <- ggplot(Spec.error.phylo, aes(x = Parasitoid, y = mean.Phylo, colour = Temp)) + 
   geom_point(position  = position_dodge(width = 0.8)) + 
   geom_errorbar(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 0.8), width = 0) +
   scale_color_discrete_diverging(palette = "Berlin",) + 
   labs(x = "parasitoids", 
        y = "phylogenetic diversity",
        fill = "temperature",
        colour = "temperature") +
   coord_cartesian(ylim = c(0.5, 3.5))
ggsave("Specialization_error.svg", plot = Spec.error.plot, width = 6, height = 4)

# Fig. S4 Heatmap counts####
P_estimate <- d |> 
   select(P, Host, Parasitoid, Temp) |> 
   mutate(Temp = factor(Temp, levels = c(20,24,28))) |> 
   mutate(Host = factor(Host, levels = c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))) |> 
   group_by(Host, Parasitoid, Temp) |> 
   summarise(prob = sum(P)/n()) |> 
   mutate(prob = round(prob, 1))

Heatmap_P <- Heatmap_PR %+% P_estimate + 
   geom_text(aes(label = prob)) +
   labs(fill = "mean parasitoids") +
   scale_y_discrete(labels = c(BIR = "*D. birchii*", 
                               PST = "*D. pseudotakahashii*" , 
                               PSA = "*D. pseudoananassae*", 
                               BIP = "*D. bipectinata*", 
                               PAL = "*D. pallidifrons*", 
                               SUL = "*D. sulfurigaster*", 
                               RUB = "*D. rubida*"),
                    expand = expansion()) +
   scale_fill_continuous_sequential(limits = c(0,50), palette = "Reds 2",
                                    breaks = c(0, 25, 50)) +
   theme(axis.text.y = ggtext::element_markdown())

ggsave("S4_Heatmap_counts.svg", plot = Heatmap_P, height = 5.5, width = 8.5)

# Fig. S5 Contrasts without phylogeny ####
S5_PR <- ggplot(data = PR_contrasts_nophylo, mapping = aes(x = Parasitoid, y = odds.ratio, colour = Host)) + 
   geom_point(position = position_dodge(width = 0.9))+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), position = position_dodge(width = 0.9)) +
   geom_hline(yintercept = 1, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   labs(y = "odds ratio", 
        title = "Parasitism success") +
   scale_colour_discrete(labels = c(BIR = "*D. birchii*", 
                                    PST = "*D. pseudotakahashii*" , 
                                    PSA = "*D. pseudoananassae*", 
                                    BIP = "*D. bipectinata*", 
                                    PAL = "*D. pallidifrons*", 
                                    SUL = "*D. sulfurigaster*", 
                                    RUB = "*D. rubida*")) +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = ggtext::element_markdown(size = 14),
         axis.title = element_text(size = 16))

S5_DI <- S5_PR %+% 
   DI_contrasts.nophylo + 
   labs(y = "odds ratio", 
        title = "Degree of infestation")

S5 <- S5_PR / S5_DI 
ggsave("S5_contrasts_nophylo.svg", S5, width = 9, height = 9)

# Fig. S6 Control contrasts without phylogeny ####
Contrast_Control_nophylo <- summary(contrast(Control_EMM, method = 'trt.vs.ctrl')) |> 
   mutate(Host = factor(Host, levels = rev(c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))))
write.table(Contrast_Control_nophylo, "Control survival contrasts.txt")
S6_Ctrl_nophylo <- ggplot(data = Contrast_Control_nophylo, mapping = aes(x = Host, y = odds.ratio)) + 
   geom_point()+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD)) +
   geom_hline(yintercept = 1, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   scale_x_discrete(labels = c(BIR = "*D. birchii*", 
                               PST = "*D. pseudotakahashii*" , 
                               PSA = "*D. pseudoananassae*", 
                               BIP = "*D. bipectinata*", 
                               PAL = "*D. pallidifrons*", 
                               SUL = "*D. sulfurigaster*", 
                               RUB = "*D. rubida*")) +
   labs(y = "odds ratio") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         axis.text.x =  ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1))
ggsave("S6_Control_nophylo.svg", S6_Ctrl_nophylo, width = 9, height = 4.5)

# Fig. S7 Specialization without phylogeny ####

p.Shannon.spec <- ggplot(Spec.plot, aes(x = Parasitoid, y = Shannon, fill = Temp)) + 
   geom_dots(side = "both", position  = position_dodge(width = 0.8), 
             linewidth = 0, binwidth = 0.09, alpha = 0.5) + 
   geom_boxplot(outliers = F, coef = 0, alpha = 0.5, position  = position_dodge(width = 0.8), 
                aes(colour = Temp)) +
   scale_fill_discrete_diverging(palette = "Berlin",) + 
   scale_color_discrete_diverging(palette = "Berlin",) + 
   labs(x = "parasitoids", 
        y = "Shannon diversity",
        fill = "temp",
        colour = "temp",
        title = "A: host diversity")

pAvail.nophylo <- ggplot(dC.Spec.nophylo, aes(x = Temp, y = Spec, fill = Temp)) + 
   geom_dots(side = "both", linewidth = 0, alpha = 0.5) + 
   geom_boxplot(outliers = F, coef = 0, alpha = 0.5, aes(colour = Temp)) +
   scale_fill_discrete_diverging(palette = "Berlin",) + 
   scale_color_discrete_diverging(palette = "Berlin",) + 
   labs(x = "temperature", 
        y = "Shannon diversity",
        title = "B: available diversity") +
   theme(legend.position="none") +
   coord_cartesian(ylim = c(1, 6.9))

DSI.plot.nophylo <- DSI |> 
   mutate(Temp = fct_relevel(Temp, "20", "24", "28")) |> 
   mutate(Parasitoid = fct_relevel(Parasitoid, "A", "G", "L", "T"))
pDSI.nophylo <- ggplot(DSI.plot.nophylo, mapping = aes(x = Parasitoid, y = DSI.noPhyl, colour = Temp)) +
   geom_point(position = position_dodge(width = 0.8)) +
   scale_color_discrete_diverging(palette = "Berlin",) +
   labs(y = "SSI*",
        title = "C: SSI*",
        x = "parasitoids") + 
   coord_cartesian(ylim = c(-1,1)) +
   theme(legend.position="none")

FigS7.Specnophylo <- p.Shannon.spec+pAvail.nophylo+pDSI.nophylo + 
   plot_layout(guides = "collect", widths = c(2,1,1))
ggsave("FigS7_Specializationnophylo.svg", plot = FigS7.Specnophylo, height = 4, width = 8)

# Fig. S8 Multiple infections ####
S8_multi.A <- ggplot(data = PRPhylo.multi.A_EMM, mapping = aes(x = Host, y = odds.ratio, colour = Temp)) + 
   geom_point(position = position_dodge(width = 0.9))+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), position = position_dodge(width = 0.9)) +
   geom_hline(yintercept = 1, linetype = 3) +
   labs(y = "odds ratio") +
   scale_x_discrete(labels = c(BIR = "*D. birchii*", 
                               PST = "*D. pseudotakahashii*" , 
                               PSA = "*D. pseudoananassae*", 
                               BIP = "*D. bipectinata*", 
                               PAL = "*D. pallidifrons*", 
                               SUL = "*D. sulfurigaster*", 
                               RUB = "*D. rubida*")) +
   scale_color_discrete_diverging(palette = "Berlin") +
   labs(title = "A: *Asobara sp.*", colour = "Temperature") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = ggtext::element_markdown(size = 14),
         axis.text.x = ggtext::element_markdown(size = 10, angle = 45, vjust = 1, hjust = 1),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         title = ggtext::element_markdown())
S8_multi.G <- S8_multi.A %+% PRPhylo.multi.G_EMM +
   labs(title = "B: *Ganapsis sp.*")
S8_multi.L <- S8_multi.A %+% PRPhylo.multi.L_EMM +
   labs(title = "C: *Leptopilina sp.*")

S8_multi <- S8_multi.A / S8_multi.G / S8_multi.L + plot_layout(guides = "collect") & 
   theme(legend.position='bottom')
ggsave("S8_multi.svg", plot = S8_multi, height = 10, width = 4.5)


# Fig. S9 Control fly weight ####

Ctrl_weight_Heatmap <- Fly.Wt |> 
   dplyr::filter(Parasitoid == "C") |> 
   mutate(Host = factor(Host, levels = c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))) |> 
   group_by(Host, Temp) |> 
   summarise(body_weight = mean(body_weight/10)) |> 
   mutate(Temp = factor(Temp, levels = c(20,24,28))) |> 
   ungroup() |> 
   ggplot(mapping = aes(x = Temp, y = Host, fill = body_weight))+
   scale_x_discrete(expand = c(0, 0)) + 
   scale_y_discrete(NULL, expand = c(0, 0)) + 
   geom_tile(color = "black", linewidth = 0.5) +
   coord_fixed() + 
   scale_fill_continuous_sequential(palette = "Purples 2") +
   labs(x = "temperature", fill = "body mass (mg)") +
   theme(legend.position = "right",
         axis.text.x = element_text(size = 14),
         axis.text.y = element_blank(),
         axis.title = element_text(size = 16))
ggsave("FigS9_Control_weightHeatmap.svg", plot = Ctrl_weight_Heatmap)

#Fig. S10 Control fly weight contrasts ####
Contrast_Control_weight <- summary(contrast(EMM.Fly.Wt.temp, method = 'trt.vs.ctrl',, ref = 2)) |> 
   mutate(Host = factor(Host, levels = rev(c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))),
          estimate = estimate/10,
          lower.HPD = lower.HPD/10,
          upper.HPD = upper.HPD/10) |> 
   dplyr::filter(Parasitoid == "C") |> 
   dplyr::filter(!(Host == "PST" & contrast == "Temp28 - Temp24"))

S10_Ctrl_weight <- ggplot(data = Contrast_Control_weight, mapping = aes(x = Host, y = estimate)) + 
   geom_point()+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD)) +
   geom_hline(yintercept = 0, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   scale_x_discrete(labels = c(BIR = "*D. birchii*", 
                               PST = "*D. pseudotakahashii*" , 
                               PSA = "*D. pseudoananassae*", 
                               BIP = "*D. bipectinata*", 
                               PAL = "*D. pallidifrons*", 
                               SUL = "*D. sulfurigaster*", 
                               RUB = "*D. rubida*")) +
   labs(y = "body mass difference (mg)") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16),
         axis.text.x =  ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1))
ggsave("FigS10_Control_weight.svg", S10_Ctrl_weight, width = 9, height = 4.5)

# Fig. S11 Fly-Parasitoid mass correlation ####

S11_ParFlyWt <- ParFly.Wt.pred |> 
   mutate(Parasitoid = factor(Parasitoid, levels = c("A", "G", "L", "D"), labels = c("A", "G", "L", "T")),
          .epred = .epred/10) |> 
   ggplot(mapping = aes(x = fly_weight/10, y = body_weight/10, color = Parasitoid, fill = Parasitoid)) +
   stat_lineribbon(aes(y = .epred), .width = c(.95, .68), alpha = 1/5, linewidth = 0.7) +
   geom_point(data = ParFly.Wt, position = position_jitter(height = 0)) + 
   labs(x = "host body mass (mg)", y = "parasitoid body mass (mg)") +
   scale_color_discrete() +
   scale_fill_discrete() +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16))

ggsave("S11_parfly_mass.svg", S11_ParFlyWt)


# Fig. S12 Temperature differences parasitoids ####
Par.counts <- Par.Wt |> 
   count(Host, Parasitoid, Temp, .drop = F) |> 
   pivot_wider(names_from = Temp, values_from = n) |> 
   mutate(Temp20 = `20` * `24`,
          Temp28 = `20` * `28`) |> 
   select(!c(`20`, `24`, `28`)) |> 
   pivot_longer(cols = Temp20:Temp28, names_to = "contrast", values_to = "present") |> 
   mutate(present = present > 0,
          contrast = paste(contrast, "- Temp24"),
          Parasitoid = fct_recode(Parasitoid, T = "D"))
S12_Pars_tempcontrasts <- summary(Contrast.Par.Wt.temp) |> 
   mutate(Parasitoid = fct_recode(Parasitoid, T = "D"),
          estimate = estimate/10,
          lower.HPD = lower.HPD/10,
          upper.HPD = upper.HPD/10) |> 
   left_join(Par.counts, by = join_by(Host, Parasitoid, contrast)) |> 
   dplyr::filter(present == T) |> 
   mutate(Host = factor(Host, levels = rev(c("BIR", "PST", "PSA", "BIP", "PAL", "SUL", "RUB"))),
          Parasitoid = fct_relevel(Parasitoid, "A", "G", "L", "T")) |> 
   ggplot(mapping = aes(x = Parasitoid, y = estimate, colour = Host)) + 
   geom_point(position = position_dodge(width = 0.9)) +
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), position = position_dodge(width = 0.9)) +
   geom_hline(yintercept = 0, linetype = 3) +
   facet_wrap(~contrast, axes = "all_y") + 
   labs(y = "body mass difference (mg)") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16))

ggsave("S12_masschange_temperature.svg", S12_Pars_tempcontrasts)

# Fig. S13 Encapsulation ####
emm.encap.plot <- as.data.frame(emm.encap) |> 
   dplyr::filter(Parasitoid != "C") |> 
   mutate(Temp = factor(Temp, levels = c("20", "24", "28")),
          Parasitoid = factor(Parasitoid, levels = c("A", "G", "L", "D")))

encap.plot <- ggplot(emm.encap.plot, aes(x = Parasitoid, y = prob, colour = Temp)) + 
   geom_point(position  = position_dodge(width = 0.8)) + 
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD),
                 position = position_dodge(width = 0.8), width = 0) +
   scale_color_discrete_diverging(palette = "Berlin", nmax = 3) + 
   scale_y_continuous(labels = scales::label_percent()) +
   labs(x = "parasitoids", 
        y = "melanization",
        fill = "temperature",
        colour = "temperature") +
   coord_cartesian(ylim = c(0, 1))
ggsave("Encapsulation.svg", plot = encap.plot, width = 6, height = 4)

# Fig. S14 Parasitoid-Control contrasts ####
S14_fly.wt.par <- ggplot(data = Contrast.Fly.Wt.par, mapping = aes(x = Parasitoid, y = estimate, colour = Host)) + 
   geom_point(position = position_dodge(width = 0.9))+
   geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), position = position_dodge(width = 0.9)) +
   geom_hline(yintercept = 0, linetype = 3) +
   facet_wrap(~Temp, axes = "all_y") + 
   labs(y = "body mass difference (mg)") +
   theme(strip.background = element_blank(),
         strip.text = element_text(size = 14),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 16))
ggsave("S14_fly_wt_par.svg", S14_fly.wt.par, width = 12, height = 4.5)
