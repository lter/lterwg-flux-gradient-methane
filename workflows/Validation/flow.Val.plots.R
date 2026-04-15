# Validation Results and Figures
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(paletteer)

load(file = fs::path(localdir,paste0("Val_PARMS_DIEL_Q10_RSHPc_wide.Rdata")))

plot.diel <- ENSEMBLE_DIELSc  %>%  ggplot()+ 
  geom_line(aes(x = Hour , y = DIEL, col=season, linetype = RSHP), size=0.5) +
  scale_colour_discrete_qualitative(palette = "Harmonic") +
  theme_bw() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "transparent", linewidth = 0.5),
        legend.title = element_blank()) + 
  ylab(expression(paste( "CH"[4], " (g m"^-2, ")"))) + 
  facet_wrap(~site+ season, scales = "free_y")


plot.totaldaily <- ENSEMBLE_DIELSc_wide  %>% 
  ggplot( aes( x= FG_DIEL*1000, y=EC_DIEL*1000, col=site)) +
  geom_point(alpha=0.5, size =2) + 
  geom_smooth(method = "lm", se =, size=0.5, formula = 'y ~ x + 0' ) +
  stat_regline_equation( aes(label = ..eq.label..), 
                         size =3.5,
                         formula = 'y ~ x + 0', fontface=2,
                         label.x.npc = "left",
                         label.y.npc = 1) +
  stat_regline_equation(aes(label = ..rr.label..), size =3.5, 
                        fontface=2,
                        #formula = 'y ~ x + 0',
                        label.x.npc = "center",
                        label.y.npc = 0.35) +
  geom_abline( slope=1, col="grey40", linetype="dashed", size =1) +
  theme_bw() +
  scale_color_discrete_qualitative( palette = "Dark2") + ylab(expression(paste( " EC CH"[4], " ( mg m"^-2, " day" ^-1, ")")))  +
  xlab(expression(paste( " GF CH"[4], " (mg m"^-2, " day" ^-1, ")"))) +
  facet_wrap( ~site, ncol=1) 





plot.base <- ENSEMBLE_Q10_eq4c_wide %>% 
  ggplot( aes( x= Rref.mean_FG, y=Rref.mean_EC, col=SITE_ID)) +
  geom_point(alpha=0.5, size =2) + 
  geom_smooth(method = "lm", se =, size=0.5, formula = 'y ~ x + 0' ) +
  stat_regline_equation( aes(label = ..eq.label..), 
                        size =3.5,
                        formula = 'y ~ x + 0', fontface=2,
                        label.x.npc = "left",
                        label.y.npc = "top") +
  stat_regline_equation(aes(label = ..rr.label..), size =3.5, 
                        fontface=2,
                       # formula = 'y ~ x + 0',
                        label.x.npc =0.6,
                        label.y.npc= 0.5) +
  geom_abline( slope=1, col="grey40", linetype="dashed", size =1) +
  theme_bw() +  
  scale_color_discrete_qualitative( palette = "Dark2") + ylab('Base Respiration (EC)') +
  xlab('Base Respiration (GF)')+ facet_wrap( ~SITE_ID, ncol=1) +
  xlim(0,2) + ylim(0,2)

plot.Q10 <- ENSEMBLE_Q10_eq4c_wide %>% 
  ggplot( aes( x= Q10.mean_FG, y= Q10.mean_EC, col=SITE_ID)) +
  geom_point(alpha=0.5, size =2) + 
  geom_smooth(method = "lm", se =, size=0.5, formula = 'y ~ x + 0' ) +
  stat_regline_equation( aes(label = ..eq.label..), 
                         size =3.5,
                         formula = 'y ~ x + 0', fontface=2,
                         label.x.npc = "left",
                         label.y =3)+
  stat_regline_equation(aes(label = ..rr.label..), size =3.5, 
                        fontface=2,
                        #formula = 'y ~ x + 0',
                        label.x.npc = "center",
                        label.y =1.1) +
  geom_abline( slope=1, col="grey40", linetype="dashed", size =1) +
  theme_bw()  +  
  scale_color_discrete_qualitative( palette = "Dark2") + 
  ylab('Q10 (EC)') +
  xlab('Q10(GF)') + facet_wrap( ~SITE_ID, ncol=1) +
  xlim(1,3) + ylim(1,3)


plot.daily.boxplot <- ENSEMBLE_DIELSc_wide  %>% ggplot( ) +
  # geom_boxplot( aes( x= site, y = EC_DIEL), fill='transparent', col="grey" ) + 
  geom_boxplot(aes( x= site, y = FG_DIEL*1000), fill='transparent') +
  geom_jitter(aes( x= site, y = FG_DIEL*1000, col= Season), alpha = 0.5, size = 4) +
  # geom_jitter( aes( x= site, y = EC_DIEL, col= season), alpha = 0.25, size = 4) + 
  theme_bw() +
  scale_color_paletteer_d("ggthemes::Seattle_Grays")  +
  ylab(expression(paste( " GF CH"[4], " (mg m"^-2, " day" ^-1, ")"))) +
  labs(color = "") + xlab("")


p1 <- ggarrange( plot.daily.boxplot, labels="A")

p2 <- ggarrange( plot.totaldaily, plot.Q10, plot.base , 
           labels = c( "B", "C", "D"), common.legend =T, nrow=1)


final.plot <- ggarrange( p1 , p2, ncol=1)

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane/FIGURES/Validation_Flux_Figure_p2.png", plot = p2 , width = 6, height = 6, units = "in")

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/FluxGradient/lterwg-flux-gradient-methane/FIGURES/Validation_Flux_Figure_p1.png", plot = p1 , width = 6, height = 3, units = "in")

# Results:
ENSEMBLE_DIELSc_wide$FG_DIEL %>% summary
ENSEMBLE_DIELSc_wide$EC_DIEL %>% summary

