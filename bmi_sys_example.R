rm(list=ls())
library(foreign); library(ggplot2); library(tidyr)
setwd("your working directory")

theme_set(
  theme_classic() +
  theme(legend.position = "right")
)

beta_true <- 1.5270575

# Load in the coefficient data
DATA <- read.dta("vary_cutoff_long.dta")
DATA <- DATA[!(DATA$group==1),]

plot <- ggplot(DATA, aes(cutoff, beta)) +
        geom_errorbar(
        aes(ymin = beta - 1.96*std, 
            ymax = beta + 1.96*std, 
            color = as.factor(group)),
        position = position_dodge(0.3), width = 0.3
    ) +
    xlab("BMI threshold") + ylab("Comparison of estimated effects with true effect") +
    geom_point(aes(color = as.factor(group)), position = position_dodge(0.3)) +
    geom_hline(aes(yintercept=beta_true, color = "#575863"), linetype="dashed") +
    scale_x_continuous(breaks = seq(22.5, 35, by = 2.5)) +
    scale_y_continuous(breaks = c(0,round(beta_true,2),seq(4,18,4))) +
    labs(colour="") + 
    scale_color_manual(name="Estimators",
                     labels=c("True","Latent","Naive"),
                     values=c("#575863","#F8766D","#00BFC4"))

plot
ggsave(filename=paste("real_data_sim",type=".pdf",sep=""),plot=plot)