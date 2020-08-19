rm(list=ls())
library(foreign); library(TwoSampleMR); library(MendelianRandomization); library(ggplot2)
setwd("your working directory")

theme_set(
  theme_classic() +
    theme(legend.position = "right",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
)

# Import two-sample data
snp_outcome <- read.csv("pasman_raw.csv")

# Generate IVW MR estimates
IVWoutput <- mr_ivw(mr_input(bx = snp_outcome$ze_beta, bxse = snp_outcome$ze_se,
                             by = snp_outcome$zy_beta, byse = snp_outcome$zy_se))
beta_ivw <- IVWoutput$Estimate
se_ivw <- IVWoutput$StdError 

# Theta parameter (genetic variance)
theta <- 0.034

# Compute genetic variance
var_z <- 2*snp_outcome$freq_allele_1*(1-snp_outcome$freq_allele_1)

# Compute estimate of genetic variance
sigma_G <- sqrt(sum((snp_outcome$ze_beta^2)*var_z))

# Compute estimator
beta_mrlat <- sigma_G*beta_ivw/sqrt(theta)

# Compute variance
se_mrlat <- sigma_G*se_ivw/sqrt(theta)

# Create figure
input <- c(0.02,0.034,0.05)
ninput <- NROW(input)

output <- data.frame(matrix(nrow=ninput,ncol=3))
output <- cbind(rep(1,ninput),input,output)
colnames(output) <- c("pos","theta","beta","ci_lower","ci_upper")

for (i in 1:ninput) {
  theta <- output$theta[i]
  if (theta == 0) { next }
  
  beta_mrlat <- sigma_G*beta_ivw/sqrt(theta)
  se_mrlat <- sigma_G*se_ivw/sqrt(theta)
  
  output[i,"beta"] <- exp(beta_mrlat)
  output[i,"ci_lower"] <- exp(beta_mrlat - 1.96*se_mrlat)
  output[i,"ci_upper"] <- exp(beta_mrlat + 1.96*se_mrlat)
}

plot <- ggplot(output, aes(pos, beta)) +
  geom_errorbar(
    aes(ymin = ci_lower, 
        ymax = ci_upper, 
        color = as.factor(theta)),
    position = position_dodge(0.3), width = 0.05
  ) +
  xlab("") + ylab("Effect estimate (OR scale)") +
  geom_point(aes(color = as.factor(theta)), position = position_dodge(0.3)) +
  geom_hline(yintercept=1, linetype="dashed", color = "#575863") +
  scale_y_continuous(breaks = seq(0.8,1.6,0.1)) +
  labs(colour="") + 
  scale_color_manual(name="Estimates",
                     labels=c(expression(paste(theta^2," = 0.02",sep="")),
                              expression(paste(theta^2," = 0.034",sep="")),
                              expression(paste(theta^2," = 0.05",sep=""))),
                     values=c("#F8766D","#00BFC4", "#7CAE00")) +
  coord_flip()

plot
ggsave(filename=paste("fig_pasman",type=".pdf",sep=""),plot=plot)