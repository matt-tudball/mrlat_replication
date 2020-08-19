rm(list=ls())
library(foreign); library(TwoSampleMR); library(MendelianRandomization); library(ggplot2)
setwd("your working directory")

theme_set(
  theme_classic() +
    theme(legend.position = "right")
)

# MR-Base ID list:
# Coronary artery disease: 7
# Type 2 diabetes: 24
# Breast cancer: 1126

latent_beta <- function(snp_exposure,id,theta) {
  snp_outcome <- extract_outcome_data(snps=snp_exposure$SNP, outcomes=paste("ieu-a-",id,sep=""))
  
  ts <- merge(snp_exposure, snp_outcome, by="SNP")
  
  # Generate IVW MR estimates
  IVWoutput <- mr_ivw(mr_input(bx = ts$beta.exposure, bxse = ts$se.exposure, 
                               by = ts$beta.outcome, byse = ts$se.outcome))
  beta_ivw <- IVWoutput$Estimate
  se_ivw <- IVWoutput$StdError 
  
  # Compute genetic variance
  var_z <- 2*snp_exposure$eaf.exposure*(1-snp_exposure$eaf.exposure)
  
  # Compute estimate of genetic variance
  sigma_G <- sqrt(sum((snp_exposure$beta.exposure^2)*var_z))
  
  # Compute estimator
  beta_lat <- sigma_G*beta_ivw/sqrt(theta)
  
  # Compute variance
  se_lat <- sigma_G*se_ivw/sqrt(theta)
  return(list(beta_lat=beta_lat,se_lat=se_lat))
}

input <- expand.grid(c(7,24,1126),c(0.01,0.02,0.05))
ninput <- NROW(input)

output <- data.frame(matrix(nrow=ninput,ncol=3))
output <- cbind(rep(c(0.5,1,1.5),3),input,output)
colnames(output) <- c("pos","id","theta","beta","ci_lower","ci_upper")

for (i in 1:ninput) {
  id <- input[i,1]; theta <- input[i,2]
  
  if (id %in% c(7,24)) { snp_exposure <- read_exposure_data("richardson_twosample.csv", sep=",") }
  if (id %in% c(1126)) { snp_exposure <- read_exposure_data("richardson_twosample_females.csv", sep=",") }
  if (theta == 0) { next }
  
  estimates <- latent_beta(snp_exposure,id,theta)
  
  output[i,"beta"] <- exp(estimates$beta_lat)
  output[i,"ci_lower"] <- exp(estimates$beta_lat - 1.96*estimates$se_lat)
  output[i,"ci_upper"] <- exp(estimates$beta_lat + 1.96*estimates$se_lat)
}

plot <- ggplot(output, aes(pos, beta)) +
  geom_errorbar(
    aes(ymin = ci_lower, 
        ymax = ci_upper, 
        color = as.factor(theta)),
    position = position_dodge(0.3), width = 0.1
  ) +
  xlab("") + ylab("Effect estimate (OR scale)") +
  geom_point(aes(color = as.factor(theta)), position = position_dodge(0.3)) +
  geom_hline(yintercept=1, linetype="dashed", color = "#575863") +
  scale_x_continuous(breaks=c(0.5,1,1.5),labels=c("Coronary artery disease","Type 2 diabetes","Breast cancer")) +
  scale_y_continuous(breaks = seq(0.6,3.1,0.2)) +
  labs(colour="") + 
  scale_color_manual(name="Estimates",
                     labels=c(expression(paste(theta^2," = 0.01",sep="")),
                              expression(paste(theta^2," = 0.02",sep="")),
                              expression(paste(theta^2," = 0.05",sep=""))),
                     values=c("#F8766D","#00BFC4", "#7CAE00")) +
  coord_flip()

plot
ggsave(filename=paste("fig_richardson",type=".pdf",sep=""),plot=plot)