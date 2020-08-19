# Replication kit for "A simple framework for Mendelian randomisation of latent continuous exposures with discrete measurements"

This is the replication kit to reproduce the results in this paper.
Contact me at matt.tudball@bristol.ac.uk if you have any questions.

## Description of Files

### Data
* `pasman_raw.csv`: Summary statistics available from the supplementary data in Pasman et al (2018).

### Data preparation
* `genetic_all.do`: Formatting genetic instruments for early life adiposity in males and females (Richardson et al, 2020).
* `genetic_female.do`: Formatting genetic instruments for early life adiposity in females (Richardson et al, 2020).
* `bmi_sys_example.do`: Generating estimates from various cut-offs in the BMI and systolic blood pressure example.

### Figures and tables
* `bmi_sys_example.R`: Generates Figure 3.
* `example_pasman.R`: Generates Figure 4.
* `example_richardson.R`: Generates Figure 5.
* `assumptions_sim.R`: Generates Tables 1 and 2.
