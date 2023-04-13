# AdoptionComplianceRisk

Version 1.0.0

This MATLAB script calculates the investment thresholds, along with optimal emission and violation levels from the simulations reported in Arguedas, C., Peinado, F. and Zofío, J.L. (2023) "Incentives for Green Technology Adoption, Imperfect Compliance, and Risk Aversion".

The script has been run in version R2022a of MATLAB and uses the function 'vpasolve' from the Symbolic Math Toolboox to determine optimal emissions, see equation (2) in the article.     

## Usage

To use this script add the m-file to the MATLAB path.

The simulations can be perfomed with the different parametrizations of the inspection probabilities and fines for non-compliance presented in Table 1 of the article. The parameters correspoding to the tax vakue (tau), degree of risk aversion (rho) and uncertainity about the efficiency of the new technbology (alpha) can be also adjusted to replicate the different results of the article.

The script guides the user step-by-tep on the solution of the optimization problems reported in the article. 

## Authors

Carmen Arguedas
Fernando Peinado 
José L. Zofío


