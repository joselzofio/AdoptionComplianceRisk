# Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty

Version 2.0.0

These MATLAB scripts calculate the investment thresholds, along with optimal emission and violation levels from the simulations reported in Arguedas, C., Peinado, F., and Zofío, J.L. (2023) "Incentives for Green Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty". Universidad Autónoma de Madrid. Manuscript submitted for publication.

The first script (Technology_Adoption_Compliance_Risk_without_Tech_Uncertainty.m) implements the simulations reported in section 3.3 of the article without technological uncertainty about the new technology.

The second script (Technology_Adoption_Compliance_Risk_with_Tech_Uncertainty.m) implements the simulations reported in section 4 of the article enhancing the baseline model with technological uncertainty 
about the efficiency of the new technology, which may have either high or low abatement costs. The simulations correspond to what is termed in the article as 'partial' and 'full' technological uncertainty. 

Both scripts have been run in version R2023a of MATLAB and use the function 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions. For the case of technological certainty and partial technological uncertainty see equation (2) in the article, minimizing firms' expected disutility in terms of compliance and non-compliance (equation (1)). For full uncertainty it solves equations (17) and (18).  

## Usage

To run the scripts, add the m-file to the MATLAB path.

The simulations can be performed with the different parametrizations of the inspection probabilities and fines for non-compliance presented in Table 1 of the article. The parameters corresponding to the tax value (tau), degree of risk aversion (rho) and uncertainty about the efficiency of the new technology (alpha) can be also adjusted.  At the end of the script, we present the parameters' values necessary to replicate all the simulations reported in the Figures of the article.

The scripts guide the user step-by-step on the solution of the optimization problems reported in the article. 

## Authors

Carmen Arguedas <br>
Fernando Peinado <br>
José L. Zofío


