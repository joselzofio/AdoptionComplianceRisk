# Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty

Version 3.0.0

This set of three MATLAB scripts calculates the investment thresholds, along with optimal emission and violation levels from the simulations reported in Arguedas, C., Peinado, F., and Zofío, J.L. (2023) "Incentives for Green Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty", Dptm. of Economics, Universidad Autónoma de Madrid.

The first script (1.Technology_Adoption_Compliance_Risk_without_Tech_Uncertainty.m) implements the simulations reported in section 3.3 of the article without technological uncertainty about the new technology.

The second script (2.Technology_Adoption_Compliance_Risk_with_Partial_Tech_Uncertainty.m) implements the simulations reported in section 4.1 of the article enhancing the baseline model with technological uncertainty about the efficiency of the new technology, which may have either high or low abatement costs. The simulations correspond to what is termed in the article as 'partial' technological uncertainty, considering that firms know the efficiency of the technology after investing in it.   

The third script (3.Technology_Adoption_Compliance_Risk_with_Full_Tech_Uncertainty.m) implements the simulations reported in section 4.2 of the article considering full technological uncertainty 
about the efficiency of the new technology. As in the scenario with partial technological uncertainty (script no. 2), the least efficient technology presents high abatement costs (h), while the most efficient technology presents low abatement costs (l). But now investment and abatement decisions are taken without knowing the efficiency of the technology even after its adoption.

Both scripts have been run in version R2023b of MATLAB and use the function 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions. For the case of technological certainty and partial technological uncertainty see equation (2) in the article, minimizing firms' expected disutility in terms of compliance and non-compliance (equation (1)). For full uncertainty it solves equations (16) and (18) depending on whether the firm complies with the regulation or not; i.e., under perfect and imperfect compliance.  

## Usage

To run the scripts, add the m-file to the MATLAB path.

The simulations can be performed with the different parametrizations of the inspection probabilities and fines for non-compliance presented in Table 1 of the article. The parameters corresponding to the tax value (tau), degree of risk aversion (rho) and uncertainty about the efficiency of the new technology (alpha) can be also adjusted. The parametrizations for the three alternative simulations considered in the paper (Table 1) are
presented in their corresponding sections of the script. 

The scripts guide the user step-by-step on the solution of the optimization problems reported in the article. 

## Authors

Carmen Arguedas <br>
Fernando Peinado <br>
José L. Zofío


