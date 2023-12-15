% Date: 31th of October, 2023

% This set of three MATLAB scripts calculates the investment thresholds, along with 
% optimal emission and violation levels from the simulations reported in 
% Arguedas, C., Peinado, F., and Zofío, J.L. (2023) "Incentives for Green 
% Technology Adoption and Compliance under Risk Aversion and Technological Uncertainty",
% Dptm. of Economics, Universidad Autónoma de Madrid.
% Optimal emissions are calculated by solving the condition presented in 
% equation (2) of the article, minimizing firms' expected disutility in terms of 
% compliance and non-compliance (see equation (1)). 
% The scripts has been run in version R2023b of MATLAB and uses the function
% 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions.

% This first script implements the simulations reported in section 3.3 
% of the article without technological uncertainty about the new
% technology.

% The simulations can be performed with the different parametrizations of 
% the inspection probabilities and fines for non-compliance presented in Table 1 
% of the article. The parameters corresponding to the tax value (tau) and degree 
% of risk aversion (rho) can be also adjusted. The parametrizations for the 
% three alternative simulations considered in the paper (Table 1) are
% presented in their corresponding sections of the script.

% Follow steps [1], [2], ... to ensure correct outputs.

% Notation:
% pi=monitoring probability
% i=fixed investment cost of intalling the cleaner/new technology
% ei0=actual emissions with the old technology
% ei1=actual emissions with the new technology
% ei0c=optimal actual emissions with the old technology
% ei1c=optimal actual emissions with the new technology
% ri0=reported (declared) emissions with old technology
% ri1=reported (declared) emissions with new technology
% tau=tax on declared emissions
% rho=degree of risk aversion
% ff=fixed part of the fine

clear all;
global pi rho i tau ei0c ei1c ff
syms ei0 ei1 ri0 ri1

% [1] Assigns values for parameters

pi=0.5
tau=20
rho=1
ff=0

% [2] Defines the abatement costs functions for the old technology (0) and the new technology (1)

ci0=(100-ei0)*ei0;
ci1=(50-ei1)*ei1;

% Calculates the optimal emissions depending on the technology (0 or 1)

opei0eqn=diff(ci0,ei0)+tau==0;
opei0 = vpasolve(opei0eqn, ei0, [-Inf Inf]);
ei0c=opei0;
opei1eqn=diff(ci1,ei1)+tau==0;
opei1 = vpasolve(opei1eqn, ei1, [-Inf Inf]);
ei1c=opei1;

% [3] Assigns abatement costs funtions like in [2] but substituting ei0 and ei1 
% by the optimal values ei0c and ei1c, respectively. Assigns sanctioning functions (fines) with ei0c and ei1c. 
% These fucntions are the same to the case with technological uncertainty (second script) 
% to allow comparability. Notice that fi0d corresponds to the first derivative of fi0 evaluated at 
% the violation level; that is, at ei0c-ri0 (which is not the same as de first derivative 
% of fi0 evaluated at ri0). fi0d is entered manually, and the same applies to fi1d

ci0=(100-ei0c)*ei0c;
ci1 = (50-ei1c)*ei1c;
fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2; % Parametrizations 1 and 3: fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2; Parametrization 2: fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2;
fi0d=ff+2*ei0c-2*ri0; % Parametrizations 1 and 3: fi0d=ff+2*ei0c-2*ri0; Parametrization 2: fi0d=ff*5+2*5*ei0c-2*5*ri0;
fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2; % Parametrizations 1 and 3: fi1h=ff*(ei1c-ri1)+(ei1c-ri1)^2; Parametrization 2: fi1=ff*5*(ei1c-ri1)+5*(ei1c-ri1)^2;
fi1d=ff+2*ei1c-2*ri1; % Parametrizations 1 and 3: fi1d=ff+2*ei1c-2*ri1; Parametrization 2: fi1d=ff*5+5*2*ei1c-5*2*ri1;

% [4] Backwards loop for several i values. Do not directly specify i values
% in i, instead, set counter=maximum value of i to be studied

for counter=2375:-1:0;

i=counter;

% Old Technology

% [5] Disutility function specified as a power function. Notice that dia0d and dib0d correspond to the first
% derivatives with respect to ai0=(ci0+tau*ri0) and bi0=(ci0+tau*ri0+fi0)
% of dia0 and dib0 respectively.? Again, they are entered manually at this step

dia0 = (ci0+tau*ri0)^(rho+1);
dia0d = (rho+1)*(ci0+tau*ri0)^(rho);
dib0 = (ci0+tau*ri0+fi0)^(rho+1);
dib0d = (rho+1)*(ci0+tau*ri0+fi0)^(rho);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dib0d*pi*fi0d)/((1-pi)*dia0d+pi*dib0d) == tau;

% [6] The function 'vpasolve' solves the equation numerically. With the range
% [0 ei0c] we consider only real solutions beteween 0 (which consitutes the
% lower bound because firms will not declare negative emissions) and the
% optimal emissions level with the old tech (which constitutes the upper
% bound because firms will not declare more emissions than actual
% emissions). Be careful to keep Fo the same as in the santioning function in [3]
% but substituting ri0 by O (capital letter 'O', the amount of declared emissions 
% with the old technology)

O = vpasolve(eqn, ri0,[0 ei0c]); % Amount of declared emissions with the old technology

Vo=ei0c-O; % Violation level with the old technology

Co=ci0; % Abatement costs with the old technology

To=O*tau; %Total taxes paid for declared emissions with the old technology

Fo=ff*(ei0c-O)+(ei0c-O)^2; % Fine for the violation level with the old technology, Parametrizations 1 and 3: Fo=ff*(ei0c-O)+(ei0c-O)^2; Parametrization 2: Fo=ff*5*(ei0c-O)+5*(ei0c-O)^2;

% Shows total costs to get an idea about how large they are numerically (being the independent variable)

Costso=Co+To+Fo;

% [7] Disutility assuming the old technology. Make sure that the elements multiplying (1-pi) 
% and pi are the same disutility functions written for [5] but substituting ci0, tau*ri and fi0
% by Co, To and Fo respectively.

Do=(1-pi)*((Co+To)^(rho+1))+pi*((Co+To+Fo)^(rho+1));

% New Technology

% [8] Follows the same reasoning as in [5] when setting these disutility
% functions. It keeps their definitions unchanged but adds the investment
% cost i to the total costs and replaces the numeral of the  
% old technology (0) with that of the new technology (1)

dia1 = (ci1+tau*ri1+i)^(rho+1);
dia1d = (rho+1)*(ci1+tau*ri1+i)^(rho);
dib1 = (ci1+tau*ri1+fi1+i)^(rho+1);
dib1d = (rho+1)*(ci1+tau*ri1+fi1+i)^(rho);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dib1d*pi*fi1d)/((1-pi)*dia1d+pi*dib1d) == tau;

% [9] The function 'vpasolve' solves the equation numerically. With the range
% [0 ei1c] we consider only real solutions beteween 0 (which consitutes the
% lower bound because firms will not declare negative emissions) and the
% optimal emissions level with the new tech (which constitutes the upper
% bound because firms will not declare more emissions than actual
% emissions). Be careful to keep Fn the same as in the santioning function in [3]
% but substituting ri1 by N (the amount of declared emissions with the new technology)

N = vpasolve(eqn, ri1,[0 ei1c]); % Amount of declared emissions with the least efficient new technology

Vn=ei1c-N; % Violation level with the least efficient new technology

Cn=ci1; % Abatement costs with the least efficient new technology

Tn=N*tau; % Total taxes paid for declared emissions with the least efficient new technology

In=i; % Fixed investment cost in the new technology

Fn=ff*(ei1c-N)+(ei1c-N)^2; % Fine for the violation level with the least efficient new technology, Parametrizations 1 and 3: Fn=ff*(ei1c-N)+(ei1c-N)^2; Parametrization 2: Fn=ff*5*(ei1c-N)+5*(ei1c-N)^2;

% Shows total costs with the new technology to get an idea about how large they are numerically (the independent variable)

Costsn=Cn+Tn+In+Fn;

% [10] Defines disutility assuming the new technology. Make sure that the elements multiplying 
% (1-pi) and pi are the same disutility functions written for [8] but substituting ci1, tau*ri1, i and 
% fi1 by Cn, Tn, In and Fn, respectively

Dn=(1-pi)*((Cn+Tn+In)^(rho+1))+pi*((Cn+Tn+In+Fn)^(rho+1));

% Now we check if the disutility of producing with the old technology is larger than the
% disutility of producing with the new technology. When this is the case, investing in
% the new technology becomes the best option. This indicates the program
% that an approximate Ii threshold has been found and therefore the loop ends

if Do-Dn>0
    break
end
    
end

X=['Indifferent if the investment cost of the new technology (Ii) is aproximately between [',num2str(i),',',num2str(i+1),')'];
disp(X)
fprintf('Optimal emissions with the old technology = %s\n',opei0);
fprintf('Optimal emissions with the new technology = %s\n',opei1);
fprintf('Violation level with the old technology = %s\n',Vo);
fprintf('Violation level with the new technology = %s\n',Vn);

