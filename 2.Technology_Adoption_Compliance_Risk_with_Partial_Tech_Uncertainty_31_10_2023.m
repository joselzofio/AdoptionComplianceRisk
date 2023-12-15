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

% This second script implements the simulations reported in section 4.1 
% of the article enhancing the baseline model with partial technological uncertainty 
% about the efficiency of the new technology. The least efficient technology 
% presents high abatement costs (h) while the most efficient technology presents low 
% abatement costs (l).

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
% ei0=actual emissions with old technology
% ei1h=actual emissions with the least efficient new technology
% ei1l=actual emissions with the most efficient new technology
% ei0c=optimal actual emissions with old technology
% ei1hc=optimal actual emissions with the least efficient new technology
% ei1lc=optimal actual emissions with the most efficient new technology
% ri0=reported (declared) emissions with old technology
% ri1=reported (declared) emissions with new technologies
% tau=tax on declared emissions
% rho=level of risk aversion
% ff=fixed part of the fine
% alpha=likelihood of least efficient abatement new technology

clear all;
global pi rho i tau ei0c ei1hc ei1lc alpha ff
syms ei0 ei1h ei1l ri0 ri1

% [1] Assigns values for parameters

pi=0.5
tau=20
rho=1
alpha=0.5
ff=0

% [2] Defines the abatement costs functions for the old technology (0), the new technology (1) with
% high (h) abatement costs--denoted by 1h--and the new technology (1) with low (l) abatement costs
% denoted by (1l)

ci0=(100-ei0)*ei0;
ci1h=(75-ei1h)*ei1h;
ci1l=(25-ei1l)*ei1l;

% Calculates the optimal emissions depending on the technology (0, 1h or 1l)
opei0eqn=diff(ci0,ei0)+tau==0;
opei0 = vpasolve(opei0eqn, ei0, [-Inf Inf]);
ei0c=opei0;
opei1heqn=diff(ci1h,ei1h)+tau==0;
opei1h = vpasolve(opei1heqn, ei1h, [-Inf Inf]);
ei1hc=opei1h;
opei1leqn=diff(ci1l,ei1l)+tau==0;
opei1l = vpasolve(opei1leqn, ei1l, [-Inf Inf]);
ei1lc=opei1l;

% [3] Assigns abatement costs funtions like in [2] but substituting ei0, ei1h and ei1l 
% by the optimal values ei0c, ei1hc and ei1lc, respectively. Assigns sanctioning functions (fines) with ei0c, ei1hc 
% and ei1lc. These fucntions  are the same to the case without technology uncertainty (first script) 
% to allow comparability. Notice that fi0d corresponds to the first derivative of fi0 evaluated at 
% the violation level, this is at ei0c-ri0 (which is not the same as de first derivative 
% of fi0 evaluated at ri0). fi0d is entered manually, and the same applies to fi1hd and fi1ld

ci0=(100-ei0c)*ei0c;
ci1h = (75-ei1hc)*ei1hc;
ci1l = (25-ei1lc)*ei1lc;
fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2; % Parametrizations 1 and 3: fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2; Parametrization 2: fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2;
fi0d=ff+2*ei0c-2*ri0; % Parametrizations 1 and 3: fi0d=ff+2*ei0c-2*ri0; Parametrization 2: fi0d=ff*5+2*5*ei0c-2*5*ri0;
fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2; % Parametrizations 1 and 3: fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2; Parametrization 2: fi1=ff*5*(ei1hc-ri1)+5*(ei1hc-ri1)^2;
fi1hd=ff+2*ei1hc-2*ri1; % Parametrizations 1 and 3: fi1d=ff+2*ei1hc-2*ri1; Parametrization 2: fi1d=ff*5+5*2*ei1hc-5*2*ri1;
fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2; % Parametrizations 1 and 3: fi1h=ff*(ei1lc-ri1)+(ei1lc-ri1)^2; Parametrization 2: fi1=ff*5*(ei1lc-ri1)+5*(ei1lc-ri1)^2;
fi1ld=ff+2*ei1lc-2*ri1; % Parametrizations 1 and 3: fi1d=ff+2*ei1hc-2*ri1; Parametrization 2: fi1d=ff*5+5*2*ei1hc-5*2*ri1;

% [4] Backwards loop for several i values. Do not directly specify i values
% in i, instead, set counter=maximum value of i to be studied

for counter=2105:-1:0;

i=counter;

% Old Technology

% [5] Disutility function specified as a power function. Notice that dia0d and dib0d 
% correspond to the first derivatives of the disutility functions dia0 and dia1 with respect to 
% ai0=(ci0+tau*ri0) and bi0=(ci0+tau*ri0+fi0), respectively. Again they must be entered manually at this step

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
% emissions). Be careful to keep Fo the same as in the sanctioning function in [3]
% substituting ri0 by O (capital letter 'O', the amount of declared emissions with 
% the old technology)

O = vpasolve(eqn, ri0,[0 ei0c]); % Amount of declared emissions with the old technology

Vo=ei0c-O; % Violation level with the old technology
 
Co=ci0; % Abatement costs with the old technology

To=O*tau; % Total taxes paid for declared emissions with the old technology

Fo=ff*(ei0c-O)+(ei0c-O)^2; % Fine for the violation level with the old technology, Parametrizations 1 and 3: Fo=ff*(ei0c-O)+(ei0c-O)^2; Parametrization 2: Fo=ff*5*(ei0c-O)+5*(ei0c-O)^2;

% Shows total costs to get an idea about how large they are numerically (the independent variable)

Costso=Co+To+Fo;

% [7] Disutility assuming the old technology. Make sure that the elements multiplying (1-pi) and pi are the same
% disutility functions written for [5] but substituting ci0, tau*ri and fi0
% by Co, To and Fo respectively.

Do=(1-pi)*((Co+To)^(rho+1))+pi*((Co+To+Fo)^(rho+1));

% Least efficient new technology (1) with high (h) abatement costs--denoted by (1h)

% [8] Follows the same reasoning as in [5] when setting these disutility
% functions. It keeps their definitions unchanged but adds the investment
% cost i to the total costs and replaces the numeral of the  
% old technology (0) with that of the least effienct new technology (1h)

dia1h = (ci1h+tau*ri1+i)^(rho+1);
dia1hd = (rho+1)*(ci1h+tau*ri1+i)^(rho);
dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1);
dib1hd = (rho+1)*(ci1h+tau*ri1+fi1h+i)^(rho);

% Proposition 1: Equation to implicitly obtain declared emissions
eqn = (dib1hd*pi*fi1hd)/((1-pi)*dia1hd+pi*dib1hd) == tau;

% [9] The function 'vpasolve' solves the equation numerically. With the range
% [0 ei1hc] we consider only real solutions beteween 0 (which consitutes the
% lower bound because firms will not declare negative emissions) and the
% optimal emissions level with the least efficient new tech (which constitutes the upper
% bound because firms will not declare more emissions than actual
% emissions). Be careful to keep Fnh the same as the santioning function in [3]
% but substituting ri1 by Nh (the amount of declared emissions with the least efficient new technology)

Nh = vpasolve(eqn, ri1,[0 ei1hc]); %Amount of declared emissions with the least efficient new technology

Vnh=ei1hc-Nh; % Violation level with the least efficient new technology

Cnh=ci1h; % Abatement costs with the least efficient new technology

Tnh=Nh*tau; % Total taxes paid for declared emissions with the least efficient new technology

Inh=i; % Fixed investment cost in the new technology

Fnh=ff*(ei1hc-Nh)+(ei1hc-Nh)^2; % Fine for the violation level with the least efficient new technology, Parametrizations 1 and 3: Fn=ff*(ei1hc-Nh)+(ei1hc-Nh)^2; Parametrization 2: Fn=ff*5*(ei1hc-Nh)+5*(ei1hc-Nh)^2;

% Shows total costs with the least effienct new technology to get an idea about how large they 
% are numerically (the independent variable)

Costsnh=Cnh+Tnh+Inh+Fnh;

% [10] Disutility assuming the least efficient new technology. Make sure that the elements multiplying (1-pi) 
% and pi are the same disutility functions written for [8] but substituting ci1h, tau*ri1, i and 
% fi1h by Cnh, Tnh, Inh and Fnh, respectively

Dnh=(1-pi)*((Cnh+Tnh+Inh)^(rho+1))+pi*((Cnh+Tnh+Inh+Fnh)^(rho+1));

% Most efficient new technology (1) with low (l) abatement cost-- denoted by (1l) 

% [11] Follows the same reasoning as in [5] when setting these disutility
% functions. It keeps their definitions unchanged but adds the investment
% cost i to the total costs and replaces the numeral of the  
% old technology (0) with that of the most effienct new technology (1l)


dia1l = (ci1l+tau*ri1+i)^(rho+1);
dia1ld = (rho+1)*(ci1l+tau*ri1+i)^(rho);
dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1);
dib1ld = (rho+1)*(ci1l+tau*ri1+fi1l+i)^(rho);

% Proposition 1: Equation to implicitly obtain declared emissions

eqn = (dib1ld*pi*fi1ld)/((1-pi)*dia1ld+pi*dib1ld) == tau;

% [12] The function 'vpasolve' solves the equation numerically. With the range
% [0 ei1lc] we consider only real solutions beteween 0 (which consitutes the
% lower bound because firms will not declare negative emissions) and the
% optimal emissions level with the most efficient new tech (which constitutes the upper
% bound because firms will not declare more emissions than actual
% emissions). Be careful to keep Fnl the same as in the santioning function in [3]
% but substituting ri1 by Nl (the amount of declared emissions with the most efficient new technology)

Nl = vpasolve(eqn, ri1,[0 ei1lc]); % Amount of declared emissions with the most efficient new technology

Vnl=ei1lc-Nl; % Violation level with the most efficient new technology

Cnl=ci1l; % Abatement costs with the most efficient new technology

Tnl=Nl*tau; % Total taxes paid for declared emissions with the most efficient new technology

Inl=i; % Fixed investment cost in the new technology

Fnl=ff*(ei1lc-Nl)+(ei1lc-Nl)^2; % Fine for the violation level with the most efficient new technology, Parametrizations 1 and 3: Fn=ff*(ei1lc-Nl)+(ei1lc-Nl)^2; Parametrization 2: Fn=ff*5*(ei1lc-Nl)+5*(ei1lc-Nl)^2;

% Shows total costs with the most effienct new technology to get an idea about how large they 
% are numerically (the independent variable)

Costsnl=Cnl+Tnl+Inl+Fnl;

% [13] Disutility assuming the most efficient new technology. Make sure that the elements multiplying 
% (1-pi) and pi are the same disutility functions written for [11] but substituting ci1l, tau*ri1, i and 
% fi1l by Cnl, Tnl, Inl and Fnl, respectively.

Dnl=(1-pi)*((Cnl+Tnl+Inl)^(rho+1))+pi*((Cnl+Tnl+Inl+Fnl)^(rho+1));

% Now we calculate the expected disutility with the possible new technologies

EDn=alpha*Dnh+(1-alpha)*Dnl;

% Now we check if the disutility of producing with the old technology is larger than the
% expected disutility of producing with the new technology. When this is the case, investing in
% the new technology becomes the best option. This indicates the program
% that an aproximate Ii threshold has been found and therefore the loop ends

if Do-EDn>0
    break
end
    
end

X=['Indifferent if the investiment cost of the new technology (Ii) is aproximately between [',num2str(i),',',num2str(i+1),')'];
disp(X)
fprintf('Optimal emissions with the old technology = %s\n',opei0);
eopei1=alpha*opei1h+(1-alpha)*opei1l;
fprintf('Expected optimal emissions with the new technology = %s\n',eopei1);
fprintf('Optimal emissions with the least efficient new technology = %s\n',opei1h);
fprintf('Optimal emissions with the most efficient new technology = %s\n',opei1l);
fprintf('Violation level with the old technology = %s\n',Vo);
eVn=alpha*Vnh+(1-alpha)*Vnl;
fprintf('Expected violation level with the new technology = %s\n',eVn);
fprintf('Violation level with the least efficient new technology = %s\n',Vnh);
fprintf('Violation level with the most efficient new technology = %s\n',Vnl);

