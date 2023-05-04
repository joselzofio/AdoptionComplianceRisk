% This MATLAB script calculates the investment thresholds, along with 
% optimal emission and violation levels from the simulations reported in 
% Arguedas, C., Peinado, F., and Zofio, J.L. (2023) "Incentives for Green 
% Technology Adoption, Imperfect Compliance, and Risk Aversion". 
% Optimal emissions are calculated by solving the condition presented in 
% equation (2) of the article, minimizing firms' expected disutility in terms of 
% compliance and non-compliance (see equation (1)). 
% The script has been run in version R2022a of MATLAB and uses the function
% 'vpasolve' from the Symbolic Math Toolbox to determine optimal
% emissions.

%The second script implements the simulations reported in section 4 
% of the article enhancing the baseline model with technological uncertainty 
% about the efficiency of the new technology. The least efficient technology 
% presents high abatement costs (h) while the most efficient technology presents low 
% abatement costs (l)

% Follow steps [1], [2], ... to ensure correct outputs.

% The simulations can be performed with the different parametrizations of 
% the inspection probabilities and fines for non-compliance presented in Table 1 
% of the article. The parameters corresponding to the tax value (tau), degree 
% of risk aversion (rho) and uncertainty about the efficiency of the new 
% technology (alpha) can be also adjusted. At the end of the script we present 
% the parameters' values necessary to replicate all the simulations reported 
% in the Figures of the article.

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
fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2;
fi0d=ff+2*ei0c-2*ri0;
fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2;
fi1hd=ff+2*ei1hc-2*ri1;
fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2;
fi1ld=ff+2*ei1lc-2*ri1;

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

Fo=ff*(ei0c-O)+(ei0c-O)^2; % Fine for the violation level with the old technology

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

Fnh=ff*(ei1hc-Nh)+(ei1hc-Nh)^2; % Fince for the violation level with the least efficient new technology

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

Fnl=ff*(ei1lc-Nl)+(ei1lc-Nl)^2; % Fine for the violation level with the most efficient new technology

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

% Data to replicate the Figures in the article:

% Figures 4, 6 and 7: Parametrization 1 only

% Alpha = 0

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 19.301772360941364490984426550317
% Vnl = 18.99245315813295213351210770804
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0
% ff=0
% Vnh = 18.726365051122342057586281780631
% Vnl = 18.224413258401657712256236912388
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0
% ff=0
% Vnh = 18.238198867726156737637539671894
% Vnl = 17.606764301703377699993234865649
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0
% ff=0
% Vnh = 17.815232599171785605063150328151
% Vnl = 17.092374312977055638457145806162
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0
% ff=0
% Vnh = 17.442841583754874384821014096979
% Vnl = 16.653229578803857460064442621404
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0
% ff=0
% Vnh = 17.110816158051062354885644555716
% Vnl = 16.271309186571179832008658998362
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0
% ff=0
% Vnh = 16.811743381243509168948245570282
% Vnl = 15.934334528935537899263003044392
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0
% ff=0
% Vnh = 16.540069775096241559687201618301
% Vnl = 15.633564799571109717441496672587
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0
% ff=0
% Vnh = 16.29152712615606868547056009844
% Vnl = 15.362558723709445960239793454227
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0
% ff=0
% Vnh = 16.062764395340996214003301452459
% Vnl = 15.116434690950190056731574439201
% Indifferent if the investment cost in the new technology is Ii=[3093,3094)

%-------

% Alpha = 0.25

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0.25
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 19.231328099509869514251662843963
% Vnl = 18.838136095015301653263893634592
% Indifferent if the investment cost in the new technology is Ii=[2570,2571)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0.25
% ff=0
% Vnh = 18.584804006866309543766869068096
% Vnl = 17.932776731674537945807077458526
% Indifferent if the investment cost in the new technology is Ii=[2472,2473)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0.25
% ff=0
% Vnh = 18.022405639992635553499389095446
% Vnl = 17.180501324395489419999467909022
% Indifferent if the investment cost in the new technology is Ii=[2366,2367)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0.25
% ff=0
% Vnh = 17.522280055162601913590753498986
% Vnl = 16.531032386159403066356037539271
% Indifferent if the investment cost in the new technology is Ii=[2256,2257)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0.25
% ff=0
% Vnh = 17.074020741942462380876569009478
% Vnl = 15.962885752802952334074762879226
% Indifferent if the investment cost in the new technology is Ii=[2152,2153)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0.25
% ff=0
% Vnh = 16.670497425950675428790422267695
% Vnl = 15.462642451290225430628689559743
% Indifferent if the investment cost in the new technology is Ii=[2058,2059)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0.25
% ff=0
% Vnh = 16.307078683583530951049671620744
% Vnl = 15.022358528867827162827976064618
% Indifferent if the investment cost in the new technology is Ii=[1977,1978)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0.25
% ff=0
% Vnh = 15.979304734498648688600920985499
% Vnl = 14.634576592865709037616310710018
% Indifferent if the investment cost in the new technology is Ii=[1909,1910)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0.25
% ff=0
% Vnh = 15.682491692799259675605846913483
% Vnl = 14.291469354802551957195937284737
% Indifferent if the investment cost in the new technology is Ii=[1852,1853)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0.25
% ff=0
% Vnh = 15.413068243172342998254077643078
% Vnl = 13.987336060537199326844380923794
% Indifferent if the investment cost in the new technology is Ii=[1805,1806)

% For 
% pi=0.5
% tau=20
% rho=11
% alpha=0.25
% ff=0
% Vnh = 15.166702341140724826759154521836
% Vnl = 13.715013443843147475849960618832
% Indifferent if the investment cost in the new technology is Ii=[1765,1766)

% For 
% pi=0.5
% tau=20
% rho=12
% alpha=0.25
% ff=0
% Vnh = 14.940581426245840790832356233845
% Vnl = 13.470164246675583855441180102428
% Indifferent if the investment cost in the new technology is Ii=[1731,1732)

% For 
% pi=0.5
% tau=20
% rho=13
% alpha=0.25
% ff=0
% Vnh = 14.732844978899559390832444007987
% Vnl = 13.250130154789668336748729518081
% Indifferent if the investment cost in the new technology is Ii=[1703,1704)

% For 
% pi=0.5
% tau=20
% rho=14
% alpha=0.25
% ff=0
% Vnh = 14.540003336174248710742674307056
% Vnl = 13.049223417686407143730063070562
% Indifferent if the investment cost in the new technology is Ii=[1678,1679)

% For 
% pi=0.5
% tau=20
% rho=15
% alpha=0.25
% ff=0
% Vnh = 14.360663871089782923863810948082
% Vnl = 12.865554650576847874274198963045
% Indifferent if the investment cost in the new technology is Ii=[1656,1657)

% For 
% pi=0.5
% tau=20
% rho=16
% alpha=0.25
% ff=0
% Vnh = 14.19363944145994320639150666981
% Vnl = 12.697528232570800704893363805991
% Indifferent if the investment cost in the new technology is Ii=[1637,1638)

% For 
% pi=0.5
% tau=20
% rho=17
% alpha=0.25
% ff=0
% Vnh = 14.037258770251788246427955959929
% Vnl = 12.542640638667378713812404897821
% Indifferent if the investment cost in the new technology is Ii=[1620,1621)

% For 
% pi=0.5
% tau=20
% rho=18
% alpha=0.25
% ff=0
% Vnh = 13.89062983397960058389175101618
% Vnl = 12.399732507992629704521902601038
% Indifferent if the investment cost in the new technology is Ii=[1605,1606)

% For 
% pi=0.5
% tau=20
% rho=19
% alpha=0.25
% ff=0
% Vnh = 13.752975961463175117879336164825
% Vnl = 12.267798759469801856500507716845
% Indifferent if the investment cost in the new technology is Ii=[1592,1593)

% For 
% pi=0.5
% tau=20
% rho=20
% alpha=0.25
% ff=0
% Vnh = 13.62296527703654245834577287994
% Vnl = 12.144852607301888145723155765472
% Indifferent if the investment cost in the new technology is Ii=[1580,1581)

%--------

% Alpha = 0.5

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0.5
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 19.15533827507332185926176043157
% Vnl = 18.653960336192927052901029425648
% Indifferent if the investment cost in the new technology is Ii=[2104,2105)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0.5
% ff=0
% Vnh = 18.452048466601192312778124874051
% Vnl = 17.631047556528398439982461081871
% Indifferent if the investment cost in the new technology is Ii=[1994,1995)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0.5
% ff=0
% Vnh = 17.851474551066332415196634534046
% Vnl = 16.80801930241896094473433897082
% Indifferent if the investment cost in the new technology is Ii=[1896,1897)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0.5
% ff=0
% Vnh = 17.331984738750575114487745118846
% Vnl = 16.128510972986086014088743529409
% Indifferent if the investment cost in the new technology is Ii=[1814,1815)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0.5
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% Vnh = 16.878301689624960163027400990489
% Vnl = 15.558219621483619790563864808335
% Indifferent if the investment cost in the new technology is Ii=[1747,1748)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0.5
% ff=0
% Vnh = 16.47940422889132184248792973263
% Vnl = 15.074874219929352759041405779093
% Indifferent if the investment cost in the new technology is Ii=[1694,1695)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0.5
% ff=0
% Vnh = 16.125173517683213066461404648824
% Vnl = 14.659181512159272075106540842903
% Indifferent if the investment cost in the new technology is Ii=[1651,1652)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0.5
% ff=0
% Vnh = 15.809419917529482587152319395947
% Vnl = 14.300413215033858692989406980426
% Indifferent if the investment cost in the new technology is Ii=[1618,1619)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0.5
% ff=0
% Vnh = 15.524206522163911127357867119026
% Vnl = 13.98436191124448927059558411127
% Indifferent if the investment cost in the new technology is Ii=[1590,1591)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0.5
% ff=0
% Vnh = 15.265520108025544335080504485313
% Vnl = 13.704763695040111463242070989103
% Indifferent if the investment cost in the new technology is Ii=[1567,1568)

% For 
% pi=0.5
% tau=20
% rho=11
% alpha=0.5
% ff=0
% Vnh = 15.029509563743148646663063589039
% Vnl = 13.455513872522624089467339297735
% Indifferent if the investment cost in the new technology is Ii=[1548,1549)

% For 
% pi=0.5
% tau=20
% rho=12
% alpha=0.5
% ff=0
% Vnh = 14.81289411794680131575979660488
% Vnl = 13.231513104318970124937153311802
% Indifferent if the investment cost in the new technology is Ii=[1532,1533)

% For 
% pi=0.5
% tau=20
% rho=13
% alpha=0.5
% ff=0
% Vnh = 14.613503365498141843693460665524
% Vnl = 13.029675780363758631650353899221
% Indifferent if the investment cost in the new technology is Ii=[1519,1520)

% For 
% pi=0.5
% tau=20
% rho=14
% alpha=0.5
% ff=0
% Vnh = 14.428196345484752846279052679555
% Vnl = 12.845003176869997560350314405612
% Indifferent if the investment cost in the new technology is Ii=[1507,1508)

% For 
% pi=0.5
% tau=20
% rho=15
% alpha=0.5
% ff=0
% Vnh = 14.255444488876372558304577726479
% Vnl = 12.675467122152169829028044212453
% Indifferent if the investment cost in the new technology is Ii=[1496,1497)

% For 
% pi=0.5
% tau=20
% rho=16
% alpha=0.5
% ff=0
% Vnh = 14.094634405715486120344601397726
% Vnl = 12.520604789654250507515610697761
% Indifferent if the investment cost in the new technology is Ii=[1487,1488)

% For 
% pi=0.5
% tau=20
% rho=17
% alpha=0.5
% ff=0
% Vnh = 13.943999135277899207313658990138
% Vnl = 12.3777646974859630915849308504
% Indifferent if the investment cost in the new technology is Ii=[1479,1480)

% For 
% pi=0.5
% tau=20
% rho=18
% alpha=0.5
% ff=0
% Vnh = 13.802597908070256297754267859645
% Vnl = 12.245750059409344807838014816612
% Indifferent if the investment cost in the new technology is Ii=[1472,1473)

% For 
% pi=0.5
% tau=20
% rho=19
% alpha=0.5
% ff=0
% Vnh = 13.668945124104500350738010521146
% Vnl = 12.122363388707191221505473680309
% Indifferent if the investment cost in the new technology is Ii=[1465,1466)

% For 
% pi=0.5
% tau=20
% rho=20
% alpha=0.5
% ff=0
% Vnh = 13.543008900012680794653408292684
% Vnl = 12.007924525186449141120002784613
% Indifferent if the investment cost in the new technology is Ii=[1459,1460)

%------

% Alpha = 0.75

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0.75
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 19.075006735921200488440562231025
% Vnl = 18.435521943549791302120340037362
% Indifferent if the investment cost in the new technology is Ii=[1695,1696)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0.75
% ff=0
% Vnh = 18.331943564584676832653056357411
% Vnl = 17.330117240600504214526008859744
% Indifferent if the investment cost in the new technology is Ii=[1628,1629)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0.75
% ff=0
% Vnh = 17.716558928903125391999174571203
% Vnl = 16.488030692754325301782472382093
% Indifferent if the investment cost in the new technology is Ii=[1576,1577)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0.75
% ff=0
% Vnh = 17.196330782109051921420280984633
% Vnl = 15.819630837417059855039958616687
% Indifferent if the investment cost in the new technology is Ii=[1537,1538)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0.75
% ff=0
% Vnh = 16.749201064103246670420840055013
% Vnl = 15.273657371646158235599972084715
% Indifferent if the investment cost in the new technology is Ii=[1508,1509)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0.75
% ff=0
% Vnh = 16.359216811322282247263421224495
% Vnl = 14.817078319175461550301859996972
% Indifferent if the investment cost in the new technology is Ii=[1486,1487)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0.75
% ff=0
% Vnh = 16.014842120192349320697500777345
% Vnl = 14.428026293043104712235611481487
% Indifferent if the investment cost in the new technology is Ii=[1469,1470)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0.75
% ff=0
% Vnh = 15.707128704902304303572024458215
% Vnl = 14.090547887541945205643292179026
% Indifferent if the investment cost in the new technology is Ii=[1455,1456)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0.75
% ff=0
% Vnh = 15.430220243484280113617827075512
% Vnl = 13.795137055605088276517168064734
% Indifferent if the investment cost in the new technology is Ii=[1444,1445)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0.75
% ff=0
% Vnh = 15.178942348177842853588390301187
% Vnl = 13.533446337689497855918542420376
% Indifferent if the investment cost in the new technology is Ii=[1435,1436)

% For 
% pi=0.5
% tau=20
% rho=11
% alpha=0.75
% ff=0
% Vnh = 14.949051128096895979150245755782
% Vnl = 13.298841323177043037204129654724
% Indifferent if the investment cost in the new technology is Ii=[1427,1428)

% For 
% pi=0.5
% tau=20
% rho=12
% alpha=0.75
% ff=0
% Vnh = 14.737680422917386837802077066634
% Vnl = 13.087244030879495145361765500965
% Indifferent if the investment cost in the new technology is Ii=[1420,1421)

% For 
% pi=0.5
% tau=20
% rho=13
% alpha=0.75
% ff=0
% Vnh = 14.543196918694797432185313353139
% Vnl = 12.896752348394795960306133597935
% Indifferent if the investment cost in the new technology is Ii=[1415,1416)

% For 
% pi=0.5
% tau=20
% rho=14
% alpha=0.75
% ff=0
% Vnh = 14.362325469217224942053001059896
% Vnl = 12.722156273915450798752218046572
% Indifferent if the investment cost in the new technology is Ii=[1410,1411)

% For 
% pi=0.5
% tau=20
% rho=15
% alpha=0.75
% ff=0
% Vnh = 14.19416627733708638710125604033
% Vnl = 12.562679018556190966961862045891
% Indifferent if the investment cost in the new technology is Ii=[1406,1407)

% For 
% pi=0.5
% tau=20
% rho=16
% alpha=0.75
% ff=0
% Vnh = 14.036699810568800863430617223368
% Vnl = 12.415324360992664192157991898818
% Indifferent if the investment cost in the new technology is Ii=[1402,1403)

% For 
% pi=0.5
% tau=20
% rho=17
% alpha=0.75
% ff=0
% Vnh = 13.889488741456850347253101910513
% Vnl = 12.279937103403074197887836634465
% Indifferent if the investment cost in the new technology is Ii=[1399,1400)

% For 
% pi=0.5
% tau=20
% rho=18
% alpha=0.75
% ff=0
% Vnh = 13.750881593554330622972391743162
% Vnl = 12.154065593136474211883247153171
% Indifferent if the investment cost in the new technology is Ii=[1396,1397)

% For 
% pi=0.5
% tau=20
% rho=19
% alpha=0.75
% ff=0
% Vnh = 13.620052349835393098570592960417
% Vnl = 12.036717564459065646313088620768
% Indifferent if the investment cost in the new technology is Ii=[1393,1394)

% For 
% pi=0.5
% tau=20
% rho=20
% alpha=0.75
% ff=0
% Vnh = 13.496969709943635320851541942486
% Vnl = 11.928233834273685956825671863341
% Indifferent if the investment cost in the new technology is Ii=[1391,1392)

%-------

% Alpha = 1

% For 
% pi=0.5
% tau=20
% rho=1
% g=1
% alpha=1
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 18.99245315813295213351210770804
% Vnl = 18.180135613620714131572189603023
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=1
% ff=0
% Vnh = 18.224413258401657712256236912388
% Vnl = 17.033974854255577887120397428573
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=1
% ff=0
% Vnh = 17.606764301703377699993234865649
% Vnl = 16.207851572058070005204813897061
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=1
% ff=0
% Vnh = 17.092374312977055638457145806162
% Vnl = 15.568675866440863723945779948547
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=1
% ff=0
% Vnh = 16.653229578803857460064442621404
% Vnl = 15.051780790353014824685110550267
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=1
% ff=0
% Vnh = 16.271309186571179832008658998362
% Vnl = 14.620879943226023683528831990857
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=1
% ff=0
% Vnh = 15.934334528935537899263003044392
% Vnl = 14.253602973474097071736282460404
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=1
% ff=0
% Vnh = 15.633564799571109717441496672584
% Vnl = 13.935205364978083113784952794178
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=1
% ff=0
% Vnh = 15.362558723709445960239793454197
% Vnl = 13.655469227062237983705902688943
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=1
% ff=0
% Vnh = 15.116434690950190056731574439302
% Vnl = 13.407031317241854432177351983922
% Indifferent if the investment cost in the new technology is Ii=[1343,1344)
%-------


%------------------------------------


% Figures 5 and 8: Alpha = 0.5 only

%-- Parametrization 1:

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0.5
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 19.15533827507332185926176043157
% Vnl = 18.653960336192927052901029425648
% Vo = 18.992644823671112477022217124303
% Indifferent if the investment cost in the new technology is Ii=[2104,2105)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0.5
% ff=0
% Vnh = 18.452048466601192312778124874051
% Vnl = 17.631047556528398439982461081871
% Vo = 18.224714777552266690538782247751
% Indifferent if the investment cost in the new technology is Ii=[1994,1995)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0.5
% ff=0
% Vnh = 17.851474551066332415196634534046
% Vnl = 16.80801930241896094473433897082
% Vo = 17.607135239022531277902034331177
% Indifferent if the investment cost in the new technology is Ii=[1896,1897)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0.5
% ff=0
% Vnh = 17.331984738750575114487745118846
% Vnl = 16.128510972986086014088743529409
% Vo = 17.092791798466836818848648827557
% Indifferent if the investment cost in the new technology is Ii=[1814,1815)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0.5
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% Vnh = 16.878301689624960163027400990489
% Vnl = 15.558219621483619790563864808335
% Vo = 16.653679485560627852964273847932
% Indifferent if the investment cost in the new technology is Ii=[1747,1748)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0.5
% ff=0
% Vnh = 16.47940422889132184248792973263
% Vnl = 15.074874219929352759041405779093
% Vo = 16.271782226905795703729384429426
% Indifferent if the investment cost in the new technology is Ii=[1694,1695)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0.5
% ff=0
% Vnh = 16.125173517683213066461404648824
% Vnl = 14.659181512159272075106540842903
% Vo = 15.934824306078698540919644438243
% Indifferent if the investment cost in the new technology is Ii=[1651,1652)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0.5
% ff=0
% Vnh = 15.809419917529482587152319395947
% Vnl = 14.300413215033858692989406980426
% Vo = 15.634066746355999108438619906361
% Indifferent if the investment cost in the new technology is Ii=[1618,1619)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0.5
% ff=0
% Vnh = 15.524206522163911127357867119026
% Vnl = 13.98436191124448927059558411127
% Vo = 15.363069484444746565174800288555
% Indifferent if the investment cost in the new technology is Ii=[1590,1591)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0.5
% ff=0
% Vnh = 15.265520108025544335080504485313
% Vnl = 13.704763695040111463242070989103
% Vo = 15.116951741873366231112337562693
% Indifferent if the investment cost in the new technology is Ii=[1567,1568)

% For 
% pi=0.5
% tau=20
% rho=11
% alpha=0.5
% ff=0
% Vnh = 15.029509563743148646663063589039
% Vnl = 13.455513872522624089467339297735
% Indifferent if the investment cost in the new technology is Ii=[1548,1549)

% For 
% pi=0.5
% tau=20
% rho=12
% alpha=0.5
% ff=0
% Vnh = 14.81289411794680131575979660488
% Vnl = 13.231513104318970124937153311802
% Indifferent if the investment cost in the new technology is Ii=[1532,1533)

% For 
% pi=0.5
% tau=20
% rho=13
% alpha=0.5
% ff=0
% Vnh = 14.613503365498141843693460665524
% Vnl = 13.029675780363758631650353899221
% Indifferent if the investment cost in the new technology is Ii=[1519,1520)

% For 
% pi=0.5
% tau=20
% rho=14
% alpha=0.5
% ff=0
% Vnh = 14.428196345484752846279052679555
% Vnl = 12.845003176869997560350314405612
% Indifferent if the investment cost in the new technology is Ii=[1507,1508)

% For 
% pi=0.5
% tau=20
% rho=15
% alpha=0.5
% ff=0
% Vnh = 14.255444488876372558304577726479
% Vnl = 12.675467122152169829028044212453
% Indifferent if the investment cost in the new technology is Ii=[1496,1497)

% For 
% pi=0.5
% tau=20
% rho=16
% alpha=0.5
% ff=0
% Vnh = 14.094634405715486120344601397726
% Vnl = 12.520604789654250507515610697761
% Indifferent if the investment cost in the new technology is Ii=[1487,1488)

% For 
% pi=0.5
% tau=20
% rho=17
% alpha=0.5
% ff=0
% Vnh = 13.943999135277899207313658990138
% Vnl = 12.3777646974859630915849308504
% Indifferent if the investment cost in the new technology is Ii=[1479,1480)

% For 
% pi=0.5
% tau=20
% rho=18
% alpha=0.5
% ff=0
% Vnh = 13.802597908070256297754267859645
% Vnl = 12.245750059409344807838014816612
% Indifferent if the investment cost in the new technology is Ii=[1472,1473)

% For 
% pi=0.5
% tau=20
% rho=19
% alpha=0.5
% ff=0
% Vnh = 13.668945124104500350738010521146
% Vnl = 12.122363388707191221505473680309
% Indifferent if the investment cost in the new technology is Ii=[1465,1466)

% For 
% pi=0.5
% tau=20
% rho=20
% alpha=0.5
% ff=0
% Vnh = 13.543008900012680794653408292684
% Vnl = 12.007924525186449141120002784613
% Indifferent if the investment cost in the new technology is Ii=[1459,1460)

%--------

%-- Parametrization 2:

% For 
% pi=0.1
% tau=20
% rho=1
% alpha=0.5
% ff=0
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2
% fi0d=ff*5+5*2*ei0c-5*2*ri0
% fi1h=ff*5*(ei1hc-ri1)+5*(ei1hc-ri1)^2
% fi1hd=ff*5+5*2*ei1hc-5*2*ri1
% fi1l=ff*5*(ei1lc-ri1)+5*(ei1lc-ri1)^2
% fi1ld=ff*5+5*2*ei1lc-5*2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vo = 15.287075111100040593793791554515
% Vnh = 15.768865417138167673737365176067
% Vnl = 14.429237563493457976951481025045
% Indifferent if the investment cost in the new technology is Ii=[2105,2106)

% For 
% pi=0.1
% tau=20
% rho=2
% alpha=0.5
% ff=0
% Vo = 13.267426232629783694500684664382
% Vnh = 13.750964791928428322652299093359
% Vnl = 12.191978626436358891332688811099
% Indifferent if the investment cost in the new technology is Ii=[1997,1998)

% For 
% pi=0.1
% tau=20
% rho=3
% alpha=0.5
% ff=0
% Vo = 12.02779868264938653202151770644
% Vnh = 12.455921887265317105798130727015
% Vnl = 10.832524954404182196409256904624
% Indifferent if the investment cost in the new technology is Ii=[1901,1902)

% For 
% pi=0.1
% tau=20
% rho=4
% alpha=0.5
% ff=0
% Vo = 11.153470492629668165594926399275
% Vnh = 11.520685661636847037479970510185
% Vnl = 9.879844058655064467636307686135
% Indifferent if the investment cost in the new technology is Ii=[1819,1820)

% For 
% pi=0.1
% tau=20
% rho=5
% alpha=0.5
% ff=0
% Vo = 10.488100123732625638154487046621
% Vnh = 10.802271728120054946627509551384
% Vnl = 9.1648152622589687387996941898462
% Indifferent if the investment cost in the new technology is Ii=[1753,1754)

% For 
% pi=0.1
% tau=20
% rho=6
% alpha=0.5
% ff=0
% Vo = 9.9566855101389639122553877421692
% Vnh = 10.226207991470341814545240565229
% Vnl = 8.6018184722700146338474197561431
% Indifferent if the investment cost in the new technology is Ii=[1699,1700)

% For 
% pi=0.1
% tau=20
% rho=7
% alpha=0.5
% ff=0
% Vo = 9.5177807131044272739745370336105
% Vnh = 9.7520149745731727818486483376326
% Vnl = 8.1467710304117987016601864550676
% Indifferent if the investment cost in the new technology is Ii=[1657,1658)

% For 
% pi=0.1
% tau=20
% rho=8
% alpha=0.5
% ff=0
% Vo = 9.1462184763573729653028049467131
% Vnh = 9.351074991914020720974055985709
% Vnl = 7.7671505878652787001000808497372
% Indifferent if the investment cost in the new technology is Ii=[1622,1623)

% For 
% pi=0.1
% tau=20
% rho=9
% alpha=0.5
% ff=0
% Vo = 8.8256450013978363083265322427283
% Vnh = 9.0076620012671486360271935891511
% Vnl = 7.4470613757903087308126149980299
% Indifferent if the investment cost in the new technology is Ii=[1595,1596)

% For 
% pi=0.1
% tau=20
% rho=10
% alpha=0.5
% ff=0
% Vo = 8.5448750685094805256352938844373
% Vnh = 8.7069015173704995849947195140691
% Vnl = 7.1691147762477807266017891366632
% Indifferent if the investment cost in the new technology is Ii=[1571,1572)

% For 
% pi=0.1
% tau=20
% rho=11
% alpha=0.5
% ff=0
% Vo = 8.2959451469937330023093078915527
% Vnh = 8.4420880071548823320802567391689
% Vnl = 6.9274714122354384420252281332098
% Indifferent if the investment cost in the new technology is Ii=[1552,1553)

% For 
% pi=0.1
% tau=20
% rho=12
% alpha=0.5
% ff=0
% Vo = 8.0729990155086776383469205830941
% Vnh = 8.2058349381004635223588382128603
% Vnl = 6.7139012852047526332454856399058
% Indifferent if the investment cost in the new technology is Ii=[1536,1537)

% For 
% pi=0.1
% tau=20
% rho=13
% alpha=0.5
% ff=0
% Vo = 7.8716138175826786711378375724624
% Vnh = 7.9929326039135861647872045700945
% Vnl = 6.5228610532888154002651426259865
% Indifferent if the investment cost in the new technology is Ii=[1522,1523)

% For 
% pi=0.1
% tau=20
% rho=14
% alpha=0.5
% ff=0
% Vo = 7.6883737569252083709760834067361
% Vnh = 7.7998785961871236864394492305205
% Vnl = 6.3509585484608592288967332197919
% Indifferent if the investment cost in the new technology is Ii=[1510,1511)

% For 
% pi=0.1
% tau=20
% rho=15
% alpha=0.5
% ff=0
% Vo = 7.520590162880793954727333430416
% Vnh = 7.6232534953049955346991724010317
% Vnl = 6.1944760969665698606343402435495
% Indifferent if the investment cost in the new technology is Ii=[1499,1500)

% For 
% pi=0.1
% tau=20
% rho=16
% alpha=0.5
% ff=0
% Vo = 7.3661117601527039671958012298725
% Vnh = 7.4614987916055568783563907777592
% Vnl = 6.0523644619730473730641829248237
% Indifferent if the investment cost in the new technology is Ii=[1490,1491)

% For 
% pi=0.1
% tau=20
% rho=17
% alpha=0.5
% ff=0
% Vo = 7.2231925374891658792864294319815
% Vnh = 7.3121654950140198197320256646413
% Vnl = 5.9218689861215842822647941626081
% Indifferent if the investment cost in the new technology is Ii=[1482,1483)

% For 
% pi=0.1
% tau=20
% rho=18
% alpha=0.5
% ff=0
% Vo = 7.0903975362497809562926151518278
% Vnh = 7.1731518863661697909782507066055
% Vnl = 5.8005982159037261169177566157308
% Indifferent if the investment cost in the new technology is Ii=[1474,1475)

% For 
% pi=0.1
% tau=20
% rho=19
% alpha=0.5
% ff=0
% Vo = 6.9665342635667011519500104124607
% Vnh = 7.0444896834543846961419238470733
% Vnl = 5.6894332341499332038298678225614
% Indifferent if the investment cost in the new technology is Ii=[1468,1469)

% For 
% pi=0.1
% tau=20
% rho=20
% alpha=0.5
% ff=0
% Vo = 6.8506018309841645833863789019694
% Vnh = 6.9238982140079540756768618016326
% Vnl = 5.5854075059148331499121876241693
% Indifferent if the investment cost in the new technology is Ii=[1462,1463)

%--------

%-- Parametrization 3:

% For 
% pi=0.5
% tau=20
% rho=1
% alpha=0.5
% ff=40 (39.99)
% ci1h = (75-ei1hc)*ei1hc;
% ci1l = (25-ei1lc)*ei1lc;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1h=ff*(ei1hc-ri1)+(ei1hc-ri1)^2
% fi1hd=ff+2*ei1hc-2*ri1
% fi1l=ff*(ei1lc-ri1)+(ei1lc-ri1)^2
% fi1ld=ff+2*ei1lc-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1h = (ci1h+tau*ri1+i)^(rho+1)
% dia1hd = rho*(ci1h+tau*ri1+i)^(rho)
% dib1h = (ci1h+tau*ri1+fi1h+i)^(rho+1)
% dib1hd = rho*(ci1h+tau*ri1+fi1h+i)^(rho)
% dia1l = (ci1l+tau*ri1+i)^(rho+1)
% dia1ld = rho*(ci1l+tau*ri1+i)^(rho)
% dib1l = (ci1l+tau*ri1+fi1l+i)^(rho+1)
% dib1ld = rho*(ci1l+tau*ri1+fi1l+i)^(rho)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[2110,2111)

% For 
% pi=0.5
% tau=20
% rho=2
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[2006,2007)

% For 
% pi=0.5
% tau=20
% rho=3
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1911,1912)

% For 
% pi=0.5
% tau=20
% rho=4
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1830,1831)

% For 
% pi=0.5
% tau=20
% rho=5
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1763,1764)

% For 
% pi=0.5
% tau=20
% rho=6
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1708,1709)

% For 
% pi=0.5
% tau=20
% rho=7
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1665,1666)

% For 
% pi=0.5
% tau=20
% rho=8
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1629,1630)

% For 
% pi=0.5
% tau=20
% rho=9
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1601,1602)

% For 
% pi=0.5
% tau=20
% rho=10
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1577,1578)

% For 
% pi=0.5
% tau=20
% rho=11
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1557,1558)

% For 
% pi=0.5
% tau=20
% rho=12
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1540,1541)

% For 
% pi=0.5
% tau=20
% rho=13
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1526,1527)

% For 
% pi=0.5
% tau=20
% rho=14
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1513,1514)

% For 
% pi=0.5
% tau=20
% rho=15
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1503,1504)

% For 
% pi=0.5
% tau=20
% rho=16
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1493,1494)

% For 
% pi=0.5
% tau=20
% rho=17
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1485,1486)

% For 
% pi=0.5
% tau=20
% rho=18
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1477,1478)

% For 
% pi=0.5
% tau=20
% rho=19
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1470,1471)

% For 
% pi=0.5
% tau=20
% rho=20
% alpha=0.5
% ff=40 (39.99)
% Vnh = 0
% Vnl = 0
% Indifferent if the investment cost in the new technology is Ii=[1464,1465)










