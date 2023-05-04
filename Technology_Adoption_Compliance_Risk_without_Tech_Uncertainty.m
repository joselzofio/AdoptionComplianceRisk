% This MATLAB script calculates the investment thresholds, along with 
% optimal emission and violation levels from the simulations reported in 
% Arguedas, C., Peinado, F., and Zofio, J.L. (2023) "Incentives for Green 
% Technology Adoption, Imperfect Compliance, and Risk Aversion". 
% Optimal emissions are calculated by solving the condition presented in 
% equation (2) of the article, minimizing firms' expected disutility in terms of 
% compliance and non-compliance (see equation (1)). 
% The script has been run in version R2022a of MATLAB and uses the function
% 'vpasolve' from the Symbolic Math Toolbox to determine optimal emissions.

% This first script implements the simulations reported in section 3.3 
% of the article without technological uncertainty about the new
% technology.

% Follow steps [1], [2], ... to ensure correct outputs.

% The simulations can be performed with the different parametrizations of 
% the inspection probabilities and fines for non-compliance presented in Table 1 
% of the article. The parameters corresponding to the tax value (tau) and degree 
% of risk aversion (rho) can be also adjusted. At the end of the script we present 
% the parameters' values necessary to replicate all the simulations reported 
% in the Figures of the article.

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
fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2;
fi0d=ff+2*ei0c-2*ri0;
fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2;
fi1d=ff+2*ei1c-2*ri1;

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

Fo=ff*(ei0c-O)+(ei0c-O)^2; % Fine for the violation level with the old technology

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

Fn=ff*(ei1c-N)+(ei1c-N)^2; % Fine for the violation level with the least efficient new technology

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

% Data to replicate the Figures in the article:

% Figures 1 and 2: Ii = Ii^*

%-- Parametrization 1:

% For 
% pi=0.5
% tau=20
% rho=1
% ff=0
% ci0=(100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2;
% fi0d=ff+2*ei0c-2*ri0;
% fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2;
% fi1d=ff+2*ei1c-2*ri1;
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 18.992644823671112477022217124303
% Vn = 18.992389253299664315722833887169
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=2
% ff=0
% Vo = 18.224714777552266690538782247751
% Vnh = 18.224312728814620915840597852519
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=3
% ff=0
% Vo = 17.607135239022531277902034331177
% Vnh = 17.606640629491407504760480571295
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=4
% ff=0
% Vo = 17.092791798466836818848648827557
% Vnh = 17.092235123156473405865961300384
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=5
% ff=0
% Vo = 16.653679485560627852964273847932
% Vnh = 16.653079581210972428179542096421
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=6
% % ff=0
% Vo = 16.271782226905795703729384429426
% Vnh = 16.271151477584728519339717218874
% Indifferent if investing cost of new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=7
% ff=0
% Vo = 15.934824306078698540919644438243
% Vnh = 15.934171241093750733511752778254
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=8
% ff=0
% Vo = 15.634066746355999108438619906361
% Vnh = 15.633397455433187564448172413402
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=9
% ff=0
% Vo = 15.363069484444746565174800288555
% Vnh = 15.362388441944105726631583802825
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=10
% ff=0
% Vo = 15.116951741873366231112337562693
% Vnh = 15.116262312876746654388060707205
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

%-------

%-- Parametrization 2:

% For 
% pi=0.1
% tau=20
% rho=1
% ff=0
% ci0=(100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2
% fi0d=ff*5+5*2*ei0c-5*2*ri0
% fi1=ff*5*(ei1c-ri1)+5*(ei1c-ri1)^2
% fi1d=ff*5+5*2*ei1c-5*2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 15.287075111100040593793791554515
% Vnh = 15.286359628018830585533203135384
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=2
% ff=0
% Vo = 13.267426232629783694500684664382
% Vnh = 13.267426232629783694500684664382
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=3
% ff=0
% Vo = 12.02779868264938653202151770644
% Vnh = 12.026971385016763025375409921001
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=4
% ff=0
% Vo = 11.153470492629668165594926399275
% Vnh = 11.152649706156618601295263424074
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=5
% ff=0
% Vo = 10.488100123732625638154487046621
% Vnh = 10.488100123732625638154487046621
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=6
% ff=0
% Vo = 9.9566855101389639122553877421692
% Vnh = 9.9558937276118008592638615913962
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=7
% ff=0
% Vo = 9.5177807131044272739745370336105
% Vnh = 9.5170048309416955635362611886562
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=8
% ff=0
% Vo = 9.1462184763573729653028049467131
% Vnh = 9.1454580943735112747023335134339
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=9
% ff=0
% Vo = 8.8256450013978363083265322427283
% Vnh = 8.8248994311692952188444651853791
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.1
% tau=20
% rho=10
% ff=0
% Vo = 8.5448750685094805256352938844373
% Vnh = 8.5448750685094805256352938844373
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

%--------

%-- Parametrization 3:

% For 
% pi=0.5
% tau=20
% rho=1
% ff=40
% ci0 = (100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2
% fi1d=ff+2*ei1c-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=2
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=3
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=4
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=5
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=6
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=7
% ff=40
% Vo = 0.00
% Vnh =  0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=8
% ff=40
% Vo = 0.00
% Vnh = 0.00
% IIndifferent if the investment cost in the new technology is Ii=[2374,2375)

% For  
% pi=0.5
% tau=20
% rho=9
% ff=40
% Vo = 0.00
% Vnh =  0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)

% For 
% pi=0.5
% tau=20
% rho=10
% ff=40
% Vo = 0.00
% Vnh = 0.00
% Indifferent if the investment cost in the new technology is Ii=[2374,2375)


%------------------------------------


% Figure 3: Ii = 1000 < Ii^*

%-- Parametrization 1:

% For 
% pi=0.5
% tau=20
% rho=1
% ff=0
% ci0 = (100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2
% fi1d=ff+2*ei1c-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 18.992644823671112477022217124303
% Vnh = 18.450157086167419520251790257272

% For 
% pi=0.5
% tau=20
% rho=2
% ff=0
% Vo = 18.224714777552266690538782247751
% Vnh = 17.411917088124442386594245410838

% For 
% pi=0.5
% tau=20
% rho=3
% ff=0
% Vo = 17.607135239022531277902034331177
% Vnh = 16.638859077313608480352102024707

% For 
% pi=0.5
% tau=20
% rho=4
% ff=0
% Vo = 17.092791798466836818848648827557
% Vnh = 16.028061568153501095997637568113

% For 
% pi=0.5
% tau=20
% rho=5
% ff=0
% Vo = 16.653679485560627852964273847932
% Vnh = 15.526588979476975095326678036994

% For 
% pi=0.5
% tau=20
% rho=6
% ff=0
% Vo = 16.271782226905795703729384429426
% Vnh = 15.103619556607416914666639066134

% For 
% pi=0.5
% tau=20
% rho=7
% ff=0
% Vo = 15.934824306078698540919644438243
% Vnh = 14.739644784033047409004986551451

% For 
% pi=0.5
% tau=20
% rho=8
% ff=0
% Vo = 15.634066746355999108438619906361
% Vnh = 14.421547875627903538002373154487

% For 
% pi=0.5
% tau=20
% rho=9
% ff=0
% Vo = 15.363069484444746565174800288555
% Vnh = 14.140095341747797759865430423882

% For 
% pi=0.5
% tau=20
% rho=10
% ff=0
% Vo = 15.116951741873366231112337562693
% Vnh = 13.888548131583097097729310765745

%--------

%-- Parametrization 2:

% For 
% pi=0.1
% tau=20
% rho=1
% ff=0
% ci0 = (100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2
% fi0d=ff*5+5*2*ei0c-5*2*ri0
% fi1=ff*5*(ei1c-ri1)+5*(ei1c-ri1)^2
% fi1d=ff*5+5*2*ei1c-5*2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 15.287075111100040593793791554515
% Vnh = 13.980938546795589217810439512991

% For 
% pi=0.1
% tau=20
% rho=2
% ff=0
% Vo = 13.267426232629783694500684664382
% Vnh = 11.839699918778410792009159299261

% For 
% pi=0.1
% tau=20
% rho=3
% ff=0
% Vo = 12.02779868264938653202151770644
% Vnh = 10.600700045173216717180975460378

% For 
% pi=0.1
% tau=20
% rho=4
% ff=0
% Vo = 11.153470492629668165594926399275
% Vnh = 9.7539982815644691842396193545668

% For 
% pi=0.1
% tau=20
% rho=5
% ff=0
% Vo = 10.488100123732625638154487046621
% Vnh = 9.1226374882394051593953634269781

% For 
% pi=0.1
% tau=20
% rho=6
% ff=0
% Vo = 9.9566855101389639122553877421692
% Vnh = 8.6256343944303408587943607270615

% For 
% pi=0.1
% tau=20
% rho=7
% ff=0
% Vo = 9.5177807131044272739745370336105
% Vnh = 8.2196264443296564964567961430436

% For 
% pi=0.1
% tau=20
% rho=8
% ff=0
% Vo = 9.1462184763573729653028049467131
% Vnh = 7.8788804180834660734215533160684

% For 
% pi=0.1
% tau=20
% rho=9
% ff=0
% Vo = 8.8256450013978363083265322427283
% Vnh = 7.5869676607955885129689307170346

% For 
% pi=0.1
% tau=20
% rho=10
% ff=0
% Vo = 8.5448750685094805256352938844373
% Vnh = 7.332809347414298703596383569944

%--------

%-- Parametrization 3:

% For 
% pi=0.5
% tau=20
% rho=1
% ff=40
% ci0 = (100-ei0c)*ei0c;
% ci1 = (50-ei1c)*ei1c;
% fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2
% fi0d=ff+2*ei0c-2*ri0
% fi1=ff*(ei1c-ri1)+(ei1c-ri1)^2
% fi1d=ff+2*ei1c-2*ri1
% dia0 = (ci0+tau*ri0)^(rho+1)
% dia0d = rho*(ci0+tau*ri0)^(rho)
% dib0 = (ci0+tau*ri0+fi0)^(rho+1)
% dib0d = rho*(ci0+tau*ri0+fi0)^(rho)
% dia1 = (ci1+tau*ri1+i)^(rho+1)
% dia1d = rho*(ci1+tau*ri1+i)^(rho)
% dib1 = (ci1+tau*ri1+fi1+i)^(rho+1)
% dib1d = rho*(ci1+tau*ri1+fi1+i)^(rho)
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=2
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=3
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=4
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=5
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=6
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=7
% ff=40
% Vo = 0.00
% Vnh =  0.00

% For 
% pi=0.5
% tau=20
% rho=8
% ff=40
% Vo = 0.00
% Vnh = 0.00

% For 
% pi=0.5
% tau=20
% rho=9
% ff=40
% Vo = 0.00
% Vnh =  0.00

% For 
% pi=0.5
% tau=20
% rho=10
% ff=40
% Vo = 0.00
% Vnh =  0.00








