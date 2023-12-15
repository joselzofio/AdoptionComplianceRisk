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

% This third script implements the simulations reported in section 4.2 
% of the article enhancing the baseline model with full technological uncertainty 
% about the efficiency of the new technology. As in the scenario with partial 
% technological uncertainty (script no. 2), the least efficient technology 
% presents high abatement costs (h) while the most efficient technology presents low 
% abatement costs (l). Now, as presented in eqs. (14), (16) and (18), the investment 
% and abatement decisions are taken without knowing the effiency of the technology 
% once it has been deployed.

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
% ri0=declared (reported) emissions with old technology
% ri1=declared (reported) emissions with new technology
% tau=tax on declared emissions
% rho=level of risk aversion
% ff=fixed part of the fine
% alpha=likelihood of least efficient abatement new technology

clear all;
global pi rho i tau ei0c ff alpha
syms ei0 ei1 ri0 ri1

% [1] Assign values for parameters
pi=0.5
tau=20
rho=1
ff=0
alpha=0.50

% [2] Assign abatement costs functions for technologies old (0), new with
% high abatement costs (1h) and new with low abatement costs (1l). Also,
% derivatives of abatement costs functions for new tech (1), ci1hd and
% ci1ld respectively
ci0=(100-ei0)*ei0;
ci1h=(75-ei1)*ei1;
ci1hd = 75 - 2 * ei1;
ci1l=(25-ei1)*ei1;
ci1ld = 25 - 2 * ei1;

% Optimal emissions with technology 0
opei0eqn=diff(ci0,ei0)+tau==0;
opei0 = vpasolve(opei0eqn, ei0, [-Inf Inf]);
ei0c=opei0;

% [3] Assign abatement cost funtion with tech 0 like in [2] substituting ei0 
% by ei0c.
% Assign santioning functions with ei0c for old tech keeping them equivalent 
% to the case with late technology uncertainty to allow comparability.
% Notice that fi0d corresponds to the first derivative of fi0 evaluated at 
% the violation level, this is at ei0c-ri0 (which is not the same as de 
% first derivative of fi0 evaluated at ri0). It must be updated manually.
% Also, describe the sanctioning function for the new technology fi1 in 
% terms of the unknown ei1.
% Finally, define the derivative of the sanctioning function with respect 
% to the violation level (ei-ri) evaluated at (ei-ri)=0
ci0=(100-ei0c)*ei0c;
fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2;  % Parametrizations 1 and 3: fi0=ff*(ei0c-ri0)+(ei0c-ri0)^2; Parametrization 2: fi0=ff*5*(ei0c-ri0)+5*(ei0c-ri0)^2;
fi0d=ff+2*ei0c-2*ri0; % Parametrizations 1 and 3: fi0d=ff+2*ei0c-2*ri0; Parametrization 2: fi0d=ff*5+2*5*ei0c-2*5*ri0;
fi1=ff*(ei1-ri1)+(ei1-ri1)^2;  % Parametrizations 1 and 3: fi1=ff*(ei1-ri1)+(ei1-ri1)^2; Parametrization 2: fi1=ff*5*(ei1-ri1)+5*(ei1-ri1)^2;
fi1d=ff+2*ei1-2*ri1; % Parametrizations 1 and 3: fi1d=ff+2*ei1-2*ri1; Parametrization 2: fi1d=ff*5+5*2*ei1-5*2*ri1;
fid_0 = ff + 2 * 0; % Parametrizations 1 and 3: fid_0 = ff + 2 * 0; Parametrization 2: fid_0 = ff * 5 + 5 * 2 * 0; 

% [4] Backwards loop for several i values. Do not directly specify i values
% in i, instead, set counter=maximum value of i to be studied

for counter=1343:-1:0

i=counter

% [5] Disutility function shaped as a power function. Notice that dia0d and dib0d correspond to the first
% derivatives with respect to ai0=(ci0+tau*ri0) and bi0=(ci0+tau*ri0+fi0)
% of dia0 and dib0 respectively (old tech), and they must be updated manually at this
% step. The same happens with dia1hd, dib1hd, dia1ld and dib1ld (new tech):
dia0 = (ci0+tau*ri0)^(rho+1);
dia0d = (rho+1)*(ci0+tau*ri0)^(rho);
dib0 = (ci0+tau*ri0+fi0)^(rho+1);
dib0d = (rho+1)*(ci0+tau*ri0+fi0)^(rho);

dia1h = (ci1h+tau*ri1+i)^(rho+1);
dia1hd = (rho+1)*(ci1h+tau*ri1+i)^(rho);
dib1h = (ci1h+tau*ri1+fi1+i)^(rho+1);
dib1hd = (rho+1)*(ci1h+tau*ri1+fi1+i)^(rho);
dia1l = (ci1l+tau*ri1+i)^(rho+1);
dia1ld = (rho+1)*(ci1l+tau*ri1+i)^(rho);
dib1l = (ci1l+tau*ri1+fi1+i)^(rho+1);
dib1ld = (rho+1)*(ci1l+tau*ri1+fi1+i)^(rho);


% Old Technology

% Proposition 1: Equation to implicitly obtain declared emissions with the
% old technology
eqn = (dib0d*pi*fi0d)/((1-pi)*dia0d+pi*dib0d) == tau;

% [6] The function 'vpasolve' solves the equation numerically. With the range
% [0 ei0c] we consider only real solutions beteween 0 (which consitutes the
% lower bound because firms won't declare negative emissions) and the
% optimal emissions level with the old tech (which constitutes the upper
% bound because firms won't declaring more emissions than actual
% emissions).
% Be careful to keep Fo the same as the sanctioning function in [3]
% substituting ri0 by Old (the amount of declared emissions with the old
% technology we find with vpasolve)

solution = vpasolve(eqn, ri0,[0 ei0c]); %Amount of declared emissions with the old technology

Old = solution;

Vo=ei0c-Old; %Violation level with the old technology

Co=ci0; %Abatement costs with the old technology

To=Old*tau; %Total taxes paid for declared emissions with the old technology

Fo=ff*(ei0c-Old)+(ei0c-Old)^2; %Sanction for the violation level with the old technology, Parametrizations 1 and 3: Fo=ff*(ei0c-Old)+(ei0c-Old)^2; Parametrization 2: Fo=ff*5*(ei0c-Old)+5*(ei0c-Old)^2;

% To get an idea about how large are total costs (the independent variable)
Costso=Co+To+Fo;

% [7] Disutility assuming the old technology. Make sure that the elements 
% multiplying (1-pi) and pi are the same disutility functions written for [5] 
% but substituting ci0, tau*ri and fi0 by Co, To and Fo respectively.
Do=(1-pi)*((Co+To)^(rho+1))+pi*((Co+To+Fo)^(rho+1));


% New Technology

% Define scenario 3 equation to find optimal pollution level and the two
% first order conditions with lambda=0:
scenario3_eq = alpha * ((rho+1)*(ci1h+tau*ei1+i)^(rho)) * (ci1hd + tau) + (1 - alpha) * ((rho+1)*(ci1l+tau*ei1+i)^(rho)) * (ci1ld + tau);
foc_ei = alpha * ((1 - pi) * dia1hd + pi * dib1hd) * ci1hd + (1 - alpha) * ((1 - pi) * dia1ld + pi * dib1ld) * ci1ld + (alpha * dib1hd + (1 - alpha) * dib1ld) * pi * fi1d;
foc_ri = alpha * ((1 - pi) * dia1hd + pi * dib1hd) * tau + (1 - alpha) * ((1 - pi) * dia1ld + pi * dib1ld) * tau - (alpha * dib1hd + (1 - alpha) * dib1ld) * pi * fi1d;
focs = [foc_ei foc_ri];

% Condition for full compliance:
condition = pi * fid_0 - tau;

% Check the condition. 
% Notice we might get two sets of solutions, one for emisions (ei1) and 
% another one for reported emisions (ri1). We will use only real positive 
% ei1 and ri1 solutions, that we store in real_positive_esol and in real_positive_rsol
% respectively.
if condition >= 0
    % If condition is satisfied, solve scenario 3
    solution = vpasolve(scenario3_eq, ei1, [0 75]);

    esol=solution;
    
    % Initialize arrays to store real positive solutions
    real_positive_esol = [];
    real_positive_rsol = [];

    % Loop through the solutions
    for j = 1:length(esol)
        if isreal(esol(j)) && esol(j) > 0 && esol(j) <= 75
            real_positive_esol = [real_positive_esol; esol(j)];
            real_positive_rsol = [real_positive_rsol; esol(j)];
        end
    end
    disp('Real and positive esol values:');
    disp(real_positive_esol);
    
    disp('Real and positive rsol values:');
    disp(real_positive_rsol);
    
    if real_positive_esol <= 25

        % [8a] Be careful to update Cnh and Cnl to account for the abatement cost functions
        % defined in [2] substituting ei1 by En.
        % Be careful to keep Fo the same as the santioning function in [3]
        % substituting ei1 and ri1 by En and Rn respectively.
        En = real_positive_esol; %Amount of emissions with the new technology
            
        Rn = real_positive_rsol; %Amount of declared emissions with the new technology
            
        Vn=En-Rn; %Violation level with the new technology
        
        Cnh = (75 - En) * En; %Abatement cost with the high abatement cost new technology
        
        Cnl = (25 - En) * En; %Abatement cost with the low abatement cost new technology
        
        Tn=Rn*tau; %Total taxes paid for declared emissions with the new technology
        
        In=i; %Fixed investment cost in the new technology
        
        Fn=ff*(En-Rn)+(En-Rn)^2; %Sanction for the violation level with the new technology, Parametrizations 1 and 3: Fn=ff*(En-Rn)+(En-Rn)^2; Parametrization 2: Fn=ff*5*(En-Rn)+5*(En-Rn)^2;
        
        % To get an idea about how large total costs are (the independent variable)
        Costsnh=Cnh+Tn+In+Fn;
        Costsnl=Cnl+Tn+In+Fn;
        
        % [9a] Expected disutility assuming the new technology. Make sure that the elements multiplying (1-pi) and pi are the same
        % disutility functions written for [5] but substituting ci1h, ci1l, tau*ri1, i and 
        % fi1 by Cnh, Cnl, Tn, In and Fn, respectively
        Dn = alpha * ((1-pi)*((Cnh+Tn+In)^(rho+1))+pi*((Cnh+Tn+In+Fn)^(rho+1))) + (1 - alpha) * ((1-pi)*((Cnl+Tn+In)^(rho+1))+pi*((Cnl+Tn+In+Fn)^(rho+1)));
        
        % Now we check if the disutility with the old technology is larger than the
        % disutility with the new technology under full uncertainty. When this is 
        % the case, investing in the new technology becomes the best option. This 
        % indicates the program that an aproximate Ii threshold has been found and 
        % therefore it can end the loop.
        if Do-Dn>0
            break
        end
    else
        scenario3_eq_l0 = alpha * ((rho+1)*(ci1h+tau*ei1+i)^(rho)) * (ci1hd + tau) + (1 - alpha) * ((rho+1)*(0+tau*ei1+i)^(rho)) * (0 + tau);
        solution = vpasolve(scenario3_eq_l0, ei1, [0 75]);
    
        esol=solution;
        
        % Initialize arrays to store real positive solutions
        real_positive_esol = [];
        real_positive_rsol = [];
    
        % Loop through the solutions
        for j = 1:length(esol)
            if isreal(esol(j)) && esol(j) > 0 && esol(j) <= 75
                real_positive_esol = [real_positive_esol; esol(j)];
                real_positive_rsol = [real_positive_rsol; esol(j)];
            end
        end
        disp('Real and positive esol values:');
        disp(real_positive_esol);
        
        disp('Real and positive rsol values:');
        disp(real_positive_rsol);
        % [8b] Be careful to update Cnh and Cnl to account for the abatement cost functions
        % defined in [2] substituting ei1 by En.
        % Be careful to keep Fo the same as the santioning function in [3]
        % substituting ei1 and ri1 by En and Rn respectively.
        En = real_positive_esol; %Amount of emissions with the new technology
            
        Rn = real_positive_rsol; %Amount of declared emissions with the new technology
            
        Vn=En-Rn; %Violation level with the new technology
        
        Cnh = (75 - En) * En; %Abatement cost with the high abatement cost new technology
        
        Cnl = 0; %Abatement cost with the low abatement cost new technology
        
        Tn=Rn*tau; %Total taxes paid for declared emissions with the new technology
        
        In=i; %Fixed investment cost in the new technology
        
        Fn=ff*(En-Rn)+(En-Rn)^2; %Sanction for the violation level with the new technology, Parametrizations 1 and 3: Fn=ff*(En-Rn)+(En-Rn)^2; Parametrization 2: Fn=ff*5*(En-Rn)+5*(En-Rn)^2;
        
        % To get an idea about how large total costs are (the independent variable)
        Costsnh=Cnh+Tn+In+Fn;
        Costsnl=Cnl+Tn+In+Fn;
        
        % [9b] Expected disutility assuming the new technology. Make sure that the elements multiplying (1-pi) and pi are the same
        % disutility functions written for [5] but substituting ci1h, ci1l, tau*ri1, i and 
        % fi1 by Cnh, Cnl, Tn, In and Fn, respectively
        Dn = alpha * ((1-pi)*((Cnh+Tn+In)^(rho+1))+pi*((Cnh+Tn+In+Fn)^(rho+1))) + (1 - alpha) * ((1-pi)*((Cnl+Tn+In)^(rho+1))+pi*((Cnl+Tn+In+Fn)^(rho+1)));
        
        % Now we check if the disutility with the old technology is larger than the
        % disutility with the new technology under full uncertainty. When this is 
        % the case, investing in the new technology becomes the best option. This 
        % indicates the program that an aproximate Ii threshold has been found and 
        % therefore it can end the loop.
        if Do-Dn>0
            break
        end
    end
else
    % If condition is not satisfied, solve the system of equations focs = [foc_ei foc_ri]
    vars = [ei1 ri1];
    range = [0 75; 0 75];
    solution = vpasolve(focs,vars,range);

    esol=solution.ei1;
    rsol=solution.ri1;
    
    % Initialize arrays to store real positive solutions
    real_positive_esol = [];
    real_positive_rsol = [];

    % Loop through the solutions
    for j = 1:length(esol)
        if isreal(esol(j)) && isreal(rsol(j)) && esol(j) > 0 && rsol(j) > 0 && esol(j) <= 75 && rsol(j) <= Old
            real_positive_esol = [real_positive_esol; esol(j)];
            real_positive_rsol = [real_positive_rsol; rsol(j)];
        end
    end
    disp('Real and positive esol values:');
    disp(real_positive_esol);
    
    disp('Real and positive rsol values:');
    disp(real_positive_rsol);
    
    if real_positive_esol <= 25

        % [8c] Be careful to update Cnh and Cnl to account for the abatement cost functions
        % defined in [2] substituting ei1 by En.
        % Be careful to keep Fo the same as the santioning function in [3]
        % substituting ei1 and ri1 by En and Rn respectively.
        En = real_positive_esol; %Amount of emissions with the new technology
            
        Rn = real_positive_rsol; %Amount of declared emissions with the new technology
            
        Vn=En-Rn; %Violation level with the new technology
        
        Cnh = (75 - En) * En; %Abatement cost with the high abatement cost new technology
        
        Cnl = (25 - En) * En; %Abatement cost with the low abatement cost new technology
        
        Tn=Rn*tau; %Total taxes paid for declared emissions with the new technology
        
        In=i; %Fixed investment cost in the new technology
        
        Fn=ff*(En-Rn)+(En-Rn)^2; %Sanction for the violation level with the new technology, Parametrizations 1 and 3: Fn=ff*(En-Rn)+(En-Rn)^2; Parametrization 2: Fn=ff*5*(En-Rn)+5*(En-Rn)^2;
        
        % To get an idea about how large total costs are (the independent variable)
        Costsnh=Cnh+Tn+In+Fn;
        Costsnl=Cnl+Tn+In+Fn;
        
        % [9c] Expected disutility assuming the new technology. Make sure that the elements multiplying (1-pi) and pi are the same
        % disutility functions written for [5] but substituting ci1h, ci1l, tau*ri1, i and 
        % fi1 by Cnh, Cnl, Tn, In and Fn, respectively
        Dn = alpha * ((1-pi)*((Cnh+Tn+In)^(rho+1))+pi*((Cnh+Tn+In+Fn)^(rho+1))) + (1 - alpha) * ((1-pi)*((Cnl+Tn+In)^(rho+1))+pi*((Cnl+Tn+In+Fn)^(rho+1)));
        
        % Now we check if the disutility with the old technology is larger than the
        % disutility with the new technology under full uncertainty. When this is 
        % the case, investing in the new technology becomes the best option. This 
        % indicates the program that an aproximate Ii threshold has been found and 
        % therefore it can end the loop.
        if Do-Dn>0
            break
        end
    else
        foc_ei_l0 = alpha * ((1 - pi) * dia1hd + pi * dib1hd) * ci1hd + 0 + (alpha * dib1hd + (1 - alpha) * (rho+1)*(0+tau*ri1+fi1+i)^(rho)) * pi * fi1d;
        foc_ri_l0 = alpha * ((1 - pi) * dia1hd + pi * dib1hd) * tau + (1 - alpha) * ((1 - pi) * (rho+1)*(0+tau*ri1+i)^(rho) + pi * (rho+1)*(0+tau*ri1+fi1+i)^(rho)) * tau - (alpha * dib1hd + (1 - alpha) * (rho+1)*(0+tau*ri1+fi1+i)^(rho)) * pi * fi1d;
        focs = [foc_ei_l0 foc_ri_l0];

        vars = [ei1 ri1];
        range = [0 75; 0 75];
        solution = vpasolve(focs,vars,range);
    
        esol=solution.ei1
        rsol=solution.ri1
        
        % Initialize arrays to store real positive solutions
        real_positive_esol = [];
        real_positive_rsol = [];
    
        % Loop through the solutions
        for j = 1:length(esol)
            if isreal(esol(j)) && isreal(rsol(j)) && esol(j) > 0 && rsol(j) > 0 && esol(j) <= 75 && rsol(j) <= Old
                real_positive_esol = [real_positive_esol; esol(j)];
                real_positive_rsol = [real_positive_rsol; rsol(j)];
            end
        end
        disp('Real and positive esol values:');
        disp(real_positive_esol);
        
        disp('Real and positive rsol values:');
        disp(real_positive_rsol);
        
        % [8d] Be careful to update Cnh and Cnl to account for the abatement cost functions
        % defined in [2] substituting ei1 by En.
        % Be careful to keep Fo the same as the santioning function in [3]
        % substituting ei1 and ri1 by En and Rn respectively.
        En = real_positive_esol; %Amount of emissions with the new technology
            
        Rn = real_positive_rsol; %Amount of declared emissions with the new technology
            
        Vn=En-Rn; %Violation level with the new technology
        
        Cnh = (75 - En) * En; %Abatement cost with the high abatement cost new technology
        
        Cnl = 0; %Abatement cost with the low abatement cost new technology
        
        Tn=Rn*tau; %Total taxes paid for declared emissions with the new technology
        
        In=i; %Fixed investment cost in the new technology
        
        Fn=ff*(En-Rn)+(En-Rn)^2; %Sanction for the violation level with the new technology, Parametrizations 1 and 3: Fn=ff*(En-Rn)+(En-Rn)^2; Parametrization 2: Fn=ff*5*(En-Rn)+5*(En-Rn)^2;
        
        % To get an idea about how large total costs are (the independent variable)
        Costsnh=Cnh+Tn+In+Fn;
        Costsnl=Cnl+Tn+In+Fn;
        
        % [9d] Expected disutility assuming the new technology. Make sure that the elements multiplying (1-pi) and pi are the same
        % disutility functions written for [5] but substituting ci1h, ci1l, tau*ri1, i and 
        % fi1 by Cnh, Cnl, Tn, In and Fn, respectively
        Dn = alpha * ((1-pi)*((Cnh+Tn+In)^(rho+1))+pi*((Cnh+Tn+In+Fn)^(rho+1))) + (1 - alpha) * ((1-pi)*((Cnl+Tn+In)^(rho+1))+pi*((Cnl+Tn+In+Fn)^(rho+1)));
        
        % Now we check if the disutility with the old technology is larger than the
        % disutility with the new technology under full uncertainty. When this is 
        % the case, investing in the new technology becomes the best option. This 
        % indicates the program that an aproximate Ii threshold has been found and 
        % therefore it can end the loop.
        if Do-Dn>0
            break
        end
    end
end

end

X=['Indifferent if the installing cost of the new technology (Ii) is aproximately between [',num2str(i),',',num2str(i+1),')'];
disp(X)
fprintf('Optimal emissions with the old technology=%s\n',opei0);
fprintf('Optimal emissions with the new technology=%s\n',En);
fprintf('Violation level with the old technology=%s\n',Vo);
fprintf('Violation level with the new technology=%s\n',Vn);
