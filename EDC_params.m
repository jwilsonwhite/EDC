function EDC_params(f_fit_xy, m_fit_xx,alphaBH)

% Parameter file for the age-structured EDC model
% These model runs designed for White, Cole, Cherr, Connon, & Brander
% manuscript (Ecol. Appl.)
% Will White - July 2016 
% whitejw@uncw.edu

% Inputs:
% f_fit_xy: fitness of XY females
% m_fit_xx: fitness of XX males
% alphaBH: slope of Beverton-Holt S-R curve (specified as the CRT)
if ~exist('L_reduction','var')
    L_reduction = 1;
end

% General parameters common to all runs:
T = 20; % number of years per simulation

% c2: probability of female success at skewed sex ratios
% range from 0 (males harass only at very skewed sex ratios) to 1 (male
% harassment is intense & linear increasing function of sex ratio)
c2 = 0; % mate harassment (smaller values = less harrassment)

%b_alpha = 1; % first shape parameter of Beta cdf mating function. Value of 1 forces it to always have negative concavity

% relative fitness of aberrant indivs
f_fit_xy = f_fit_xy; % fitness of XY females
f_fit_yy = 0; % fitness of YY females
m_fit_xx = m_fit_xx; % fitness of XX males
m_fit_yy = 0; % fitness of YY males


% fish life history parameters
amax = 2; % number of age classes, in years
amat = 1; % mature after 1st year
n_init = 0.1; % initial population size
M = 2.10; % baseline natural mortality rate (mean of seasonal rates if there is seasonal variation)
         % From Fishbase
         
Linf = 13; %  (cm) Fishbase (M. menidia)
k = 1.39; % (1/year) Fishbase (M. menidia, female)
t0 = -0.31; % (years) Fishbase (M. menidia)

% No length-fecundity relationship, so assume cubic
fec_cons = 1;
fec_exp = 3;

% Length & Fec vectors
Age = 0:(amax-1);
Len = Linf*(1-exp(-k.*(Age-t0)));
Fec = fec_cons.*Len.^fec_exp;
Fec(Age<amat) = 0; % no fecundity below age of maturity

% Survival to age
L = exp(-M*(0:(amax-1))); 

% Creata a transition matrix to advance each age class with mortality
% Could make this have changing mortality rates with age
A = diag(ones(1,amax-1).*exp(-M)); % vector of constant mortality with age. 
A(:,end+1) = 0;
A = [zeros(1,size(A,2));A]; 


% Beverton-Holt
%alphaBH = 0.31; 
betaBH = 1; % Beverton-Holt max (arbitrary)

save EDC_params.mat