function [Ntotal, Ph_rat] = EDC_popmodel_as(varargin)


% Age-structured EDC model
% As described in White, Cole, Cherr, Connon, & Brander manuscript (Ecol. Appl.)
% Will White - July 2016 
% whitejw@uncw.edu

% EDC_popmodel_as.m
% Implements modified version of genotype-specific model from Cotton &
% Wedekind (2009, Cons Biol)
% This version: age-structured discrete-time model with Beverton-Holt DD,
% and adjustable spawning functional response

% This code is designed to be called by EDC_runs.m

load EDC_params.mat

% Input arguments:
% c1: probability of finding a male successfully.
% range from 0 (only need a few males - group spawning) to 1 (linear
% increase - monogamous pair spawning)
c1 = varargin{1};

b_alpha = varargin{2}; % first shape parameter of Beta cdf mating function. Value of 1 forces it to always have negative concavity

% Vectors of masculinization/feminization rates
P = varargin{3}; % masculinization rate
Q = varargin{4}; % feminization rate

% relative fitness of all females due to EDC effects
f_fit_EDC = varargin{5};

% length reduction due to EDCs
L_reduction = varargin{6};

% State variables: (may need to divvy these up using structures if
% dimensionality gets too high
% Nm: males.  Tx3xPxQ matrix: [XX XY YY]
% Nf: females: same dimensions as Nm
Nm= nan(T,3,length(P),length(Q));
Nf = Nm;
NN = nan(T,length(P),length(Q)); % total pop size
XX = NN; % sex ratio

%%% First: get no-EDC solution

% Probability of mate-finding (first argument is 1 bc the function maxes at
% sex ratio of 0.5)
p_male = min(1,betacdf(0.5*2, b_alpha, c1)); % 
% Probability of female success (due to male harassment at high male ratio)
p_female = (1-0.5)./(1- 0.5 + c2);

% include these in fecundity
Fec_tmp = Fec.*p_male.*p_female;

R0 = sum(0.5*L(:).*Fec_tmp(:)); % total reproductive output per female recruit (hence factor 0.5)
alphaBH = 1./(alphaBH*R0); % translate CRT into number of fish

Req = betaBH.*(1-1/(alphaBH*R0)); % this is the solution to the bev-holt equation (number of recruits at equilibrium)
Neq = sum(Req.*L); % total population size at equilibrium

% If there is an EDC effect on length, apply that here (after R0 has been
% calculated for unimpacted population):
Len = Len*L_reduction;
Fec = fec_cons.*Len.^fec_exp;
Fec(Age<amat) = 0;

% Loop over EDC effect rates
for pp = 1:length(P)
    for qq = 1:length(Q)
    
    p = P(pp);
    q = Q(qq);

    % Define age-structured vectors (temporary, do not store)
    % A x 3 x T matrix: [XX XY YY]
    Nmt = nan(amax,3,T); % males
    Nft = Nmt;         % females
    N = nan(amax,T); % total pop size
    X = nan(T,1); % sex ratio
    
% initial conditions (natural sex ratio)
Nmt(:,:,1) = [0*L(:) n_init.*L(:) 0.*L(:)];
Nft(:,:,1) = [n_init.*L(:) 0*L(:)  0.*L(:)];

%tally total pop size & sex ratio
Nmtt = sum(sum(Nmt(:,:,1)));
Nftt = sum(sum(Nft(:,:,1)));
N(:,1) = sum( Nmt(:,:,1) + Nft(:,:,1), 2);
X(1) = Nmtt./sum(N(:,1)); 


% Iterations
for t = 2:T
 
Nmt(:,:,t) = A*Nmt(:,:,t-1);
Nft(:,:,t) = A*Nft(:,:,t-1);

% Probability of a female finding a male, based on Rankin & Kokko (2007,
% Oikos). Multiply by 2 so it reaches 1 at sex ratio of 0.5.
p_male = min(1,betacdf(2*X(t-1),b_alpha,c1)); 

% Probability of female success (due to male harassment at high male ratio)
if X(t-1) == 1 % all males
    p_female = 0;
else
p_female = (1-X(t-1))./(1- X(t-1) + c2);
end

% Female fitness reduction due to EDC effects:
p_female = p_female .* f_fit_EDC;

% Get fecundity at age
Fec_tmp = Fec.*p_male.*p_female; % include probability of finding & not being harassed

% calculate reproductive outcomes.  
[Nmt(1,:,t), Nft(1,:,t)] = m_comb(Nmt(:,:,t-1),Nft(:,:,t-1),p,q,Fec_tmp,f_fit_xy,f_fit_yy,m_fit_xx,m_fit_yy);

% Bev-Holt spawning limitation
Rtot = sum(Nmt(1,:,t) + Nft(1,:,t));
Nmt(1,:,t) = BH(alphaBH,betaBH,Rtot).*Nmt(1,:,t);
Nft(1,:,t) = BH(alphaBH,betaBH,Rtot).*Nft(1,:,t);

%tally total pop size & sex ratio
Nmtt = sum(sum(Nmt(:,:,t)));
Nftt = sum(sum(Nft(:,:,t)));
N(:,t) = sum( Nmt(:,:,t) + Nft(:,:,t), 2);
X(t) = Nmtt./sum(N(:,t)); 
if isnan(X(t)); 
    X(t) = 0;   
end;

end % end t = 1:T

Nm(:,:,pp,qq) = squeeze(sum(Nmt,1))';
Nf(:,:,pp,qq) = squeeze(sum(Nft,1))';

end % end loop over P
end % end loop over Q

% Summary stats
Nm_end = squeeze(Nm(end,:,:,:));
Nf_end = squeeze(Nf(end,:,:,:));

% Phenotypic sex ratio:
Ph_rat = squeeze(sum(Nm_end,1)./(sum(Nf_end,1) + sum(Nm_end,1)));

% Genotypic sex ratio
Ge_rat = squeeze((sum(Nf_end([2 3],:,:),1) + sum(Nm_end([2 3],:,:),1))./(sum(Nf_end,1) + sum(Nm_end,1)));

% Population size (as a fraction of unfished)
Ntotal = squeeze((sum(Nm_end,1) + sum(Nf_end,1)))./Neq;
Ntotal_un = squeeze((sum(Nm_end,1) + sum(Nf_end,1)));
Ntotal_un = Ntotal_un(1,1);



