function kxy = mkkern(x,y,F,fixparm,T)

% adapted from Will's code adapted from Easterling's code, used for PISCO
% rockfish data set

%DEFINE WHO CAN GET FISHED
%isjuv = 1 - normcdf(x,fixparm(4),diff(x(1,1:2))/2); 
isjuv = 1 - normcdf(x,fixparm(5),fixparm(7)); 

%define pr(reproductive) using a maturity ogive
%not needed now because open population
%ismat = 1 - normcdf(x,fixparm(5),diff(x(1,1:2))/2);

%SURVIVAL PART OF KERNEL
%add in stochasticity for process error, but keep constant across size 
%Assume no stochasticity in fishing
%mortality.
M = fixparm(4); %normrnd(fixparm(3), varparm(1)); %CHECK VALUE FOR STD
m = ones(size(x)).*M + (1-isjuv).*F; %this is a matrix size x, 
                                              %mortality for each size 
                                              
p1 = exp(-m*T); % convert mortality rate to survivorship, iterate over time steps

%GROWTH PART OF KERNEL
%add variability in growth to k
Linf = fixparm(1);
k = fixparm(2);
%k = max(realmin,normrnd(fixparm(2),varparm(2)));  
%growth
pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
%add variability around von Bertalanffy growth
psig1 = pmean1*fixparm(7); % make this a parameter; %.*varparm(3); 

%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

kxy = p1.*p2;
%add fecundity to kernel if a closed population

