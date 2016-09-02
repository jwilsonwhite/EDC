function EDC_runs

% Make assorted runs of the age-structured EDC model
% These model runs designed for White, Cole, Cherr, Connon, & Brander
% manuscript (Ecol. Appl.)
% Will White - July 2016 whitejw@uncw.edu

% First set of runs: generic model with multiple alternative mating
% functions
% Second set of runs: Menidia-specific model
do1 = false;
do2 = true;
doRuns = false; % Need to actually do the Menidia model runs, or are they already complete?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do1
% First set of results: baseline with several different mating functions,
% no overall reduction in fecundity due to EDC, but 2 different values of
% the penalty on sex-changed individuals
P = 0:0.01:0.99; % feminization rate
Q = 0:0.01:0.99; % masculinization rate

% Parameters of the mating function:
Mf = [40, 1, 15]; % polygyny, monogamy, minimum-male-threshold. % mate-finding probability.  > 20 = group spawning, < 10 ~ linear increase
Mf2 = [1, 1, 6];
Ff = [1, 0.5]; % fertility reduction due to sex change
ColM = {'r','b','k'};
LsF = {'-','--'};
alphaBH = 0.37; % this is 1/slope of BH function.
                % values used: 0.14 (baseline), 0.05, 0.37 (range from
                % Myers et al 1999); 0.31 (for Menidia, based on clupeids)

% Setup figures:
figure(1)
clf
set(gcf,'units','cent','position',[10,40,17,18])

figure(2)
clf
set(gcf,'units','cent','position',[30,40,8,14])

SP = [1,2;3,4;5,6]; % subplot indices

for m = 1:length(Mf)
    for f = 1:length(Ff)
        
% setup parameters
EDC_params(Ff(f),Ff(f),alphaBH)

% mating function (Beta CDF parameters)
c1 = Mf(m);
c2 = Mf2(m);

% fitness reduction due to EDCs.  (not used in these runs)
f_fit_reduction = 1; 

[Ntotal, SexRat] = EDC_popmodel_as(c1,c2, P, Q,f_fit_reduction,1);

%-------------------------------------------------------
% contourplot
figure(1)
subplot(3,2,SP(m,f))
colormap(1-gray)
contourf(Ntotal)
xlabel('Probability of masculinization');
ylabel('Probability of feminization');
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',[1:20:length(P) 100],'xticklabel',[P(1:20:end),1])
set(gca,'ytick',[1:20:length(P) 100],'yticklabel',[P(1:20:end),1])
caxis([0 1])
ch = colorbar;
set(ch,'tickdir','out','ticklength',[0.02])
ylabel(ch,'Population size (proportion of Nmax)')
axis equal
%-------------------------------------------------------

%-------------------------------------------------------
% Line plot
figure(2)
 
% for increasing masculinization
subplot(2,1,1)
hold on
    plot(Q,Ntotal(1,:),'color',ColM{m},'linestyle',LsF{f})
    xlabel('Probability of masculinization')
    ylabel('Population size (proportion of Nmax')
    ylim([-0.05 1.1])
    set(gca,'tickdir','out','ticklength',[0.02 0.02])
    
    % for increasing feminization
subplot(2,1,2)
    hold on
    plot(P,Ntotal(:,1),'color',ColM{m},'linestyle',LsF{f})
    xlabel('Probability of feminization')
    ylabel('Population size (proportion of Nmax')
    ylim([-0.05 1.1])
    set(gca,'tickdir','out','ticklength',[0.02 0.02])

    end
end
% End section 1
end % end if do1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 2: Menidia-specific parameters & mating function
if do2

P = 0:0.01:0.95; % feminization rate
Q = 0:0.01:0.95; % masculinization rate

c1 = 1.002;  % estimated from spawning experiment
c2 = 1.2823; % estimated from spawning experiment

% Fertility reduction due to sex change.
Ff = 0.5; % fertility reduction due to sex change (based on ~50% reduction in milt production by intersex Rutilus, Jobling et al. 2002)

% Site-specific data [Georgiana, Napa, Sacramento, Suisun, Denverton]
% See data in Summary_data_from_BCole.xlsx
MeanSR = [0.39, 0.26, 0.16, 0.69, 0.65];
LowerSR = [0.23, 0.12, 0.07, 0.54, 0.51]; % lower 95%CI
UpperSR = [0.62, 0.23, 0.39, 1.24, 1.16];

alphaBH = 0.31;

% setup parameters
EDC_params(Ff,Ff,alphaBH) 

% fitness reduction due to EDCs.
% Based on choriogenin expression in Brander et al.(2016) Aquatic Tox,
% we see the following pattern due to bifenthrin [0, 0.5, 5, 50] (ng/L):
Bif_reduction = [10^-4.1, 10^-4.3, 10^-4.65, 10^-4.8];
f_fit_reduction = Bif_reduction/Bif_reduction(1);

% The first value (63%) corresponds with effects seen in spawning trials with 0.5 ng/L bifenthrin (egg production is 65% of normal)

% Length reduction due to EDC
% Based on difference between urban & ranch sites in Brander et al. (2013)
% PLoS ONE
L_reduction = 0.91; % 9% reduction

% Setup figures:
figure(3)
clf
set(gcf,'units','cent','position',[10,40,16,24])
subplots = reshape(1:8,[2,4])';

figure(4)
clf
set(gcf,'units','cent','position',[10,40,16,24])

figure(5)
clf
set(gcf,'units','cent','position',[10,40,15,17])
Col = [0.1 0.1 0.9;...
       0.6 0.1 0.8;...
       0.9 0.1 0.1;...
       0.5 0.05 0.05]; % colormap for plotting


lnames = {'','lred'};
lreds = [1,L_reduction];
for f = 1:length(f_fit_reduction)
    for l = 1:2 % length reduction
fname = strcat('EDC_Menidia_runs_f_fit_reduct',num2str(f),'_',lnames{l});
if doRuns   
[Ntotal, SexRat] = EDC_popmodel_as(c1,c2, P, Q,f_fit_reduction(f),lreds(l));
save(fname)
end
%else % load previous runs
load(fname)


figure(3)
subplot(4,2,subplots(f,l))
colormap(1-gray)
contourf(Ntotal,0:0.1:1)
xlabel('Probability of masculinization');
ylabel('Probability of feminization');
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',[1:20:length(P) 100],'xticklabel',[P(1:20:end),1])
set(gca,'ytick',[1:20:length(P) 100],'yticklabel',[P(1:20:end),1])
caxis([0 1])
ch = colorbar;
set(ch,'tickdir','out','ticklength',[0.015])
ylabel(ch,'Population size (proportion of Nmax)')

figure(4)
subplot(4,2,subplots(f,l))
colormap(1-gray)
contourf(SexRat)
xlabel('Probability of masculinization');
ylabel('Probability of feminization');
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',[1:20:length(P) 100],'xticklabel',[P(1:20:end),1])
set(gca,'ytick',[1:20:length(P) 100],'yticklabel',[P(1:20:end),1])
caxis([0 1])
ch = colorbar;
set(ch,'tickdir','out','ticklength',[0.015])
ylabel(ch,'Sex ratio (proportion male)')


% For figure 5, make a 6-panel figure, one for each site + 1 schematic
figure(5)
Subplots = [2,4,6,3,5]; % put Cole sites on the right, SB sites on the left 

for i = 1:5
    subplot(3,2,Subplots(i))
hold on

% Create a mesh
Smesh = 0:0.01:1;
Sv = SexRat(:);
Nv = Ntotal(:);
Zmax = zeros(size(Smesh));
Zmin = Zmax;
for s = 1:(length(Smesh)-1);
    OK = Sv >= Smesh(s) & Sv < Smesh(s+1);
    if any(OK)
    Zmax(s) = max(Nv(OK));
    Zmin(s) = min(Nv(OK));
    end
end
%contourf(Smesh,Nmesh,Z>0);
plot(Smesh,Zmax,'color',Col(f,:))
plot(Smesh,Zmin,'color',Col(f,:))

% Now plot the intersection of the appropriate site with the curves
do_intersection(MeanSR(i),Smesh,Zmin,'k');
do_intersection(MeanSR(i),Smesh,Zmax,'k');

ylabel('Equilibrium population size')
xlabel('Sex ratio (proportion male')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
end % end loop over sites
end % end loop over length reduction
end % end loop over f
%-------------------------------------------------------

end % end if do2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function do_intersection(X,Mesh,Y,Col)

% Plot the intersection lines
Diff = abs(Mesh-X);
OK = find(Diff == min(Diff),1);

% Plot vertical line:
plot([Mesh(OK) Mesh(OK)],[0,Y(OK)],'k-','color',Col);
% Plot horizontal lines:
plot([0 Mesh(OK)],[Y(OK),Y(OK)],'k-','color',Col);
