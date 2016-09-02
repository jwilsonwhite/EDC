function spawning_fr

% Generate mating function from Menidia spawning trials

% For White, Cole, Cherr, Connon, & Brander
% manuscript (Ecol. Appl.)
% Will White - whitejw@uncw.edu
% April 2013

Trials = [1  3   5:7 8 9:10];
Tanks = 3:4;


% Read in fish IDs
ID = importdata('Data_Oct2013/FishIDs_F2013.csv');
Code = ID.textdata(2:end,2);
Sex = ID.textdata(2:end,1);
ID.data(ID.data==-9999) = NaN; % code for missing data
SL = ID.data(:,1); % add NaN b/c last line of data is lopped off
Mass = ID.data(:,2);
Gonad = ID.data(:,3);

% Read in trial info
TR = importdata('Data_Oct2013/TrialIDs_F2013.csv');
TR_ID = TR.textdata(2:end,1);
TRnum = TR.data(:,1);
TRtank = TR.data(:,2);

% Read in Egg data
E = importdata('Data_Oct2013/Eggs_F2013.csv');
E = E.data; % % columns are trial, tank, #fert eggs, #unfert eggs, day of trial

D = [];
% Loop over trials
for tr = 1:length(Trials)
    Tr = Trials(tr);
    
% Loop over tanks
for ta = 1:length(Tanks)
    Ta = Tanks(ta);
    
    % Get Fish IDs in trial
    TR_IDtmp = TR_ID(TRnum==Tr & TRtank == Ta);
    
    % Get sex of those fish
    Sex_tmp = cell(length(TR_IDtmp),1);
    SL_tmp = nan(length(TR_IDtmp),1);
    Mass_tmp = nan(length(TR_IDtmp),1);
    Gonad_tmp = nan(length(TR_IDtmp),1);
    
    for tt = 1:length(TR_IDtmp)
        
        Ind = strcmp(TR_IDtmp{tt},Code);
        if any(Ind)
            ID_tmp = find(Ind,1);
            Sex_tmp{tt} = Sex{ID_tmp};
            SL_tmp(tt) = SL(ID_tmp);
            Mass_tmp(tt) = Mass(ID_tmp);
            Gonad_tmp(tt) = Gonad(ID_tmp);
        else
            Sex_tmp{tt} = 'NaN';
        end
            
    end % end tt
    
    
    M = sum(strcmp('m',Sex_tmp));
    F = sum(strcmp('f',Sex_tmp));
    Ratio = M./(M+F);
    Unk = sum(strcmp('NaN',Sex_tmp));
    
    Mass_F = nansum(Mass_tmp(strcmp('f',Sex_tmp)));
    Mass_Fmean = nanmean(Mass_tmp(strcmp('f',Sex_tmp)));
    Mass_M = nansum(Mass_tmp(strcmp('m',Sex_tmp)));
    Mass_ratio = Mass_M./(Mass_M+Mass_F);
    Gonad_F = nansum(Gonad_tmp(strcmp('f',Sex_tmp)));
    Gonad_Fmean = nanmean(Gonad_tmp(strcmp('f',Sex_tmp)));
    
    E_ok = E(:,1) == Tr & E(:,2) == Ta; % & E(:,5) == 1; % only take 1st day of trial
    if sum(E_ok)>=1
    E_tmp = nansum(E(E_ok,3)); % Total fertilized eggs produced
    else E_tmp = NaN;
    end

    % Scale by mass, assuming missing fish had a 50:50 sex ratio
    E_tmp = E_tmp./(Gonad_F + Gonad_Fmean.*Unk.*(1-0.5));  
  
    Fe =  (E(E_ok,3)); % fertilized eggs
    Ufe = (E(E_ok,4)); % unfertilized eggs
    F_tmp = nanmean(Fe./(Fe+Ufe)); % proportion fertilized
    
    D = [D; Tr, Ta, Ratio, Mass_ratio, E_tmp, F_tmp, Unk, Mass_F, Gonad_F];
    
    
    
end % end loop over tanks
end % end loop over trials


D = D(D(:,7)<5,:);
badrow = D(:,1) ==8 & D(:,2) ==3;
D = D(~badrow,:);

% Egg production vs. sex ratio
Dtmp = D(D(:,5)<20,:);
Dtmp = D;
OKrows = true(size(Dtmp,1),1);
OKrows(6)=false;
%Dtmp = Dtmp(OKrows,:);
[Bt, ~, Rt, ~, Stats] = regress(Dtmp(:,5),[ones(size(Dtmp,1),1),Dtmp(:,4)]); % Eggs vs. biomass sex ratio

% Make a plot (not used for publication)
figure(1)
clf
sh(1) = subplot(2,2,1);
hold on
plot(D(:,3),D(:,5),'ko','markersize',8)
xlabel('Sex ratio (proportion male)')
ylabel('Egg production (g female^-^1 d^-^1')

sh(2) = subplot(2,2,3);
hold on
plot(D(:,4),D(:,5),'ko','markersize',8)
xlabel('Biomass sex ratio (proportion male)')
ylabel('Egg production (g female^-^1 d^-^1')

sh(3) = subplot(2,2,2);
hold on
plot(D(:,3),D(:,6),'ko','markersize',8)
xlabel('Sex ratio (proportion male)')
ylabel('Fertilization rate')

sh(4) = subplot(2,2,4);
hold on
plot(D(:,4),D(:,6),'ko','markersize',8)
xlabel('Biomass sex ratio (proportion male)')
ylabel('Fertilization rate')


set(sh(:),'tickdir','out','ticklength',[0.02 0.02])
set(sh(:),'xtick',0:0.1:1)

xlabel('Biomass sex ratio (proportion male)')
ylabel('Fertilization rate')

csvwrite('fr_results.csv',Dtmp)

% Fit Beta cdf using maximum likelihood:
% Procedure for 2-D likelihood surface & profile confidence intervals
% See Bolker's book, Ch 6 for details & justification

% Profile likelihood:
A = linspace(1e-10,7,1e3);
B = linspace(1e-10,7,1e3);
AA = repmat(A(:),[1,length(B)]);
BB = repmat(B(:)',[length(A),1]);

% fit a Beta cdf function to each combination of a & b
for i = 1:length(D)
L(:,:,i) =  betacdf(D(i,3),AA,BB)-D(i,5)./max(D(:,5)) ;
end

Sig = std(L,[],3); % estimate standard deviation from the data
N = normpdf(L,zeros(size(L)),repmat(Sig,[1,1,size(L,3)]));
LL = -1*sum( log(N), 3);  % negative log-likelihood

minLL = min(min(LL));

MLE = [AA(LL==minLL), BB(LL==minLL)];

Conf = chi2inv(0.9,1)/2;  % LRT bivariate confidence interval value

% Profile CIs for each parameter
Aprof = min(LL); % minimum NLL for each value of B
Aprof_red = A(Aprof<(minLL+Conf));
Aprofile = [Aprof_red(1), Aprof_red(end)];

Bprof = min(LL'); % minimum NLL for each value of B
Bprof_red = B(Bprof<(minLL+Conf));
Bprofile = [Bprof_red(1), Bprof_red(end)];

% Likelihoods that fall within the confidence region:
LL_ok = LL <= minLL + Conf;
Aconf = AA(LL_ok);
Bconf = BB(LL_ok);
Aconf = Aconf(:);
Bconf = Bconf(:);

% Plot likelihood surface
figure(2)
set(gcf,'units','cent','position',[20 20 12 10])
clf
hold on
[~,h]=contourf(B, A, LL-minLL, 0:10);
[~,h2]=contour(B,A,LL-minLL,Conf);
set(h,'linestyle','none')
set(h2,'linecolor','k')
set(gca,'tickdir','out')
xlabel('b parameter')
ylabel('a parameter')
colormap(repmat((1:-0.01:0.2)',[1,3]))
cb = colorbar;
set(cb,'ytick',0:20,'tickdir','out')
ylabel(cb,'Delta negative log-likelihood')
set(gca,'clim',[0 12])
plot(MLE(2),MLE(1),'ko')



% Now here's a figure that just has the sex ratio business
figure(3)
set(gcf,'units','cent','position',[20 20 8 6])
clf
hold on
plot(D(:,3),D(:,5),'ko','markersize',8)
set(gca,'TickDir','out','Ticklength',[0.015 0.015])
set(gca,'xtick',0.2:0.1:0.6,'ytick',0:100:1000,'fontsize',10)
set(gca,'ylim',[50 450],'xlim',[0.2 0.6])
xlabel('Sex ratio (proportion male)','fontsize',12)
ylabel('Egg production (g female^-^1 d^-^1','fontsize',12)

Xd = linspace(0.25,0.6,1e2);
Bm = betacdf(Xd,MLE(1),MLE(2)).*max(D(:,5));
plot(Xd,Bm,'k')

% Confidence intervals:
for i = 1:length(Aconf)
    Bc(i,:) = betacdf(Xd,Aconf(i),Bconf(i)).*max(D(:,5));
end

Bl = quantile(Bc,0.025);
Bu = quantile(Bc,0.975);


%Bl = betacdf(Xd,Aprofile(1),Bprofile(1)).*max(D(:,5));
plot(Xd,Bl,'k:')

%Bu = betacdf(Xd,Aprofile(2),Bprofile(2)).*max(D(:,5));
plot(Xd,Bu,'k:')




