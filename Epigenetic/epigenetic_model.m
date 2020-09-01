function [Nf_total, Nm_total, dx] =epigenetic_model(exposure_type,epiF1,epiF2,paternal_meth_only,chemical,RandSim,RangeSim,Index)

% Model of epigenetic population dynamics
% For EPA STAR project led by SMB
% Parameterization based on SMB & DeCourten experiments w silversides

% Notes for paper - could do this using IBM, but...fish are so fecund with
% very high initial mortality that actually doing those calculations is even
% more impractical than what we are doing here

% OK Only track paternal metholome, possibly. Make this a switching option
% Add temp-dependent mortality to reproduce advantage of female size

if ~exist('epiF1','var');epiF1 = true;end
if ~exist('epiF2','var');epiF2 = true;end
if ~exist('paternal_meth_only','var');paternal_meth_only = false;end
if ~exist('chemical','var'); chemical = 'EE2';end

% Move some of this to a params file?
T = 16*6; % length of simulation
Amax = 12; % maximum age of any individual
Age = 1:Amax;
MatF = [50 1]; % mean & SD of maturity ogive 
MatM = [60 1]; % mean & SD of maturity ogive
Sex = {'f','m'};
% [Linf k CV_L M ismat repro_c repro_e mean_r_size sd_r_size]
fixparmM = [130, 0.24, 0.1, 0.35, 20, 1e-3, 3, 9.5, 1];
fixparmF = [130, 0.24, 0.1, 0.35, 20, 1e-3, 3, 9.5, 1];


% Bev-Holt
alfa = 0.31;
alfa = 0.03;
beta = 1;

% IPM parameters:
max_x = fixparmF(1)*2;
n_grid = 50;
x = linspace(0,max_x,n_grid);
dx = diff(x(1:2));
Sy = makeSimpVec(dx,n_grid); % for integration
Symat = repmat(Sy(:)',[length(Sy),1]);
Symat3 = repmat(Sy(:),[1,Amax,Amax]);

% Maturity ogives
isMatF = normcdf(x,MatF(1),MatF(2));
isMatM = normcdf(x,MatM(1),MatM(2));


% Baseline:
Exposure = zeros(1,T); % exposure concentration in each time step

switch exposure_type
    
case 'Baseline'
% nothing happens
    
% Winter acute exposure
case 'Single_winter'
Exposure(71:73) = 1;

case 'Single_summer'
% Summer acute exposure
Exposure(74:76) = 1;

case 'Chronic'
% Chronic exposure:
Exposure(71:end) = 1;
end % end switch exposure type


Temp = 25*ones(1,T); % water temperature in each time step

% Temp data: From CDMO (SF Bay NERR), 1st Mallard Station. Approximate
% range is 8 to 23 ºC annually.
% But currently the model has a breakpoint temp of 25º so vary accordingly
Temp = (sin(((1:T)-1)/6*2*pi)/2+0.5)*(30-20)+20;

% Have to initialize first T0 years to build up backlog of epigenetic
% history
T0 = Amax*3; % both parents and grandparents

% State variable: N
% Indexed by size, age, sex, parental age, grandparental age, time
% Initialize at stable size & age distribution, then apportion parental
% ages.
% Move this to a separate function: 
K_gm_f = kernmatSimp_epi(x,fixparmF,0,0,0,0,1,1,epiF1,epiF2,paternal_meth_only,'growth_mort',chemical,NaN,NaN,NaN,RandSim,RangeSim,Index);
K_gm_m = kernmatSimp_epi(x,fixparmM,0,0,0,0,1,1,epiF1,epiF2,paternal_meth_only,'growth_mort',chemical,NaN,NaN,NaN,RandSim,RangeSim,Index);
K_r = kernmatSimp_epi(x,fixparmF,0,0,0,0,1,1,epiF1,epiF2,paternal_meth_only,'repro',chemical,15,0.5,isMatF,RandSim,RangeSim,Index);

N_init = zeros(n_grid,2,T0);
% recruit distribution: 

% Get eigs to adjust alpha
Eig = max(eig(K_r+K_gm_f));
%alfa = 1/(alfa*Eig);
alfa = alfa/Eig;


R0 = normpdf(x,fixparmF(8),fixparmF(9));
N_init(:,:,1) = repmat(R0(:),[1,2]);

%  integrate
K_gm_f = Symat.*K_gm_f;
K_gm_m = Symat.*K_gm_m;
K_r = Symat.*K_r;

for t = 2:100
    N_init(:,1,t) = K_gm_f*N_init(:,1,t-1);
    N_init(:,2,t) = K_gm_m*N_init(:,2,t-1);
    R = K_r*N_init(:,1,t-1);
    DD_surv = alfa/(1+alfa/beta*sum(R));
    R = R*DD_surv;
    N_init(:,1,t) = N_init(:,1,t) + 0.5*R;
    N_init(:,2,t) = N_init(:,2,t) + 0.5*R;
end % end initial time loop


Nf0 = N_init(:,1,T0)/100;
Nm0 = N_init(:,2,T0)/100;

% apportion initial size distribution & ages
N = zeros(n_grid,Amax,2,Amax,Amax,T);
N(:,:,1,:,:,1:T0) = repmat(Nf0,[1,Amax,1,Amax,Amax,T0])./(Amax.^3);
N(:,:,2,:,:,1:T0) = repmat(Nm0,[1,Amax,1,Amax,Amax,T0])./(Amax.^3);

SR = 0.5*ones(1,T);

for t = (T0+1):T
    
    if mod(t,20)==0; disp(t); end
    % calculate the current sex ratio
    % move to a separate function. 
    MatF = squeeze(N(:,:,1,:,:,t-1));
    MatM = squeeze(N(:,:,2,:,:,t-1));
    MatF = Sy(:)'*(isMatF.*sum(sum(sum(MatF,4),3),2));
    MatM = Sy(:)'*(isMatM.*sum(sum(sum(MatM,4),3),2));
    SR(t) = MatM/(MatM+MatF);
    
    if SR(t) == 0; keyboard; end
    
    % calculate the current distribution of maternal ages
    F_agesum = Sy(:)'*sum(sum(squeeze(N(:,:,1,:,:,t-1)),4),3); % sum of numbers in each age
    F_agedist = F_agesum/sum(F_agesum);
    % calculate the current distribution of grandparental ages (mother's
    % side)
    % More complicated bc have to account for current age to count
    % backwards
    % This iss
    GF_ages = squeeze(sum(squeeze(sum(N(:,:,:,:,:,t-1),5)),3));
    GF_agemat = Age(:)+Age(:)'; % matrix of possible GF ages + F ages
    GF_agesum = squeeze(sum(Symat3.*GF_ages));
    GF_agedist = GF_agesum./sum(GF_agesum(:));
    GF_ageprop = GF_agedist.*GF_agemat;
    meanGFage = sum(GF_ageprop(:));
    sdGFage = sqrt(mean( (GF_ageprop(:) - meanGFage).^2 ));
    GF_agevec = 1:(Amax*2); % MAYBE THIS SHOULD BE Amax*3???
    GF_agedist = normcdf(GF_agevec+0.5,meanGFage,sdGFage)-normcdf(GF_agevec-0.5,meanGFage,sdGFage);
    % TO DO: 
    % Also back-calculate paternal-side grandfather age
    % distribution). This will be based on F_agedist_m, weighted by current
    % age distribution.
    % WEIGHT THESE DISTRIBUTIONS BY CURRENT AGE DISTRIBUTION 
    
    %Track parental & grandparental (father) age distributions
    F_agedist_m(:,t) = F_agedist;
    GF_agedist_m(:,t) = GF_agedist;
    
  
    
    % separate kernel for each age in the population:
    for a = 1:(Amax-1)
        
        % separate kernel for each parental age
        for pa = 1:Amax
        
        % separate kernel for each grandparental age
        for gpa = 1:Amax
        
            % generate growth & mortality kernel as function of sex, age, pa, gpa, & current
            % exposure. Also now includes paternal age distribution
            Exp_agevec = fliplr(Exposure( t-a-Amax*2+1 : t-a)); % fliplr so that age counts backwards properly
            %if t == 80; keyboard; end
           % keyboard
            K_gm_f = kernmatSimp_epi(x,fixparmF,Exposure(t),Exposure(t-pa),Exposure(t-(pa+gpa)),...
                                     Exp_agevec,F_agedist_m(:,t-a),GF_agedist_m(:,t-a),epiF1,epiF2,...
                                     paternal_meth_only,'growth_mort',chemical,NaN,NaN,NaN,RandSim,RangeSim,Index);
            K_gm_m = kernmatSimp_epi(x,fixparmM,Exposure(t),Exposure(t-pa),Exposure(t-(pa+gpa)),...
                                     Exp_agevec,F_agedist_m(:,t-a),GF_agedist_m(:,t-a),epiF1,epiF2,...
                                     paternal_meth_only,'growth_mort',chemical,NaN,NaN,NaN,RandSim,RangeSim,Index);

            % advance age & growth for this age, pa, gpa
            N(:,a+1,1,pa,gpa,t)= (Symat.*K_gm_f)*N(:,a,1,pa,gpa,t-1);
            N(:,a+1,2,pa,gpa,t)= (Symat.*K_gm_m)*N(:,a,2,pa,gpa,t-1);

            
            % generate reproductive kernel (incl. penalty on egg production or survival) based on pa, gpa,
            % current exposure and parental & grandparental age distribution. 
            [K_r,SR_offspring] = kernmatSimp_epi(x,fixparmF,Exposure(t),Exposure(t-pa),Exposure(t-gpa),...
                                                 Exp_agevec,F_agedist,GF_agedist,epiF1,epiF2,...
                                                 paternal_meth_only,'repro',chemical,Temp(t-1),SR(t-1),isMatF,RandSim,RangeSim,Index);
            
            % Add new recruits to the appropriate age, pa, gpa vector
            R = (Symat.*K_r)*N(:,a,1,pa,gpa,t-1);
            Rect_tmp_tmp(:,1,a,pa,gpa) = R*(1-SR_offspring); % female recruits
            Rect_tmp_tmp(:,2,a,pa,gpa) = R*SR_offspring; % male recruits
            
            
        end % end loop over grandparental age
        
        Rect_tmp = sum(Rect_tmp_tmp,5);
        
        end % end loop over parental age
    end % end loop over age
    
    
    %--------------------------------------
    % TO DO: 
    % reweight recruit parental age distribution to account for paternal age?
    
    %--------------------------
    % Apply density dependence:
    % Total recruits in all classes:
    Total_recruits = sum(sum(sum(Rect_tmp,4),3),2);
    % Integrate to get numbers
    Total_recruits = Sy(:)'*Total_recruits(:);
    % Apply BH
    BH_surv = alfa/(1+alfa/beta*Total_recruits);
    %disp(Total_recruits)
   % if Total_recruits<1; keyboard; end
    %keyboard
    % Now apply survival to recruits:
    Rect_DD = BH_surv*Rect_tmp;
    Rect_DD(:,:,end+1,:) = 0; % nobody has the oldest-age fathers
    % Now add back in to the first age class:
    N(:,1,1,:,:,t) = reshape(Rect_DD(:,1,:,:),[length(x),1,1,Amax,Amax]); % females
    N(:,1,2,:,:,t) = reshape(Rect_DD(:,2,:,:),[length(x),1,1,Amax,Amax]); % males
    
    
 %   Sy(:)'*sum(sum(sum(Rect_DD,4),3),2)
    
   % if t == 58; keyboard; end
   % keyboard
    %--------------------------
    
    
           % if Exposure(t); keyboard; end
    
end % end loop over time

Nf_total = squeeze(N(:,:,1,:,:,:));
Nm_total = squeeze(N(:,:,2,:,:,:));
Nf_total = squeeze(sum(sum(sum(Nf_total,4),3),2));
Nm_total = squeeze(sum(sum(sum(Nm_total,4),3),2));

%keyboard




