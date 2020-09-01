function [kmat, SR_offspring] = kernmatSimp_epi(x,fixparm,exposure,exposure_p,exposure_gp,exposure_vec,F_dist,GF_dist,epiF1,epiF2,paternal_meth_only,type,chemical,temp,SR,isMatF,RandSim,RangeSim,Index,T)

% Set up the integration mesh kernel using Simpsons
% For epigenetic EDC model

%adapted from Will from Easterling. Evenly spaced grid, now add weights so
%use Simpson's rule using Marissa's code.

if ~exist('T','var'); T = 1; end
if ~exist('temp','var'); temp = NaN; end
if ~exist('SR','var'); SR = 0.5; end
if ~exist('isMatF','var'); isMatF = 1; end

% turn of epigenetics if desired
if ~epiF1; exposure_p = exposure_p*0; end
if ~epiF2; exposure_gp = exposure_gp*0; exposure_vec = exposure_vec*0;  end
if paternal_meth_only;  exposure_vec = exposure_vec*0;  end
y = x;
%this creates a vector (y) that is equal to x

[x,y] = meshgrid(x,y); % Matlab built in function
%x is an array (original x by original x) with each row equal to original
%vector x
%y is an array (original y by original y) with each column equal to
%original vector y
%X corresponds to size at time t
%Y corresponds to size at time t+1

% Modify fixparm to account for exposure effects:

% Contents of fixparm are
% [Linf k CV_L M ismat repro_c repro_e mean_r_size sd_r_size]

% Make the kernel from the grid
kmat = mkkern_epi(x,y,fixparm,exposure,exposure_p,exposure_gp,exposure_vec,...
    F_dist,GF_dist,type,chemical,temp,SR,isMatF,RandSim,RangeSim,Index,T);

kmat = max(0,kmat); 

% Make sex ratio temperature, exposure, and sex-ratio dependent
% Assume probability distribution with 50:50 at 25ºC, as in DeCourten %
% Brander (2017)
if exposure==0
SR_offspring = normcdf(temp,25,8); % this approximately matches dC&B result
else
    switch chemical
        case 'EE2'
% assuming EE2 exposure like in dC&B
SR_offspring = 1-normcdf(temp,20,8); % (this does not match exactly)

        case 'Bif'
          SR_offspring = normcdf(temp,23,8);  % (this does not match exactly)
        otherwise
        SR_offspring = normcdf(temp,25,8);    
          
    end % end switch
    
    end


