function kxy = mkkern_epi(x,y,fixparm,exposure,exposure_p,exposure_gp,exposure_vec,F_dist,GF_dist,type,chemical,temp,SR,isMatF,RandSim,RangeSim,Index,T)

% adapted from Will's code adapted from Easterling's code, 
% used for epigenetic model
Coef = epigenetic_coefs;

switch type
    case 'growth_mort'
%SURVIVAL PART OF KERNEL
%add in stochasticity for process error, but keep constant across size 
%Assume no stochasticity in fishing
%mortality.
M = fixparm(4); %normrnd(fixparm(3), varparm(1)); %CHECK VALUE FOR STD


% Add deformity mortality. Assume deformities represent a proportional
% increase in mortality. 
%if exposure_p
P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;

c = Coef(1).DeformF1(1).(chemical).c;
if RandSim
Exp = normrnd(c(1),c(2)*c(1)); 
elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
else Exp = c(1); 
end

M = M*(1-P_exp) + M*P_exp*Exp;
%end % end if exposure


% Add in reduced deformity mortality in F0
% [currently turned off bc WTF would this happen]
if exposure
    switch chemical
    case 'EE2'
    case 'Bif'
   %     M = M*0.5;
    case 'TB' 
   %    M = M*0.5;
    case 'Levo'
    end % end switch
end % end if exposure

% Add in deformity mortality observed in DeCourten & Brander 2017:
%if exposure
c = Coef(1).DeformF0_dCB(1).(chemical).c;
if RandSim
Exp = normrnd(c(1),c(2)*c(1)); 
elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
else
    Exp = c(1);
end        
       % M = M*Exp;
        M = M*(1-exposure) + M*Exp*exposure;
%end % end if exposure

% Add in deformity mortality observed in F1 in dC&B 2017: 
c = Coef(1).DeformF1_dCB(1).(chemical).c;
if RandSim
Exp = normrnd(c(1),c(2)*c(1)); 
elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
else Exp = c(1);
end % emd of RandSim  
    switch chemical
    case 'Bif'
        if temp < 25 % THIS IS NOT QUITE RIGHT BC IT SHOULD BE DEPENDENT ON TIME AT BIRTH, OI!!!    
        P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;
        else
            P_exp = 0;
        end
    case 'EE2'

        if temp > 25 % THIS IS NOT QUITE RIGHT BC OF sIGNIF in dC&B 2017
        P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;
        else
            P_exp = 0;
        end
        
        otherwise
            P_exp=0;
    end % end switch
    M = M*(1-P_exp) + M*P_exp*Exp;  


    % Add in deformity mortality in F2:
    c = Coef(1).DeformF2(1).(chemical).c;
if RandSim
Exp = normrnd(c(1),c(2)*c(1)); 
elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
else Exp = c(1);
end % end if RandSim

    % Get combined grandparental exposure
        GP_exp = (exposure_gp + sum(exposure_vec(:).*GF_dist(:)))*0.5;
        M = M * (1-GP_exp) + M * (Exp * GP_exp);

        

m = ones(size(x)).*M; %this is a matrix size x, 
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
psig1 = pmean1*fixparm(3); % make this a parameter; %.*varparm(3); 

%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

kxy = p1.*p2;
%add fecundity to kernel if a closed population

    case 'repro'
                
        %ismat = normcdf(x,fixparm(5),diff(x(1,1:2))/2);
        % Get combined grandparental exposure
        GP_exp = (exposure_gp + sum(exposure_vec(:).*GF_dist(:)))*0.5;
        
        % LENGTH EFFECTS:
        % Fish were ~ 1 mm longer at hatch when exposed to EE2
        %if exposure
            
            c = Coef(1).HatchSize(1).(chemical).c + GP_exp*Coef(1).HatchSizeF2(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            
            elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
            else
                Exp = c(1);
            end % end if RandSim 
            
            R_len_fact = Exp*exposure;
            
            c = Coef(1).HatchSizeF2(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            
            elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
            else
                Exp = c(1);
            end % end if RandSim 
            
            R_len_fact = R_len_fact+ Exp*GP_exp;
            
            
        %end % end if exposure
                
        
        Rvec = normpdf(y,fixparm(8)+R_len_fact,fixparm(9));
   
        Fec = fixparm(6)*(isMatF.*x).^fixparm(7); 
        
        
        % ATRESIA EFFECTS + TOTAL HATCHING + Egg production
        % penalty for exposure. baseline = 20% incidence. 
        %if exposure
            
            % Atresia
            c = Coef(1).Atresia(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
                if Index == 1; Exp = c(1);
                elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
                elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
                end
            else
                Exp = c(1);
            end % end if RandSim
             
           
           Fec = Fec*Exp*exposure + Fec*(1-exposure);
            
           
           % Egg hatching
            c = Coef(1).HatchSuccessF0(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
                if Index == 1; Exp = c(1);
                elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
                elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
                end
            else Exp = c(1);
            end
              
        Fec = Fec*Exp*exposure + Fec*(1-exposure);
        
        % Egg production
            c = Coef(1).EggProdF0(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
                if Index == 1; Exp = c(1);
                elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
                elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
                end
            else Exp = c(1);
            end
              
        Fec = Fec*Exp*exposure + Fec*(1-exposure);
        
        
        %end
        
        % penalty/bonus for HATCHING SUCCESS. 
       % if exposure
            
            c = Coef(1).HatchSuccessF1(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
            else Exp = c(1);
            end 
            
            Fec = Fec*Exp*exposure_p + Fec*Exp*(1-exposure_p);
        %end
        
        % F1 bonus for hatching success: (this is wrong bc effects would be
        % felt by F0 as they are producing F1
      %  P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;
        
        
%            c = Coef(1).HatchSuccessF1(1).(chemical).c;
%            if RandSim
%            Exp = normrnd(c(1),c(2)*c(1)); 
%            elseif RangeSim
%if Index == 1; Exp = c(1);
%elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
%elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
%end
 %           else Exp = c(1);
 %           end 
        
  %          HatchFact = Exp;
            
   %     Fec = Fec*(1-P_exp) + Fec*P_exp*HatchFact;
        
        
        % penalty for LARVAL SURVIVAL:
      %  if exposure
            
            % This bit is deprecated bc F0 larval survival now handled
            % outside kernel
%            c = Coef(1).LarvSurvF0(1).(chemical).c;
%            if RandSim
%            Exp = normrnd(c(1),c(2)*c(1)); 
%            elseif RangeSim
%if Index == 1; Exp = c(1);
%elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
%elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
%end
%            else Exp = c(1);
%            end 
            
 %           Fec = Fec*Exp;

      %  end
        
        % F2 effects:
        c = Coef(1).LarvSurvF2(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
            else
                Exp = c(1);
            end 
            
            
            
            
        % Effect on F2 larval success would be felt in F1 kernel
        P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;
        Fec = Fec * (1-P_exp) + Fec * (Exp * P_exp);
        
        %%%%%%%
        % F1 effects of egg production 
       c = Coef(1).EggProdF1(1).(chemical).c;
            if RandSim
            Exp = normrnd(c(1),c(2)*c(1)); 
            elseif RangeSim
if Index == 1; Exp = c(1);
elseif Index == 2; Exp = c(1)+1.96*c(2)*c(1);
elseif Index == 3; Exp = c(1)-1.96*c(2)*c(1);
end
            else
                Exp = c(1);
            end 
            
            
        % Effect on F2 larval success would be felt in F1 kenrel
        P_exp = (exposure_p + sum(exposure_vec(1:length(F_dist)).*F_dist(:)'))*0.5;
        Fec = Fec * (1-P_exp) + Fec * (Exp * P_exp);
%%%%%%


        % Ensure it's positive
        Fec = max(Fec,realmin);
        
        % Now also apply mating function: (based on White et al. 2017)
        mating_func = min(1,betacdf(SR*2,1,1.3));
        
        kxy = mating_func*Fec.*Rvec;
        
end % end switch case

