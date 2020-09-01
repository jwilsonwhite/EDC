function Coef = epigenetic_coefs

% Function to store means & SDs of all the epigenetic effects
% Will's attempt to keep things more organized

% Each vector is mean, CV


% At this time, larval length effects are not included

% Add deformity mortality. Assume deformities represent a proportional
% increase in mortality. 

% F0 effects (reduced deformites...not included because WTF)
Coef(1).DeformF0(1).EE2.c = [1,0];
Coef(1).DeformF0(1).Bif.c = [NaN,NaN];
Coef(1).DeformF0(1).TB.c = [NaN,NaN];
Coef(1).DeformF0(1).Levo.c = [1,0];


%For EE2, 0.35 vs control 0.1. 250% increase. SD = 0.15, so CV = 0.15/0.35
Coef(1).DeformF1(1).EE2.c = [2.5, 0.15/0.35];
Coef(1).DeformF1(1).Bif.c = [1 0];
% 0.4% vs. control 0.1%. 300% increase
Coef(1).DeformF1(1).TB.c = [3, 0.1/0.3];
Coef(1).DeformF1(1).Levo.c = [1 0];

% Add in deformity mortality observed in DeCourten & Brander 2017 (using 22 deg C):
Coef(1).DeformF0_dCB(1).EE2.c = [1,0];
Coef(1).DeformF0_dCB(1).Bif.c = [0.15/0.1,0.04/0.15];
Coef(1).DeformF0_dCB(1).TB.c = [1,0];
Coef(1).DeformF0_dCB(1).Levo.c = [1,0];

% Add in deformity mortality observed in F1 in dC&B 2017:
Coef(1).DeformF1_dCB(1).EE2.c = [4,0.05/0.07];
Coef(1).DeformF1_dCB(1).Bif.c = [2,0.07/0.05];
Coef(1).DeformF1_dCB(1).TB.c = [1,0];
Coef(1).DeformF1_dCB(1).Levo.c = [1,0];

% F2 deformity effects (NEED TO ADD THIS TO CODE)
Coef(1).DeformF2(1).EE2.c = [1,0];
Coef(1).DeformF2(1).Bif.c = [1 0];
Coef(1).DeformF2(1).TB.c = [1,0];
Coef(1).DeformF2(1).Levo.c = [2.5 0.1/0.35];


% effect on hatch size
Coef(1).HatchSize(1).EE2.c = [1, 0.1 ];
Coef(1).HatchSize(1).Bif.c = [0, 0 ];
Coef(1).HatchSize(1).TB.c = [1, 0.1 ];
Coef(1).HatchSize(1).Levo.c = [1, 0.1 ];

% Atresia
% EE2 exposure = 65%
        % EE2; So 35% good vs 80% good, reduction is 35/80 = 44%
Coef(1).Atresia(1).EE2.c = [0.44, 0.5/2.5];
           % 50% good vs. 80% good. Reduction is 50/80 = 0.625
Coef(1).Atresia(1).Bif.c = [0.625, 0.5/2.5];
Coef(1).Atresia(1).TB.c = [1,0];
Coef(1).Atresia(1).Levo.c = [1,0];

% Hatch success
Coef(1).HatchSuccessF0(1).EE2.c = [1,0];
Coef(1).HatchSuccessF0(1).Bif.c = [1,0];
% F0 reduction
Coef(1).HatchSuccessF0(1).TB.c = [70/82, 10/70];
Coef(1).HatchSuccessF0(1).Levo.c = [1,0];

% Hatch success (F1 effect) (% may need to do these as binomial)
Coef(1).HatchSuccessF1(1).EE2.c = [1+82/65, 7/82];
Coef(1).HatchSuccessF1(1).Bif.c = [1+85/65, 5/85];
Coef(1).HatchSuccessF1(1).TB.c = [1+80/65, 10/80 ];
Coef(1).HatchSuccessF1(1).Levo.c = [1+86/65, 12/86];

% Larval Survival F0
Coef(1).LarvSurvF0(1).EE2.c = [1,0 ];
Coef(1).LarvSurvF0(1).Bif.c = [1+(90-78)/78, 7/90 ];
Coef(1).LarvSurvF0(1).TB.c = [1,0];
Coef(1).LarvSurvF0(1).Levo.c = [1,0];

% Larval Survival F2
        % In F2, reduced by 0.75/0.85= 0.88
Coef(1).LarvSurvF2(1).EE2.c = [0.8824, 15/75];
Coef(1).LarvSurvF2(1).Bif.c = [1,0];
Coef(1).LarvSurvF2(1).TB.c = [1,0 ];
Coef(1).LarvSurvF2(1).Levo.c = [65/85, 7/85 ];


