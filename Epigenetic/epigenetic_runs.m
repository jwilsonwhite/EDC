function epigenetic_runs

% Wrapper function to generate plots for PRIMO 2019
% Update Oct 2019: Add in parameter uncertainty for SETAC 2019

% Read in the variability in each parameter here...




doBase = false;
doRuns = [0 0 0 0];
RunNum = 1e2; % how many random sims to do?
if RunNum < 2; RandSim = false; else RandSim = true; end
RandSim = true;
RangeSim = false;

Chems = {'EE2','Bif','TB','Levo'};
Types = {'Chronic','Single_summer'}; %,'Single_winter'};
Types2 = {'Full_epi','No_epi','F1_only','paternal_only'}; % 'F2_only
TT = [1 1 0; 0 0 0; 1 0 0; 0 1 0; 1 1 1];
Cols = [0.9 0.1 0.1; 0.1 0.1 0.9; 0.7 0 0.7; 0.5 0.5 0.5];
Lty = {'-','-','-','-'};

% Labels:
Labs = {'No exposure','F0+F1+F2','F0 only','F0+F1','Paternal only'};
% Take out F2 only
        

Dir = 'runs_RangeSim_Oct2019/';

% Baseline run:
if doBase
    
    for i = 1:RunNum
    [Nf_tmp,Nm_tmp, dx]=epigenetic_model('Baseline',1,1,0,'EE2',RandSim,RangeSim,i);
    Nf(:,:,i) = Nf_tmp;
    Nm(:,:,i) = Nm_tmp;
         end % end loop over RunNum
             
         fname = strcat(Dir,'Baseline','.mat');
         save(fname,'Nf','Nm','dx')
    
%[Nf,Nm,dx]=epigenetic_model('Baseline',1,1,0,'EE2',RandSim,RangeSim,1);
%fname = strcat(Dir,'Baseline','.mat');
%save(fname,'Nf','Nm','dx')
end

for c = 1:length(Chems)
    
    if doRuns(c)
    for t  = 1:length(Types)
     for tt  = 1:length(Types2)  
         for i = 1:RunNum
    [Nf_tmp,Nm_tmp, dx]=epigenetic_model(Types{t},TT(tt,1),TT(tt,2),TT(tt,3),Chems{c},RandSim,RangeSim,i);
    Nf(:,:,i) = Nf_tmp;
    Nm(:,:,i) = Nm_tmp;
         end % end loop over RunNum
             
         fname = strcat(Dir,Chems{c},'_',Types{t},'_',Types2{tt},'.mat');
         save(fname,'Nf','Nm','dx')
         end % end loop over Types2 
    end % end loop over Types
    
    end % end if doRuns
end % end loop over Chems



% Now make plots
for c = 4%1:length(Chems)
    
    
    for t  = 1%:length(Types)
        
    figure
    clf
    set(gcf,'units','cent','position',[10 10 12 13])
    
    hold on
    % Baseline:
    fname = strcat(Dir,'Baseline.mat');
    load(fname,'Nf','Nm','dx')
    Nt = squeeze(sum(Nf,1))*dx + squeeze(sum(Nm,1))*dx; % now these are time x rep
    Nt = Nt./repmat(Nt(70,:),[size(Nt,1),1]); % scale to t = 70
    T = size(Nt,1);
    Time = (1:T)/6 - 70/6;
    Nt_med = quantile(Nt,0.5,2);
    Nt_up = quantile(Nt,0.95,2);
    Nt_low = quantile(Nt,0.05,2);
    plot(Time,Nt_med,'k-','linewidth',2)
    plot(Time,Nt_low,'k-','linewidth',2)
    plot(Time,Nt_up,'k-','linewidth',2)
    text(4.1,Nt_med(94),Labs{1},'color','k','fontsize',12)

    
       
    for tt  = 1:4%length(Types2)   
    fname = strcat(Dir,Chems{c},'_',Types{t},'_',Types2{tt},'.mat');
    load(fname,'Nf','Nm','dx')
    Nt = squeeze(sum(Nf,1))*dx + squeeze(sum(Nm,1))*dx; % now these are time x rep
    Nt = Nt./repmat(Nt(70,:),[size(Nt,1),1]); % scale to t = 70
    
    Nt_med = quantile(Nt,0.5,2);
    Nt_up = quantile(Nt,0.95,2);
    Nt_low = quantile(Nt,0.05,2);

    
    plot(Time,Nt_med,'k-','linewidth',2,'color',Cols(tt,:),'linestyle',Lty{tt});
    plot(Time,Nt_up,'k-','linewidth',2,'color',Cols(tt,:),'linestyle',Lty{tt});
    plot(Time,Nt_low,'k-','linewidth',2,'color',Cols(tt,:),'linestyle',Lty{tt});

    text(4.1,Nt_med(94),Labs{tt+1},'color',Cols(tt,:),'fontsize',12)
    
     end % end loop over Types2 
    
    
    set(gca,'tickdir','out','ticklength',[0.015 0.015])
    set(gca,'xcolor','k','ycolor','k')
    set(gca,'ylim',[0.0 1.1])
    xlabel('Time (y)','fontsize',14)
    ylabel('Relative population size','fontsize',14)
    axis square
    switch Types{t}
        case 'Chronic'
            
            set(gca,'xlim',[0,4])
            
        case {'Single_summer','Single_winter'}
            
            set(gca,'xlim',[0,4])
        
    end % end switch Types
    
    
    fname2 = strcat(Dir,Chems{c},'_',Types{t},'.eps');
    
    print(fname2,'-depsc2')
    

    
    
    
        end % end loop over Types
end % end loop over Chems

