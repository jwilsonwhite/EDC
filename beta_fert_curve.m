function beta_fert_curve

% Figure to illustrate some different possible shapes of the fertilization
% curve

figure(10)
clf
set(gcf,'position',[200 500 400 300])

hold on

R = linspace(0,0.5,100);
a = [1 1];
b = [1, 40];
LS = {'-','--'};

for bb = 1:length(b)
    
    plot(R,betacdf(R*2,a(bb),b(bb)),'linestyle',LS{bb},'color','k')
    
end

set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
set(gca,'fontsize',12)
ylabel('Fertilization rate','fontsize',14)
xlabel('Sex ratio (proportion male)','fontsize',14)
axis square

% Now one more
a = 6;
b = 15;
plot(R,betacdf(R*2,a,b),'linestyle',':','color','k')



