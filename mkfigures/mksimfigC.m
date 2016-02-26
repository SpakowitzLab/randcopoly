f1=figure;hold;set(gca,'fontsize',30)
set(f1,'position',[0,0,800,600])

LAM=-0.75;
plotsim2(0.01,LAM);
plotsim2(0.10,LAM);
plotsim2(1.00,LAM);
savename = sprintf('../../results/randcopoly-results/random-simulation/Figure7A');
saveas(gcf,savename,'epsc')

f1=figure;hold;set(gca,'fontsize',30)
set(f1,'position',[0,0,800,600])

LAM=0.00;
plotsim2(0.01,LAM);
plotsim2(0.10,LAM);
plotsim2(1.00,LAM);
%legend('N_M=0.05','N_M=0.50','N_M=5.00')
savename = sprintf('../../results/randcopoly-results/random-simulation/Figure7B');
saveas(gcf,savename,'epsc')
