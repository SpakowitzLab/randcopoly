% Plot R/L probability distribution from ideal chains
% and Monte-Carlo simulations
clear;close all

f1=figure;hold;set(gca,'fontsize',30)
set(f1,'position',[0,0,800,600])

LAM=-0.75;
plotsim2(0.01,LAM);
plotsim2(0.10,LAM);
plotsim2(1.00,LAM);
text(0.075,7,'N_M=5.00','FontSize',20)
text(0.4,3,'N_M=0.50','FontSize',20)
text(0.8,8.8,'N_M=0.05','FontSize',20)
savename = sprintf('../../results/randcopoly-results/random-simulation/Figure7A');
saveas(gcf,savename,'epsc')

f1=figure;hold;set(gca,'fontsize',30)
set(f1,'position',[0,0,800,600])

LAM=0.00;
plotsim2(0.01,LAM);
plotsim2(0.10,LAM);
plotsim2(1.00,LAM);
text(0.075,7,'N_M=5.00','FontSize',20)
text(0.4,3,'N_M=0.50','FontSize',20)
text(0.8,8.8,'N_M=0.05','FontSize',20)
savename = sprintf('../../results/randcopoly-results/random-simulation/Figure7B');
saveas(gcf,savename,'epsc')