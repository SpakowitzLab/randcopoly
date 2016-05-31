%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;close all
FILENAME = 'PEG30undopedMW1500';

% Load experiment data.
addpath('functions/')
data = load(strcat('exp-data/',FILENAME,'.csv'));  % SAXS data with q in A^(-1)
                                           % 30wt ~= 16mol%
% preset parameters in theory
N=100;  % total of 100 monomers
[FA,rm]=calcmol(1500,0.3);
q = data(:,1);
s = data(:,2:end);
load(strcat('savedata/',FILENAME,'.mat'));

%% PLOT 1: Structure factors
NFIT = length(TK);
sfit = zeros(length(qf),NFIT);
for IT = 1:NFIT
    sfit(:,IT) = SCALE./(-2*CHI(IT)+s2invwlc(N,NM,FA,LAM,qf*rm));
end

figure;hold
for IT = 1:NFIT
    col = (IT-1)/(NFIT-1);
    plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
    plot(qf*rm,sfit(:,IT),'--','LineWidth',2,'color',[col 0 1-col]) % theoretical fit
end
xlabel('R_Mq');ylabel('S(q)');
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([0.02,0.6]*rm);ylim([7,1e2]);
box on

%% PLOT 2: CHI vs T, and 1/S(q*) vs T
figure;
subplot(2,1,1);
plot(-1./TK,CHI,'ko')
xlabel('-1/T (K^{-1})');ylabel('\chi')

iq = find(q==qf(1));
subplot(2,1,2);
plot(-1./TK,1./s(iq,:),'ko')
xlabel('-1/T (K^{-1})');ylabel('S^{-1}(q=0)')

%% PLOT 3: Phase diagram
FAV = linspace(0.1,0.9,101);
CHISV = zeros(length(FAV),1);
for ii = 1:length(FAV)
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FAV(ii),LAM);
    CHISV(ii)=0.5*sval;  % spinodal
end

figure;hold;set(gca,'fontsize',18)
plot(FAV,CHISV*NM,'k-','linewidth',2)
for IT=1:NFIT
    col = (IT-1)/(NFIT-1);
    plot(FA,CHI(IT)*NM,'.','color',[col 0 1-col],...
        'MarkerSize',20);
    plot(0.2382,CHI(IT)*NM,'.','color',[col 0 1-col],...
        'MarkerSize',20);
end
xlabel('f_A');ylabel('\chi_S v N_M')
ylim([0.1,1.2]);xlim([0.1,0.9]);box on