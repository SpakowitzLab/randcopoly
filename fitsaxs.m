%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;close all
% FILENAME = 'PEG30undopedMW1500';
FILENAME = 'PEG30undopedMW900';

global sf qf N rm FA NFIT

% Load experiment data, SAXS data with q in A^(-1)
addpath('functions/')
data = load(strcat('exp-data/',FILENAME,'.csv'));  

% preset parameters in theory
N=100;  % total of 100 monomers
[FA,rm]=calcmol(1500,0.3);

if strcmp(FILENAME,'PEG30undopedMW1500')
    TV = [22,40:20:180];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
elseif strcmp(FILENAME,'PEG30undopedMW900')
    TV = [22,40:20:160];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
end
q = data(:,1);
s = data(:,2:end);
%% 
% Define a small q cutoff and estimate peak intensity S(q~=0).

if strcmp(FILENAME,'PEG30undopedMW1500')
    % if PEG wt% = 0.30
    qmin0 = 0.03803;
    qmin0f = 0.04094;
    qfinal = 0.22138;
elseif strcmp(FILENAME,'PEG30undopedMW900')
    qmin0 = 0.0577;
    qmin0f = 0.08312;
    qfinal = 0.31367;
end

iqmin0 = find(q==qmin0);
iqmin0f = find(q==qmin0f);
iqfinal = find(q==qfinal);
smin = mean(s(iqmin0:iqmin0f,:));
sminstd = std(1./s(iqmin0:iqmin0f,:));

%% 
% Start fitting analytical S(q) to experimental data.

% number of fitting temperatures
NFIT = length(TV);

% define fitting range
qf = q(iqmin0:iqfinal);
sf = zeros(length(qf),NFIT);
for IT = 1:NFIT
    sf(:,IT) = s(iqmin0:iqfinal,IT);
end

%% Fit with random copolymer model
% initial fit guess
% x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT
x0 = [100, 0.55,   0.1, ones(1,NFIT)*10];
lb = [ 50,  -1,   0.01, ones(1,NFIT)*0.1];
ub = [500,   1,   10.0, ones(1,NFIT)*100.0];

% start fitting
options = optimset('MaxFunEvals',1e2,'MaxIter',1e2,'Display','iter');
fit = lsqnonlin(@saxsfitfunc,x0,lb,ub,options);

%% % reconstruct fitted function
sfit = zeros(length(qf),NFIT);
for IT = 1:NFIT
    % mdl = @(x,qf) x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm)*x(3))';
    sfit(:,IT) = fit(1)./(-2*fit(3+IT)+s2invwlc(N,fit(3),FA,fit(2),qf*rm));
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

% Save fitted results
SCALE = fit(1);
LAM = fit(2);
NM = fit(3);
CHI = fit(4:end);
save(strcat('savedata/',FILENAME,'.mat'),'SCALE','LAM','NM','CHI','TK','qf');




%% % Fit with Ornstein-Zernike function (Lorentzian)
% % initial fit guess
% % x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT
% x0 = [1,ones(1,NFIT)*10];
% 
% % start fitting
% % options = optimset('TolX',1e-2,'TolFun',1e0,'MaxFunEvals',1e2,'MaxIter',1e2,...
% %     'Display','iter');
% % fit = lsqnonlin(@ozfit,x0,lb,ub,options);
% fit = lsqnonlin(@ozfit,x0,[],[]);
% 
% %% % reconstruct fitted function
% sfit = zeros(length(qf),NFIT);
% for IT = 1:NFIT
%     % mdl = @(x,qf) x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm)*x(3))';
%     sfit(:,IT) = fit(1)*fit(1+IT).^2./(1+qf.^2*fit(1+IT)^2);
% end
% 
% figure;hold
% for IT = 1:NFIT
%     col = (IT-1)/(NFIT-1);
%     plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
%     plot(qf*rm,sfit(:,IT),'--','LineWidth',2,'color',[col 0 1-col]) % theoretical fit
% end
% xlabel('R_Mq');ylabel('S(q)');
% set(gca,'xscale','log');set(gca,'yscale','log');
% xlim([0.02,0.6]*rm);ylim([7,1e2]);
% box on