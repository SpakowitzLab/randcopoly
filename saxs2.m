%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;close all
global sf qf N rm FA NFIT

% Load experiment data.
addpath('functions/')
data = load('exp-data/PEG30undopedMW1500.csv');  % SAXS data with q in A^(-1)
                                           % 30wt ~= 16mol%
% preset parameters in theory
N=100;  % total of 100 monomers
rm = 31.84; % estimate end-to-end distance of a "chemical monomer" in unit Angstrom
FA=0.126;   % equal chemical composition

TV = [22,40:20:180];  % temperature in degree C
TK = TV+273.15;       % temperature in Kelvin
q = data(:,1);
s = data(:,2:end);
%% 
% Define a small q cutoff and estimate peak intensity S(q~=0).

% if PEG wt% = 0.30
qmin0 = 0.03803;
qmin0f = 0.04094;
qfinal = 0.22138;

iqmin0 = find(q==qmin0);
iqmin0f = find(q==qmin0f);
iqfinal = find(q==qfinal);
smin = mean(s(iqmin0:iqmin0f,:));
sminstd = std(1./s(iqmin0:iqmin0f,:));

%% 
% Start fitting analytical S(q) to experimental data.

% number of fitting temperatures
NFIT = 9;

% define fitting range
qf = q(iqmin0:iqfinal);
sf = zeros(length(qf),NFIT);
for IT = 1:NFIT
    sf(:,IT) = s(iqmin0:iqfinal,IT);
end

% initial fit guess
% x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT
x0 = [25,0.6,   0.1,ones(1,NFIT)];
lb = [1  ,0.001,0.01, ones(1,NFIT)*0.001];
ub = [1e3,1.0  ,10.0, ones(1,NFIT)*10.0 ];

% start fitting
options = optimset('TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',1e2,'MaxIter',1e2);
fit = lsqnonlin(@saxsfit2,x0,lb,ub,options);

% Save fitted results
SCALE = fit(1);
LAM = fit(2);
NM = fit(3);
CHI = fit(4:end);
save('PEG30undopedMW1500.mat','SCALE','LAM','NM','CHI','TK','qf');