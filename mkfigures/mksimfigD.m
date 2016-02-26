% plot radius of gyration of single chains
clear all;

folder = '../../results/randcopoly-results/scalcbatch-12-15-15';
folder2 = '../../results/randcopoly-results/rdata-12-15-15';

% Plot Parameters
EPSV = [0.01,0.10,1.00];
LAM = -0.75;

% col = ['r','b','k'];
figure;hold;set(gca,'fontsize',20)

for ii = 1:3
EPS = EPSV(ii);
col = (ii-1)/(length(EPSV)-1);

% Simulation Parameters
N = 8;G = 5;Ree = 2;

% Calculate lengths of worm-like chains
L0 = Ree*EPS*((-0.5+0.5*exp(-2*EPS*G)+EPS*G)^(-0.5));  % Submonomer contour length
Lc = (N*G-1)*L0;  % Polymer contour length
Lp = L0/2/EPS;  % Persistent length
Rgid = power((1/3)*Lp*Lc-Lp^2+2*Lp^3/Lc*(1-(Lp/Lc)*(1-exp(-Lc/Lp))),1/2);

% Load data
simparam=load([folder,'/chivals']);
% find corresponding simulation at EPS and LAM
SIMNUM = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
% finding simulation index(indices)
chemparam=load([folder,'/chemind']);
CHEMNUM=chemparam(chemparam(:,1)==SIMNUM,2);

% Read in data
filename = sprintf('/RMC_SIM%dCHEM%d',SIMNUM,CHEMNUM);
R = load([folder2,filename]);
% CHIV = linspace(0,4.5,length(R));

% Plot radius of gyration
plot(R(:,1),R(:,2)/Rgid,'.-','linewidth',3,'color',[col 0 1-col],'markersize',20)

% Plot monomer length
% Lm = (G-1)*L0;  % Monomer contour length
%R2id = (1/3)*Lp*Lm-Lp^2+2*Lp^3/Lm*(1-(Lp/Lm)*(1-exp(-Lm/Lp)));
%R2id = (2*Lp)^2*(Lm/2/Lp-0.5*(1-exp(-Lm/Lp)));
%R2id = sqrt(R2id)
%plot(CHIV,R(:,3)/sqrt(R2id),'s--','linewidth',3,'color',col(ii))
end

xlabel('\chivG');ylabel('R_g/R_{id}')
axis([0 20 .98 1.15])