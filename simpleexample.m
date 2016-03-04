% Some simple examples using functions to evaluate
% phase behavior of random copolymers of wormlike chains
% kmaxwlc.m and s2invwlc
%   -> use kmaxgc.m and s2invgc.m for Gaussian chains
%   -> use kmaxrr.m and s2invrr.m for perfectly rigid rods
%      (expect rigid rod functions to be relatively slow with
%       large number of monomers)
addpath('functions')
close all
clear

% Simple Example 1: find spinodal, critical wavelength, and peak sharpness
N=100;  % total of 100 monomers
NM=0.1; % each monomer has 0.1 Kuhn steps
FA=0.5; % equal chemical composition
LAM=-0.75;  % statistically anti-correlated random copolymer
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

[KS,sval,D2GAM2]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*sval;  % spinodal

fprintf('Results :\nSpinodal chivN_M= %.2f\n',CHIS*NM)
fprintf('Critical wavemode R_Mq* = %.2f\n',KS*RM)
fprintf('Second der. of structure factor at q* = %.2f\n',D2GAM2/(sval^2*RM^2*NM))

% Simple Example 2: plot density-density correlations
CHI=0.5*CHIS;  % Flory-Huggins parameter
K0=1e-2;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=21;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

S2INV=s2invwlc(N,NM,FA,LAM,K);
figure;set(gca,'fontsize',20)
loglog(RM*K,1./(-2*CHI+S2INV),'k-','linewidth',2);
xlabel('R_Mq');ylabel('S(q)')
axis([K0 KF 1e-4 1e1])
