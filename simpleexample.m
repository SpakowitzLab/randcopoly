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
NM=100; % each monomer has 100 Kuhn steps
FA=0.5; % equal chemical composition
LAM=0;  % ideal random copolymer

[KS,sval,D2GAM2]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*sval;  % spinodal

sprintf('Spinodal chivN_M= %.2f',CHIS*NM)
sprintf('Critical wavemode q* = %.2f',KS)
sprintf('Second der. of structure factor at q* = %.2f',D2GAM2)

% Simple Example 2: plot density-density correlations
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
LAM=-0.75; % anti-correlated random copolymer
FA=0.5;    % equal chemical composition
CHI=0.1/NM;  % Flory-Huggins parameter

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-2;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=201;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

% evaluate s2inv
S2INV=s2invwlc(N,NM,FA,LAM,K);
figure;loglog(RM*K,1./(-2*CHI+S2INV));
xlabel('R_Mq');ylabel('S(q)')
axis([K0 KF 1e-1 1e2])