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

% Example 1: plot density-density correlations
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
LAM=-0.75; % anti-correlated random copolymer
FA=0.5;    % equal chemical composition
CHI=0.1/NM;  % Flory-Huggins parameter

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-5;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=201;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

% evaluate s2inv
S2INV=s2invwlc(N,NM,FA,LAM,K);
figure;loglog(RM*K,1./(-2*CHI+S2INV));
xlabel('R_Mq');ylabel('S(q)')
axis([K0 KF 1e-2 1e1])

% Example 2: find spinodal
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
LAM=0; % ideal random copolymer

FAV = linspace(0.1,0.9,38);
CHIS = zeros(length(FAV),1);
for ii = 1:length(FAV)
    FA = FAV(ii);
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
    CHIS(ii)=0.5*sval;  % spinodal
end
figure;plot(FAV,CHIS*NM)
xlabel('f_A');ylabel('\chi_S v N_M')

% Example 3: find critical wavemode
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
FA=0.5;    % equal chemical composition
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

LAMV = linspace(-1,.99,84);
KS = zeros(length(LAMV),1);
for ii = 1:length(LAMV)
    LAM = LAMV(ii);
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
    KS(ii)=kval;
end
figure;plot(LAMV,RM*KS)
xlabel('f_A');ylabel('R_Mq^*')

% Example 4: find peak sharpness
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
FA=0.5;    % equal chemical composition
CHI=0.1/NM;  % Flory-Huggins parameter

LAMV = linspace(-1,.99,84);
D2GAM2 = zeros(length(LAMV),1);
for ii = 1:length(LAMV)
    LAM = LAMV(ii);
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
    D2GAM2(ii)=d2gam2;
end
figure;plot(LAMV,log(-D2GAM2))
xlabel('f_A');ylabel('log|sharpness|')