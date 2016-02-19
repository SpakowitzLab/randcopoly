% Some simple examples using functions to evaluate
% phase behavior of random copolymers of wormlike chains
% kmaxwlc.m and s2invwlc
%   -> use kmaxgc.m and s2invgc.m for Gaussian chains
%   -> use kmaxrr.m and s2invrr.m for perfectly rigid rods
%      (expect rigid rod functions to be relatively slow with
%       large number of monomers)
addpath('functions')
%close all
clear

% Example 1: plot density-density correlations vs wavevector at different CHI
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
LAM=0.99; % anti-correlated random copolymer
FA=0.5;    % equal chemical composition

% find spinodal CHIS
[kval,sval]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*sval;
CHI=CHIS*[0 0.2 0.4 0.6 0.8];

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-2;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=201;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

% evaluate s2inv
[SINV]=s2invwlc(N,NM,FA,LAM,K);
figure;hold
for I=1:length(CHI)
    COL=(I-1)/(length(CHI)-1);
    loglog(RM*K,1./(-2*CHI(I)+SINV),'-','LineWidth',2,'Color',[COL 0 1-COL])
end
xlabel('R_Mq');ylabel('S(q)')
set(gca,'xscale','log');set(gca,'yscale','log');
axis([K0 KF 1e-1 1e2])

% Example 2: find spinodal vs. fraction of A monomers
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

% Example 3: find critical wavemode and spinodal vs. chemical correlation
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
FA=0.5;    % equal chemical composition
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

LAMV = linspace(-1,.99,501);
KS = zeros(length(LAMV),1);
CHIS = zeros(length(LAMV),1);
D2S = zeros(length(LAMV),1);
for ii = 1:length(LAMV)
    LAM = LAMV(ii);
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
    KS(ii)=kval;
    CHIS(ii)=0.5*sval;  % spinodal
    D2S(ii)=1/(sval^2*RM^2)*d2gam2;
end
figure;plot(LAMV,RM*KS)
xlabel('f_A');ylabel('R_Mq^*')

figure;plot(LAMV,CHIS*NM)
xlabel('f_A');ylabel('\chi_Sv')

figure;plot(LAMV,D2S/NM);

% Example 4: find peak sharpness vs. chemical correlation
N=100;  % total of 100 monomers
NM=1; % each monomer has 100 Kuhn steps
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