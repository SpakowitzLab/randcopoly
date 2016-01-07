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

% Simple example: wormlike chain code
N=100;  % total of 100 monomers
NM=100; % each monomer has 100 Kuhn steps
FA=0.5; % equal chemical composition
LAM=0;  % ideal random copolymer

[KS,sval,D2GAM2]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*sval;  % spinodal

sprintf('Spinodal chivN_M= %.2f',CHIS*NM)
sprintf('Critical wavemode q* = %.2f',KS)
sprintf('Second der. of structure factor at q* = %.2f',D2GAM2)