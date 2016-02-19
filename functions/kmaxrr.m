function [kval,sval,d2gam2]=kmaxrr(N,NM,FA,LAM)
% Evaluate structure factor of random copolymer melt near peak, at
% zero Flory-Huggins parameter. Chains are modeled as perfectly rigid rods
% Usage :: [kval,sval,d2gam2]=kmaxrr(N,NM,FA,LAM)
% Outputs ::
%   kval = critical wavemode of instability
%   sval = inverse of structure factor (s2invwlc) at kval
%   d2gam2 = second derivative of sval at kval
% Inputs ::
%   N = number of monomers
%   NM = number of Kuhn steps per monomer
%   FA = fraction of A monomers
%   LAM = degree of chemical correlation
% Shifan Mao (1/6/16)

% Check conditions on chemical composition and correlation
if ~( FA>=0 && 1-FA>=0 && LAM>=-1 && 1-LAM>=0 && ...
     FA*(1-LAM)+LAM>=0 && FA*(LAM-1)+1>=0 )
 error('chemical composition and correlation constraints not satisfied')
end

R2=NM^2;

NK=50;
KV=transpose(logspace(-2,2,NK))/sqrt(R2);

S=s2invrr(N,NM,FA,LAM,KV);
[~,IND]=min(S);

if IND ==1
    kval=1e-2/sqrt(R2);
    sval=s2invrr(N,NM,FA,LAM,kval);
    kval=0;
else

    KV2=transpose(linspace(KV(IND-1),KV(IND+1),NK));
    DK=KV2(2)-KV2(1);
    S=s2invrr(N,NM,FA,LAM,KV2);
    [~,IND]=min(S);
    K=KV2(IND);
    
    A=s2invrr(N,NM,FA,LAM,K-DK);
    B=s2invrr(N,NM,FA,LAM,K);
    C=s2invrr(N,NM,FA,LAM,K+DK);
    KAP=(A+C-2*B)/DK^2;
    kval=(A-C)/(2*KAP*DK)+K;

    sval=s2invrr(N,NM,FA,LAM,kval);

end

% find peak sharpness
G = @(k) s2invrr(N,NM,FA,LAM,k);

if (IND==1 || kval*sqrt(R2)<=1e-1)  % forward differences
    ks = 1/sqrt(R2)*5e-2;
    dks = 1/sqrt(R2)*5e-2;
    d2gam2 = (G(ks+2*dks)-2*G(ks+dks)+G(ks))/(dks^2);
else  % central differences
    ks = kval;
    dks = 1/sqrt(R2)*5e-2;
    d2gam2 = (G(ks+dks)-2*G(ks)+G(ks-dks))/(dks^2);
end
end