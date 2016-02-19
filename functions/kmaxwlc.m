function [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer)
% Evaluate structure factor of random copolymer melt near peak, at
% zero Flory-Huggins parameter. Chains are modeled as wormlike chains.
% Usage :: [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer)
% Outputs ::
%   kval = critical wavemode of instability
%   sval = inverse of structure factor (s2invwlc) at kval
%   d2gam2 = second derivative of sval at kval
% Inputs ::
%   N = number of monomers
%   NM = number of Kuhn steps per monomer
%   FA = fraction of A monomers
%   LAM = degree of chemical correlation
%   d = number of dimensions, default 3
%   ORDmax = maximum number of eigenvalues, default 20
%   ORD = number of eigenvalues, default 20
%   ResLayer = number of residual layers, default 500
% Shifan Mao (1/6/16)

% Check conditions on chemical composition and correlation
if ~( FA>=0 && 1-FA>=0 && LAM>=-1 && 1-LAM>=0 && ...
     FA*(1-LAM)+LAM>=0 && FA*(LAM-1)+1>=0 )
 error('chemical composition and correlation constraints not satisfied')
end

% Fill in unset optional values

switch nargin
    case 4
        d=3;
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 5
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 6
        ORD=20;
        ResLayer=500;
    case 7
        ResLayer=500;        
end

R2=r2wlc(NM);

NK=50;

KV=transpose(logspace(-2,2,NK))/sqrt(R2);

S=s2invwlc(N,NM,FA,LAM,KV,d,ORDmax,ORD,ResLayer);
[~,IND]=min(S);

if IND ==1
    kval=1e-2/sqrt(R2);
    sval=s2invwlc(N,NM,FA,LAM,kval,d,ORDmax,ORD,ResLayer);
    kval=0;
else
    
    KV2=transpose(linspace(KV(IND-1),KV(IND+1),NK));
    DK=KV2(2)-KV2(1);
    S=s2invwlc(N,NM,FA,LAM,KV2,d,ORDmax,ORD,ResLayer);
    [~,IND]=min(S);
    K=KV2(IND);
    
    A=s2invwlc(N,NM,FA,LAM,K-DK,d,ORDmax,ORD,ResLayer);
    B=s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer);
    C=s2invwlc(N,NM,FA,LAM,K+DK,d,ORDmax,ORD,ResLayer);
    KAP=(A+C-2*B)/DK^2;
    kval=(A-C)/(2*KAP*DK)+K;
    
    sval=s2invwlc(N,NM,FA,LAM,kval,d,ORDmax,ORD,ResLayer);
    
end

% find peak sharpness
G = @(k) s2invwlc(N,NM,FA,LAM,k,d,ORDmax,ORD,ResLayer);

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