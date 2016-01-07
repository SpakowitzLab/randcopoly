function [kval,sval,d2gam2]=kmaxrr(N,NM,FA,LAM)

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
ks = kval;
dks = 1/sqrt(R2)*5e-2;
G = @(k) s2invrr(N,NM,FA,LAM,k);

if ks>1e-1  % central differences
    d2gam2 = (G(ks+dks)-2*G(ks)+G(ks-dks))/(dks^2);
else  % forward differences
    d2gam2 = (G(ks+2*dks)-2*G(ks+dks)+G(ks))/(dks^2);
end
d2gam2 = -1/(NM*R2)*d2gam2./power(G(ks),2);
end