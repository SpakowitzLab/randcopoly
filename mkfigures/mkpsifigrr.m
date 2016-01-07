addpath('../functions')
close all
clear

N=100;
NM=0.01;
LAM=-0.9;
FA=0.5;
d=3;

RM=sqrt(r2wlc(NM));
K0=1e-3;
KF=1e3;
NK=2001;
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

ORD=25;
ORDmax=25;
ResLayer=500;

[kval,sval]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer);

CHIS=0.5*sval;

CHI=CHIS*[0];

[kval,sval]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer);
KS=kval;
SS=sval;

[SINV]=s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer);
[SINVRR]=s2invrr(N,NM,FA,LAM,K);
RMRR=NM;

for I=1:length(CHI)
%    COL=(I-1)/(length(CHI)-1);
    COL=1;
    
    figure(1)
    semilogx(RM*K,1./(-2*CHI(I)+SINV),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    semilogx(RM*KS,1/(-2*CHI(I)+SS),'o','LineWidth',2,'MarkerSize',4,'Color',[COL 0 1-COL])

    figure(2)
    loglog(RM*K,1./(-2*CHI(I)+SINV),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    loglog(RM*KS,1/(-2*CHI(I)+SS),'o','LineWidth',2,'MarkerSize',4,'Color',[COL 0 1-COL])
    loglog(RMRR*K,1./(-2*CHI(I)+SINVRR),'--','LineWidth',2,'Color',[COL 0 1-COL])
    
end

figure(1)
ylabel('<\psi(q)\psi(-q)>','FontSize',18)
xlabel('R_{m}q','FontSize',18)
set(gca,'FontSize',14)
axis([K0 KF 0 1e-1])

SC=3/3;
PSIMIN=min(1./SINV);
IND0=1600;
INDF=length(K);
figure(2)
loglog(RM*K(IND0:INDF),2*PSIMIN*power(K(IND0:INDF)/(KF/RM),-SC),'k--')

figure(2)
ylabel('<\psi(q)\psi(-q)>','FontSize',18)
xlabel('R_{m}q','FontSize',18)
set(gca,'FontSize',14)
text(0.13*KF,70*PSIMIN,['~q^{-1}'],'FontSize',14)
axis([K0 KF 1e-5 1e-1])

