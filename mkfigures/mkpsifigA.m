% This is an example calculating structure factor
% of rigid, anti-correlated random copolymer melts,
% correpsonding to Figure 2A in "Impact of ..." paper
addpath('../functions')
close all
clear

SAVEON=1;

N=100;
NM=0.1;
LAM=-0.75;
FA=0.5;
d=3;

RM=sqrt(r2wlc(NM));
K0=1e-2;
KF=1e2;
NK=2001;
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

ORD=20;
ORDmax=20;
ResLayer=500;

[kval,sval]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer);

CHIS=0.5*sval;

CHI=CHIS*[0 0.2 0.4 0.6 0.8];

[kval,sval]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer);
KS=kval;
SS=sval;

[SINV]=s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer);

for I=1:length(CHI)
    COL=(I-1)/(length(CHI)-1);

    figure(1)
    loglog(RM*K,1./(-2*CHI(I)+SINV),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    loglog(RM*KS,1/(-2*CHI(I)+SS),'o','LineWidth',2,'MarkerSize',4,'Color',[COL 0 1-COL])
    
end
SC=1;
PSIMIN=min(1./SINV);
IND0=1600;
INDF=length(K);

figure(1)
loglog(RM*K(IND0:INDF),2*PSIMIN*power(K(IND0:INDF)/(KF/RM),-SC),'k:','LineWidth',2)
ylabel('<\psi(q)\psi(-q)>','FontSize',18)
xlabel('R_{m}q','FontSize',18)
set(gca,'FontSize',14)
text(0.13*KF,20*PSIMIN,['~q^{-1}'],'FontSize',14)
axis([K0 KF 1e-4 1e-1])
set(gca,'Xtick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
set(gca,'Ytick',[1e-4 1e-3 1e-2 1e-1])

if SAVEON==1
    figure(1)
    saveas(gcf,'fig2A.eps','epsc')    
end
