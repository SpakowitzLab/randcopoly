close all
clear

file='lamdata.mat';
CALC=0;
SAVEON=1;

NM0=0.01;
NMF=100;
NNM=101;
NM=transpose(logspace(log10(NM0),log10(NMF),NNM));
SKIP=10;

FA=0.5;
N=100;

ORD=25;
ORDmax=25;
ResLayer=500;

d=3;

NLAM=501;
LAM=transpose(linspace(-1,1,NLAM));

KS=zeros(NNM,NLAM);
CHIS=zeros(NNM,NLAM);
RM=zeros(NNM,1);

KSGC=zeros(NLAM,1);
CHISGC=zeros(NLAM,1);
NMGC=NM(NNM);
RMGC=sqrt(NMGC);

KSRR=zeros(NLAM,1);
CHISRR=zeros(NLAM,1);
NMRR=NM(1);
RMRR=NMRR;


if CALC==1
    
    for INM=1:NNM
        RM(INM)=sqrt(r2wlc(NM(INM)));
        for IL=1:NLAM
            [kval,sval]=kmaxwlc(N,NM(INM),FA,LAM(IL),d,ORDmax,ORD,ResLayer);
            CHIS(INM,IL)=0.5*sval;
            KS(INM,IL)=kval*RM(INM);
            if INM==1
                [kval,sval]=kmaxrr(N,NMRR,FA,LAM(IL));
                CHISRR(IL)=0.5*sval;
                KSRR(IL)=kval*RMRR;
            end
            
            if INM==NNM
                [kval,sval]=kmaxgc(N,NMGC,FA,LAM(IL));
                CHISGC(IL)=0.5*sval;
                KSGC(IL)=kval*RMGC;            
            end
            
        end
        INM
                
    end
    
    save(['data/',file],'N','NM','FA','LAM','CHIS','KS','RM',...
        'CHISGC','KSGC','RMGC','NMGC','CHISRR','KSRR','RMRR','NMRR','ORD','ORDmax','ResLayer','N','FA')
    
else
    load(['data/',file])
end

for INM=1:length(NM)
    KS(INM,:)=KS(INM,:).*heaviside(KS(INM,:)-0.022);
end

for INM=1:SKIP:NNM
    
    COL=(INM-1)/(NNM-1);

    figure(1)
    plot(LAM(2:501),KS(INM,2:501),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    
    figure(2)
    plot(LAM(2:501),CHIS(INM,2:501)*NM(INM),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on

    figure(3)
    plot(LAM(2:501),2*pi./KS(INM,2:501),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    
    
end

figure(1)
plot(LAM,KSGC,'k-','LineWidth',2)
hold on

figure(2)
plot(LAM,CHISGC*NMGC,'k-','LineWidth',2)
hold on

figure(3)
plot(LAM,2*pi./KSGC,'k-','LineWidth',2)
hold on

figure(1)
plot(LAM,KSRR,'k--','LineWidth',2)
hold on

figure(2)
plot(LAM,CHISRR*NMRR,'k--','LineWidth',2)
hold on

figure(3)
plot(LAM,2*pi./KSRR,'k--','LineWidth',2)
hold on

figure(1)
set(gca,'FontSize',14)
xlabel('\lambda','FontSize',18)
ylabel('R_{m}q^{*}','FontSize',18)

figure(2)
set(gca,'FontSize',14)
xlabel('\lambda','FontSize',18)
ylabel('4f_{A}f_{B}v\chiN_{m}','FontSize',18)

figure(3)
axis([-1 0 0 6])
set(gca,'FontSize',14)
xlabel('\lambda','FontSize',18)
ylabel('D=2 \pi/(R_{m}q^{*})','FontSize',18)


LAML=zeros(NNM,1);
DEL=4;
POW=4;
for INM=1:NNM
    CRS=0;
    IND=3;
    while CRS==0
        if KS(INM,IND)==0
            P=polyfit(transpose(KS(INM,(IND-1-DEL):(IND-1))),LAM((IND-1-DEL):(IND-1)),POW);
            
            LAML(INM)=polyval(P,0);
            CRS=1;
        else
            IND=IND+1;
        end
    end
    
end


CRS=0;
IND=3;
while CRS==0
    if KSRR(IND)==0
        P=polyfit(KSRR((IND-1-DEL):(IND-1)),LAM((IND-1-DEL):(IND-1)),POW);
        
        LAMLRR=polyval(P,0);
        CRS=1;
    else
        IND=IND+1;
    end
end
CRS=0;
IND=3;
while CRS==0
    if KSGC(IND)==0
        P=polyfit(KSGC((IND-1-DEL):(IND-1)),LAM((IND-1-DEL):(IND-1)),POW);
        
        LAMLGC=polyval(P,0);
        CRS=1;
    else
        IND=IND+1;
    end
end

%LAMLGC=-(2-sqrt(3));

figure(4)
semilogx(NM,LAML,'b-','LineWidth',2)
hold on
semilogx(NM,LAMLGC+0*NM,'k-','LineWidth',2)
semilogx(NM,LAMLRR+0*NM,'k--','LineWidth',2)
axis([0.01 100 -0.28 -0.10])
set(gca,'FontSize',14)
xlabel('N_{M}','FontSize',18)
ylabel('\lambda_{L}','FontSize',18)

figure(5)

for INM=1:SKIP:NNM
    
    COL=(INM-1)/(NNM-1);
    
    loglog((LAML(INM)-LAM(2:501)),KS(INM,2:501),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
    
end

figure(5)
loglog((LAMLRR-LAM(2:501)),KSRR(2:501),'k--','LineWidth',2)
loglog((LAMLGC-LAM(2:501)),KSGC(2:501),'k-','LineWidth',2)

figure(5)

loglog((LAML(INM)-LAM(150:201)),1.8*power((LAML(INM)-LAM(150:201))/0.1,0.5),'k:','LineWidth',2)
text(0.06,1.3,'~(\lambda_{L}-\lambda)^{1/2}','FontSize',14)

axis([1e-2 1e0 1e0 5e0])
set(gca,'FontSize',14)
xlabel('\lambda_{L}-\lambda','FontSize',18)
ylabel('R_{m}q^{*}','FontSize',18)

IND=63;
NFA=1000;
FA=transpose(linspace(0,1,NFA));

figure(6)
for INM=1:SKIP:NNM
    COL=(INM-1)/(NNM-1);
    plot(FA,NM(INM)*CHIS(INM,IND)./(FA.*(1-FA)*4),'-','LineWidth',2,'Color',[COL 0 1-COL])
    hold on
end
plot(FA,NMRR*CHISRR(IND)./(FA.*(1-FA)*4),'k--','LineWidth',2)
plot(FA,NMGC*CHISGC(IND)./(FA.*(1-FA)*4),'k-','LineWidth',2)


axis([0 1 0 15])
set(gca,'FontSize',14)
xlabel('f_{A}','FontSize',18)
ylabel('v \chi_{S} N_{M}','FontSize',18)

IND=251;
NFA=1000;
FA=transpose(linspace(0,1,NFA));

figure(7)
%for INM=1:SKIP:NNM
%    COL=(INM-1)/(NNM-1);
%    plot(FA,NM(INM)*CHIS(INM,IND)./(FA.*(1-FA)*4),'-','LineWidth',2,'Color',[COL 0 1-COL])
%    hold on
%end
%plot(FA,NMRR*CHISRR(IND)./(FA.*(1-FA)*4),'k--','LineWidth',2)
plot(FA,NMGC*CHISGC(IND)./(FA.*(1-FA)*4),'k-','LineWidth',2)

axis([0 1 0 15])
set(gca,'FontSize',14)
xlabel('f_{A}','FontSize',18)
ylabel('v \chi_{S} N_{M}','FontSize',18)


if SAVEON==1
    figure(1)
    saveas(gcf,'qfig.eps','epsc')
    
    figure(2)
    saveas(gcf,'chifig.eps','epsc')
    
    figure(3)
    saveas(gcf,'dfig.eps','epsc')
    
    figure(4)
    saveas(gcf,'lamfig.eps','epsc')

    figure(5)
    saveas(gcf,'qlogfig.eps','epsc')

    figure(6)
    saveas(gcf,'phase75fig.eps','epsc')

    figure(7)
    saveas(gcf,'phase0fig.eps','epsc') 
end