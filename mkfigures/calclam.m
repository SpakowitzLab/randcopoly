addpath('../functions')
close all
clear

file='lamdata.mat';
CALC=1;

NM0=0.01;
NMF=100;
NNM=101;
NM=transpose(logspace(log10(NM0),log10(NMF),NNM));

FA=0.5;
N=100;

ORD=20;
ORDmax=20;
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
                
    
        save(['data/',file],'N','NM','FA','LAM','CHIS','KS','RM',...
        'CHISGC','KSGC','RMGC','NMGC','CHISRR','KSRR','RMRR','NMRR','ORD','ORDmax','ResLayer','N','FA')
    end
    
else
    load(['data/',file])
end
