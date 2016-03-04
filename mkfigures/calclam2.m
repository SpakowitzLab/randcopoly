addpath('../functions')

N=100;  % total of 100 monomers
FA=0.5;    % equal chemical composition

% % Gaussian chain
% NM=100; % each monomer has 100 Kuhn steps
% RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
% 
% LAMV = linspace(-.99,.99,501);
% KS = zeros(length(LAMV),1);
% CHIS = zeros(length(LAMV),1);
% D2S = zeros(length(LAMV),1);
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii);
%     [kval,sval,d2gam2]=kmaxgc(N,NM,FA,LAM);
%     KS(ii)=kval;
%     CHIS(ii)=0.5*sval;  % spinodal
%     D2S(ii)=d2gam2./(sval^2*RM^2);
% end
% 
% data = [LAMV',RM*KS,CHIS*NM,D2S/NM];
% dlmwrite(sprintf('data/GC'),data)
% 
% % WLC
% NMV = logspace(-2,2,11);
% LAML_WLC = zeros(length(NMV),1);
% cnt = 1;
% for NM = NMV
%     col = (cnt-1)/(length(NMV)-1);
%     NM
%     RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
% 
%     LAMV = linspace(-1,.99,501);
%     KS = zeros(length(LAMV),1);
%     CHIS = zeros(length(LAMV),1);
%     D2S = zeros(length(LAMV),1);
%     for ii = 1:length(LAMV)
%         LAM = LAMV(ii)
%         [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
%         KS(ii)=kval;
%         CHIS(ii)=0.5*sval;  % spinodal
%         D2S(ii)=d2gam2./(sval^2*RM^2);
%     end
%     
%     data = [LAMV',RM*KS,CHIS*NM,D2S/NM];
%     dlmwrite(sprintf('data/WLC_NM%.2f',NM),data)
%     cnt = cnt+1;
% end

% Rigid rod
NM=0.01; % each monomer has 100 Kuhn steps
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

LAMV = linspace(-.99,.99,501);
KS = zeros(length(LAMV),1);
CHIS = zeros(length(LAMV),1);
D2S = zeros(length(LAMV),1);
for ii = 1:length(LAMV)
    LAM = LAMV(ii)
    [kval,sval,d2gam2]=kmaxrr(N,NM,FA,LAM);
    KS(ii)=kval;
    CHIS(ii)=0.5*sval;  % spinodal
    D2S(ii)=d2gam2./(sval^2*RM^2);
    data = [LAMV',RM*KS,CHIS*NM,D2S/NM];
    dlmwrite(sprintf('data/RR'),data)
end