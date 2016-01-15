% Example 2: find spinodal vs. fraction of A monomers
addpath('../functions')
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

% figure(6)
% for INM=1:SKIP:NNM
%     COL=(INM-1)/(NNM-1);
%     plot(FA,NM(INM)*CHIS(INM,IND)./(FA.*(1-FA)*4),'-','LineWidth',2,'Color',[COL 0 1-COL])
%     hold on
% end
% plot(FA,NMRR*CHISRR(IND)./(FA.*(1-FA)*4),'k--','LineWidth',2)
% plot(FA,NMGC*CHISGC(IND)./(FA.*(1-FA)*4),'k-','LineWidth',2)
% 
% axis([0 1 0 15])
% set(gca,'FontSize',14)
% xlabel('f_{A}','FontSize',18)
% ylabel('v \chi_{S} N_{M}','FontSize',18)
% 
% IND=251;
% NFA=1000;
% FA=transpose(linspace(0,1,NFA));
% 
% figure(7)
% %for INM=1:SKIP:NNM
% %    COL=(INM-1)/(NNM-1);
% %    plot(FA,NM(INM)*CHIS(INM,IND)./(FA.*(1-FA)*4),'-','LineWidth',2,'Color',[COL 0 1-COL])
% %    hold on
% %end
% %plot(FA,NMRR*CHISRR(IND)./(FA.*(1-FA)*4),'k--','LineWidth',2)
% plot(FA,NMGC*CHISGC(IND)./(FA.*(1-FA)*4),'k-','LineWidth',2)
% 
% axis([0 1 0 15])
% set(gca,'FontSize',14)
% xlabel('f_{A}','FontSize',18)
% ylabel('v \chi_{S} N_{M}','FontSize',18)
