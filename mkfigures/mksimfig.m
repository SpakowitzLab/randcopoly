clear;close all
cd ../

% start code
LAMV = 0.5;
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
RM=2;
Lbox=20;
Lbin=1;

figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',20)
color = ['r','b','k'];

cnt=1;p1=[];p2=[];
for LAM = LAMV
    for EPS = EPSV
        [chis,ks,chiv,ksim,sinv_theory,sinv_sim]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;

        figure(1);
        p1(cnt)=plot(chiv*G,sinv_sim,'o','color',color(cnt));
        plot(chiv*G,sinv_theory,'-','linewidth',3,'color',color(cnt))

        figure(2);
        p2(cnt)=plot(chiv*G,ksim,'o-','color',color(cnt));
        plot(chiv*G,ones(length(ksim),1)*ks*sqrt(R2),'--','color',color(cnt),'linewidth',3)
    %     plot(chiv*G,ones(length(ksim),1)*2*pi/Lbox*RM,'k--')
    %     plot(chiv*G,ones(length(ksim),1)*2*pi/Lbin*RM,'k--')
        cnt = cnt+1;
    end
end

figure(1);
xlabel('\chivG');ylabel('S^{-1}(q^*)')
% legend(p1,'N_M=0.05','N_M=0.50','N_M=5.00');box on
% if LAM==0
%     ylim([0,1]);xlim([0,8]);
%     set(gca,'Ytick',0:0.2:1.0)
%     set(gca,'YtickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})
% elseif LAM==-0.75
%     ylim([0,2]);xlim([0,20]);
%     set(gca,'Ytick',0:0.4:2.0)
%     set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0'})
% else
%     ylim([0,2]);xlim([0,20]);
%     set(gca,'Ytick',0:0.4:2.0)
%     set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0'})
% end

figure(2);
xlabel('\chivG');ylabel('q^*R_M')
legend(p2,'N_M=0.05','N_M=0.50','N_M=5.00');box on
% if LAM==0
%     ylim([0,4]);xlim([0,8]);
%     set(gca,'Ytick',[0,1,2,3,4])
%     set(gca,'YtickLabel',{'0','1','2','3','4'})
% elseif LAM==-0.75
%     ylim([2,5]);xlim([0,20]);
%     set(gca,'Ytick',2:0.5:5)
%     set(gca,'YtickLabel',{'2.0','2.5','3.0','3.5','4.0','4.5','5.0'})
% else
%     ylim([0,5]);xlim([0,20]);
%     set(gca,'Ytick',0:1.0:5)
%     set(gca,'YtickLabel',{'0','1','2','3','4','5'})
% end

figure(1)
savename = sprintf('../results/random-simulation/ssim-lam%.2f.eps',LAM);
saveas(gcf,savename,'epsc')

figure(2)
savename = sprintf('../results/random-simulation/qsim-lam%.2f.eps',LAM);
saveas(gcf,savename,'epsc')

% end code

cd mkfigures/