clear;close all
cd ../

% start code
LAMV = -0.75;
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
RM=2;
Lbox=20;
Lbin=1;

% plot range
ind = 1:3:41;

figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',30)

% plot mean-field theory results
for LAM = LAMV
    ieps = 1;
    for EPS = EPSV
        col = (ieps-1)/(length(EPSV)-1);
        
        [chis,chiv,ks,ksim,SINV_MF,SINV_SIM]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        
        % plot range
        chiv = chiv(ind);
        SINV_SIM = SINV_SIM(ind);
        SINV_MF = SINV_MF(ind);

        figure(1);
        plot(chiv*G,SINV_MF,'--','linewidth',3,'color',[col 0 1-col])

        figure(2);
        plot(chiv*G,ones(length(chiv),1)*ks*sqrt(R2),'--','color',[col 0 1-col],'linewidth',3)
        ieps = ieps+1;
    end
end

% plot simulation results
cnt=1;p1=[];p2=[];
for LAM = LAMV
    ieps = 1;
    for EPS = EPSV
        col = (ieps-1)/(length(EPSV)-1);
        
        [chis,chiv,ks,ksim,SINV_MF,SINV_SIM]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        
        % plot range
        chiv = chiv(ind);
        SINV_SIM = SINV_SIM(ind);
        ksim = ksim(ind);

        figure(1);
        p1(cnt)=plot(chiv*G,SINV_SIM,'o',...
            'MarkerFaceColor',[col 0 1-col],'MarkerEdgeColor',[col 0 1-col],...
            'linewidth',3,'markersize',10);

        figure(2);
        p2(cnt)=plot(chiv*G,ksim,'.-','color',[col 0 1-col],...
            'linewidth',3,'markersize',20);
        ieps = ieps+1;
        cnt = cnt+1;
    end
end

figure(1);
xlabel('\chivG');ylabel('S^{-1}(q^*)');box on
legend(p1,{'N_M=0.05','N_M=0.50','N_M=5.00'},'location','northeast')
    %'Position',[0.72,0.25,0.1,0.1])
if LAM==0
    ylim([0,1]);xlim([0,8]);
    set(gca,'Ytick',0:0.2:1.0)
    set(gca,'YtickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})
elseif LAM==-0.75
    ylim([0,2]);xlim([0,20]);
    set(gca,'Ytick',0:0.4:2.0)
    set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0'})
else
    ylim([0,2]);xlim([0,20]);
    set(gca,'Ytick',0:0.4:2.0)
    set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0'})
end

figure(2);
xlabel('\chivG');ylabel('q^*R_M');box on
% legend(p2,'N_M=0.05','N_M=0.50','N_M=5.00')
if LAM==0
    ylim([0,4]);xlim([0,8]);
    set(gca,'Ytick',[0,1,2,3,4])
    set(gca,'YtickLabel',{'0','1','2','3','4'})
elseif LAM==-0.75
    ylim([2,5]);xlim([0,20]);
    set(gca,'Ytick',2:0.5:5)
    set(gca,'YtickLabel',{'2.0','2.5','3.0','3.5','4.0','4.5','5.0'})
else
    ylim([0,5]);xlim([0,20]);
    set(gca,'Ytick',0:1.0:5)
    set(gca,'YtickLabel',{'0','1','2','3','4','5'})
end

% end code

cd mkfigures/


figure(1)
savename = sprintf('../../results/randcopoly-results/random-simulation/ssim-lam%.2f.eps',LAM);
saveas(gcf,savename,'epsc')

figure(2)
savename = sprintf('../../results/randcopoly-results/random-simulation/qsim-lam%.2f.eps',LAM);
saveas(gcf,savename,'epsc')