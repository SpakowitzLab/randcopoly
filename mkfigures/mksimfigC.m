clear;close all

% Plot SINV

% start code
LAMV = [-0.75,0.00];
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
RM=2;
Lbox=20;
Lbin=1;

% plot range
ind = 1:2:41;

% plot mean-field theory results
for LAM = LAMV
    figure;hold;set(gca,'fontsize',20)

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
        plot(chiv*G,SINV_MF,'--','linewidth',3,'color',[col 0 1-col])
        ieps = ieps+1;
    end

% plot simulation results
    cnt=1;p1=[];
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
        p1(cnt)=plot(chiv*G,SINV_SIM,'o',...
            'MarkerEdgeColor',[col 0 1-col],'linewidth',3,'markersize',10,'linewidth',2);
        ieps = ieps+1;
        cnt = cnt+1;
    end

    xlabel('\chivG');ylabel('S^{-1}(q^*)');box on
    legend(p1,{'N_M=0.05','N_M=0.50','N_M=5.00'},'location','northeast')
        %'Position',[0.72,0.25,0.1,0.1])
    if LAM==0
        ylim([0,1.]);xlim([0,8]);
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

    % end code
    savename = sprintf('../../results/randcopoly-results/random-simulation/ssim-lam%.2f.eps',LAM);
    saveas(gcf,savename,'epsc')

end