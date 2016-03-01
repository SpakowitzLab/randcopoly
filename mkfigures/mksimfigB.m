close all
clear;

% start code
LAMV_SIM = [-0.75:0.25:0.25];
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
Lbox=20;
Lbin=1;

figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',20)

ieps = 1;
for EPS = EPSV
    col = (ieps-1)/(length(EPSV)-1);

    % plot mean-field theory results
    NM=G*EPS;  % number of Kuhn steps per monomer
    data = load(sprintf('data/WLC_NM%.2f',NM));
    LAMV_MF = data(:,1);
    KS_MF = data(:,2);
    D2S_MF = data(:,4);
    
    figure(1);plot(LAMV_MF,KS_MF,'--','linewidth',3,'color',[col 0 1-col])
    figure(2);plot(LAMV_MF,D2S_MF,'--','linewidth',3,'color',[col 0 1-col])
    ieps = ieps+1;
end


ieps = 1;
for EPS = EPSV
    col = (ieps-1)/(length(EPSV)-1);
    
    % plot simulation results
    KSV_SIM = zeros(length(LAMV_SIM),46);
    D2SV_SIM = zeros(length(LAMV_SIM),46);
    ilam=0;
    for LAM = LAMV_SIM
        ilam = ilam+1;
        [~,~,~,KS_SIM,~,~,~,D2S_SIM]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        KSV_SIM(ilam,:) = KS_SIM;
        D2SV_SIM(ilam,:) = D2S_SIM;
    end
    
    figure(1);plot(LAMV_SIM,KSV_SIM(:,42),'o',...
        'MarkerEdgeColor',[col 0 1-col],'markersize',12,'linewidth',2)
    figure(2);plot(LAMV_SIM,D2SV_SIM(:,42),'o',...
        'MarkerEdgeColor',[col 0 1-col],'markersize',12,'linewidth',2)
    ieps = ieps+1;
end

figure(1);
xlim([-1,.5]);ylim([0,5])
xlabel('\lambda');ylabel('R_Mq^*');box on
set(gca,'Xtick',-1:0.25:0.5)
set(gca,'XtickLabel',{'-1','-0.75','-0.5','-0.25','0','0.25'})

figure(2);
xlabel('\lambda');
ylabel('Peak sharpness \Delta_\psi');
box on
xlim([-1,.5]);
ylim([0,10]);
set(gca,'yscale','linear')
% set(gca,'yscale','log');
% ylim([1e-2,1e2])
set(gca,'Xtick',-1:0.25:0.5)
set(gca,'XtickLabel',{'-1','-0.75','-0.5','-0.25','0','0.25'})

% end code

figure(1);
savename = sprintf('../../results/randcopoly-results/random-simulation/qstar.eps');
saveas(gcf,savename,'epsc')

figure(2);
savename = sprintf('../../results/randcopoly-results/random-simulation/d2s.eps');
saveas(gcf,savename,'epsc')