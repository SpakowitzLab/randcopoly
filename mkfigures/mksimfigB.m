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
    
    if EPS ==1.00
        plotrange = 2;
        figure(1);plot(LAMV_SIM(1:2),KSV_SIM(1:2,plotrange),'o',...
            'MarkerEdgeColor',[col 0 1-col],...
            'markersize',15,'linewidth',2)
        figure(2);plot(LAMV_SIM(1:2),-D2SV_SIM(1:2,plotrange)/G,'o',...
            'MarkerEdgeColor',[col 0 1-col],...
            'markersize',15,'linewidth',2)
    else
        plotrange = 2;
        figure(1);plot(LAMV_SIM(1:3),KSV_SIM(1:3,plotrange),'o',...
            'MarkerEdgeColor',[col 0 1-col],...
            'markersize',15,'linewidth',2)
        figure(2);plot(LAMV_SIM(1:3),-D2SV_SIM(1:3,plotrange)'/G,'o',...
            'MarkerEdgeColor',[col 0 1-col],...
            'markersize',15,'linewidth',2)
    end
    
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
    
    plotrange = 41;
    figure(1);plot(LAMV_SIM,KSV_SIM(:,plotrange),'o',...
        'MarkerFaceColor',[col 0 1-col],'MarkerEdgeColor',[col 0 1-col],...
        'markersize',15)
    figure(2);plot(LAMV_SIM,-D2SV_SIM(:,plotrange)/G,'o',...
        'MarkerFaceColor',[col 0 1-col],'MarkerEdgeColor',[col 0 1-col],...
        'markersize',15)
    
    ieps = ieps+1;
end

figure(1);
xlim([-1,.5])
xlabel('\lambda');ylabel('R_Mq^*');box on

figure(2);
xlabel('\lambda');
ylabel('Peak sharpness \Delta_\psi');
set(gca,'yscale','log');box on
xlim([-1,.5]);ylim([1e-3,1e4])

% end code

figure(1);
savename = sprintf('../../results/randcopoly-results/random-simulation/qstar.eps');
saveas(gcf,savename,'epsc')

figure(2);
savename = sprintf('../../results/randcopoly-results/random-simulation/d2s.eps');
saveas(gcf,savename,'epsc')