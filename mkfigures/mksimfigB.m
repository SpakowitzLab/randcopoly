close all
clear;
cd ../

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
    FA=0.5;    % equal chemical composition
    RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

    LAMV_MF = linspace(-1,.5,51);
    KS_MF = zeros(length(LAMV_MF),1);
    D2S_MF = zeros(length(LAMV_MF),1);
    for ii = 1:length(LAMV_MF)
    	LAM = LAMV_MF(ii);
        [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
        KS_MF(ii)=kval;
        D2S_MF(ii) = -1/(sval^2*RM^2*EPS)*d2gam2;
    end
    figure(1);plot(LAMV_MF,RM*KS_MF,'--','linewidth',3,'color',[col 0 1-col])
    figure(2);plot(LAMV_MF,-D2S_MF/G,'--','linewidth',3,'color',[col 0 1-col])
    
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

    plotrange = 42;
    figure(1);plot(LAMV_SIM,KSV_SIM(:,plotrange),'o',...
        'MarkerFaceColor',[col 0 1-col],'MarkerEdgeColor',[col 0 1-col],...
        'markersize',15)
    
    plotrange = 1;
    figure(2);plot(LAMV_SIM,-D2SV_SIM(:,plotrange)/G,'o',...
        'MarkerFaceColor',[col 0 1-col],'MarkerEdgeColor',[col 0 1-col],...
        'markersize',15)
    
    ieps = ieps+1;
end
figure(1);xlabel('\lambda');ylabel('R_Mq^*')
figure(2);
xlabel('\lambda');ylabel('\Delta_\psi');set(gca,'yscale','log')
ylim([1e-2,1e5])
    
% end code

cd mkfigures/

figure(1);
savename = sprintf('../../results/randcopoly-results/random-simulation/qstar.eps');
saveas(gcf,savename,'epsc')

figure(2);
savename = sprintf('../../results/randcopoly-results/random-simulation/d2s.eps');
saveas(gcf,savename,'epsc')