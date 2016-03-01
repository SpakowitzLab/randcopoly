figure(1);set(gca,'fontsize',20);hold
figure(2);set(gca,'fontsize',20);hold
figure(3);set(gca,'fontsize',20);hold

% WLC
NMV = logspace(-2,2,11);
cnt = 1;
for NM = NMV
    col = (cnt-1)/(length(NMV)-1);
    
    data = load(sprintf('data/WLC_NM%.2f',NM));
    LAMV = data(:,1);
    KS = data(:,2);
    CHIS = data(:,3);
    D2S = data(:,4);
    
    figure(1);plot(LAMV,KS,'linewidth',3,'color',[col 0 1-col])
    figure(2);plot(LAMV,CHIS,'linewidth',3,'color',[col 0 1-col])
    
    % Find Lifshitz point
    IND = find(KS<=1e-1);
    LAML_WLC = LAMV(IND(1));
    figure(3);plot(LAML_WLC-LAMV,D2S,'linewidth',2,'color',[col 0 1-col])
    
    cnt = cnt+1;
end

% Gaussian chain
data = load(sprintf('data/GC'));
LAMV = data(:,1);
KS = data(:,2);
CHIS = data(:,3);
D2S = data(:,4);

figure(1);plot(LAMV,KS,'linewidth',3,'color',[0 0 0])
figure(2);plot(LAMV,CHIS,'linewidth',3,'color',[0 0 0])

% Find Lifshitz point
IND = find(KS<=1e-1);
LAML_GC = LAMV(IND(1));
figure(3);plot(LAML_GC-LAMV,D2S,'linewidth',3,'color',[0 0 0])

% Rigid rod
data = load(sprintf('data/RR'));
LAMV = data(:,1);
KS = data(:,2);
CHIS = data(:,3);
D2S = data(:,4);

figure(1);plot(LAMV,KS,'--','linewidth',3,'color',[0 0 0])
figure(2);plot(LAMV,CHIS,'--','linewidth',3,'color',[0 0 0])

% Find Lifshitz point
IND = find(KS<=1e-1);
LAML_RR = -0.1;
figure(3);plot(LAML_RR-LAMV,D2S,'--','linewidth',3,'color',[0 0 0])

figure(1);xlabel('\lambda');ylabel('R_Mq^*');box on
figure(2);xlabel('\lambda');ylabel('\chi_Sv');box on
figure(3);set(gca,'yscale','log');box on
xlabel('\lambda_L-\lambda');ylabel('Peak sharpness \Delta_\psi')
ylim([1e-2,1e3])

figure(3);saveas(gcf,'peaksharpness','epsc')