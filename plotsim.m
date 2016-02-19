function [chis,chiv,KS_MF,KS_SIM,SINV_MF,SINV_SIM,D2S_MF,D2S_SIM]=plotsim(EPS,LAM,PLOTON)
%% Plots density-density correlation of sim. and theory
% INPUTS::
%   EPS, number of Kuhn steps per monomer
%   LAM, degree of chemical correlation

% simulation folder
folder = '../results/randcopoly-results/scalcbatch-12-15-15';
addpath('functions/')

% simulation constants
FA=0.5;  % fraction of A blocks
M=8;        % number of blocks
G=5;  % number of discrete monomers
NM=G*EPS;  % number of Kuhn steps per monomer

% range of chi params.
chiind=fliplr([1:2:11,21,31,41]);
plotind=chiind(1:end);

% load simulation parameters
simparam=load([folder,'/chivals']);
% find corresponding simulation at EPS and LAM
ind = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
% finding simulation index(indices)
chemparam=load([folder,'/chemind']);
SIMNUM = ind;
CHEMNUM=chemparam(chemparam(:,1)==ind,2);
% Flory-Huggins parameter
chiv = load(sprintf([folder,'/sdata-%d-%d/Sdata/chilist'],SIMNUM,CHEMNUM));
chiv = chiv/G;

% Find spinodal CHI and critical k
[KS_MF,sval,d2gam2]=kmaxwlc(M,NM,FA,LAM);
chis=0.5*sval;

R2=-0.5+0.5*exp(-2*NM)+NM;
k=logspace(log10(0.6),log10(12),100)./sqrt(R2); % wavevectors

% Plot: density-density correlation
if PLOTON==1
    f1=figure;hold;set(gca,'fontsize',30)
    set(f1,'position',[0,0,800,600])
    cnt=find(chiind==plotind(1));
    for ii = plotind
        CHI = chiv(ii);
        col = cnt/length(chiind);

        if CHI<chis*EPS*0.85
            % Plot analytical theory (note conversion of monomer volume v, factor of EPS)
            val = s2invwlc(M,NM,FA,LAM,k);
            plot(k.*sqrt(R2),1./(-2*CHI+EPS*val),'--','color',[1-col 0 col],'linewidth',3)
        end

        % Plot simulation results
        filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
            SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
        S = load(filename);
        plot(S(:,1),S(:,2),'.-','color',[1-col 0 col],'linewidth',3,'markersize',20)
        cnt = cnt+1;
    end

    xlabel('R_Mq');
    ylabel('<\psi(q)\psi(-q)>')
    xlim([.6,12]);ylim([1e-1,2e2]);
    set(gca,'Ytick',[1e-1 1e0 1e1 1e2])
    set(gca,'Xtick',[1,10])
    set(gca,'xscale','log');set(gca,'yscale','log')
    box on
    
    savename = sprintf('../results/randcopoly-results/structure-figures/sfig-eps%.2f-lam%.2f.eps',EPS,LAM);
    saveas(gcf,savename,'epsc')
end

% Find peak of structure factors
KS_SIM = zeros(length(chiv),1);
SINV_MF = zeros(length(chiv),1);
SINV_SIM = zeros(length(chiv),1);
D2S_SIM = zeros(length(chiv),1);
D2S_MF = zeros(length(chiv),1);

for ii = 1:length(chiv)
    CHI = chiv(ii);
    
    % Find peak position
    SINV_MF(ii)=-2*CHI+EPS*sval;

    % Plot simulation results
    filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
        SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
    S = load(filename);
    
    ind = find(S(:,2)==max(S(:,2)));IND = ind(1);
    KS_SIM(ii) = S(IND,1);
    SINV_SIM(ii) = 1./S(IND,2);
    
    if IND>2  % central differences
%         DK1 = S(IND,1)-S(IND-1,1);
%         DK2 = S(IND+1,1)-S(IND,1);
%         SP1 = S(IND-1,2);
%         SP2 = S(IND,2);
%         SP3 = S(IND+1,2);
%         D2S_SIM(ii) =
%         (2/(DK1+DK2)/DK2^2)*(DK1*SP3-(DK1+DK2)*SP2+DK2*SP1);
        Kfit = S(IND-2:IND+2,1);
        Sfit = S(IND-2:IND+2,2);
    else  % forward differences
        Kfit = S(IND+1:IND+3,1);
        Sfit = S(IND+1:IND+3,2);
    end
    fit = polyfit(Kfit,Sfit,2);
    SINV_SIM(ii) = 1./polyval(fit,-fit(2)/(2*fit(1)));
    KS_SIM(ii) = -fit(2)/(2*fit(1));
    D2S_SIM(ii) = 2*fit(1);
    D2S_MF(ii) = -1/(sval^2*R2*EPS)*d2gam2;
end