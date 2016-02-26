function [dum]=plotsim2(EPS,LAM)
% Plots single chain correlations in Monte-Carlo simulation
% INPUTS::
%   EPS, number of Kuhn steps per monomer
%   LAM, degree of chemical correlation

% simulation folder
folder = '../../results/randcopoly-results/scalcbatch-12-15-15';
addpath('../functions/')

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

R2=-0.5+0.5*exp(-2*NM)+NM;

folder = '../../results/randcopoly-results/Pdata-12-15-15';
cnt=find(chiind==plotind(1));
for ii = plotind
    CHI = chiv(ii);
    col = cnt/length(chiind);

    % Plot simulation results
    filename = sprintf([folder,'/PMC_SIM%dCHEM%dCHI%.8f'],SIMNUM,CHEMNUM,CHI*G);
    R = load(filename);
    plot(R(:,1),R(:,2),'.-','color',[1-col 0 col],'linewidth',3,'markersize',20)
    cnt = cnt+1;
end

xlabel('R_Mq');
ylabel('<\psi(q)\psi(-q)>')
axis([0,1,0,10])
box on

%     savename = sprintf('../results/randcopoly-results/structure-figures/sfig-eps%.2f-lam%.2f.eps',EPS,LAM);
%     saveas(gcf,savename,'epsc')

L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);
L=(M*G-1)*L0;
LP=L0/(2*EPS);
N=M*G-1;
FNUM=round(100*L/(2*LP));

XG=transpose(linspace(0,1,5001));
if FNUM<=2000
    G=load(['../../results/data//out',int2str(FNUM),'.txt']);
else
    DXG=XG(2)-XG(1);
    G=exp(-1.5*N*XG.^2).*XG.*XG;
    G=G./(sum(G).*DXG);
end

if FNUM<=2000
    plot(XG,G.*power(XG,2)*(4*pi),'k-','LineWidth',2)
else
    plot(XG,G,'k-','LineWidth',2)
end