function [x,sfit,resnorm]=saxsfit(qf,sf,x0,N,NM,FA,RM,LAM)
% % define fitting function
% % param(1) = scattering scale
% % param(2) = CHI
% % param(3) = LAM
% Sfun = @(param,qfit) param(1)./(-2*param(2)*CHIS+s2invwlc(N,NM,FA,param(3),qfit*RM)*NM);
% 
% % start fitting
% lb = [1  ,0.001,0];
% ub = [1e3,1.0  ,1];
% 
% options = optimset('Display','on',...
%     'TolX',1e-3,'TolFun',1e0,'MaxFunEvals',1e2,'MaxIter',1e2);
% x = lsqcurvefit(Sfun,x0,qf,sf,lb,ub,options);
% 
% % evaluate fitted structure factor intensity
% sfit = Sfun(x,qf);
% 
% % evaluate s2inv
% [~,SS]=kmaxwlc(N,NM,FA,LAM);
% CHIS=0.5*SS*NM;
% 
% Sfun = @(param,qfit) param(1)./(-2*param(2)*CHIS+s2invwlc(N,NM,FA,LAM,qfit*RM)*NM);
% 
% % start fitting
% lb = [1  ,0.001];
% ub = [1e3,1.0  ];
% 
% options = optimset('Display','on',...
%     'TolX',1e-3,'TolFun',1e0,'MaxFunEvals',1e2,'MaxIter',1e2);
% [x,resnorm] = lsqcurvefit(Sfun,x0,qf,sf,lb,ub,options);
% 
% % evaluate fitted structure factor intensity
% sfit = Sfun(x,qf);

% define fitting function
% param(1) = scattering scale
% param(2) = CHI(T)
% param(3) = LAM

% evaluate s2inv
[~,SS]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*SS*NM;

Sfun = @(param,qfit) param(1)./(-2*param(2)*CHIS+s2invwlc(N,NM,FA,LAM,qfit*RM)*NM);

% start fitting
lb = [1  ,0.001];
ub = [1e3,1.0  ];

options = optimset('Display','on',...
    'TolX',1e-3,'TolFun',1e0,'MaxFunEvals',1e2,'MaxIter',1e2);
[x,resnorm] = lsqcurvefit(Sfun,x0,qf,sf,lb,ub,options);

% evaluate fitted structure factor intensity
sfit = Sfun(x,qf);