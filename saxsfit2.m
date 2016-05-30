function F = saxsfit2(x)

global sf qf N rm FA NFIT

% function to minimize
% mdl = @(x,qf) x(1)./(-2*x(3+IT)+s2invwlc(N,NM,FA,x(2),qf*rm)*x(3))';
NQ = length(qf);
F = zeros(NFIT*NQ,1);

for IT = 1:NFIT
    % x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT
    fun = x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm))-sf(:,IT);
    F(1+(IT-1)*NQ:IT*NQ) = fun;
end

end