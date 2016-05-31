function F = ozfit(x)

global sf qf NFIT

NQ = length(qf);
F = zeros(NFIT*NQ,1);

for IT = 1:NFIT
    % x(1) = intensity scale, x(2) = xi
    
    % function to minimize
    % mdl = @(x,qf) x(1)./(1+qf.^2*x(2)^2);
    fun = x(1)*x(1+IT)^2./(1+qf.^2*x(1+IT)^2)-sf(:,IT);
    F(1+(IT-1)*NQ:IT*NQ) = fun;
end

end