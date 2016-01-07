function [val]=s2invgc(N,NM,FA,LAM,k)
% Calculate the s matrix

[GAMQ]=gammaq2(N,NM,LAM,k);

val=real(GAMQ/(2*FA*(1-FA)*NM));

function [val]=gammaq2(N,NM,LAM,k)

% Reset N to a row vector if entered as a column

if iscolumn(N)==1
    N=transpose(N);
end

val=zeros(length(k),length(N));

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

for j=1:length(k)
    if k(j)*sqrt(NM)<1e-2
        % zero wavemode limit
        GAMQ0=2*power(1+2/(N*(1-1/LAM))*((LAM-LAM^N)/(1-LAM)-N+1),-1);
        val(j,:)=1/GAMQ0*NM^2;
    else
        R=-k(j)^2/6;
        Z0=exp(R*NM);
        Z1=Z0*LAM;

        valeq=(Z0/R^2-1/R^2-NM/R);
        valne=(1/N)*2/R^2*Z1*(Z1.^N-N*Z1+N-1)/(1-Z1)^2*(cosh(R*NM)-1);

        valeq(isnan(valeq))=0;
        valne(isnan(valne))=0;

        valeq(isinf(valeq))=0;
        valne(isinf(valne))=0;

        val(j,:)=val(j,:)+(valeq+valne);
    end
end

val=NM^2*power(val,-1);