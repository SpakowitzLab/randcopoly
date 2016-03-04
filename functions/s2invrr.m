function val=s2invrr(N,NM,FA,LAM,k)
% Evaluate structure factor of random copolymer melt, at
% zero Flory-Huggins parameter. Chains are modeled as perfectly rigid
% rods.
% Usage :: val=s2invrr(N,NM,FA,LAM,k)
% Output :: val = inverse of structure factor
% Inputs ::
%   N = number of monomers
%   NM = number of Kuhn steps per monomer
%   FA = fraction of A monomers
%   LAM = degree of chemical correlation
%   k = wavevector, Fourier variable
% Andrew Spakowitz (4/14/15)
% Shifan Mao (1/6/16)

% Check conditions on chemical composition and correlation
if ~( FA>=0 && 1-FA>=0 && LAM>=-1 && 1-LAM>=0 && ...
     FA*(1-LAM)+LAM>=0 && FA*(LAM-1)+1>=0 )
 error('chemical composition and correlation constraints not satisfied')
end

[GAMQ]=gammaq2(N,NM,LAM,k);

val=real(GAMQ/(2*FA*(1-FA)*NM));

function [val]=gammaq2(N,NM,LAM,k)

% Reset N to a row vector if entered as a column

% if iscolumn(N)==1
%     N=transpose(N);
% end

val=zeros(length(k),length(N));

for j=1:length(k)
    if k(j)*NM<1e-2
        % zero wavemode limit
        GAMQ0=2*power(1+2/(N*(1-1/LAM))*((LAM-LAM^N)/(1-LAM)-N+1),-1);
        val(j,:)=1/GAMQ0*NM^2;
    else
        valeq=power(k(j),-2).*(-1+cos(k(j)*NM)+NM*k(j).*sinint(k(j)*NM));
        val(j,:)=val(j,:)+valeq;

        for J=2:N
            A=NM*(J-1)*k(j);
            valne=power(k(j),-2).*(-2*cos(A)+2*cos(k(j)*NM).*cos(A)...
                -2*A.*sinint(A)+(A-k(j)*NM).*sinint(A-k(j)*NM)+(A+k(j)*NM).*sinint(A+k(j)*NM));
            val(j,:)=val(j,:)+(N-J+1)/N*LAM^(J-1)*valne;
        end
    end
end
val=NM^2*power(val,-1);
