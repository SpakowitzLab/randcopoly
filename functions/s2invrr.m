function [val]=s2invrr(N,NM,FA,LAM,k)

% % Calculate the Fourier transform of the Green function
% % for the wormlike chain in d-dimensions
% %
% % Andrew Spakowitz (4/14/15)
% 
% % Calculate the s matrix
% 
% [SAA,SAB,SBA,SBB]=s2rr(N,NM,FA,LAM,k);
% 
% DET=SAA.*SBB-SAB.*SBA;
% 
% val=real(N*NM*(SAA+SBB+2*SAB)./DET);

[GAMQ]=gammaq2(N,NM,LAM,k);

val=real(GAMQ/(2*FA*(1-FA)*NM));

function [val]=gammaq2(N,NM,LAM,k)

% Reset N to a row vector if entered as a column

if iscolumn(N)==1
    N=transpose(N);
end

val=zeros(length(k),length(N));

for j=1:length(k)
    if k(j)*sqrt(NM)<1e-2
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