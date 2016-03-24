function val=s2invwlc(N,NM,FA,LAM,k,d,ORDmax,ORD,ResLayer)
% Evaluate structure factor of random copolymer melt, at
% zero Flory-Huggins parameter. Chains are modeled as wormlike chains.
% Usage :: val=s2invwlc(N,NM,FA,LAM,k,d,ORDmax,ORD,ResLayer)
% Output :: val = inverse of structure factor
% Inputs ::
%   N = number of monomers
%   NM = number of Kuhn steps per monomer
%   FA = fraction of A monomers
%   LAM = degree of chemical correlation
%   k = wavevector, Fourier variable
%   d = number of dimensions, default 3
%   ORDmax = maximum number of eigenvalues, default 20
%   ORD = number of eigenvalues, default 20
%   ResLayer = number of residual layers, default 500
% Andrew Spakowitz (4/14/15)
% Shifan Mao (1/6/16)

% Check conditions on chemical composition and correlation
if ~( FA>=0 && 1-FA>=0 && LAM>=-1 && 1-LAM>=0 && ...
     FA*(1-LAM)+LAM>=0 && FA*(LAM-1)+1>=0 )
 error('chemical composition and correlation constraints not satisfied')
end

% Fill in unset optional values

switch nargin
    case 5
        d=3;
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 6
        ORDmax=20;
        ORD=20;
        ResLayer=500;        
    case 7
        ORD=20;
        ResLayer=500;        
    case 8
        ResLayer=500;        
end

% If dimensions is 2, reset value to small perturbation above 2

if d==2
    d=2+1e-10;
end

[GAMQ]=gammaq2(N,NM,LAM,k,d,ORDmax,ORD,ResLayer);

val=real(GAMQ/(2*FA*(1-FA)*NM));

function [val]=gammaq2(N,NM,LAM,k,d,ORDmax,ORD,ResLayer)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Fill in unset optional values

switch nargin
    case 4
        d=3;
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 5
        ORDmax=20;
        ORD=20;
        ResLayer=500;        
    case 6
        ORD=20;
        ResLayer=500;        
    case 7
        ResLayer=500;        
end

% If dimensions is 2, reset value to small perturbation above 2

if d==2
    d=2+1e-10;
end

% Reset N to a row vector if entered as a column

if iscolumn(N)==1
    N=transpose(N);
end

val=zeros(length(k),length(N));

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

for j=1:length(k)
    if k(j)*sqrt(r2wlc(NM))<=1e-2
        % zero wavemode limit
        GAMQ0=2*power(1+2/(N*(1-1/LAM))*((LAM-LAM^N)/(1-LAM)-N+1),-1);
        val(j,:)=1/GAMQ0*NM^2;
    else
        % calculate the eigenvalues
        R=MatRoots(k(j),d,ORD);
        NR=ORDmax;

        [j0,dj0]=CalRes0(k(j),d,ResLayer);
        valeq=(NM/j0-dj0/j0^2);

        val(j,:)=val(j,:)+valeq;

        % get the residues for all roots of each k(j)
        Residue=CalRes(R(1:NR),k(j),d,ResLayer);

        for I=1:NR
            Z0=exp(R(I)*NM);
            Z1=Z0*LAM;
            valeq=Z0/R(I)^2;
            valne=(2*Z1./N).*(Z1.^N-N*Z1+N-1)/(1-Z1)^2*(cosh(R(I)*NM)-1)/R(I)^2;

            valeq(isnan(valeq))=0;
            valne(isnan(valne))=0;

            valeq(isinf(valeq))=0;
            valne(isinf(valne))=0;

            val(j,:)=val(j,:)+Residue(I)*(valeq+valne);
        end
    end
end

val=NM^2*power(val,-1);

function Eig=MatRoots(k,d,ORD)

% find roots of denominator (eigenvalues) by solving eigenvalue problem

Eig=zeros(ORD,1);

if k>8000

    % use large k asmyptotic expansion for large k
    
    for I=1:floor(ORD/2)
        l=I-1;
        alpha=1/sqrt(8*k);
        Eig(2*l+1)=1i*k-Epsilon(l,d,alpha);
        Eig(2*l+2)=conj(Eig(2*l+1));
    end
        
else

    % use matrix method for intermediate and small k regime
    n=4*ORD;
    E=zeros(n,n);
    for m=1:n
        if k<=1
            a=complex(0,-k*sqrt(m*(m+d-3)/(2*m+d-2)/(2*m+d-4)));            
            if m>1 
                b=complex(0,-k*sqrt((m-1)*((m-1)+d-3)/(2*(m-1)+d-2)/(2*(m-1)+d-4)));
            end
            if m==1
                E(m,1:2)=[(m-1)*(m+d-3),a];
            elseif m==n
                E(m,n-1:n)=[b,(m-1)*(m+d-3)];
            else
                E(m,m-1:m+1)=[b,(m-1)*(m+d-3),a];
            end
        else
            a=complex(0,-sqrt(m*(m+d-3)/(2*m+d-2)/(2*m+d-4)));            
            if m>1 
                b=complex(0,-sqrt((m-1)*((m-1)+d-3)/(2*(m-1)+d-2)/(2*(m-1)+d-4)));
            end
            if m==1
                E(m,1:2)=[(m-1)*(m+d-3)/k,a];
            elseif m==n
                E(m,n-1:n)=[b,(m-1)*(m+d-3)/k];
            else
                E(m,m-1:m+1)=[b,(m-1)*(m+d-3)/k,a];
            end
        end
    end
    TempMat=eig(E);
    [~,index]=sort(real(TempMat));
    TempMat=TempMat(index);
    if k<=1
        Eig=-TempMat(1:ORD);
    else
        Eig=-TempMat(1:ORD)*k;
    end
end

function value=Epsilon(l,d,alpha)

% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=(d-3)/2;
n=2*l+m+1;

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n^2+3-3*m^2)-m*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n*(n^2+3-9*m^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n^4+34*n^2+9)-(102*n^2+42)*m^2+33*m^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n*(33*n^4+410*n^2+405)-(1230*n^2+1722)*m^2+813*m^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n^6+140*n^4+327*n^2+54-(420*n^4+1350*n^2+286)*m^2+(495*n^2+314)*m^4-82*m^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;
  
function Res=CalRes(R,k,d,ResLayer)

% calculate the residual of all eigenvalues

ResThreshold=1e-12;     % threshold to go from small k asymptot to matrix method
ImagThreshold=1e-8;     % threshold to go from small k asymptot to matrix method
NR=length(R);           % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative
Res=zeros(NR,1);    % residual vector, each row corresponds to each eigenvalue

% find the residual
for n=1:NR

    % use asymptotic residual to figure out whether calculate the residual
    % using continued fraction or stay with the asymptotic form for small k
    % limit

    Res(n)=SmallAsympRes(k,n,d);
    if abs(Res(n))>ResThreshold
        if k<=1

            % residual using continued fraction for small k
            p=R(n);
            W=p+(ResLayer+d-2)*ResLayer;
            Wprime=1;
            for L=ResLayer:-1:1
                AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
                Wprime=1-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm+AL^2/W;
            end
            Res(n)=1/Wprime;
        else

            % residual using continued fraction for large k
            p=R(n);
            W=(p+(ResLayer+d-2)*ResLayer)/k;
            Wprime=1/k;
            for L=ResLayer:-1:1
                AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
                Wprime=1/k-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm/k+AL^2/W;
            end
            Res(n)=1/(k*Wprime);
        end
    end
    if abs(imag(R(n)))<ImagThreshold
        Res(n)=real(Res(n));
    end
end

function Res=SmallAsympRes(K,n,d)

% calculate the residual using small k asymptot

l=n-1;
Res=1;
Wl=l*(l+d-2);

for j=0:(l-1)
    Wj=j*(j+d-2);
    ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
    Res=Res*ajp1^2/(Wl-Wj)^2;
end

Res=Res*K^(2*l)*(-1)^l;

function [j0,dj0]=CalRes0(k,d,ResLayer)

if k<=1
    
    % residual using continued fraction for small k
    p=0;
    W=p+(ResLayer+d-2)*ResLayer;
    Wprime=1;
    for L=ResLayer:-1:1
        AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
        Wprime=1-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm+AL^2/W;
    end
    dj0=Wprime;
    j0=W;
else
    
    % residual using continued fraction for large k
    p=0;
    W=(p+(ResLayer+d-2)*ResLayer)/k;
    Wprime=1/k;
    for L=ResLayer:-1:1
        AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
        Wprime=1/k-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm/k+AL^2/W;
    end
    dj0=k*Wprime;
    j0=k*W;
end
