function S = sequence(FA,LAM,N)
% generate chemical sequence of random copolymers
% or user-defined sequence
% note the sequence output is not deterministic
% 0 == A
% 1 == B

PAA=FA*(1.-LAM)+LAM;
PBB=FA*(LAM-1.)+1.;
PAB=1.-PBB;

S = zeros(1,N);
TEST = rand();
if (TEST<FA) 
    S(1)=0;
else
    S(1)=1;
end

for ii = 2:N 
    TEST = rand();
    if (S(ii-1)==0)
        if (TEST<PAA)
              S(ii)=0;
        else
              S(ii)=1;
        end
    else
        if (TEST<PAB)
              S(ii)=0;
           else
              S(ii)=1;
        end
    end
end