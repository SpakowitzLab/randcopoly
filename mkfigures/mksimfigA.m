% plot structure factors in MF theory and simulations

LAMV = -0.75:0.25:0.25;
EPSV = [0.01,0.10,1.00];
PLOTON = 1;

for LAM = LAMV
    for EPS = EPSV
        plotsim(EPS,LAM,PLOTON);
    end
end