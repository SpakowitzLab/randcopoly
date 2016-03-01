


%test indbin --> ix,iy,iz transformation
nbinx = 20;
nbiny = 20;
nbinz = 20;

nbin = nbinx*nbiny*nbinz;

indbin = 1:nbin;  %bin indices
indxyz = [];

for ii = 1:length(indbin)
    ix = mod(indbin(ii),nbinx);
    if ix == 0
        ix = nbinx;
    end

    iz = ceil(indbin(ii)/nbinx/nbiny);
    iy = (indbin(ii)-ix-(iz-1).*nbinx.*nbiny)./nbinx + 1;
    indxyz = [indxyz;ix,iy,iz];
end