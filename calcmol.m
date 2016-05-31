function [mol,rm] = calcmol(MWPEG,WTPERC)
% calculates mol percent of hydrophilic segments in
% ODPA-AP6F-PEG random copolymer system
% MWPEG = molecular weight of PEG in g/mol
% WTPERC = weight percent of PEG

mol = 608*WTPERC/(MWPEG - (MWPEG+306)*WTPERC+ 608*WTPERC);
if MWPEG == 1500
    RMB = 44.6;
elseif MWPEG == 900
    RMB = 39;
end

rm  = 30*(1-mol)+RMB*mol;