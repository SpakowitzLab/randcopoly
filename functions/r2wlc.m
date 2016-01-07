function [r2]=r2wlc(N)
% Mean-squared end-to-end distance of wormlike chain
% Usage :: [r2]=r2wlc(N)
% Output :: r2 = mean-squared end-to-end distance of wormlike chains
% Input :: N = number of Kuhn steps of wormlike chain

r2=-0.5+0.5*exp(-2*N)+N;
