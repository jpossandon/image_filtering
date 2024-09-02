function [U, V] = dftuv(M,N)
%DFTUV Computes meshgrid frequency matrices
%   [U,V] = DFTUV(M,N) computes meshgrid frequency matrices U and V. U and
%   V are useful computing frequency domain filter functions that can be
%   uses with DFTFILT. U and V are both M by N%   
%   23/07/07 copied from GONZALEZ RC. et al, Digital Image Processing Using
%   MATLAB (2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up range of variables
u = 0:(M-1);
v = 0:(N-1);

% Compute the indice for use in meshgrid
idx = find(u>M/2);
u(idx) = u(idx) - M;
idy = find(v>N/2);
v(idy) = v(idy) - N;

% Compute the meshgrid arrays
[V,U] = meshgrid(v,u);