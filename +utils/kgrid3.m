function [kx,ky,kz,w] = kgrid3(Nk, shift)
% Uniform 3D Monkhorst–Pack in [-π,π)
if nargin<2, shift=[0 0 0]; end
base = (0:Nk-1)/Nk - 0.5 + shift(:)';
[kx3,ky3,kz3] = ndgrid(2*pi*base, 2*pi*base, 2*pi*base);
kx = kx3(:); ky = ky3(:); kz = kz3(:);
w = ones(numel(kx),1);
end