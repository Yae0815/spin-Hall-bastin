function [Sx,Sy,Sz] = spinops()
% Spin-1/2 operators embedded in 4×4 basis (orbital ⊗ spin)
% Orbital Pauli on left (2×2), spin Pauli on right (2×2).
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1]; s0 = eye(2);
Sx = kron(eye(2), sx);
Sy = kron(eye(2), sy);
Sz = kron(eye(2), sz);
end