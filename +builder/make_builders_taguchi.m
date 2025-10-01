function model = make_builders_taguchi(p)
% make_builders_taguchi  Taguchi (2020) DSM/TDSM lattice model (4x4)
% Inputs (struct p):
%   eta, txy, tz, M, beta, gamma, a (optional lattice const, default 1)
%
% Returns struct model with fields:
%   Norb=4, H(k), dHdkx(k), dHdky(k), dHdkz(k), Sx,Sy,Sz, a

    eta   = p.eta;
    txy   = p.txy;
    tz    = p.tz;
    M     = p.M;
    beta  = p.beta;
    gamma = p.gamma;
    if isfield(p,'a'), a = p.a; else, a = 1; end

    % Pauli (orbital on left, spin on right)
    sx = [0 1; 1 0];
    sy = [0 -1i; 1i 0];
    sz = [1 0; 0 -1];
    s0 = eye(2);

    % Spin operators embedded (I_orbital ⊗ s_i)
    Sx = kron(eye(2), sx);
    Sy = kron(eye(2), sy);
    Sz = kron(eye(2), sz);

    K  = @(A,B) kron(A,B);

    % ----- a_i(k) scalars (all ASCII) -----
    a1 = @(k)  eta * sin(k(1));                                   % sigma_x ⊗ s_z
    a2 = @(k) -eta * sin(k(2));                                   % sigma_y ⊗ s_0
    a3 = @(k)  M - txy*(cos(k(1))+cos(k(2))) - tz*cos(k(3));      % sigma_z ⊗ s_0
    a4 = @(k) (beta + gamma) * sin(k(3)) * (cos(k(2)) - cos(k(1)));         % sigma_x ⊗ s_x
    a5 = @(k) -(beta - gamma) * sin(k(3)) * sin(k(2)) * sin(k(1));         % sigma_x ⊗ s_y

    % Hamiltonian
    H = @(k) ( a1(k) * K(sx,sz) + ...
               a2(k) * K(sy,s0) + ...
               a3(k) * K(sz,s0) + ...
               a4(k) * K(sx,sx) + ...
               a5(k) * K(sx,sy) );

    % ----- analytic derivatives of a_i(k) -----
    da1dkx = @(k)  eta * cos(k(1));
    da2dky = @(k) -eta * cos(k(2));
    da3dkx = @(k)  txy * sin(k(1));
    da3dky = @(k)  txy * sin(k(2));
    da3dkz = @(k)  tz  * sin(k(3));

    % a4 = (beta+gamma) sin kz (cos ky - cos kx)
    da4dkx = @(k)  (beta + gamma) * sin(k(3)) * sin(k(1));
    da4dky = @(k) -(beta + gamma) * sin(k(3)) * sin(k(2));
    da4dkz = @(k)  (beta + gamma) * cos(k(3)) * (cos(k(2)) - cos(k(1)));

    % a5 = -(beta-gamma) sin kz sin ky sin kx
    da5dkx = @(k) -(beta - gamma) * sin(k(3)) * sin(k(2)) * cos(k(1));
    da5dky = @(k) -(beta - gamma) * sin(k(3)) * cos(k(2)) * sin(k(1));
    da5dkz = @(k) -(beta - gamma) * cos(k(3)) * sin(k(2)) * sin(k(1));

    % dH/dk_i
    Hkx = @(k) ( da1dkx(k) * K(sx,sz) + ...
                 da3dkx(k) * K(sz,s0) + ...
                 da4dkx(k) * K(sx,sx) + ...
                 da5dkx(k) * K(sx,sy) );

    Hky = @(k) ( da2dky(k) * K(sy,s0) + ...
                 da3dky(k) * K(sz,s0) + ...
                 da4dky(k) * K(sx,sx) + ...
                 da5dky(k) * K(sx,sy) );

    Hkz = @(k) ( da3dkz(k) * K(sz,s0) + ...
                 da4dkz(k) * K(sx,sx) + ...
                 da5dkz(k) * K(sx,sy) );

    % package
    model.Norb  = 4;
    model.H     = @(k) H(k);
    model.dHdkx = @(k) Hkx(k);
    model.dHdky = @(k) Hky(k);
    model.dHdkz = @(k) Hkz(k);
    model.Sx    = Sx;  model.Sy = Sy;  model.Sz = Sz;
    model.a     = a;
end
