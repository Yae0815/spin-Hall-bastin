function out = bastin_main(params)
% σ^{s_γ}_{αβ} via Bastin (f_n - f_m) on a 3D k-grid

    % ---------- unpack ----------
    model = params.model;
    Nk    = params.Nk;      eta  = params.eta;   Ef = params.Ef;
    mu    = params.mu;      T    = params.T;     e  = params.e;
    hbar  = params.hbar;

    alpha = lower(params.alpha);
    beta  = lower(params.beta);
    gamma = lower(params.gamma);

    if isfield(params,'shift'),      shift = params.shift;      else, shift = [0 0 0]; end
    if isfield(params,'units'),  unit_mode = params.units;      else, unit_mode = 'lattice'; end
    if isfield(params,'use_parfor'), use_parfor = params.use_parfor; else, use_parfor = true; end

    assert(model.Norb==4,'This template assumes Norb=4.');

    % ---------- operators ----------
    Sg = pick_spin(model, gamma);
    [vifun, vjfun] = pick_velocity_funs(model, alpha, beta, hbar); % function handles v_i(k), v_j(k)

    % ---------- k grid ----------
    [kxs, kys, kzs, wk] = utils.kgrid3(Nk, shift);
    Ktot = numel(kxs);

    % ---------- Bastin sum over k ----------
    if use_parfor
        parts = zeros(Ktot,1);
        parfor t = 1:Ktot
            k  = [kxs(t),kys(t),kzs(t)];
            parts(t) = local_k_contrib_static( ...
                model, Sg, vifun, vjfun, k, Ef, mu, T, eta, wk(t) );
        end
        sigma_sum = sum(parts);
    else
        sigma_sum = 0.0;
        for t = 1:Ktot
            k  = [kxs(t),kys(t),kzs(t)];
            sigma_sum = sigma_sum + local_k_contrib_static( ...
                model, Sg, vifun, vjfun, k, Ef, mu, T, eta, wk(t) );
        end
    end

    % ---------- prefactor & units ----------
    % Uniform grid in [-π,π): Δk=2π/Nk ⇒ ∫ d^3k/(2π)^3 ≈ (Δk^3/(2π)^3)Σ = (1/Nk^3)Σ
    BZfac = 1/(Nk^3);
    val = (e/hbar) * BZfac * sigma_sum;

    if strcmpi(unit_mode,'SI')
        if isfield(model,'a'), a = model.a; else, a = 1; end
        val = bastin.units.to_SI(val, a, e, hbar); % 目前 pass-through，之後可在 units.m 實作
    end

    out.sigma = val;
    out.meta  = struct('Nk',Nk,'eta',eta,'Ef',Ef,'mu',mu,'T',T, ...
                       'alpha',alpha,'beta',beta,'gamma',gamma,'units',unit_mode);
end

% ================== subfunctions (非巢狀，可被 parfor 呼叫) ==================

function val = local_k_contrib_static(model, Sg, vifun, vjfun, k, Ef, mu, T, eta, w)
    % H(k) & eigen
    H  = model.H(k);
    [V,D] = eig(H,'vector');
    E = D(:);

    % Fermi-Dirac at μ_eff = Ef+mu
    fn = bastin.fermi(E - (Ef+mu), T);

    % k-dependent velocities & spin current
    vi  = vifun(k);
    vj  = vjfun(k);
    Jsi = 0.5*(vi*Sg + Sg*vi);

    % band basis
    Jb = V' * Jsi * V;
    Vb = V' * vj  * V;

    % Bastin sum over n≠m
    N = numel(E);
    acc = 0.0;
    for n = 1:N
        for m = 1:N
            if m==n, continue; end
            dEmn  = E(n) - E(m);
            denom = dEmn^2 + eta^2;
            acc = acc + (fn(n) - fn(m)) * imag( Jb(n,m) * Vb(m,n) ) / denom;
        end
    end
    val = acc * w;
end

function S = pick_spin(model, comp)
    switch comp
        case 'x', S = model.Sx;
        case 'y', S = model.Sy;
        case 'z', S = model.Sz;
        otherwise, error('gamma must be x/y/z');
    end
end

function [vifun, vjfun] = pick_velocity_funs(model, alpha, beta, hbar)
    switch alpha
        case 'x', vifun = @(k) (1/hbar)*model.dHdkx(k);
        case 'y', vifun = @(k) (1/hbar)*model.dHdky(k);
        case 'z', vifun = @(k) (1/hbar)*model.dHdkz(k);
        otherwise, error('alpha must be x/y/z');
    end
    switch beta
        case 'x', vjfun = @(k) (1/hbar)*model.dHdkx(k);
        case 'y', vjfun = @(k) (1/hbar)*model.dHdky(k);
        case 'z', vjfun = @(k) (1/hbar)*model.dHdkz(k);
        otherwise, error('beta must be x/y/z');
    end
end
