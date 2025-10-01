% Scan M/tz vs σ^{s_z}_{xy} at fixed txy/tz (Bastin)
clear; clc;

addpath(pwd); addpath(genpath(pwd));
%% model knobs (Taguchi lattice model)
eta = 0.89; % eV
beta = 0.67; gamma = 0.335; % example SOC (β/tz=2γ/tz=0.67 in paper → tune with tz)


Nk = 21; % odd preferred
eta_b = 5e-3; % broadening (eV)
Ef = 0.0; % reference
mu = 0.0; T = 0; % T=0 clean limit


alpha='x'; beta_dir='y'; gamma_pol='z'; 


% parameter grids
ratios_txy_over_tz = [0.5, 1.0, 1.5];
M_over_tz_grid = linspace(0.2, 1.8, 25);


colors = lines(numel(ratios_txy_over_tz));
figure; hold on;
legend_str={};


for it = 1:numel(ratios_txy_over_tz)
txy_over_tz = ratios_txy_over_tz(it);
vals = zeros(size(M_over_tz_grid));


for iM = 1:numel(M_over_tz_grid)
M_over_tz = M_over_tz_grid(iM);
tz = -3.4*eta; % example scale (tune to your baseline)
txy = txy_over_tz * tz;
Mval = M_over_tz * tz;


p = struct('eta',eta,'txy',txy,'tz',tz,'M',Mval,'beta',beta,'gamma',gamma);
model = builder.make_builders_taguchi(p);


params = struct('model',model,'Nk',Nk,'eta',eta_b,'Ef',Ef,'mu',mu,'T',T, ...
'e',1,'hbar',1,'alpha','x','beta','y','gamma','z', ...
'units','lattice','shift',[0 0 0],'use_parfor',true);
out = bastin_main(params);
vals(iM) = out.sigma; % (ħ/e) lattice units if e=ħ=1
end


plot(M_over_tz_grid, vals, '-o', 'DisplayName', sprintf('t_{xy}/t_z = %.2f', txy_over_tz), 'Color', colors(it,:));
legend_str{end+1}=sprintf('t_{xy}/t_z=%.2f',txy_over_tz); %#ok<SAGROW>
end


xlabel('M/t_z'); ylabel('\sigma^{s_z}_{xy} (lattice units)'); grid on; legend show; title('Bastin SHC (Taguchi lattice)');