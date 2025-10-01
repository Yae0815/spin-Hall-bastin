function f = fermi(E, T)
% Fermi-Dirac with T in Kelvin, E in eV
kB = 8.617333262e-5; % eV/K
if T<=0
f = double(E<0);
else
f = 1 ./ (1 + exp(E/(kB*T)));
end
end