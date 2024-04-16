function S = solar_flux_R01(d)
% Outputs solar flux in W/m2 given distance from sun to target in meters.
% https://www.ess.uci.edu/~yu/class/ess200a/lecture.2.global.pdf
% Version 1.0 completed 1/5/2024

% Solar Luminosity
L = 3.83*10^26; % W

% Solar flux density
S = L./(4.*pi.*d.^2); % W/m2

end