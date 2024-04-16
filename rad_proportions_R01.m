function RPM = rad_proportions_R01(eps,F,cutoff)
% Calculates the radiation proportionality matrix for a set of elements for
% reflections under the cutoff value.

% Version 1.0 completed 2/20/2024

[N,~] = size(F);

% Reflectivities
rho = 1 - eps;
drF = diag(rho)*F;

temp = eps .* F;
RPM = zeros(N);
diff = Inf;
while diff > cutoff
    diff = max(temp,[],"all");
    RPM = RPM + temp;
    temp = drF * temp;
end
end