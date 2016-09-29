% func: compute net force/acceleration on the system caused by coupling of
%       tidal force and d_rho, manifested by a 1st order degree "child"
%       mode
% input args:
%       fP,fP: Plm and Blm components of f_tide = -d_rho*grad(V_td)
%       r: vector of nodal radius
%       mass: mass of the planet
% outputs:
%       a: P1m and B1m components of the nonzero acceleration, aP == aB
%       f_tot: nonzero net force

function [a,f_tot] = net_force(fP,fB,r,mass)
    I = 0;
    n = length(r);
    % lc = 1, mc = -1,0,1 have the same form in computing total force
    f = (fP + 2*fB).*r.^2;
    % integrate f throughout the mantle for total force
    for i = 1:n-1
        I = I + (f(i)+f(i+1))*(r(i+1)-r(i))/2;
    end
    % see note for detail...
    f_tot = sqrt(4*pi/3)*I;
    % total acceleration
    a_tot = f_tot/(mass*4*pi);
    % P1m and B1m components, aP == aB
    a = sqrt(4*pi/3)*a_tot;
end