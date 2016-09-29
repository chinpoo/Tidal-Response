% func: compute phi-dependent function of a real form SH and its first and
%       second order derivatives, for purpose of numerical integration.
% input args: 
%       m: harmonic order
%       phi: row vector of longitudinal angle "phi"
% outputs:
%       f: vector of values of phi-dependent function evalauted at phi
%       fp: vector of values of 1st order derivative evalauted at phi
%       fpp: vector of values of 2nd order derivative evalauted at phi

function [f,fp,fpp] = phi_func(m,phi)
    n = length(phi);
    if m == 0
        f = ones(1,n);
        fp = zeros(1,n);
        fpp = zeros(1,n);
    elseif m > 0
        f = cos(m*phi);
        fp = -m*sin(m*phi);
        fpp = -m^2*cos(m*phi);
    else
        mp = abs(m);
        f = sin(mp*phi);
        fp = mp*cos(mp*phi);
        fpp = -mp^2*sin(mp*phi);
    end
end     % end of func