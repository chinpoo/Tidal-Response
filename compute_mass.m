% func: compute the total non-dim mass below certain radius for a given density profile. 
% input args: 
%     r:     column vector of radius, from small to large radius
%     rho_0: column vector of density profile, from bottom to surface
%     R_0:   column vector of radius of layer boundary, from bottom to
%                     surface; has the same size as rho_0.
% outputs:
%     mass:  column vector of mass below r

function mass = compute_mass(r,rho_0,R_0)
    if size(rho_0) ~= size(R_0)
        error('incorrect size of 1-D profile...\n');
    else
        N0 = length(rho_0);
        if r ~= sort(r)
            error('[1] incorrect r vector...\n');
        else
            if r(1) < 0 || r(end) > R_0(end)      % R_0(end) should be R_SURF
                error('[2] incorrect r vector...\n');
            else
                n = length(r);
                % c = 4/3*pi;         % dimensional coef
                c = 1/3;            % non-dim coef
                mass = zeros(n,1);  % vector of mass
                r_tmp = unique([r;R_0]);
                n1 = length(r_tmp);
                rho_tmp = zeros(n1,1);
                for i = N0:-1:1
                    ind = (r_tmp <= R_0(i));
                    rho_tmp(ind) = rho_0(i);
                end
                vol_tmp = c*r_tmp.^3;       % total volume below radius r_tmp
                m = zeros(n1,1);
                % compute mass
                for i = 1:n1
                    if i == 1
                        m(i) = rho_tmp(i)*vol_tmp(i);
                    else
                        m(i) = m(i-1) + rho_tmp(i)*(vol_tmp(i)-vol_tmp(i-1));
                    end
                end
                % output mass at r
                for i = 1:n
                    mass(i) = m(r_tmp == r(i));
                end
            end
        end
    end
end
    
    