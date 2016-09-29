% func: compute prefactorized associated legendre polynomial (theta 
%       dependent)and its derivatives, for purpose of numerical integration.
% input args:
%       l,m: harmonic degree and order
%       theta: row vector of co-latitude angle "theta"
% outputs:
%       f: vector of values of theta-dependent function evalauted at theta
%       ft: vector of values of 1st order derivative evalauted at theta
%       ftt: vector of values of 2nd order derivative evalauted at theta
%       fds: vector of values of P_lm(theta)/sin(theta) evaluated at theta
% ** add fds to deal with sigularity in integration in vsh_expan.

function [f,ft,ftt,fds] = theta_func(l,m,theta)
    mp = abs(m);
    n = length(theta);
    x = cos(theta);
    y = sin(theta);

    if mp > l      % incorrect (l,m)
        f = zeros(1,n);
        ft = zeros(1,n);
        ftt = zeros(1,n);
        fds = zeros(1,n);
        fprintf('Warning: m > l, Ylm = 0...\n');
        return;
    else
        % prefactor...
        nlm = sqrt((2*l+1)/4/pi*factorial(l-mp)/factorial(l+mp));
        % 1-by-n
        f = assc_plm(l,mp);
        if l == 0
            ft = zeros(1,n);
            ftt = zeros(1,n);
            fds = zeros(1,n);
        else
            ft = -((l+mp)*(l-mp+1)*assc_plm(l,mp-1) - assc_plm(l,mp+1))/2;
            ftt = ((l+mp)*(l-mp+1)*(l+mp-1)*(l-mp+2)*assc_plm(l,mp-2) -...
                  2*(l^2-mp^2+l)*f + assc_plm(l,mp+2))/4;
            if mp==0
                fds = zeros(1,n);
            else
                fds = -1/2/mp*(assc_plm(l+1,mp+1) + (l-mp+1)*(l-mp+2)*assc_plm(l+1,mp-1));
            end
        end
        
        % prefactor for real-form SH
        if m == 0
            f = nlm*f;
            ft = nlm*ft;
            ftt = nlm*ftt;
        else
            f = nlm*sqrt(2)*f;
            ft = nlm*sqrt(2)*ft;
            ftt = nlm*sqrt(2)*ftt;
            fds = nlm*sqrt(2)*fds;
        end      
    end

    % nested func: build associated Legendre polynomial recursively (m>0)
    function [plm] = assc_plm(l,m)
        if abs(m) > l
            plm = zeros(1,n);
            return;
        end
        if m < 0
            plm = (-1)^m*factorial(l+m)/factorial(l-m)*assc_plm(l,-m);
            return;
        end
        if l == 0 && m == 0
            plm = ones(1,n);
            return;
        else
            if l == m
                plm = (-1)^l*y.^l;
                for i = 1:2:2*l-1
                    plm = plm*i;
                end
                return;
            elseif l == m+1
                plm = (2*l-1)*x.*assc_plm(m,m);
            else
                mm = m;
                p1 = assc_plm(mm,mm);
                p2 = (2*mm+1)*x.*p1;
                while mm+2 <= l
                    plm = ((2*mm+3)*x.*p2 - (mm+m+1)*p1)/(mm-m+2);
                    mm = mm + 1;
                    p1 = p2;
                    p2 = plm;
                end
            end
        end
    end     % end of assc_plm

end     % end of theta_func
