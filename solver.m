% func: slove linear equation for each mode, write solutions to file and
%       output response
% input args:
%       RK: structure of Runge-Kutta method ingredients
%       MC: structure of info of a specific mode, including "parent" and coupling coefficients
%       MD: structure of 1-D and 3-D model
% outputs:
%       resp: vector of tidal response
%       out: solution output matrix

function [resp,out] = solver(RK,MC,MD)

    % tidal Love numbers
    global H0 K0 L0;
    
    % arbirary magnitude of tide, Q(a)...
    Q_sf = 1;
    % Different definitions of potential results in different signs of h and l response.
    % sign_pttl = 1: potential; sign_pttl = -1: potential energy (in Citcom).
    sign_pttl = -1;

    % propagator matrix from CMB to surface
    prop_tot = RK.prop_cum(:,:,end);
    % "forcing" and traction jump...
    if MC.order == 0
        if ~isempty(RK.f) || ~isempty(RK.f_cum) || ~isempty(RK.trac_cum) || ~isempty(RK.trac_sf)
            error('no forcing vector or traction jump for 0th order equation...\n');
        end
    else
        % sum of "forcing" and traction jump modification at surface
        f_tot = RK.f_cum(:,:,end) - RK.trac_sf;
    end
    
    % build the three starting solution vectors at CMB for spheroidal mode...
    if MC.type == 1
        if MC.l == 0 && MC.m == 0 % (0,0)
            sol_v = zeros(3,2);
            sol_v(:,1) = [1,MD.eta*MD.rho(MD.l_oc)*MD.g0_cmb,0]';
            sol_v(:,2) = [0,MD.eta*MD.rho(MD.l_oc),1]';
            
            tmp = prop_tot*sol_v;
            % solving for coefficients for each starting solution
            coef = linsolve(tmp([2,3],:),f_tot([2,3]));
            sol_sf = tmp*coef - f_tot;      % X(a)
            sol_cmb = sol_v*coef;           % X(r_cmb)
            
            % tidal response
            h = (2*MC.L0+1)*MD.g0_sf*sol_sf(1)/Q_sf*sign_pttl;
            k = (2*MC.L0+1)*sol_sf(3)/Q_sf;
            h_rel = h/H0;
            k_rel = k/K0;
            l = 0;
            l_rel = 0;
        else
            sol_v = zeros(6,3);
            sol_v(:,1) = [1,0,MD.eta*MD.rho(MD.l_oc)*MD.g0_cmb,0,0,MD.rho(MD.l_oc)]';
            sol_v(:,2) = [0,1,0,0,0,0]';
            sol_v(:,3) = [0,0,MD.eta*MD.rho(MD.l_oc),0,1,(2*MC.l+1)/MD.r_cmb]';
            
            % degree-1 spheroidal mode has translational mode...
            if MC.l == 1
                % net force of the system: consistency relation
                CR = f_tot(3)+2*f_tot(4)-MD.eta*MD.g0_sf*f_tot(6);
                str = sprintf('*** net force for S(%d,%d): %g ***\n',MC.l,MC.m,CR);
                fprintf('%s',str);
                
                % solve for degree-1 mode...
                tmp = prop_tot*sol_v(:,1:2);
                coef = linsolve(tmp([3,4],:),f_tot([3,4]));
                sol_tmp = tmp*coef - f_tot;
                % remove translational mode to obtain solution at surface...
                sol_sf = sol_tmp + sol_tmp(5)/MD.g0_sf*[1,1,0,0,-MD.g0_sf,0]';
                % CMB solution with translational mode...
                % ******** NOTE ********
                % including translational mode or not does not affect stress and 
                % displacement derivatives in higher order calculations
                % **********************
                sol_cmb = sol_v(:,1:2)*coef;
            % non-degree-1 modes
            else
                tmp = prop_tot*sol_v;
                if MC.order == 0
                    % tidal loading boundary condition...
                    bc = [0,0,0,0,0,Q_sf]';       % set Q(a)=1
                    coef = linsolve(tmp([3,4,6],:),bc([3,4,6]));
                    f_tot = -bc;
                % higher order modes...
                else
                    coef = linsolve(tmp([3,4,6],:),f_tot([3,4,6]));
                end
                sol_sf = tmp*coef - f_tot;
                sol_cmb = sol_v*coef;
            end
            
            % tidal response...
            h = (2*MC.L0+1)*MD.g0_sf*sol_sf(1)/Q_sf*sign_pttl;
            l = (2*MC.L0+1)*MD.g0_sf*sol_sf(2)/Q_sf*sign_pttl;
            k = (2*MC.L0+1)*sol_sf(5)/Q_sf;
            if MC.order == 0
                k = k - 1;     % subtract tidal potential itself
                H0 = h;
                K0 = k;
                L0 = l;
            end
            h_rel = h/H0;
            k_rel = k/K0;
            l_rel = l/L0;   
        end
        resp = [h_rel,k_rel,l_rel,h,k,l];
        
    elseif MC.type == -1
        if MC.l == 1
            % ***** degree-1 toroidal motion will cause zero net rotation *****
            % I do not compute displacement for degree-1 toroidal mode, we
            % only need stress for its "child" mode calculation.
            % P11*Wc + P12*Ws = f_tot(1), Wc, Ws are horizontal displacements at
            % CMB and surface, and contain a net rotation R0. R0 can be
            % arbitrary and will not affect the child mode response. Let Wc=0,
            % then Ws=f_tot(1)/P12. It works in obtaining correct coupling
            % response. We try to remove R0.
            % *****
            sol_cmb = [0,0]';
            sol_sf = [-f_tot(1),0]';
        else
            tmp = [prop_tot(1,1),-1;prop_tot(2,1),0];
            coef = linsolve(tmp,f_tot);
            
            sol_cmb = [coef(1),0]';
            sol_sf = [coef(2),0]';
        end
        w = (2*MC.L0+1)*MD.g0_sf*sol_sf(1)/Q_sf*sign_pttl;
        w_rel = w/L0;
        resp = [0,0,w_rel,0,0,w];
    end     % end of output resp
    
    % -------------------------- solution output --------------------------
    if MC.type == 1
        col = 11;
    elseif MC.type == -1
        col = 4;
    end
    out = zeros(length(RK.r),col);
    
    if MC.order == 0
        f = zeros(RK.dim,1);
        f_cum = zeros(RK.dim,1);
        trac = zeros(RK.dim,1);
    end
    for i = 1:length(RK.r)
        prop = RK.prop_cum(:,:,i);
        a = RK.A(:,:,i);
        r = RK.r(i);
        if MC.order > 0
            f = RK.f(:,:,i);
            f_cum = RK.f_cum(:,:,i);
            trac = RK.trac_cum(:,:,i);
        end
        
        % solution X(r) at each node
        % note +/- sign here before f_cum and trac
        sol = prop*sol_cmb - f_cum + trac;
        % 1st order derivative of solution vector dX(r)/dr
        solp = a*sol - f;
        % if spheroidal, write Upp
        if MC.type == 1
            if MC.l == 0 && MC.m == 0
                upp = a(1,:)*solp - a(1,1)/r*sol(1);
            else
                upp = a(1,:)*solp - a(1,1)/r*sol(1) - a(1,2)/r*sol(2);
            end
            % auxillary variable
            x = solp(1) + (2*sol(1) - MC.l*(MC.l+1)*sol(2))/r;
            xp = upp + (2*solp(1) - MC.l*(MC.l+1)*solp(2))/r - (2*sol(1) - MC.l*(MC.l+1)*sol(2))/r^2;
            if MC.l == 0 && MC.m == 0
                out(i,:) = [r,sol(1),0,0,solp(1),0,upp,x,xp,sol(3),solp(3)];
            else
                out(i,:) = [r,sol(1),sol(2),sol(4),solp(1),solp(4),upp,x,xp,sol(5),solp(5)];
            end
        elseif MC.type == -1
            out(i,:) = [r,sol(1),sol(2),solp(2)];
        end
    end

end
