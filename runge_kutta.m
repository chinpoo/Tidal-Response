% func: 4th order Runge-Kutta based propagator matrix method, used to solve 
%       1st order ODE in matrix form. Main functionalities:
%       1. Construction of propagator matrix
%       2. Integration of (mode coupling) forcing vectors
%       3. Dealing with traction and Q (due to d_rho) discontinuities
% input args: 
%       nn: vector of # of nodes in each mantle layer
%       MC: structure of info of a specific mode, including "parent" and coupling coefficients
%       MD: structure of 1-D and 3-D model
%       F:  structure of flags
% outputs:
%       RK: structure of Runge-Kutta method ingredients

function [RK] = runge_kutta(ng,MC,MD,nv,F)

    global ACC;
    
    % vector of # of nodes in each mantle layer (two nodes at layer interface)
    nn = ng + 1;
    % # of nodes from CMB to surface
    N = sum(nn);
    % -------------------------- Initialization ---------------------------
    % construct a structure to store all intermediate results...
    RK = struct;
    % dimension of propagator matrix
    RK.dim = 0;
    % if it is a spheroidal mode...
    if MC.type == 1
        % in particular, a (0,0) mode in response
        if MC.l == 0 && MC.m == 0
            RK.dim = 3;
        else
            RK.dim = 6;
        end
    % if it is a toroidal mode...
    elseif MC.type == -1
        RK.dim = 2;    
    end
    % differential propagator matrix across two neighboring nodes: prop(ri->ri+1), i = 0,1...N
    RK.prop_diff = zeros(RK.dim,RK.dim,N);
    % cumulative propagator matrix from CMB to any node: prop(r0->ri), i = 0,1...N
    p_cum = eye(RK.dim);
    RK.prop_cum = zeros(RK.dim,RK.dim,N);
    % A matrix at every node: A(ri), i = 0,1...N
    RK.A = zeros(RK.dim,RK.dim,N);
    % radius of each node
    RK.r = zeros(N,1);
    % mode coupling "forcing" vector and traction jump needed for higher orders
    if MC.order > 0
        % forcing at every node, f(ri), i = 0,1...N
        RK.f = zeros(RK.dim,1,N);
        % "forcing" and traction modification integrated from CMB to any node
        % f(r0->ri), i = 0,1...N
        f_cum = zeros(RK.dim,1);
        RK.f_cum = zeros(RK.dim,1,N);
        % trac(r0->ri), i = 0,1...N
        RK.trac_cum = zeros(RK.dim,1,N);
    else
        RK.f = [];
        RK.f_cum = [];
        RK.trac_cum = [];
    end
    
    % for higher order...
    if MC.order > 0
        % determine indices of the nodes at layer interfaces...
        % number of nodes at interfaces
        N_T = 2*length(ng);
        node_T = zeros(N_T,1);
        % nodal indices...
        i_top = 0;
        for i = 1:length(ng)
            node_T(2*i-1) = i_top + 1;       % bottom node of the layer
            node_T(2*i) = node_T(2*i-1) + 2*ng(i);
            i_top = node_T(2*i);
        end
        % ---------------------------------------------------------------------
        % -------------------- Construct "forcing" vectors --------------------
        % read solutions from "parent" modes (full set of "parents" obtained from MC.mu)
        for i_p = 1:size(MC.mu,1)
            type_p = MC.mu(i_p,1);
            l_p = MC.mu(i_p,2);
            m_p = MC.mu(i_p,3);
            fname_p = file_sol(type_p,l_p,m_p,MC.order-1,MC.L0,MC.M0,MC.L1,MC.M1);
            if exist(fname_p,'file') ~= 2
                error(['parent solution file',fname_p,'does NOT exist...\n']);
            else
                % read in "parent" solution array
                tmp = load(fname_p);
                sol = tmp.sol;
                r_p = sol(:,1);            % col#1: radius of nodes at "parent" resolution
                N_p = length(r_p);
                if N_p ~= sum(2*ng + 1) || N_p ~= length(nv.mu)
                    error(['incorrect solution file',fname_p,'...\n']);
                else
                    % gravitational acceleration g0
                    g0_p = compute_mass(r_p,MD.rho,MD.r)./r_p.^2;
                    % initialization of "force" and discontinuities from sol...
                    if i_p == 1
                        % spheroidal...
                        if MC.type == 1
                            % due to delta_mu
                            fP_dm = zeros(N_p,1);       % Plm component of d_mu induced forcing
                            fB_dm = zeros(N_p,1);       % Blm component of d_mu induced forcing
                            R_dm = zeros(N_T,1);        % R component of d_mu induced traction term
                            S_dm = zeros(N_T,1);        % S component of d_mu induced traction term
                            % due to delta_lambda
                            fP_dl = zeros(N_p,1);
                            fB_dl = zeros(N_p,1);
                            R_dl = zeros(N_T,1);
                            % NO d_lambda induced S component, S_dl is always 0
                            % due to delta_rho
                            fP_dr = zeros(N_p,1);
                            fB_dr = zeros(N_p,1);
                            % delta rho causes extra "forcing" and discontinuity in Poisson's eqn
                            fQ_dr = zeros(N_p,1);
                            Q_dr = zeros(N_T,1);
                        elseif MC.type == -1
                            % due to delta_mu
                            fC_dm = zeros(N_p,1);       % Clm component of d_mu induced forcing
                            T_dm = zeros(N_T,1);        % T component of d_mu induced traction term
                            % NO toroidal mode from delta_lambda, forcing is always zero
                            % due to d_rho
                            fC_dr = zeros(N_p,1);
                        end
                    end     % "forcing" initialization
                    
                    % read in solutions...
                    if type_p == 1      % N_p-by-11
                        U_p = sol(:,2);        % radial displacement
                        V_p = sol(:,3);        % horizontal displacement
                        S_p = sol(:,4)./nv.mu;    % horizontal component of traction
                        Up_p = sol(:,5);       % 1st order derivative of U_p
                        Sp_p = sol(:,6)./nv.mu;   % 1st order derivative of S_p
                        Upp_p = sol(:,7);      % 2nd order dericative of U_p
                        % more auxillary variables for delta_lambda
                        X_p = sol(:,8);        % x = div(u)
                        Xp_p = sol(:,9);       % 1st order derivative of X_p
                        K_p = sol(:,10);       % potential increment
                        Kp_p = sol(:,11);      % 1st order derivative of K_p
                        % Vp_p = sol(:,12);      % 1st order derivative of V_p
                        
                        % Case #1: S + LH ----> S
                        if MC.type == 1
                            % ----- d_mu -----
                            if F.IS_D_MU && ~isempty(MC.mu)
                                c = MC.mu(i_p,4:9);
                                % fP = d_mu*[c1*(Upp+2*Up/r-2*U/r^2)+c2*S/r+c3*V/r^2]
                                % fB = d_mu*[c4*(Sp+3*S/r)+c5*U/r^2+c6*V/r^2]
                                fP_dm = fP_dm + nv.d_mu.*(c(1)*(Upp_p+2*Up_p./r_p-2*U_p./r_p.^2)+c(2)*S_p./r_p+c(3)*V_p./r_p.^2);
                                fB_dm = fB_dm + nv.d_mu.*(c(4)*(Sp_p+3*S_p./r_p)+c(5)*U_p./r_p.^2+c(6)*V_p./r_p.^2);
                                % R = d_mu*c1*Up;
                                % S = d_mu*c4*S;
                                R_dm = R_dm + c(1)*nv.d_mu(node_T).*Up_p(node_T);
                                % ??? note that if mode_p = (0,0), V, S, Sp are zeros ???
                                S_dm = S_dm + c(4)*nv.d_mu(node_T).*S_p(node_T);
                            end
                            
                            % ----- d_lambda -----
                            if F.IS_D_LD && ~isempty(MC.ld)
                                loc = ismember(MC.ld(:,1:3),[type_p,l_p,m_p],'rows');
                                if any(loc)
                                    c = MC.ld(loc,4:5);
                                    % fP = d_lambda*c1*Xp
                                    % fB = d_lambda*c2*X/r
                                    % ??? if mode_c = (0,0),c2 = 0 ???
                                    fP_dl = fP_dl + c(1)*nv.d_ld.*Xp_p;
                                    fB_dl = fB_dl + c(2)*nv.d_ld.*X_p./r_p;
                                    % R = d_lambda*c1*X
                                    % S = 0;
                                    R_dl = R_dl + c(1)*nv.d_ld(node_T).*X_p(node_T);
                                end
                            end
                            
                            % ----- d_rho -----
                            if F.IS_D_RHO
                                % term #1: -d_rho*grad(Phi_1) or -d_rho*grad(Phi_0 + V_td)
                                if ~isempty(MC.rho_1)
                                    loc = ismember(MC.rho_1(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        c = MC.rho_1(loc,4:5);
                                        % fP = d_rho*c1*Kp
                                        % fB = d_rho*c2*(K/r)
                                        fP_dr = fP_dr + c(1)*nv.d_rho.*Kp_p;
                                        fB_dr = fB_dr + c(2)*nv.d_rho.*K_p./r_p;
                                        % **************************************
                                        % 1st order tidal forcing term -d_rho*grad(V(t)) may induce a degree-1 mode that 
                                        % causes a net force on the system. Our goal here is to remove the net force by
                                        % introducing an inertial force onto every portion of the system
                                        % Note that K_p and Kp_p already contains tidal potential itself at 1st order, i.e. psi+V_tf 
                                        if MC.order == 1 && MC.l == 1
                                            % f_tide = -d_rho*grad(V_td(t))
                                            fP_tide = 2*c(1)*nv.d_rho.*r_p/(2*MC.L0+1);
                                            fB_tide = c(2)*nv.d_rho.*r_p/(2*MC.L0+1);
                                            % compute the total acceleration caused by f_tide
                                            [ACC,f_net] = net_force(fP_tide,fB_tide,r_p,MD.mass(end));
                                            % subtract inertial force from every portion of the system
                                            fP_dr = fP_dr - ACC*nv.rho;
                                            fB_dr = fB_dr - ACC*nv.rho;
                                            % ---- validation ----
                                            % In inertial frame, the core sense an inertial force that cancels the total
                                            % force in the mantle.
                                            f_core = -sqrt(3/4/pi)*ACC*(MD.mass(2)*4*pi);
                                            fP_tide_1 = fP_tide - ACC*nv.rho;
                                            fB_tide_1 = fB_tide - ACC*nv.rho;
                                            [tmp,f_man] = net_force(fP_tide_1,fB_tide_1,r_p,MD.mass(end));
                                            fprintf('*** net force in accelerating frame: %g ***\n',f_man+f_core);
                                            % --------------------
                                        end
                                        % **************************************
                                    end
                                end
                                % term #2: in both eqn of motion and Poisson's eqn
                                if ~isempty(MC.rho_2)
                                    loc = ismember(MC.rho_2(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        % ??? double check!! ???
                                        c = MC.rho_2(loc,4:5);
                                        % fQ = d_rho*(c1*X+c2*V/r)
                                        % eqn of motion, Plm component only
                                        fq = nv.d_rho.*(c(1)*X_p + c(2)*V_p./r_p);
                                        fP_dr = fP_dr + fq.*g0_p;
                                        % Poisson's eqn
                                        fQ_dr = fQ_dr + fq;
                                        % Q = d_rho*c1*U
                                        Q_dr = Q_dr + c(1)*nv.d_rho(node_T).*U_p(node_T);
                                    end
                                end
                                % term #3: only for 2nd order, from coupling of d_rho and ACC
                                if MC.order == 2 && ~isempty(MC.rho_3) && l_p == 1 && ACC ~= 0
                                    loc = ismember(MC.rho_3(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        c = MC.rho_3(loc,4:5);
                                        % fP = c1*d_rho*ACC
                                        % fB = c2*d_rho*ACC
                                        fP_dr = fP_dr - c(1)*ACC*nv.d_rho;
                                        fB_dr = fB_dr - c(2)*ACC*nv.d_rho;
                                    end
                                end
                            end
                            
                        % Case #2: S + LH ----> T
                        elseif MC.type == -1
                            % ----- d_mu -----
                            if F.IS_D_MU && ~isempty(MC.mu)
                                c = MC.mu(i_p,4:6);
                                % fC = d_mu*[c1*(Sp+3*S/r)+c2*U/r^2+c3*V/r^2] (c2 = 0)
                                fC_dm = fC_dm + nv.d_mu.*(c(1)*(Sp_p+3*S_p./r_p)+c(3)*V_p./r_p.^2);
                                % traction discontinuity at layer interfaces
                                % T = d_mu*c1*S;
                                T_dm = T_dm + c(1)*nv.d_mu(node_T).*S_p(node_T);
                            end
                            
                            % ----- d_rho -----
                            if F.IS_D_RHO
                                % term #1
                                if ~isempty(MC.rho_1)
                                    loc = ismember(MC.rho_1(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        c = MC.rho_1(loc,4);
                                        % fC = c1*d_rho*(K/r)
                                        fC_dr = fC_dr + c*nv.d_rho.*K_p./r_p;
                                    end
                                end
                                % term #3
                                if MC.order == 2 && ~isempty(MC.rho_3) && l_p == 1 && ACC ~= 0
                                    loc = ismember(MC.rho_3(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        c = MC.rho_3(loc,4);
                                        % fC = c*d_rho*ACC
                                        fC_dr = fC_dr - c*ACC*nv.d_rho;
                                    end
                                end
                            end
                        end   % MC.type
                        
                    elseif type_p == -1        % N_p-by-4
                        W_p = sol(:,2);        % horizontal displacement
                        T_p = sol(:,3)./nv.mu;    % horizontal component of traction
                        Tp_p = sol(:,4)./nv.mu;   % 1st order derivative of T_p
                        % Wp_p = sol(:,5);       % 1st order derivative of W_p
                        
                        % Case #3: T + LH ----> S
                        if MC.type == 1
                            % ----- d_mu -----
                            if F.IS_D_MU && ~isempty(MC.mu)
                                c = MC.mu(i_p,4:6);
                                % fP = d_mu*c1*T/r
                                % fB = d_mu*[c2*(Tp+3*T/r)+c3*W/r^2]
                                fP_dm = fP_dm + c(1)*nv.d_mu.*T_p./r_p;
                                fB_dm = fB_dm + nv.d_mu.*(c(2)*(Tp_p+3*T_p./r_p)+c(3)*W_p./r_p.^2);
                                % R = 0
                                % S = d_mu*c2*T
                                S_dm = S_dm + c(2)*nv.d_mu(node_T).*T_p(node_T);
                            end
                            
                            % ----- d_rho -----
                            if F.IS_D_RHO
                                if ~isempty(MC.rho_2)
                                    loc = ismember(MC.rho_2(:,1:3),[type_p,l_p,m_p],'rows');
                                    if any(loc)
                                        c = MC.rho_2(loc,4);
                                        % fQ = d_rho*c*W/r
                                        fq = c*nv.d_rho.*W_p./r_p;
                                        % eqn of motion, Plm component only
                                        fP_dr = fP_dr + fq.*g0_p;
                                        % Poisson's eqn
                                        fQ_dr = fQ_dr + fq;
                                        % Q = 0
                                    end
                                end
                            end
                            
                        % Case #4: T + LH ----> T
                        elseif MC.type == -1
                            % ----- d_mu -----
                            if F.IS_D_MU && ~isempty(MC.mu)
                                c = MC.mu(i_p,4:5);
                                % fC = d_mu*[c1*(Tp+3*T/r)+c2*W/r^2]
                                fC_dm = fC_dm + nv.d_mu.*(c(1)*(Tp_p+3*T_p./r_p)+c(2)*W_p./r_p.^2);
                                % T = d_mu*c1*T
                                T_dm = T_dm + c(1)*nv.d_mu(node_T).*T_p(node_T);
                            end
                        end  % MC.type
                    end % type_p
                end % N_p
            end % read in "parent" mode
        end % loop of parent modes
        % combine "forcing" from mu, lambda and rho...
        if MC.type == 1
            fP = fP_dm + fP_dl + MD.eta*fP_dr;
            fB = fB_dm + fB_dl + MD.eta*fB_dr;
            fQ = fQ_dr;
            R = R_dm + R_dl;
            S = S_dm;
            Q = Q_dr;
        elseif MC.type == -1
            fC = fC_dm + MD.eta*fC_dr;
            T = T_dm;
        end
    end % MC.order > 0
    
    % ---------------------------- Runge-Kutta ----------------------------
    for i_l = 1:length(ng)
        % # of subdivisions in each layer
        n = ng(i_l);
        % indexing of bottom layer interface
        ind = sum(nn(1:i_l)) - n;
        % indexing of bottom layer interface at finer grid
        ind_finer = 2*ind - i_l;
        % indexing of material layer
        l_mat = i_l + MD.l_oc;
        % step size in each layer
        h = (MD.r(l_mat)-MD.r(l_mat-1))/n;
        
        % dealing with the 1st node of each layer...
        RK.prop_diff(:,:,ind) = eye(RK.dim);
        if i_l == 1     % at r_cmb
            RK.prop_cum(:,:,ind) = eye(RK.dim);
        else
            RK.prop_cum(:,:,ind) = RK.prop_cum(:,:,ind-1);
            if MC.order > 0
                RK.f_cum(:,:,ind) = RK.f_cum(:,:,ind-1);
            end
        end
        % ----- traction jump at interfaces -----
        if MC.order > 0
            if MC.type == 1
                if MC.l == 0 && MC.l == 0 % (0,0)
                    if i_l == 1  % r_cmb
                        RK.trac_cum(:,:,ind) = -[0,R(i_l),0]';
                    else
                        RK.trac_cum(:,:,ind) = RK.trac_cum(:,:,ind-1) + [0,R(2*i_l-2)-R(2*i_l-1),0]';
                    end
                else
                    if i_l == 1  % r_cmb
                        RK.trac_cum(:,:,ind) = -[0,0,R(i_l),S(i_l),0,Q(i_l)]';
                    else
                        RK.trac_cum(:,:,ind) = RK.trac_cum(:,:,ind-1) + [0,0,R(2*i_l-2)-R(2*i_l-1),S(2*i_l-2)-S(2*i_l-1),0,Q(2*i_l-2)-Q(2*i_l-1)]';
                    end
                end
            elseif MC.type == -1
                if i_l == 1  % r_cmb
                    RK.trac_cum(:,:,ind) = -[0,T(i_l)]';
                else
                    RK.trac_cum(:,:,ind) = RK.trac_cum(:,:,ind-1) + [0,T(2*i_l-2)-T(2*i_l-1)]';
                end
            end
        end     % end traction at interfaces
        
        % loop over n sublayers
        for i = 1:n
            % 3 evaluation points within each sublayer
            k = [2*i-2,2*i-1,2*i]';
            % radius of 3 points 
            r = MD.r(l_mat-1) + h/2*k;
            % g0 at r...
            g0 = compute_mass(r,MD.rho,MD.r)./r.^2;
            % A matrix at r...
            if i == 1
                a0 = create_a_matrix(r(1),g0(1),l_mat,MC.l,MC.type,MD);
            else
                a0 = a2;
            end
            a1 = create_a_matrix(r(2),g0(2),l_mat,MC.l,MC.type,MD);
            a2 = create_a_matrix(r(3),g0(3),l_mat,MC.l,MC.type,MD);
            
            P0 = a0;
            P1 = a1 + h/2*a1*P0;
            P2 = a1 + h/2*a1*P1;
            P3 = a2 + h*a2*P2;
            % -------- prop[r(i-1) -> r(i)] --------
            p_diff = eye(RK.dim) + h/6*(P0+2*P1+2*P2+P3);
            RK.prop_diff(:,:,ind+i) = p_diff;
            % -------- prop[r_cmb -> r(i)] --------
            p_cum = p_diff*p_cum;
            RK.prop_cum(:,:,ind+i) = p_cum;
            % -------- r(i) & A[r(i)] --------
            if i == 1
                RK.r(ind) = r(1);
                RK.A(:,:,ind) = a0;
            end
            RK.r(ind+i) = r(3);
            RK.A(:,:,ind+i) = a2;
            
            % "forcing" vector at r...
            if MC.order > 0
                j = ind_finer+k;
                if MC.type == 1
                    if MC.l == 0 && MC.l == 0 % (0,0)
                        % 07/20/2015: no fQ_dr for (0,0)?
                        fv0 = [0,fP(j(1)),0]';
                        fv1 = [0,fP(j(2)),0]';
                        fv2 = [0,fP(j(3)),0]';
                    else
                        fv0 = [0,0,fP(j(1)),fB(j(1)),0,fQ(j(1))]';
                        fv1 = [0,0,fP(j(2)),fB(j(2)),0,fQ(j(2))]';
                        fv2 = [0,0,fP(j(3)),fB(j(3)),0,fQ(j(3))]';
                    end
                elseif MC.type == -1
                    fv0 = [0,fC(j(1))]';
                    fv1 = [0,fC(j(2))]';
                    fv2 = [0,fC(j(3))]';
                end
                F0 = fv0;
                F1 = fv1 + h/2*a1*F0;
                F2 = fv1 + h/2*a1*F1;
                F3 = fv2 + h*a2*F2;
                % -------- f[r(i-1) -> r(i)] --------
                f_diff = h/6*(F0+2*F1+2*F2+F3);                
                % -------- f[r_cmb -> r(i)] --------
                f_cum = p_diff*f_cum + f_diff;
                RK.f_cum(:,:,ind+i) = f_cum; 
                % -------- f[r(i)] --------
                if i == 1
                    RK.f(:,:,ind) = fv0;
                end
                RK.f(:,:,ind+i) = fv2;
                % -------- trac[r_cmb -> r(i)] --------
                RK.trac_cum(:,:,ind+i) = p_diff*RK.trac_cum(:,:,ind+i-1);
            end   
        end % loop over sublayers
    end % loop of i_l
    
    % record trac_sf: traction jump matching surface boundary condition
    if MC.order > 0
        if MC.type == 1
            if MC.l == 0 && MC.l == 0 % (0,0)
                RK.trac_sf = RK.trac_cum(:,:,end) + [0,R(end),0]';
            else
                RK.trac_sf = RK.trac_cum(:,:,end) + [0,0,R(end),S(end),0,Q(end)]';
            end
        elseif MC.type == -1
            RK.trac_sf = RK.trac_cum(:,:,end) + [0,T(end)]';
        end
    else
        RK.trac_sf = [];
    end

end     % end of func
