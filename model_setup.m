% func: set up the 1-D and 3-D models before tidal response calculation 
% input args: 
%     fname_1D: input file of 1-D profile of the planetary body
%               ***** radius, density, vp, vs *****
%     fname_3D: input file of 3-D structure in mu, lambda and rho
%               ***** # of layers with 3-D structure, # of harmonics in 3-D structure *****
%               ***** lower boundary, upper boundary *****
%               ***** order, l1, m1, delta_mu, delta_lambda, delta_rho (lateral variability)*****
% outputs:
%     Model: structure of 1-D and 3-D model

function [Model] = model_setup(fname_1D,fname_3D)
    
    % output variables
    Model = struct;

    % *************************** read inputs *****************************
    % read in 1-D radial profiles of density, vp and vs
    if exist(fname_1D,'file') ~= 2
        error('1D profile file does not exist...\n');
    else
        m_1d = load(fname_1D);
        rad = m_1d(:,1);                % upper bound radius of material layers
        rho = m_1d(:,2);                % material density
        vp = m_1d(:,3);                 % Vp
        vs = m_1d(:,4);                 % Vs
        mu = vs.^2.*rho;                % shear modulus
        lambda = vp.^2.*rho - 2*mu;     % lambda
        l_oc = find(vs == 0.0);         % layer of outer core
        R_CMB = rad(l_oc);              % radius of CMB
        R_SURF = rad(end);              % radius of surface
    end
    
    % read in 3-D (laterally varying) structures in density and elastic
    % moduli in form of lateral variability in different depth ranges
    if exist(fname_3D,'file') ~= 2
        error('3D structure file does not exist...\n');
    else
        m_3d = dlmread(fname_3D);
        n_depth = m_3d(1,1);            % number of depth ranges of different lateral variabilities
        n_harm = m_3d(1,2);             % number of harmonics eigenstructures
        % create new radial profiles of the mantle such that in each layer
        % 1) material properties are uniform
        % 2) lateral variabilities are uniform
        r_lower = zeros(n_depth,1);
        r_upper = zeros(n_depth,1);
        lv = cell(n_depth,1);
        size_of_lv = 6;                 % number of elements per row
        for i = 1:n_depth
            ind = 2 + (i-1)*(n_harm+1);
            r_lower(i) = m_3d(ind,1);
            r_upper(i) = m_3d(ind,2);
            lv{i} = m_3d(ind+1:ind+n_harm,1:size_of_lv);
        end
        RAD = unique([rad;r_lower;r_upper]);        % radius of new upper layer boundaries
        if RAD(l_oc) ~= R_CMB || RAD(end) ~= R_SURF
            error('Mantle is incorrectly layered...\n');  
        else
            % create new model layering
            N = size(RAD,1);            % number of layers
            N_MAN = N - 2;              % number of mantle layers (2 core layers)
            L_HET = zeros(N,1);         % indicator of lateral heterogeneity
            % mark down which layers are laterally heterogeneous
            for i = 1:n_depth
                i_lower = find(RAD == r_lower(i));
                i_upper = find(RAD == r_upper(i));
                if i_lower >= i_upper
                    error('[1] incorrect lower or upper boundary...\n');
                else
                    if any(L_HET(i_lower+1:i_upper))
                        error('[2] incorrect lower or upper boundary...\n');
                    else
                        L_HET(i_lower+1:i_upper) = i;
                    end
                end 
            end
            % new 1-D profile
            RHO = zeros(N,1);         % new density profile
            MU = zeros(N,1);          % new mu profile
            LAMBDA = zeros(N,1);      % new lambda profile
            for i = 1:size(rad)
                if i == 1
                    ind = (RAD <= rad(i));
                else
                    ind = (RAD > rad(i-1) & RAD <= rad(i));
                end
                RHO(ind) = rho(i);
                MU(ind) = mu(i);
                LAMBDA(ind) = lambda(i);
            end         
        end
    end
    
    % *********************** non-dimensionalization **********************
    % equations are non-dimensionalized before being solved
    G = 6.67e-11;           % gravitational constant
    % normalization scalings
    R_0 = R_SURF;
    MU_0 = MU(1);
    RHO_0 = RHO(1);
    
    RAD = RAD/R_0;                          % non-dim radius
    RHO = RHO/RHO_0;                        % non-dim density
    MU = MU/MU_0;                           % non-dim shear modulus
    LAMBDA = LAMBDA/MU_0;                   % non-dim lambda
    
    ETA = 4*pi*G*RHO_0^2*R_0^2/MU_0;        % non-dim coefficient eta enters in A matrix
    BETA = LAMBDA + 2*MU;                   % beta enters in A matrix
    GAMMA = MU.*(2*MU + 3*LAMBDA)./BETA;    % gamma enters in A matrix

    MASS = compute_mass(RAD,RHO,RAD);       % non-dim mass below layer boundary
    
    R_CMB = R_CMB/R_0;                      % non-dim R_CMB
    R_SURF = R_SURF/R_0;                    % non_dim R_SURF
    G_CMB = MASS(l_oc)/R_CMB^2;             % non_dim g0 at R_CMB
    G_SURF = MASS(end)/R_SURF^2;            % non_dim g0 at R_SURF

    % *************************** grid set-up *****************************
    % 1-D radial grid for implementation of propagator matrix method
    % implementation of 4th order Runge-Kutta method request:
    % 1) coarsest grid at 2nd order of perturbation
    % 2) lower order of perturbation doubles the resolution of higher order
    
    % default upperbound of step size
    H0 = 0.005;
    % vector of mantle discretization in each layer, size N_MAN
    NM = ceil((RAD(l_oc+1:N) - RAD(l_oc:N-1))/H0);
    
    % ***************************** outputs *******************************
    % constants
    Model.G = G;
    Model.eta = ETA;
    Model.R_0 = R_0;
    Model.RHO_0 = RHO_0;
    Model.MU_0 = MU_0;
    Model.r_cmb = R_CMB;
    Model.r_sf = R_SURF;
    Model.g0_cmb = G_CMB;
    Model.g0_sf = G_SURF;
    
    % model
    Model.r_origin = rad/R_0;
    Model.N = N;
    Model.l_oc = l_oc;
    Model.N_MAN = N_MAN;
    Model.NM = NM;
    Model.r = RAD;
    Model.rho = RHO;
    Model.mu = MU;
    Model.lambda = LAMBDA;
    Model.beta = BETA;
    Model.gamma = GAMMA;
    Model.mass = MASS;
    
    Model.L_HET = L_HET;
    Model.n_depth = n_depth;
    Model.n_harm = n_harm;
    Model.r_lower = r_lower/R_0;
    Model.r_upper = r_upper/R_0;
    Model.L1 = lv{1}(:,2);
    Model.M1 = lv{1}(:,3);
    Model.lv_mu = cell(n_harm,1);
    Model.lv_ld = cell(n_harm,1);
    Model.lv_rho = cell(n_harm,1);
    % constuct 3-D structure through the mantle for each (l1,m1)
    for i = 1:n_harm
        Model.lv_mu{i} = zeros(N,1);
        Model.lv_ld{i} = zeros(N,1);
        Model.lv_rho{i} = zeros(N,1);
        for j = 1:n_depth
            ind = (L_HET == j);
            Model.lv_mu{i}(ind) = lv{j}(i,size_of_lv-2);
            Model.lv_ld{i}(ind) = lv{j}(i,size_of_lv-1);
            Model.lv_rho{i}(ind) = lv{j}(i,size_of_lv);
        end
    end

end         % end of func
