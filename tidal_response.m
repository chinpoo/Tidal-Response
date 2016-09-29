% =========================================================================
% ================= Version 2.0, by CHUAN QIN, 2016-03-10 =================
% =========================================================================
% This program is used to calculate the elastic tidal response of a planet
% with 3-D elastic and density structures in any depth ranges of the manle.
% Perturbation theory is adopted for the formulation and spectral method is
% used for solving mode couplings and tidal response. The total tidal 
% response is represented by the sum of modal responses at each order of 
% perturbation: 0, 1 and 2. Each model response is determined using
% propagator matrix method, which is based on Runge-Kutta numerical scheme.
%
% ---------------------------- Implementation -----------------------------
% 1. Initialization: model setups, variables and flags
% 2. Outer loop: for each (l1,m1) in 3-D structure
%                ----- determine mode coupling hierarchy (selection rule)
%                ----- store mode coupling information (VSH expansion)
% 3. Inner loop: for each (l,m) in tidal response (from 0th to 2nd order)
%                ----- construct matrix equation (4 categories)
%                ----- apply propagator matrix method for solution
%                ----- solving tidal response for (l,m)
% 4. Output: get total response and write it to file
%
% References: 
% =========================================================================

clear; clc;
t0 = cputime;
% **************************** Initialization *****************************
% ================= Manual Settings Start =================
% harmonic (l0,m0) of body tide force component considered...
L_tide = 2;         % l0
M_tide = 0;         % m0

% type in model file names and case name...
model_1d = 'moon.dat';
model_3d = 'test_3D_4.dat';
casename = 'test_1_TD20';
% model_1d = 'ct1.dat';
% model_3d = 'ct1_spectrum.dat';
% casename = 'crust_TD22';

% control flags...
Flag = struct;
Flag.IS_D_MU = 1;       % consider lateral heterogeneity in mu (1) or not (0)
Flag.IS_D_LD = 0;       % consider lateral heterogeneity in lambda (1) or not (0)
Flag.IS_D_RHO = 0;      % consider lateral heterogeneity in rho (1) or not (0)
Flag.IS_VIS_MODEL = 1;  % plot model for visualization (1) or not (0)

% order of perturbation, starting at 0th order, truncating at most 2nd order
Order_init = 0;
Order_end = 2;
% ================= Manual Settings End =================

% output file directories and naming...
global Output;
Output = struct;
Output.prefix = casename;
Output.dir_sol = './sol/';
Output.dir_vsh = './vsh_expan/';
Output.dir_model = './model/';
Output.dir_resp = './resp/';

% if system experiences nonzero net force due to d_rho
global ACC;
ACC = 0;

% tidal Love numbers, evaluated in solver.m
global H0 K0 L0;
H0 = 0; K0 = 0; L0 = 0;

% full path of input files...
inputfile_1d = [Output.dir_model,model_1d];
inputfile_3d = [Output.dir_model,model_3d];

% model setup from inputfiles...
[M] = model_setup(inputfile_1d,inputfile_3d);

% map layered properties on to nodes for purpose of constructing "forcing" vector
% for order of perturbation D = 1 or 2.
if Order_end > 0
    NV = cell(1,Order_end+1);
    Flag.IS_1D_MAPPING = ones(1,Order_end);    % control flag
end

% RESP holds tidal response solutions for each eigenstructure (l1,m1)
% for each (l1,m1), RESP(i_h) contains
% Spheroidal: [order,1,lc,mc,h_rel,k_rel,l_rel,h,k,l] 
% Toroidal: [order,-1,lc,mc,0,0,w_rel,0,0,w] 
RESP = cell(1,M.n_harm);
col_RESP = 10;
n_resp = 0;

% visualize the model...
if Flag.IS_VIS_MODEL
    vis_model(M.r_origin,M.r_lower,M.r_upper);
end

fprintf('Initialization completed...\n');
% -------------------------------------------------------------------------

% ****************************** Outer loop *******************************
% loop over (l1,m1)...
for i_h = 1:M.n_harm         
    L1 = M.L1(i_h);
    M1 = M.M1(i_h);
    fprintf('(l1,m1): (%d,%d)...\n',L1,M1);
    
    % full mode coupling hierarchy for (L_tide,M_tide)x(l1,m1) up to order
    % of perturbation D = 2
    MODE = mode_coupling(L_tide,M_tide,L1,M1,Output.dir_vsh);
    
    % modes in the final response...
    index = MODE.modes(:,1) <= Order_end;
    n_mode = sum(index);
    % initializing RESP for (l1,m1)
    RESP{i_h} = zeros(n_mode,col_RESP);
    RESP{i_h}(:,1:4) = MODE.modes(index,:);
    n_resp = n_resp + n_mode;

    % *************************** Inner loop ******************************
    % loop over response modes...
    for i_order = Order_init:Order_end
        % Note: According to 4th order Runge-Kutta scheme, in order to propagate
        % solution at 2i-1 to 2i+1, we need nodal values at 2i-1, 2i, and
        % 2i+1. To achieve so, 1st (0th) order grid needs as double resolution
        % 2nd (1st) order (this is especially important for a degree-1 mode).
        n_grid = 2^(Order_end-i_order)*M.NM;
        
        % mapping layer properties on to nodes, for only once...
        if i_order < Order_end
            if Flag.IS_1D_MAPPING(i_order+1) == 1
                NV{i_order+2}.mu = nodal_mapping(M.mu(M.l_oc+1:end),n_grid);
                NV{i_order+2}.ld = nodal_mapping(M.lambda(M.l_oc+1:end),n_grid);
                NV{i_order+2}.rho = nodal_mapping(M.rho(M.l_oc+1:end),n_grid);
                Flag.IS_1D_MAPPING(i_order+1) = 0;
            end
            NV{i_order+2}.d_mu = nodal_mapping(M.lv_mu{i_h}(M.l_oc+1:end),n_grid).*NV{i_order+2}.mu;
            NV{i_order+2}.d_ld = nodal_mapping(M.lv_ld{i_h}(M.l_oc+1:end),n_grid).*NV{i_order+2}.ld;
            NV{i_order+2}.d_rho = nodal_mapping(M.lv_rho{i_h}(M.l_oc+1:end),n_grid).*NV{i_order+2}.rho;
        end
        
        % % modes at order of perturbation i_order...
        modes = MODE.modes(MODE.modes(:,1)==i_order,:);
        
        % loop over modes at each order of perturbation
        for i_mode = 1:size(modes,1)
            mc = struct;
            mc.L0 = L_tide;
            mc.M0 = M_tide;
            mc.L1 = L1;
            mc.M1 = M1;
            mc.order = i_order;
            mc.type = modes(i_mode,2);
            mc.l = modes(i_mode,3);
            mc.m = modes(i_mode,4);
            info = [mc.order,mc.type,mc.l,mc.m];
            % get coupling info...
            if mc.order > 0
                mc.mu = MODE.mode_mu(ismember(MODE.mode_mu(:,1:4),info,'rows'),5:end);
                if isempty(mc.mu)
                    error('incorrect child mode info mc[1]...\n')
                end
                % if T->T only, mc.rho_1 is empty
                mc.rho_1 = MODE.mode_rho{1}(ismember(MODE.mode_rho{1}(:,1:4),info,'rows'),5:end);
                if mc.type == 1
                    mc.ld = MODE.mode_ld(ismember(MODE.mode_ld(:,1:4),info,'rows'),5:end);
                    mc.rho_2 = MODE.mode_rho{2}(ismember(MODE.mode_rho{2}(:,1:4),info,'rows'),5:end);
                    if isempty(mc.ld) || isempty(mc.rho_2)
                        error('incorrect child mode info mc[2]...\n')
                    end
                elseif mc.type == -1
                    mc.ld = [];
                    mc.rho_2 = [];
                end
                if mc.order == 2 && ~isempty(MODE.mode_rho{3})
                    loc = ismember(MODE.mode_rho{3}(:,1:4),info,'rows');
                    if any(loc)
                        mc.rho_3 = MODE.mode_rho{3}(loc,5:end);
                    else
                        mc.rho_3 = [];
                    end
                else
                    mc.rho_3 = [];
                end
            % 0th order has no coupling
            elseif mc.order == 0
                mc.mu = [];
                mc.ld = [];
                mc.rho_1 = [];
                mc.rho_2 = [];
                mc.rho_3 = [];
            end
            
            % ------------------ Propagator Matrix Method -----------------
            [RK] = runge_kutta(n_grid,mc,M,NV{i_order+1},Flag);
            
            % ---------------------- solving response ---------------------
            [resp,sol] = solver(RK,mc,M);
           
            RESP{i_h}(ismember(RESP{i_h}(:,1:4),info,'rows'),5:end) = resp;
            
            out = file_sol(mc.type,mc.l,mc.m,mc.order,mc.L0,mc.M0,mc.L1,mc.M1);
            save(out,'sol');
        end     % modes of each order of pert
    end     % order of pert  
end     % (l1,m1)'s

% ---------------------------- total response -----------------------------
% sum up 1st and 2nd order modal responses and output to file...
tmp = zeros(n_resp - M.n_harm,col_RESP);    % skip 0th order in each RESP
k = 1;
for i_h = 1:M.n_harm
    n = size(RESP{i_h},1) - 1;
    tmp(k:k+n-1,:) = RESP{i_h}(2:end,:);
    k = k + n;
end
% sorted unique response modes... 
md = sortrows(unique(tmp(:,1:4),'rows'),1:4);
% 0th order...
m0 = RESP{1}(1,:);
% create final response matrix...
n = size(md,1);
RESP_TOT = zeros(n+1,col_RESP+1);
% first row is always 0th order response...
RESP_TOT(1,:) = [1,m0];    % "1" is ordering
% loop over all higher order modes...
ind_p = 0;  % if p = -1
for i = 1:n
    % find [type,lc,mc] no matter what order of perturbation
    loc = ismember(tmp(:,2:4),md(i,2:4),'rows');
    tot = sum(tmp(loc,5:end),1);
    if any(tot)
        % "mean" order of perturbation
        p = sum(tmp(loc,1))/sum(loc);
        if p > 1 && p < 2
            p = 3;      % indicating this mode mixes 1st and 2nd order resposes
        elseif p < 1 || p > 2
            error('incorrect order of perturbation...\n');
        end
    else
        % indicating total response for this mode is 0. This may result
        % from the case in which only d_lambda exist, such that there is
        % no toroidal mode
        p = -1;   
        ind_p = 1;
    end
    RESP_TOT(i+1,:) = [i+1,p,md(i,2:4),tot];
end
% remove any p = -1 mode (zero total response)
if ind_p == 1
    RESP_TOT = RESP_TOT((RESP_TOT(:,2) ~= -1),:);
end

% write to a ASCII file...
fname_resp = [Output.dir_resp,'resp.',Output.prefix,'.dat'];
fid = fopen(fname_resp,'w');
str_format = '%4d%4d%4d%4d%4d%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n';
fprintf(fid,str_format,RESP_TOT');
fclose(fid);
% -------------------------------------------------------------------------

fprintf('--- Time Used: %10.6f s ---\n',cputime - t0);

% delete old case solution files...
% delete([Output.dir_sol,Output.prefix,'*.mat']);