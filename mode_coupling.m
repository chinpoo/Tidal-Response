% func: create full mode coupling hierarchy, including "child" modes up to 
%       2nd order of perturbation and the associated coupling coefficients.
% input args: 
%     l0,m0: harmonic of the 0th order "parent" mode (spheroidal) 
%     l1,m1: harmonic of the eigenstructure
%     dir: directory where mode coupling files are placed
% outputs:
%     MODE: a structure containing all mode coupling information

function [MD] = mode_coupling(l0,m0,l1,m1,dir)
    
    order_1 = 1;
    order_2 = 2;
    % ----------------------- Data Structure --------------------------
    MD = struct;
    nm = 500;   % default max number of "child" modes from 1st and 2nd order coupling
    MD.cnt_mu = 0;
    MD.cnt_ld = 0;
    MD.cnt_rho = zeros(1,3);
    MD.modes = zeros(nm,4);     % [order,mode,l,m] of all 1st and 2nd order modes
    MD.N = 0;                   % total number of unique 1st and 2nd order modes
    % following arrays contain [order,mode_c,lc,mc,mode_p,lp,mp,expan_coeffs]
    MD.mode_mu = zeros(nm,13);
    MD.mode_ld = zeros(nm,9);
    MD.mode_rho = cell(1,3);
    MD.mode_rho{1} = zeros(nm,9);
    MD.mode_rho{2} = zeros(nm,9);
    % -----------------------------------------------------------------
    
    % 1st order of perturbation...
    mode_p = 1;     % 0th order mode is always spheroidal
    CP = get_coupling_info(l0,m0,l1,m1,mode_p,dir);     % 1st "child" modes
    read_coupling_info(CP,order_1,[mode_p,l0,m0]);      % read coupling info from CP into MD
    n1 = MD.cnt_mu;     % # of all 1st order modes
    
    % 2nd order of perturbation...
    for i = 1:n1
        mode_1 = MD.modes(i,2);
        lp_1 = MD.modes(i,3);
        mp_1 = MD.modes(i,4);
        CP = get_coupling_info(lp_1,mp_1,l1,m1,mode_1,dir);     % 2nd "child" modes
        read_coupling_info(CP,order_2,[mode_1,lp_1,mp_1]);      % read coupling info from CP into MD      
    end
    
    MD.modes = unique(MD.modes,'rows');
    if ~any(MD.modes(1,:))
        % first row stores 0th order
        MD.modes(1,:) = [0,1,l0,m0];
    end
    MD.N = size(MD.modes,1);
    MD.mode_mu = MD.mode_mu(1:MD.cnt_mu,:);
    MD.mode_ld = MD.mode_ld(1:MD.cnt_ld,:);
    MD.mode_rho{1} = MD.mode_rho{1}(1:MD.cnt_rho(1),:);
    MD.mode_rho{2} = MD.mode_rho{2}(1:MD.cnt_rho(2),:);
    % --------------------------------------------------------------------

    % nested function
    % get coupling info from existing file or vsh_expan...
    function [cp] = get_coupling_info(lp,mp,l1,m1,mode,dir)
        fname = file_vsh(lp,mp,l1,m1,mode,dir);
        if exist(fname,'file') == 2
            var = load(fname,'CP');     % CP is the name of the structure
            cp = var.CP;
        else
            cp = vsh_expan(lp,mp,l1,m1,mode);
        end
    end

    % read coupling info from CP to MD...
    % input args:
    %           cp: mode coupling info from single coupling
    %           order: order of perturbation, 1 or 2
    %           parent: parent mode info: [mode_p,lp,mp]
    function read_coupling_info(cp,order,parent)
        % -------- mu --------
        ind_start = MD.cnt_mu + 1;
        ind_end = MD.cnt_mu + cp.cnt_mu;
        MD.cnt_mu = ind_end;
        MD.mode_mu(ind_start:ind_end,[1,5:7]) = ones(cp.cnt_mu,1)*[order,parent];
        MD.mode_mu(ind_start:ind_end,[2:4,8:end]) = cp.mode_mu;
        % all "child" mode
        MD.modes(ind_start:ind_end,:) = [order*ones(cp.cnt_mu,1),cp.mode_mu(:,1:3)];
        
        % -------- lambda --------
        if parent(1) == 1       % spheroidal
            ind_start = MD.cnt_ld + 1;
            ind_end = MD.cnt_ld + cp.cnt_ld;
            MD.cnt_ld = ind_end;
            MD.mode_ld(ind_start:ind_end,[1,5:7]) = ones(cp.cnt_ld,1)*[order,parent];
            MD.mode_ld(ind_start:ind_end,[2:4,8:end]) = cp.mode_lambda;
        elseif parent(1) == -1
            if cp.cnt_ld ~= 0 || ~isempty(cp.mode_lambda)
                error('incorrect coupling info from d_lambda...\n');
            end
        end
        
        % -------- rho --------
        % term #1
        if parent(1) == 1       % spheroidal
            ind_start = MD.cnt_rho(1) + 1;
            ind_end = MD.cnt_rho(1) + cp.cnt_rho(1);
            MD.cnt_rho(1) = ind_end;
            MD.mode_rho{1}(ind_start:ind_end,[1,5:7]) = ones(cp.cnt_rho(1),1)*[order,parent];
            MD.mode_rho{1}(ind_start:ind_end,[2:4,8:end]) = cp.mode_rho{1};
        elseif parent(1) == -1
            if cp.cnt_rho(1) ~= 0 || ~isempty(cp.mode_rho{1})
                error('incorrect coupling info from d_rho term #1...\n');
            end
        end
        % term #2
        ind_start = MD.cnt_rho(2) + 1;
        ind_end = MD.cnt_rho(2) + cp.cnt_rho(2);
        MD.cnt_rho(2) = ind_end;
        MD.mode_rho{2}(ind_start:ind_end,[1,5:7]) = ones(cp.cnt_rho(2),1)*[order,parent];
        MD.mode_rho{2}(ind_start:ind_end,[2:4,8:end]) = cp.mode_rho{2};
        % term #3: 2nd order term involing a_tot
        if order == 1   % read only from 1st order 'CP'
            if cp.cnt_rho(3) && ~isempty(cp.mode_rho{3})
                p = cp.mode_rho{3}(1,1:4);
                if ismember(p,MD.mode_mu(1:MD.cnt_mu,1:4))
                    MD.cnt_rho(3) = cp.cnt_rho(3);
                    MD.mode_rho{3} = zeros(MD.cnt_rho(3),9);
                    MD.mode_rho{3}(:,1) = order+1;
                    MD.mode_rho{3}(:,5:7) = ones(cp.cnt_rho(3),1)*p(2:end);
                    MD.mode_rho{3}(:,[2:4,8:end]) = cp.mode_rho{3}(2:end,:);
                else
                    error('incorrect coupling info from d_rho term #3...\n')
                end   
            else
                MD.mode_rho{3} = [];
            end
        end
    end     % end of nested func
        
end     % end of func