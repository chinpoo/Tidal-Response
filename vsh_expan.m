% func: determine "child" modes and their VSH expansion coefficients for a
%       given "parent" mode and an eigenstructure in mu, lambda and rho
% input args: 
%     l0,m0: harmonic of the "parent" mode
%     l1,m1: harmonic of the eigenstructure
%     mode:  spheroidal (1) or toroidal (-1) of the "parent" mode
% outputs:
%     CP: a structure of "child" modes and VSH expan coefficients

function [CP] = vsh_expan(l0,m0,l1,m1,mode)

    global Output;

    if abs(m0)>l0 || abs(m1)>l1 || l0 < 0 || l1 < 0 || (mode ~= 1 && mode ~= -1)
        error('incorrect input arguments for VSH expansion...\n');
    else
        % ----------------------- Data Structure --------------------------
        % define a structure to store all the results...
        nc = 100;           % default max number of "child" modes from coupling
        CP = struct;        % CP stands for Coupling
        CP.cnt_mu = 0;      % # of "child" modes due to (l1,m1) in mu
        CP.mode_mu = zeros(nc,9);
        CP.cnt_ld = 0;      % # of "child" modes due to (l1,m1) in lambda
        if mode == 1        % only spheroidal "parent" causes "child" modes
            CP.mode_lambda = zeros(nc,5);
        else
            CP.mode_lambda = [];
        end
        CP.cnt_rho = zeros(3,1);      % # of "child" modes due to (l1,m1) in rho, for each of 3 terms
        CP.mode_rho = cell(1,3);
        CP.mode_rho{1} = zeros(nc,5);
        CP.mode_rho{2} = zeros(nc,5);
        % flag: net force is non-zero (1) or zero (0) from a specific coupling
        IS_NET_FORCE = 0;
        % -----------------------------------------------------------------
        
        % -------------------- Numerical Integration ----------------------
        % trapezoidal rule for integral of phi-dependent part, Simpson's
        % rule for integral of theta-dependent part
        N = 10000;           % # of division
        h = pi/N;           % step
        phi = 0:2*h:2*pi;   % phi vector
        theta = 0:h:pi;     % theta vector
        LL0 = l0*(l0+1);
        
        % (l0,m0)
        [f_0,fp_0,fpp_0] = phi_func(m0,phi);
        [plm_0,plmt_0,plmtt_0,plmds_0] = theta_func(l0,m0,theta);
        % perturbation (l1,m1)
        [f_1,fp_1,fpp_1] = phi_func(m1,phi);
        [plm_1,plmt_1,plmtt_1,plmds_1] = theta_func(l1,m1,theta);
        
        % define function handles...
        int_f = @trapez;
        int_t = @simpson;
        
        % loop over possible child modes, general selection rule applies...
        for lc = abs(l0-l1):(l0+l1)
            LL = lc*(lc+1);
            for mc = -lc:lc
                [f_c,fp_c,fpp_c] = phi_func(mc,phi);
                [plm_c,plmt_c,plmtt_c,plmds_c] = theta_func(lc,mc,theta);
                
                % phi components in integrands...
                F1 = f_0.*f_1.*f_c;
                F2 = fp_0.*fp_1.*f_c;
                F3 = f_0.*f_1.*fp_c;
                F4 = fp_0.*fp_1.*fp_c;
                F5 = fp_0.*f_1.*fp_c;
                F6 = fp_0.*f_1.*f_c;
                F7 = f_0.*fp_1.*fp_c;
                F8 = f_0.*fp_1.*f_c;
                F9 = fpp_0.*fp_1.*fp_c;
                F10 = fpp_0.*fp_1.*f_c;
                
                % theta components in integrands...
                T1 = plm_0.*plm_1.*plm_c.*sin(theta);
                T2 = plmt_0.*plmt_1.*plm_c.*sin(theta);
                T3 = plmds_0.*plmds_1.*plm_c.*sin(theta);                
                T4 = plmt_0.*plm_1.*plmt_c.*sin(theta);
                T5 = plm_0.*plmt_1.*plmt_c.*sin(theta);
                T6 = plmtt_0.*plmt_1.*plmt_c.*sin(theta);
                T7 = plmt_0.*plmds_1.*plmt_c;
                T8 = plmds_0.*plmds_1.*plmt_c.*cos(theta);
                T9 = plmt_0.*plm_1.*plmds_c.*sin(theta);
                T10 = plm_0.*plmt_1.*plmds_c.*sin(theta);
                T11 = plmtt_0.*plmt_1.*plmds_c.*sin(theta);
                T12 = plmt_0.*plmds_1.*plmds_c;
                T13 = plmds_0.*plmds_1.*plmds_c.*cos(theta);
                T14 = plmds_0.*plm_1.*plmds_c.*sin(theta);
                T15 = plm_0.*plmds_1.*plmds_c.*sin(theta);
                T16 = plmds_0.*plmt_1.*plmds_c.*cos(theta);
                T17 = plmt_0.*plmt_1.*plmds_c;
                T18 = plmds_0.*plmds_1.*plmds_c;
                T19 = plmt_0.*plmds_1.*plmds_c.*cos(theta);
                T20 = plmds_0.*plm_1.*plmt_c.*sin(theta);
                T21 = plm_0.*plmds_1.*plmt_c.*sin(theta);
                T22 = plmds_0.*plmt_1.*plmt_c.*cos(theta);
                T23 = plmt_0.*plmt_1.*plmt_c;
                T24 = plmds_0.*plmds_1.*plmt_c;
                T25 = plmt_0.*plmds_1.*plmt_c.*cos(theta);
                T26 = plmt_0.*plmds_1.*plm_c.*sin(theta);
                T27 = plmds_0.*plmt_1.*plm_c.*sin(theta);
                T28 = plmtt_0.*plmds_1.*plmt_c.*sin(theta);
                T29 = plmtt_0.*plmds_1.*plmds_c.*sin(theta);
                
                % if parent is spheroidal...
                if mode == 1
                    % ************** delta mu ***************
                    % Fr = (Upp+2Up/r-2U/r^2)*X1 + S/mu/r*X2 + V/r^2*X3
                    % Ft = (Sp/mu+3S/mu/r)*X4 + U/r^2*X5 + V/r^2*X6
                    % Ff = (Sp/mu+3S/mu/r)*X7 + U/r^2*X8 + V/r^2*X9
                    % VSH expan coeffs:(note that FP and FB are decoupled with FC)
                    % FP = Int{Fr*Ylm}dA = C1*(Upp+2Up/r-2U/r^2) + C2*S/mu/r + C3*V/r^2
                    % FB = Int{Ft*Ylmt+Ff*Ylmp/sint}dA/LL = C4*(Sp/mu+3S/mu/r) + C5*U/r^2 + C6*V/r^2
                    % FC = Int{-Ft*Ylmp/sint+Ff*Ylmt}dA/LL = C7*(Sp/mu+3S/mu/r) + C8*U/r^2 + C9*V/r^2     
                    % -------- FP --------
                    C = zeros(1,9);
                    C(1) = 2*int_f(F1)*int_t(T1);
                    C(2) = int_f(F1)*(int_t(T2)-LL0*int_t(T1))+int_f(F2)*int_t(T3);
                    C(3) = 2*LL0*int_f(F1)*int_t(T1);
                    % -------- FB & FC --------
                    if lc == 0
                        C(4:end) = 0;
                    else
                        % FB
                        C(4) = (int_f(F1)*int_t(T4)+int_f(F5)*int_t(T14))/LL;
                        C(5) = 2*C(4)+2*(int_f(F1)*int_t(T5)+int_f(F7)*int_t(T15))/LL;
                        c6_1 = 2*(int_f(F1)*int_t(T6)-(LL0-1)*int_f(F1)*int_t(T4)+int_f(F2)*int_t(T7)-int_f(F2)*int_t(T8))/LL;
                        c6_2 = 2*((1-LL0)*int_f(F5)*int_t(T14)-int_f(F5)*int_t(T16)+int_f(F5)*int_t(T17)+int_f(F9)*int_t(T18)+int_f(F7)*int_t(T19))/LL;
                        C(6) = c6_1+c6_2;
                        % FC
                        C(7) = (-int_f(F3)*int_t(T9)+int_f(F6)*int_t(T20))/LL;
                        C(8) = 2*C(7)+2*(-int_f(F3)*int_t(T10)+int_f(F8)*int_t(T21))/LL;
                        c9_1 = -2*(int_f(F3)*int_t(T11)-(LL0-1)*int_f(F3)*int_t(T9)+int_f(F4)*int_t(T12)-int_f(F4)*int_t(T13))/LL;
                        c9_2 = 2*((1-LL0)*int_f(F6)*int_t(T20)-int_f(F6)*int_t(T22)+int_f(F6)*int_t(T23)+int_f(F10)*int_t(T24)+int_f(F8)*int_t(T25))/LL;
                        C(9) = c9_1+c9_2;
                    end
                    C(abs(C) < 1e-8) = 0;
                    % store "child" mode info
                    if any(C(1:6))         % a spheroidal "child"
                        if any(C(7:9))     % FC must be zero if FP and FB are not
                            error('VSH expansion for d_mu is wrong...\n');
                        else
                            CP.cnt_mu = CP.cnt_mu + 1;
                            CP.mode_mu(CP.cnt_mu,1:end) = [1,lc,mc,C(1:6)];
                        end
                    else
                        if any(C(7:9))     % a toroidal "child"
                            CP.cnt_mu = CP.cnt_mu + 1;
                            CP.mode_mu(CP.cnt_mu,1:6) = [-1,lc,mc,C(7:9)];
                        end
                    end
                    
                    % ************ delta lambda **************
                    % X = Up+(2U-k^2*V)/r
                    % Fr = Xp*Y1
                    % Ft = X/r*Y2
                    % Ff = X/r*Y3
                    % VSH expan coeffs:
                    % FP = Int{Fr*Ylm}dA = D1*Xp
                    % FB = Int{Ft*Ylmt+Ff*Ylmp/sint}dA/LL = D2*X/r
                    % FC = Int{-Ft*Ylmp/sint+Ff*Ylmt}dA/LL = D3*X/r (D3 = 0)       
                    % -------- FP --------
                    D = zeros(1,3);
                    D(1) = int_f(F1)*int_t(T1);
                    % -------- FB & FC --------
                    if lc == 0
                        D(2:end) = 0;
                    else
                        d2_1 = (int_f(F1)*int_t(T4)+int_f(F5)*int_t(T14))/LL;
                        d2_2 = (int_f(F1)*int_t(T5)+int_f(F7)*int_t(T15))/LL;
                        D(2) = d2_1 + d2_2;
                        d3_1 = (-int_f(F3)*int_t(T9)+int_f(F6)*int_t(T20))/LL;
                        d3_2 = (-int_f(F3)*int_t(T10)+int_f(F8)*int_t(T21))/LL;
                        D(3) = d3_1 + d3_2;
                    end
                    D(abs(D) < 1e-8) = 0;
                    
                    if any(D)
                        if D(3) ~= 0   % FC must be zero for d_lambda
                            error('VSH expansion for d_lambda is wrong...\n');
                        else
                            CP.cnt_ld = CP.cnt_ld + 1;
                            CP.mode_lambda(CP.cnt_ld,1:end) = [1,lc,mc,D(1:2)];
                        end
                    end
                    
                    % ************ delta rho **************
                    % term #1: -d_rho*grad(psi_1) or -d_rho*grad(psi+V_td)
                    % Fr_1 = -d_rho_0*Kp*Z1_1
                    % Ft_1 = -d_rho_0*(K/r)*Z1_2
                    % Ff_1 = -d_rho_0*(K/r)*Z1_3
                    % VSH expan coeffs:
                    % FP_1 = Int{Fr_1*Ylm}dA = E1_1*d_rho_0*Kp
                    % FB_1 = Int{Ft_1*Ylmt+Ff_1*Ylmp/sint}dA/LL = E1_2*d_rho_0*(K/r)
                    % FC_1 = Int{-Ft_1*Ylmp/sint+Ff_1*Ylmt}dA/LL = E1_3*d_rho_0*(K/r)
                    % -------- FP --------
                    E1 = zeros(1,3);
                    E1(1) = -int_f(F1)*int_t(T1);
                    % -------- FB & FC --------
                    if lc == 0
                        E1(2:end) = 0;
                    else
                        % FB_1
                        E1(2) = -(int_f(F1)*int_t(T4)+int_f(F5)*int_t(T14))/LL;
                        % FC_1
                        E1(3) = (int_f(F3)*int_t(T9)-int_f(F6)*int_t(T20))/LL;
                    end
                    E1(abs(E1) < 1e-8) = 0;

                    if any(E1(1:2))
                        if E1(3)~=0     % FC must be zero if FP and FB are not
                            error('VSH expansion for term #1 of d_rho is wrong...\n');
                        else
                            CP.cnt_rho(1) = CP.cnt_rho(1) + 1;
                            CP.mode_rho{1}(CP.cnt_rho(1),1:end) = [1,lc,mc,E1(1:2)];
                            % ------------------ Net Force ------------------
                            % determine if net force can be induced (when lc == 1)
                            if lc == 1
                                if IS_NET_FORCE == 0
                                    IS_NET_FORCE = 1;
                                    CP.mode_rho{3} = zeros(nc,5);
                                    % first row stores (order_1,1,lc,mc)
                                    CP.mode_rho{3}(1,1:4) = [1,1,lc,mc];
                                    a.l = lc;
                                    a.m = mc;
                                    a.phi = {f_c,fp_c};
                                    a.theta = {plm_c,plmt_c,plmds_c};
                                else
                                    error('Incorrect: >1 degree-1 modes generated...\n');
                                end 
                            end
                            % ------------------------------------------------
                        end
                    else
                        if E1(3)~=0
                            CP.cnt_rho(1) = CP.cnt_rho(1) + 1;
                            CP.mode_rho{1}(CP.cnt_rho(1),1:4) = [-1,lc,mc,E1(3)];
                        end
                    end
                    
                    % term #2: div(d_rho*u)*g0e_r or in Poisson's eqn
                    % F_2 = d_rho_0*X*Z2_1 + d_rho_0*V/r*Z2_2
                    % SH expan
                    % G =  Int{F_2*Ylm}dA = d_rho_0*(E2_1*X + E2_2*V/r)
                    E2 = zeros(1,2);
                    E2(1) = int_f(F1)*int_t(T1);
                    E2(2) = int_f(F1)*int_t(T2)+int_f(F2)*int_t(T3);
                    E2(abs(E2) < 1e-8) = 0;
                    if any(E2)
                        CP.cnt_rho(2) = CP.cnt_rho(2) + 1;
                        CP.mode_rho{2}(CP.cnt_rho(2),1:end) = [1,lc,mc,E2(1:2)];
                    end

                % if parent is toroidal...
                elseif mode == -1
                    % ************** delta mu ***************
                    % Fr = T/mu/r*X1
                    % Ft = (Tp/mu+3T/mu/r)*X2 + W/r^2*X3
                    % Ff = (Tp/mu+3T/mu/r)*X4 + W/r^2*X5
                    % VSH expan coeffs:(note that FP and FB are decoupled with FC)
                    % FP = Int{Fr*Ylm}dA = C1*T/mu/r
                    % FB = Int{Ft*Ylmt+Ff*Ylmp/sint}dA = C2*(Tp/mu+3T/mu/r) + C3*W/r^2
                    % FC = Int{-Ft*Ylmp/sint+Ff*Ylmt}dA = C4*(Tp/mu+3T/mu/r) + C5*W/r^2
                    % -------- FP --------
                    C = zeros(1,5);
                    C(1) = int_f(F8)*int_t(T26)-int_f(F6)*int_t(T27);
                    % -------- FB & FC --------
                    if lc == 0
                        C(2:end) = 0;
                    else
                        % FB
                        C(2) = (-int_f(F6)*int_t(T20)+int_f(F3)*int_t(T9))/LL;
                        c3_1 = (int_f(F6)*((LL0-2)*int_t(T20)+2*int_t(T22)-2*int_t(T23))+int_f(F8)*(LL0*int_t(T21)+2*int_t(T28)))/LL;
                        c3_2 = (int_f(F3)*((2-LL0)*int_t(T9)+LL0*int_t(T10)+2*int_t(T11))+2*int_f(F4)*(int_t(T12)-int_t(T13)))/LL;
                        C(3) = c3_1+c3_2;
                        % FC
                        C(4) = (int_f(F5)*int_t(T14)+int_f(F1)*int_t(T4))/LL;
                        c5_1 = -(int_f(F5)*((LL0-2)*int_t(T14)+2*int_t(T16)-2*int_t(T17))+int_f(F7)*(LL0*int_t(T15)+2*int_t(T29)))/LL;
                        c5_2 = (int_f(F1)*((2-LL0)*int_t(T4)+LL0*int_t(T5)+2*int_t(T6))+2*int_f(F2)*(int_t(T7)-int_t(T8)))/LL;
                        C(5) = c5_1+c5_2;
                    end
                    C(abs(C) < 1e-8) = 0;
                
                    if any(C(1:3))
                        if any(C(4:5))
                            error('VSH expansion for d_mu is wrong...\n');
                        else
                            CP.cnt_mu = CP.cnt_mu + 1;
                            CP.mode_mu(CP.cnt_mu,1:6) = [1,lc,mc,C(1:3)];
                        end
                    else
                        if any(C(4:5))
                            CP.cnt_mu = CP.cnt_mu + 1;
                            CP.mode_mu(CP.cnt_mu,1:5) = [-1,lc,mc,C(4:5)];
                        end
                    end
                    
                    % ************ delta lambda ***************
                    % d_lambda does not induce toroidal field
                    % Fr = N/A
                    % Ft = N/A
                    % Ff = N/A
                    
                    % ************ delta rho ***************
                    % term #1 disappears for toroidal "parent"
                    
                    % term #2
                    % F_2 = d_rho_0*W/r*Z
                    % Q = Int{F_2*Ylm}dA = d_rho_0*E*W/r
                    E = int_f(F8)*int_t(T26)-int_f(F6)*int_t(T27);
                    if abs(E) > 1e-8
                        CP.cnt_rho(2) = CP.cnt_rho(2) + 1;
                        CP.mode_rho{2}(CP.cnt_rho(2),1:4) = [1,lc,mc,E];
                    end
                end
            end     % end of mc loop
        end     % end of lc loop
        
        CP.mode_mu = CP.mode_mu(1:CP.cnt_mu,:);
        CP.mode_lambda = CP.mode_lambda(1:CP.cnt_ld,:);
        CP.mode_rho{1} = CP.mode_rho{1}(1:CP.cnt_rho(1),:);
        CP.mode_rho{2} = CP.mode_rho{2}(1:CP.cnt_rho(2),:);
        
        % if net force exist, solve mode coupling for second-order coupling
        % d_rho*a_tot
        if IS_NET_FORCE == 1 && ~isempty(CP.mode_rho{3})
            for lh = abs(a.l-l1):(a.l+l1)
                llh = lh*(lh+1);
                for mh = -lh:lh
                    [f_h,fp_h,fpp_h] = phi_func(mh,phi);
                    [plm_h,plmt_h,plmtt_h,plmds_h] = theta_func(lh,mh,theta);
                       
                    F01 = a.phi{1}.*f_1.*f_h;
                    F02 = a.phi{1}.*f_1.*fp_h;
                    F03 = a.phi{2}.*f_1.*fp_h;
                    F04 = a.phi{2}.*f_1.*f_h; 
                    
                    T01 = a.theta{1}.*plm_1.*plm_h.*sin(theta);
                    T02 = a.theta{2}.*plm_1.*plmt_h.*sin(theta);
                    T03 = a.theta{2}.*plm_1.*plmds_h.*sin(theta);
                    T04 = a.theta{3}.*plm_1.*plmds_h.*sin(theta);
                    T05 = a.theta{3}.*plm_1.*plmt_h.*sin(theta);
                    
                    % Fh_r = d_rho_0*A_P*Zh_1
                    % Fh_t = d_rho_0*A_B*Zh_2
                    % Fh_f = d_rho_0*A_B*Zh_3
                    % VSH expan coeffs:
                    % FP_h = Int{Fh_r*Ylm}dA = Eh_1*d_rho_0*A_P
                    % FB_h = Int{Fh_t*Ylmt+Fh_f*Ylmp/sint}dA/LL = Eh_2*d_rho_0*A_B
                    % FC_h = Int{-Fh_t*Ylmp/sint+Fh_f*Ylmt}dA/LL = Eh_3*d_rho_0*A_B
                    % -------- FP --------
                    Eh = zeros(1,3);
                    Eh(1) = int_f(F01)*int_t(T01);
                    % -------- FB & FC --------
                    if lh == 0
                        Eh(2:end) = 0;
                    else
                        Eh(2) = (int_f(F01)*int_t(T02)+int_f(F03)*int_t(T04))/llh;
                        Eh(3) = (-int_f(F02)*int_t(T03)+int_f(F04)*int_t(T05))/llh;
                    end
                    Eh(abs(Eh) < 1e-8) = 0;
                    
                    if any(Eh(1:2))
                        if Eh(3)~=0    % FC must be zero if FP and FB are not
                            error('VSH expan of net force coupling term is wrong...\n');
                        else
                            CP.cnt_rho(3) = CP.cnt_rho(3) + 1;
                            CP.mode_rho{3}(CP.cnt_rho(3)+1,1:end) = [1,lh,mh,Eh(1:2)];
                        end
                    else
                        if Eh(3)~=0
                            CP.cnt_rho(3) = CP.cnt_rho(3) + 1;
                            CP.mode_rho{3}(CP.cnt_rho(3)+1,1:4) = [-1,lh,mh,Eh(3)];
                        end
                    end
                end     % end of mh
            end     % end of lh
            CP.mode_rho{3} = CP.mode_rho{3}(1:CP.cnt_rho(3)+1,:);
        end
        
        % -------- Write VSH to File --------
        fname = file_vsh(l0,m0,l1,m1,mode,Output.dir_vsh);
        save(fname,'CP');
        
    end     % end of main routine


    % ------------------------- Nested Functions --------------------------
    % trapezoidal rule...
    function [int] = trapez(f_phi)
        int = (2*h)/2*(f_phi(1)+f_phi(N+1)+2*sum(f_phi(2:N)));
    end
    % Simpson's rule...
    function [int] = simpson(f_theta)
        int = h/3*(f_theta(1)+f_theta(N+1)+4*sum(f_theta(2:2:N))+2*sum(f_theta(3:2:N-1)));
    end
    % ---------------------------------------------------------------------
    
end

