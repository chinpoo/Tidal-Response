% func: create A matrix in matrix equation dX/dr = A*X + F
% input args:
%           r: radius
%           g0: g0 at radius r
%           l_mat: indexing of material layer
%           L: harmonic degree of the mode
%           mode: spheroidal (1) or toroidal (-1)
%           model: structure of the model
% outputs:
%           a: squre A matrix

function [a] = create_a_matrix(r,g0,l_mat,L,mode,model)
    rho = model.rho(l_mat);
    mu = model.mu(l_mat);
    lambda = model.lambda(l_mat);
    beta = model.beta(l_mat);
    gamma = model.gamma(l_mat);
    coef = model.eta;

    % spheroidal:
    if mode == 1
        % for (0,0) mode, A matrix is 3-by-3
        if L == 0
            a = zeros(3,3);
            
            a(1,1) = -2*lambda/r/beta;
            a(1,2) = 1/beta;
            a(2,1) = 4/r*(gamma/r-coef*rho*g0);
            a(2,2) = -4*mu/r/beta;
            a(3,1) = -rho;
        % A matrix is 6-by-6
        else
            a = zeros(6,6);
            
            a(1,1) = -2*lambda/r/beta;
            a(1,2) = L*(L+1)*lambda/r/beta;
            a(1,3) = 1/beta;
            a(2,1) = -1/r;
            a(2,2) = 1/r;
            a(2,4) = 1/mu;
            a(3,1) = 4/r*(gamma/r-coef*rho*g0);
            a(3,2) = -L*(L+1)/r*(2*gamma/r-coef*rho*g0);
            a(3,3) = -4*mu/r/beta;
            a(3,4) = L*(L+1)/r;
            a(3,5) = -coef*rho*(L+1)/r;
            a(3,6) = coef*rho;
            a(4,1) = 1/r*(coef*rho*g0-2*gamma/r);
            a(4,2) = -1/r/r*(2*mu-L*(L+1)*(gamma+mu));
            a(4,3) = -lambda/r/beta;
            a(4,4) = -3/r;
            a(4,5) = coef*rho/r;
            a(5,1) = -rho;
            a(5,5) = -(L+1)/r;
            a(5,6) = 1;
            a(6,1) = -(L+1)*rho/r;
            a(6,2) = L*(L+1)*rho/r;
            a(6,6) = (L-1)/r;
        end
    % toroidal: A matrix is 2-by-2
    elseif mode == -1
        a = zeros(2,2);
        
        a(1,1) = 1/r;
        a(1,2) = 1/mu;
        a(2,1) = (L+2)*(L-1)*mu/r^2;
        a(2,2) = -3/r;
    end
    
end % function