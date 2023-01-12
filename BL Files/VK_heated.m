%% Solves Heated Von Karman flow for rotating disk
% Inputs:
% Pr - Prandtl number                  lambda - non-dim sensitivity constant
% Zmax - maximum value to solve along  N - number of discretisation points
% Outputs:
% Vel - cell array of velocities, pressure and temperature distribution
% z - vector which von karman equations are solved upon
function [Vel,z] = VK_heated(Pr,lambda,Zmax,N)
    % determine log spacing, so better precision near sphere surface
    z = [0,logspace(-2,log10(Zmax),N)];

    % solve von karman system
    options = bvpset('RelTol',1e-8,'AbsTol',1e-10,'NMax',10^4); % set tolerances
    solinit = bvpinit(z,[1 0 1 0 1 1 0]); % set initial guess
    % numerical solver
    Y = bvp4c(@(x,U) VK_system(x,U,Pr,lambda),@(Ua,Ub) VK_BC(Ua,Ub),solinit,options);
    Sol = deval(Y,z)';

    % save solution to cell array
    for i = 1:7
        Vel{i} = Sol(:,i);
    end
end

%% Set up von karman/O(1) Banks system
function dudz = VK_system(x,U,Pr,lambda)
    mu = 1./(1+lambda*U(6));
    dudz = [U(2); lambda*mu*U(7)*U(2)+(U(1)^2-U(3)^2+U(5)*U(2))/mu;
            U(4); lambda*mu*U(7)*U(4)+(U(5)*U(4)+2*U(1)*U(3))/mu;
            -2*U(1);
            U(7); Pr*U(5)*U(7)];
end

%% Set boundary conditions
function BC = VK_BC(Ua,Ub)
    BC = [Ua(1);Ub(1);Ua(3)-1;Ub(3);Ua(5);Ua(6)-1;Ub(6)];
end