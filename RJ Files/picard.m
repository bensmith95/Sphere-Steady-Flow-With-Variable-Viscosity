%% Picard iterative solver
% Solves system until correction < 0.1 
% in order to provide a close enough guess for Newton method
% Inputs:
% U - guess for U component        V - guess for V component
% W - guess for W component        T - guess for Temp.
% Pr - Prandtl number              lambda - non-dim sensitivity constant
% r - radius vector                beta - angle/height vector       
% i - current position index       alpha - under-relaxation parameter
% Outputs:
% U - improved guess for U component  V - improved guess for V component
% W - improved guess for W component  T - improved guess for Temp.

function [U,V,W,T] = picard(U,V,W,T,Pr,lambda,r,beta,i,alpha)
    %% initialise
    % Space marching parameters
    dbeta = beta(2)-beta(1); Nbeta = length(beta); 
    h1 = r(i)-r(i-1); h2 = r(i)-r(i-2);

    % set BCs
    VBC = V(1,i); TBC = T(1,i);

    % Set intitial guess
    Uo = U(:,i-1);  Vo = V(:,i-1);  Wo = W(:,i-1); To = T(:,i-1);
    
    %% iterate
    A = zeros(4*Nbeta); b = zeros(4*Nbeta,1); q = 1;
    while max(abs(q))>1e-1
        % Determine if W<0
        a = Wo>0;
        
        % Viscous terms
        mu0 = 1./(1+lambda*To);
        dT_dbeta = (To(3:end)-To(1:end-2))/(2*dbeta);
        dmu_dbeta = -lambda*mu0(2:end-1).^2.*dT_dbeta;

        % beta->Inf BCs
        A(1,1+Nbeta*0+[2 1 0]) = [-1 4 -3]; % U
        A(2,1+Nbeta*1) = 1; b(2) = VBC; % V
        A(3,1+Nbeta*2) = 1; % W
        A(4,1+Nbeta*3) = 1; b(4) = TBC; % T

        for j=2:Nbeta-1
            % Continuity
            A(1+(j-1)*4,j-1+Nbeta*0) = -1/(2*dbeta);
            A(1+(j-1)*4,j+1+Nbeta*0) = 1/(2*dbeta);
            A(1+(j-1)*4,j+Nbeta*2) = r(i)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2; 
            b(1+(j-1)*4) = r(i)*(h1^2*W(j,i-2)-h2^2*W(j,i-1))/(h1^2*h2-h1*h2^2);
            % V Momentum
            A(2+(j-1)*4,j-1+Nbeta*1) = -r(i)*Uo(j)/(2*dbeta) - mu0(j)/dbeta^2 + dmu_dbeta(j-1)/(2*dbeta);
            A(2+(j-1)*4,j+Nbeta*1) = a(j)*r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + r(i)*Wo(j) + 2*mu0(j)/dbeta^2;
            A(2+(j-1)*4,j+1+Nbeta*1) = r(i)*Uo(j)/(2*dbeta) - mu0(j)/dbeta^2 - dmu_dbeta(j-1)/(2*dbeta);
            b(2+(j-1)*4) = a(j)*r(i)^2*Wo(j)*(h1^2*V(j,i-2)-h2^2*V(j,i-1))/(h1^2*h2-h1*h2^2);
            % W Momentum
            A(3+(j-1)*4,j+Nbeta*1) = -r(i)*Vo(j);
            A(3+(j-1)*4,j-1+Nbeta*2) = -r(i)*Uo(j)/(2*dbeta) - mu0(j)/dbeta^2 + dmu_dbeta(j-1)/(2*dbeta);
            A(3+(j-1)*4,j+Nbeta*2) = a(j)*r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2*mu0(j)/dbeta^2;
            A(3+(j-1)*4,j+1+Nbeta*2) = r(i)*Uo(j)/(2*dbeta) - mu0(j)/dbeta^2 - dmu_dbeta(j-1)/(2*dbeta);
            b(3+(j-1)*4) = a(j)*r(i)^2*Wo(j)*(h1^2*W(j,i-2)-h2^2*W(j,i-1))/(h1^2*h2-h1*h2^2);
            % Heat
            A(4+(j-1)*4,j-1+Nbeta*3) = -r(i)*Uo(j)/(2*dbeta) - 1/(Pr*dbeta^2);
            A(4+(j-1)*4,j+Nbeta*3) = a(j)*r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2/(Pr*dbeta^2);
            A(4+(j-1)*4,j+1+Nbeta*3) = r(i)*Uo(j)/(2*dbeta) - 1/(Pr*dbeta^2);
            b(4+(j-1)*4) = a(j)*r(i)^2*Wo(j)*(h1^2*T(j,i-2)-h2^2*T(j,i-1))/(h1^2*h2-h1*h2^2);
        end

        % beta=0 BCs
        A(end-3,Nbeta*4-[2 1 0]) = [1 -4 3]; % T
        A(end-2,Nbeta*1) = 1; % U
        A(end-1,Nbeta*2-[2 1 0]) = [1 -4 3]; % V
        A(end,Nbeta*3-[2 1 0]) = [1 -4 3]; % W

        % Solve system
        s = A\b;

        % Under-relax
        Us = (1-alpha)*Uo+s([1:Nbeta]+Nbeta*0)*alpha;
        Vs = (1-alpha)*Vo+s([1:Nbeta]+Nbeta*1)*alpha;
        Ws = (1-alpha)*Wo+s([1:Nbeta]+Nbeta*2)*alpha;
        Ts = (1-alpha)*To+s([1:Nbeta]+Nbeta*3)*alpha;

        % Compute change in solution
        Q = abs([Us-Uo;Vs-Vo;Ws-Wo;Ts-To]);
        q = max(Q,[],'all');

        % set new guess
        Uo = Us; Vo = Vs; Wo = Ws; To = Ts;

    end
    % set new solution
    U(:,i) = Uo; V(:,i) = Vo; W(:,i) = Wo; T(:,i) = To;
end