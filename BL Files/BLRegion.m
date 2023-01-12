%% Solves the boundary layer equations with a temperature dependant viscosity 
% of the form mu(T)=1/(1+lambda*T)
% Inputs:
% Pr - Prandtl number
% lambda - non-dim sensitivity constant

function BLRegion(Pr,lambda)
%% Pre-Condition
    % Space marching parameters
    eta_max = 30; Neta = 301; Ntheta = 1501;
    eta = linspace(0,eta_max,Neta); deta = eta(2)-eta(1);
    theta = linspace(0,pi/2,Ntheta); dtheta = theta(2)-theta(1);

    % Solve Von Karman equations
    str1 = fprintf('Solving Heated Von Karman equations...\n');
    [Vel,z] = VK_heated(Pr,lambda,eta_max+1,Neta);
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solved Heated Von Karman equations.\n');

    % assign velocity profiles
    U1 = Vel{1}; V1 = Vel{3}; W1 = Vel{5}; T1 = Vel{6}; 

    % initialise field
    % field discretised as theta x eta, or M_ij = [theta(i),eta(j)]
    U = repmat(theta',1,Neta).*repmat(spline(z,U1,eta),Ntheta,1);
    V = repmat(theta',1,Neta).*repmat(spline(z,V1,eta),Ntheta,1);
    W = repmat(spline(z,W1,eta),Ntheta,1);
    T = repmat(spline(z,T1,eta),Ntheta,1);

%% NEWTON RAPHSON MULTIVARIATE SCHEME
    str2 = fprintf('Solving Heated Boundary Layer equations at theta =  %.2f\n',0);
    for i=3:Ntheta
        deg = theta(i)*180/pi; fprintf(repmat('\b',1,str2)); 
        str2 = fprintf('Solving Heated Boundary Layer equations at theta =  %.2f\n',deg);

        % Use previous latitude as an intitial guess
        Uo = U(i-1,:);  Vo = V(i-1,:);  Wo = W(i-1,:); To = T(i-1,:);

        % Iterate
        gg = 0.95; q = 1;
        A = zeros(4*Neta); B = zeros(4*Neta,1);
        while max(abs(q))>1e-5
            mu0 = 1./(1+lambda*To);
            % FINITE DIFFERENCES (excluding boundaries)
            % 2nd order centered differences
            dU_deta = (Uo(3:end)-Uo(1:end-2))/(2*deta);
            dV_deta = (Vo(3:end)-Vo(1:end-2))/(2*deta);
            dW_deta = (Wo(3:end)-Wo(1:end-2))/(2*deta); 
            dT_deta = (To(3:end)-To(1:end-2))/(2*deta);
            dmu_deta = -lambda*mu0(2:end-1).^2.*dT_deta;

            % 2nd order three point backwards differences
            dU_dtheta = (3*Uo(2:end-1)-4*U(i-1,2:end-1)+U(i-2,2:end-1))/(2*dtheta);
            dV_dtheta = (3*Vo(2:end-1)-4*V(i-1,2:end-1)+V(i-2,2:end-1))/(2*dtheta);
            dT_dtheta = (3*To(2:end-1)-4*T(i-1,2:end-1)+T(i-2,2:end-1))/(2*dtheta);

            % 2nd order centered differences
            d2V_deta2 = (Vo(3:end)-2*Vo(2:end-1)+Vo(1:end-2))/deta^2;
            d2U_deta2 = (Uo(3:end)-2*Uo(2:end-1)+Uo(1:end-2))/deta^2;
            d2T_deta2 = (To(3:end)-2*To(2:end-1)+To(1:end-2))/deta^2;

            % checking how well we satisfy the equations at the given i-th
            % location in the interior points
            RHS_cont = -(dU_dtheta + dW_deta + Uo(2:end-1)*cot(theta(i)));
            RHS_U = dmu_deta.*dU_deta + mu0(2:end-1).*d2U_deta2 - Uo(2:end-1).*dU_dtheta - Wo(2:end-1).*dU_deta + Vo(2:end-1).^2*cot(theta(i));
            RHS_V = dmu_deta.*dV_deta + mu0(2:end-1).*d2V_deta2 - Uo(2:end-1).*dV_dtheta - Wo(2:end-1).*dV_deta - Uo(2:end-1).*Vo(2:end-1)*cot(theta(i));
            RHS_T = 1/Pr*d2T_deta2 - Uo(2:end-1).*dT_dtheta - Wo(2:end-1).*dT_deta;

            % JACOBIAN CONSTRUCTION
            % wall boundary conditions
            A(1,1+Neta*0) = 1;
            B(1,1) = -Uo(1);

            A(2,1+Neta*1) = 1;
            B(2,1) = -Vo(1) + sin(theta(i));

            A(3,1+Neta*2) = 1;
            B(3,1) = -Wo(1);

            A(4,1+Neta*3) = 1;
            B(4,1) = -To(1)+1;
            
            % internal points
            for j=2:Neta-1
                % continuity
                A(1+(j-1)*4,j+Neta*0) = 3/(2*dtheta) + cot(theta(i));
                A(1+(j-1)*4,j+1+Neta*2) = 1/(2*deta);
                A(1+(j-1)*4,j-1+Neta*2) = -1/(2*deta);
                B(1+(j-1)*4,1) = RHS_cont(j-1);
                % U momentum
                A(2+(j-1)*4,j+Neta*0) = dU_dtheta(j-1) + 3*Uo(j)/(2*dtheta) - (-2*mu0(j)/deta^2);
                A(2+(j-1)*4,j+1+Neta*0) = (Wo(j)-dmu_deta(j-1))/(2*deta) - (mu0(j)/deta^2);
                A(2+(j-1)*4,j-1+Neta*0) = -(Wo(j)-dmu_deta(j-1))/(2*deta) - (mu0(j)/deta^2);
                A(2+(j-1)*4,j+Neta*1) = -2*Vo(j)*cot(theta(i));
                A(2+(j-1)*4,j+Neta*2) = dU_deta(j-1);
                A(2+(j-1)*4,j+Neta*3) = lambda*mu0(j)^2*d2U_deta2(j-1) - 2*lambda^2*mu0(j)^3*dT_deta(j-1)*dU_deta(j-1);
                A(2+(j-1)*4,j+1+Neta*3) = lambda*mu0(j)^2*dU_deta(j-1)/(2*deta);
                A(2+(j-1)*4,j-1+Neta*3) = -lambda*mu0(j)^2*dU_deta(j-1)/(2*deta);
                B(2+(j-1)*4,1) = RHS_U(j-1);
                % V momentum
                A(3+(j-1)*4,j+Neta*0) = dV_dtheta(j-1) + Vo(j)*cot(theta(i));
                A(3+(j-1)*4,j+Neta*1) = 3*Uo(j)/(2*dtheta) + Uo(j)*cot(theta(i)) - (-2*mu0(j)/deta^2);
                A(3+(j-1)*4,j+1+Neta*1) = (Wo(j)-dmu_deta(j-1))/(2*deta) - (mu0(j)/deta^2);
                A(3+(j-1)*4,j-1+Neta*1) = -(Wo(j)-dmu_deta(j-1))/(2*deta) - (mu0(j)/deta^2);
                A(3+(j-1)*4,j+Neta*2) = dV_deta(j-1);
                A(3+(j-1)*4,j+Neta*3) = lambda*mu0(j)^2*d2V_deta2(j-1) - 2*lambda^2*mu0(j)^3*dT_deta(j-1)*dV_deta(j-1);
                A(3+(j-1)*4,j+1+Neta*3) = lambda*mu0(j)^2*dV_deta(j-1)/(2*deta);
                A(3+(j-1)*4,j-1+Neta*3) = -lambda*mu0(j)^2*dV_deta(j-1)/(2*deta);       
                B(3+(j-1)*4,1) = RHS_V(j-1);  
                % T momentum
                A(4+(j-1)*4,j+Neta*0) = dT_dtheta(j-1);
                A(4+(j-1)*4,j+Neta*2) = dT_deta(j-1);
                A(4+(j-1)*4,j+Neta*3) = 3*Uo(j)/(2*dtheta) - (-2/Pr/deta^2);            
                A(4+(j-1)*4,j+1+Neta*3) = Wo(j)/(2*deta) - (1/Pr/deta^2);
                A(4+(j-1)*4,j-1+Neta*3) = -Wo(j)/(2*deta) - (1/Pr/deta^2);
                B(4+(j-1)*4,1) = RHS_T(j-1);
            end

            % top boundary condition
            A(end-3,Neta*4) = 1;
            B(end-3,1) = -To(end);

            A(end-2,Neta) = 1;
            B(end-2,1) = -Uo(end);

            A(end-1,Neta*2) = 1;
            B(end-1,1) = -Vo(end);

            A(end,3*Neta-[2 1 0]) = [1 -4 3];
            B(end,1) = -(3*Wo(end)-4*Wo(end-1)+Wo(end-2));

            % solve for correction, q
            q = A\B;
            % Newton-Raphson iteration
            Uo = Uo+q([1:Neta]+Neta*0)'*gg;
            Vo = Vo+q([1:Neta]+Neta*1)'*gg;
            Wo = Wo+q([1:Neta]+Neta*2)'*gg;
            To = To+q([1:Neta]+Neta*3)'*gg;        
        end
        % set solution
        U(i,:) = Uo; V(i,:) = Vo; W(i,:) = Wo; T(i,:) = To;
    end
    fprintf(repmat('\b',1,str2)); str2 = fprintf('All latitudes solved for.\n');
    
%% Post Process 
    str3 = fprintf('Post-Processing...\n');

    % Determine pressure
    dP = U.^2+V.^2;
    P = cumtrapz(flip(eta),flip(dP,2),2); P = flip(P,2);

    % determine 2nd order finite differences for derivatives 
    dU_dtheta = [-(3*U(1,:)-4*U(2,:)+U(3,:)); U(3:end,:)-U(1:end-2,:); (3*U(end,:)-4*U(end-1,:)+U(end-2,:))]/2/dtheta;
    dV_dtheta = [-(3*V(1,:)-4*V(2,:)+V(3,:)); V(3:end,:)-V(1:end-2,:); (3*V(end,:)-4*V(end-1,:)+V(end-2,:))]/2/dtheta;
    dW_dtheta = [-(3*W(1,:)-4*W(2,:)+W(3,:)); W(3:end,:)-W(1:end-2,:); (3*W(end,:)-4*W(end-1,:)+W(end-2,:))]/2/dtheta;
    dP_dtheta = [-(3*P(1,:)-4*P(2,:)+P(3,:)); P(3:end,:)-P(1:end-2,:); (3*P(end,:)-4*P(end-1,:)+P(end-2,:))]/2/dtheta;
    dT_dtheta = [-(3*T(1,:)-4*T(2,:)+T(3,:)); T(3:end,:)-T(1:end-2,:); (3*T(end,:)-4*T(end-1,:)+T(end-2,:))]/2/dtheta;

    dU = [-(3*U(:,1)-4*U(:,2)+U(:,3)), U(:,3:end)-U(:,1:end-2), (3*U(:,end)-4*U(:,end-1)+U(:,end-2))]/2/deta;
    dV = [-(3*V(:,1)-4*V(:,2)+V(:,3)), V(:,3:end)-V(:,1:end-2), (3*V(:,end)-4*V(:,end-1)+V(:,end-2))]/2/deta;
    dW = [-(3*W(:,1)-4*W(:,2)+W(:,3)), W(:,3:end)-W(:,1:end-2), (3*W(:,end)-4*W(:,end-1)+W(:,end-2))]/2/deta;
    dT = [-(3*T(:,1)-4*T(:,2)+T(:,3)), T(:,3:end)-T(:,1:end-2), (3*T(:,end)-4*T(:,end-1)+T(:,end-2))]/2/deta;

    % organise data into a cell array
    clear Vel
    VelBL{1} = U; VelBL{2} = V; VelBL{3} = W; VelBL{4} = P; VelBL{5} = T;
    VelBL{6} = dU; VelBL{7} = dV; VelBL{8} = dW; VelBL{9} = dP; VelBL{10} = dT; 
    VelBL{11} = dU_dtheta; VelBL{12} = dV_dtheta; VelBL{13} = dW_dtheta; VelBL{14} = dP_dtheta; VelBL{15} = dT_dtheta;

    % save data to file
    filename = ['../Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat'];
    if exist('../Flows/BL','dir')==0; mkdir ../Flows BL; end
    save(filename, 'VelBL', 'eta', 'theta')
    fprintf(repmat('\b',1,str3)); str3 = fprintf('Flow saved in %s\n', filename); pause(1)
    fprintf(repmat('\b',1,str3));fprintf(repmat('\b',1,str2));fprintf(repmat('\b',1,str1));
end