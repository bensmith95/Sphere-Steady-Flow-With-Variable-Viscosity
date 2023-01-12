%% Gauss-Seidel-Newton block iterative scheme to relax/solve SF-Vort. equations
% Inputs:
% Psi - Approximation for SF           W - Approximation for Vort.
% V - Approximation for azi. vel.      T - Approximation for Temp.        
% Sp - RHS of SF eqn                   Sw - RHS of Vort. eqn 
% Sv - RHS of V eqn                    St - RHS of T eqn
% BC - Bound. Cond. function           h - size of grid
% Re - Reynolds Number                 Pr - Prandtl number
% lambda - sensitivity constant        n - # of iterations
% a - under-relaxation parameter
% Outputs:
% Psi - Updated approx. for SF         W - Updated approx. for Vort.
% V - Updated approx for azi. vel.     T - Updated approx. for Temp.    
function [Psi,W,V,T] = Relax(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,n,a)
    %% Initalise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; 
    eps = 1/Re; hi = 1/h; [Pxxt,~] = BC(h,Re,Pr,lambda);
    
    %% Iterate
    % direction of iterations follows overall direction of flow
    for k = 1:n
        for j = Ny:-1:2
            for i = 3:Nx+1
                dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi; 
                dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
                dTdx = (T(j,i+1)-T(j,i-1))*0.5*hi; 
                dTdy = (T(j+1,i)-T(j-1,i))*0.5*hi;
                mu = 1/(1+lambda*T(j,i));
                dVdy = (V(j+1,i)-V(j-1,i))*0.5*hi; 
                if j>2
                    % Internal Points
                    if i<Nx+1
                        A = zeros(4); b = zeros(4,1);

                        % SF
                        Cx = dPsidy + 2*lambda*eps*mu^2*dTdx;
                        if Cx>=0
                            dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; X = 3*Cx*0.5*hi;
                        else
                            dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; X = -3*Cx*0.5*hi;
                        end
                        Cy = -(dPsidx - 2*lambda*eps*mu^2*dTdy);
                        if Cy>=0
                            dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; Y = 3*Cy*0.5*hi;
                        else
                            dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; Y = -3*Cy*0.5*hi;
                        end
                        LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;
                        DT = (T(j-1,i)+T(j+1,i)-T(j,i-1)-T(j,i+1))*hi^2;
                        DPsi = (Psi(j-1,i)+Psi(j+1,i)-Psi(j,i-1)-Psi(j,i+1))*hi^2;
                        visc1 = (2*lambda^2*mu^3*(dTdx^2-dTdy^2) + lambda*mu^2*DT)*DPsi;
                        d2Tdxdy = (T(j+1,i+1)-T(j-1,i+1)-T(j+1,i-1)+T(j-1,i-1))*0.25*hi^2;
                        d2Psidxdy = (Psi(j+1,i+1)-Psi(j-1,i+1)-Psi(j+1,i-1)+Psi(j-1,i-1))*0.25*hi^2;
                        visc2 = (2*lambda^2*mu^3*dTdx*dTdy - lambda*mu^2*d2Tdxdy)*d2Psidxdy;
                        % Jacobian terms
                        A(1,2) = X+Y+4*eps*mu*hi^2;
                        A(1,3) = 2*eps*dVdy;
                        A(1,4) = -4*eps*lambda^2*mu^3*(dTdx*dWdx+dTdy*dWdy) + eps*lambda*mu^2*LapW + ...
                                    2*lambda^2*eps*mu^3*(3*lambda*mu*(dTdx^2-dTdy^2)+DT)*DPsi - ...
                                        8*lambda^2*eps*mu^3*(3*lambda*mu*dTdx*dTdy-d2Tdxdy)*d2Psidxdy;
                        b(1) = -(Cx*dWdx + Cy*dWdy + 2*eps*V(j,i)*dVdy - eps*mu*LapW - eps*visc1 + 4*eps*visc2 - Sp(j,i));

                        % Vorticity
                        LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;
                        A(2,1) = -4*hi^2; A(2,2) = 1;
                        b(2) = -(W(j,i) + LapPsi - Sw(j,i));
                        
                        % V momentum
                        CVx = dPsidy+eps*lambda*mu^2*dTdx;
                        if CVx>=0
                            dV_dx = (3*V(j,i)-4*V(j,i-1)+V(j,i-2))*0.5*hi; CW = 3*CVx*0.5*hi;
                        else
                            dV_dx = -(3*V(j,i)-4*V(j,i+1)+V(j,i+2))*0.5*hi; CW = -3*CVx*0.5*hi;
                        end
                        CVy = -dPsidx+eps*lambda*mu^2*dTdy;
                        if CVy>=0
                            dV_dy = (3*V(j,i)-4*V(j-1,i)+V(j-2,i))*0.5*hi; CU = 3*CVy*0.5*hi;
                        else
                            dV_dy = -(3*V(j,i)-4*V(j+1,i)+V(j+2,i))*0.5*hi; CU = -3*CVy*0.5*hi;
                        end
                        LapV = (V(j-1,i)+V(j,i-1)-4*V(j,i)+V(j,i+1)+V(j+1,i))*hi^2;
                        A(3,3) = CW + CU + eps*dPsidy + 4*eps*mu*hi^2;
                        A(3,4) = -2*eps*lambda^2*mu^3*(dTdx*dV_dx+dTdy*dV_dy) + eps*lambda*mu^2*LapV;
                        b(3) = -(CVx*dV_dx + CVy*dV_dy + eps*V(j,i)*dPsidy - eps*mu*LapV - Sv(j,i));
                        
                        % Temperature
                        if dPsidy>=0
                            dT_dx = (3*T(j,i)-4*T(j,i-1)+T(j,i-2))*0.5*hi; Q = 3*dPsidy*0.5*hi;
                        else
                            dT_dx = -(3*T(j,i)-4*T(j,i+1)+T(j,i+2))*0.5*hi; Q = -3*dPsidy*0.5*hi;
                        end
                        if -dPsidx>=0
                            dT_dy = (3*T(j,i)-4*T(j-1,i)+T(j-2,i))*0.5*hi; Z = -3*dPsidx*0.5*hi;
                        else
                            dT_dy = -(3*T(j,i)-4*T(j+1,i)+T(j+2,i))*0.5*hi; Z = 3*dPsidx*0.5*hi;
                        end
                        LapT = (T(j-1,i)+T(j,i-1)-4*T(j,i)+T(j,i+1)+T(j+1,i))*hi^2;
                        A(4,4) = Q+Z+4*eps/Pr*hi^2;
                        b(4) = -(dPsidy*dT_dx - dPsidx*dT_dy - eps/Pr*LapT - St(j,i));

                        % solve & correct
                        s = A\b;
                        Psi(j,i) = Psi(j,i) + a*s(1); 
                        W(j,i) = W(j,i) + a*s(2);
                        V(j,i) = V(j,i) + a*s(3);
                        T(j,i) = T(j,i) + a*s(4);

                        % Update BC's 
                        if ismember(i,[3,4])
                            W(j,2) = -2*Psi(j,3)*hi^2; W(j,1) = 3*W(j,2) - 3*W(j,3) + W(j,4); 
                            Psi(j,1) = Psi(j,3); 
                            V(j,1) = 3*V(j,2) - 3*V(j,3) + V(j,4);
                            T(j,1) = 3 - 3*T(j,3) + T(j,4);
                        end
                        if ismember(i,[Nx-1,Nx])
                            if j<Ny; W(j+1,Nx+1) = (2*Psi(j+1,Nx+1)-Psi(j,Nx)-Psi(j+2,Nx))*hi^2; end
                            if j>3; W(j-1,Nx+1) = (2*Psi(j-1,Nx+1)-Psi(j-2,Nx)-Psi(j,Nx))*hi^2; end
                            W(j,Nx+2) = 3*W(j,Nx+1) - 3*W(j,Nx) + W(j,Nx-1); Psi(j,Nx+2) = Psi(j,Nx); 
                            V(j,Nx+2) = V(j,Nx); T(j,Nx+2) = T(j,Nx); 
                        end
                        if ismember(j,[3,4]) 
                            W(1,i) = -3*W(3,i) + W(4,i); Psi(1,i) = -3*Psi(3,i) + Psi(4,i); 
                            V(1,i) = V(3,i); T(1,i) = T(3,i); 
                        end
                        if ismember(j,[Ny-1,Ny])
                            W(Ny+1,i) = 2*(Psi(Ny+1,i)-Psi(Ny,i))*hi^2 - Pxxt(i-1)'; W(Ny+2,i) = 3*W(Ny+1,i) - 3*W(Ny,i) + W(Ny-1,i); 
                            Psi(Ny+2,i) = Psi(Ny,i); 
                            V(Ny+2,i) = 3*V(Ny+1,i) - 3*V(Ny,i) + V(Ny-1,i);
                            T(Ny+2,i) = 3*T(Ny+1,i) - 3*T(Ny,i) + T(Ny-1,i);
                        end
                    % Neumann BC
                    else
                        Psi(j,Nx+1) = -(Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx)-Sp(j,Nx+1))/(-2);
                        W(j,Nx+1) = (2*Psi(j,Nx+1)-Psi(j-1,Nx)-Psi(j+1,Nx))*hi^2;
                        
                        A = zeros(2); b = zeros(2,1);
                        % V mom
                        CVy = lambda*mu^2*dTdy;
                        if CVy>=0
                            dV_dy = (3*V(j,Nx+1)-4*V(j-1,Nx+1)+V(j-2,Nx+1))*0.5*hi; CU = 3*CVy*0.5*hi;
                        else
                            dV_dy = -(3*V(j,Nx+1)-4*V(j+1,Nx+1)+V(j+2,Nx+1))*0.5*hi; CU = -3*CVy*0.5*hi;
                        end
                        LapV = (V(j-1,Nx+1)+V(j,Nx)-4*V(j,Nx+1)+V(j,Nx+2)+V(j+1,Nx+1))*hi^2;
                        A(1,1) = CU + eps*dPsidy + 4*mu*hi^2;
                        A(1,2) = -2*lambda^2*mu^3*dTdy*dV_dy + lambda*mu^2*LapV;
                        b(1) = -(CVy*dV_dy + eps*V(j,Nx+1)*dPsidy - mu*LapV - Sv(j,Nx+1)); 
                        
                        % Heat
                        LapT = T(j-1,Nx+1)+T(j,Nx)-4*T(j,Nx+1)+T(j,Nx+2)+T(j+1,Nx+1);
                        A(2,2) = -4; b(2) = -(LapT - St(j,Nx+1));
                        s = A\b;
                        V(j,Nx+1) = V(j,Nx+1) + s(1);
                        T(j,Nx+1) = T(j,Nx+1) + s(2);
                                                
                        % BCs
                        W(j,Nx+2) = 3*W(j,Nx+1) - 3*W(j,Nx) + W(j,Nx-1);
                        if ismember(j,[3,4]); Psi(1,Nx+1) = -3*Psi(3,Nx+1)+Psi(4,Nx+1); V(1,Nx+1) = V(3,Nx+1); T(1,Nx+1) = T(3,Nx+1); end 
                        if ismember(j,[Ny-1,Ny]); Psi(Ny+2,Nx+1) = Psi(Ny,Nx+1); V(Ny+2,Nx+1) = 3*V(Ny+1,Nx+1) - 3*V(Ny,Nx+1) + V(Ny-1,Nx+1); 
                           T(Ny+2,Nx+1) = 3*T(Ny+1,Nx+1) - 3*T(Ny,Nx+1) + T(Ny-1,Nx+1); end 
                    end
                else
                    if i<Nx+1
                        A = zeros(2); b = zeros(2,1);
                        % V mom
                        CVx = dPsidy+eps*lambda*mu^2*dTdx;
                        if CVx>=0
                            dV_dx = (3*V(2,i)-4*V(2,i-1)+V(2,i-2))*0.5*hi; CW = 3*CVx*0.5*hi;
                        else
                            dV_dx = -(3*V(2,i)-4*V(2,i+1)+V(2,i+2))*0.5*hi; CW = -3*CVx*0.5*hi;
                        end
                        LapV = (V(1,i)+V(2,i-1)-4*V(2,i)+V(2,i+1)+V(3,i))*hi^2;
                        A(1,1) = CW + eps*dPsidy + 4*eps*mu*hi^2;
                        A(1,2) = -2*eps*lambda^2*mu^3*dTdx*dV_dx + eps*lambda*mu^2*LapV;
                        b(1) = -(CVx*dV_dx + eps*V(2,i)*dPsidy - eps*mu*LapV - Sv(2,i));
                        
                        % Heat 
                        if dPsidy>=0
                            dT_dx = (3*T(2,i)-4*T(2,i-1)+T(2,i-2))*0.5*hi; Q = 3*dPsidy*0.5*hi;
                        else
                            dT_dx = -(3*T(2,i)-4*T(2,i+1)+T(2,i+2))*0.5*hi; Q = -3*dPsidy*0.5*hi;
                        end
                        LapT = (T(1,i)+T(2,i-1)-4*T(2,i)+T(2,i+1)+T(3,i))*hi^2;
                        A(2,2) = Q+4*eps/Pr*hi^2;
                        b(2) = -(dPsidy*dT_dx - eps/Pr*LapT - St(2,i));
                        s = A\b;
                        V(2,i) = V(2,i) + s(1);
                        T(2,i) = T(2,i) + s(2);
                        
                        % BCs 
                        if ismember(i,[3,4]); V(2,1) = 3 - 3*V(2,3) + V(2,4); T(2,1) = 3 - 3*T(2,3) + T(2,4); end
                        if i==Nx; V(2,Nx+2)=V(2,Nx); T(2,Nx+2)=T(2,Nx); end
                    else
                        V(2,Nx+1) = (V(1,Nx+1)+V(2,Nx)+V(2,Nx+2)+V(3,Nx+1)+Sv(2,Nx+1)*h^2)/(4+dPsidy*h^2);
                        T(2,Nx+1) = (T(1,Nx+1)+T(2,Nx)+T(2,Nx+2)+T(3,Nx+1)-St(2,Nx+1))/4;
                    end
                end
            end
        end
        % Update BC's
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % Extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = -3*Psi(3,2:Nx+1) + Psi(4,2:Nx+1);
        Psi(Ny+2,2:Nx+1) = Psi(Ny,2:Nx+1);
        W(2:Ny+1,1) = 3*W(2:Ny+1,2) - 3*W(2:Ny+1,3) + W(2:Ny+1,4);
        W(2:Ny+1,Nx+2) = 3*W(2:Ny+1,Nx+1) - 3*W(2:Ny+1,Nx) + W(2:Ny+1,Nx-1);
        W(1,2:Nx+1) = -3*W(3,2:Nx+1) + W(4,2:Nx+1);
        W(Ny+2,2:Nx+1) = 3*W(Ny+1,2:Nx+1) - 3*W(Ny,2:Nx+1) + W(Ny-1,2:Nx+1);
        V(2:Ny+1,1) = 3*V(2:Ny+1,2) - 3*V(2:Ny+1,3) + V(2:Ny+1,4);
        V(2:Ny+1,Nx+2) = V(2:Ny+1,Nx);
        V(1,2:Nx+1) = V(3,2:Nx+1);
        V(Ny+2,2:Nx+1) = 3*V(Ny+1,2:Nx+1) - 3*V(Ny,2:Nx+1) + V(Ny-1,2:Nx+1);
        T(2:Ny+1,1) = 3 - 3*T(2:Ny+1,3) + T(2:Ny+1,4);
        T(2:Ny+1,Nx+2) = T(2:Ny+1,Nx);
        T(1,2:Nx+1) = T(3,2:Nx+1);
        T(Ny+2,2:Nx+1) = 3*T(Ny+1,2:Nx+1) - 3*T(Ny,2:Nx+1) + T(Ny-1,2:Nx+1);
    end
end