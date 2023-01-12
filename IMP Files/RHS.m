%% Determines new Right Hand Side of equations for V-Cycles
% Inputs:
% Psi - Approximation for SF           W - Approximation for Vort.
% V - Approximation for azi. vel.      T - Approximation for Temp.        
% Pres - residual of SF eqn            Wres - residual of Vort. eqn   
% Vres - residual of V eqn             Tres - residual of T eqn
% h - size of grid                     Re - Reynolds Number 
% Pr - Prandtl number                  lambda - sensitivity constant
% Outputs:
% Sp - new RHS of SF equation          Sw - new RHS of Vort. eqn
% Sv - new RHS of V eqn                St - new RHS of Temp. eqn              
function [Sph,Swh,Svh,Sth] = RHS(Psi,W,V,T,Pres,Wres,Vres,Tres,h,Re,Pr,lambda)
    %% Initialise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
    Sph = zeros(N); Swh = zeros(N); Svh = zeros(N); Sth = zeros(N);
    
    %% Determine A(v2h)
    eps = 1/Re; hi = 1/h;
    for j = Ny:-1:2
        for i = 3:Nx+1
            dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi; 
            dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
            dTdx = (T(j,i+1)-T(j,i-1))*0.5*hi; mu = 1/(1+lambda*T(j,i));
            dTdy = (T(j+1,i)-T(j-1,i))*0.5*hi;
            dVdy = (V(j+1,i)-V(j-1,i))*0.5*hi; 
            if j>2
                if i<Nx+1
                    % SF terms
                    Cx = dPsidy + 2*lambda*eps*mu^2*dTdx;
                    if Cx>=0
                        dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; 
                    else
                        dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; 
                    end
                    Cy = -(dPsidx - 2*lambda*eps*mu^2*dTdy);
                    if Cy>=0
                        dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; 
                    else
                        dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; 
                    end
                    LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;
                    DT = (T(j-1,i)+T(j+1,i)-T(j,i-1)-T(j,i+1))*hi^2;
                    DPsi = (Psi(j-1,i)+Psi(j+1,i)-Psi(j,i-1)-Psi(j,i+1))*hi^2;
                    visc1 = (2*lambda^2*mu^3*(dTdx^2-dTdy^2) + lambda*mu^2*DT)*DPsi;
                    d2Tdxdy = (T(j+1,i+1)-T(j-1,i+1)-T(j+1,i-1)+T(j-1,i-1))*0.25*hi^2;
                    d2Psidxdy = (Psi(j+1,i+1)-Psi(j-1,i+1)-Psi(j+1,i-1)+Psi(j-1,i-1))*0.25*hi^2;
                    visc2 = (2*lambda^2*mu^3*dTdx*dTdy - lambda*mu^2*d2Tdxdy)*d2Psidxdy;

                    % SF terms
                    LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;

                    % V mom. terms
                    CVx = dPsidy+eps*lambda*mu^2*dTdx;
                    if CVx>=0
                        dV_dx = (3*V(j,i)-4*V(j,i-1)+V(j,i-2))*0.5*hi; 
                    else
                        dV_dx = -(3*V(j,i)-4*V(j,i+1)+V(j,i+2))*0.5*hi; 
                    end
                    CVy = -dPsidx+eps*lambda*mu^2*dTdy;
                    if CVy>=0
                        dV_dy = (3*V(j,i)-4*V(j-1,i)+V(j-2,i))*0.5*hi; 
                    else
                        dV_dy = -(3*V(j,i)-4*V(j+1,i)+V(j+2,i))*0.5*hi; 
                    end
                    LapV = (V(j-1,i)+V(j,i-1)-4*V(j,i)+V(j,i+1)+V(j+1,i))*hi^2;

                    % Heat terms
                    if dPsidy>=0
                        dT_dx = (3*T(j,i)-4*T(j,i-1)+T(j,i-2))*0.5*hi; 
                    else
                        dT_dx = -(3*T(j,i)-4*T(j,i+1)+T(j,i+2))*0.5*hi; 
                    end
                    if -dPsidx>=0
                        dT_dy = (3*T(j,i)-4*T(j-1,i)+T(j-2,i))*0.5*hi; 
                    else
                        dT_dy = -(3*T(j,i)-4*T(j+1,i)+T(j+2,i))*0.5*hi; 
                    end
                    LapT = (T(j-1,i)+T(j,i-1)-4*T(j,i)+T(j,i+1)+T(j+1,i))*hi^2;

                    % Determine new RHS
                    Sph(j,i) = (Cx*dWdx + Cy*dWdy + 2*eps*V(j,i)*dVdy - eps*mu*LapW - eps*visc1 + 4*eps*visc2) + Pres(j,i);
                    Swh(j,i) = (W(j,i) + LapPsi) + Wres(j,i);
                    Svh(j,i) = (CVx*dV_dx + CVy*dV_dy + eps*V(j,i)*dPsidy - eps*mu*LapV) + Vres(j,i);
                    Sth(j,i) = (dPsidy*dT_dx - dPsidx*dT_dy - eps/Pr*LapT) + Tres(j,i);
                else
                    % SF
                    Sph(j,Nx+1) = (Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)-2*Psi(j,Nx+1)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx)) + Pres(j,Nx+1); 
                    
                    % V mom.
                    CVy = lambda*mu^2*dTdy;
                    if CVy>=0
                        dV_dy = (3*V(j,Nx+1)-4*V(j-1,Nx+1)+V(j-2,Nx+1))*0.5*hi; 
                    else
                        dV_dy = -(3*V(j,Nx+1)-4*V(j+1,Nx+1)+V(j+2,Nx+1))*0.5*hi; 
                    end
                    LapV = (V(j-1,Nx+1)+V(j,Nx)-4*V(j,Nx+1)+V(j,Nx+2)+V(j+1,Nx+1))*hi^2;
                    Svh(j,Nx+1) = (CVy*dV_dy + eps*V(j,Nx+1)*dPsidy - mu*LapV) + Vres(j,Nx+1);
                    
                    % Heat
                    Sth(j,Nx+1) = (T(j-1,Nx+1)+T(j,Nx)-4*T(j,Nx+1)+T(j,Nx+2)+T(j+1,Nx+1)) + Tres(j,Nx+1);                    
                end
            else
                if i<Nx+1
                    % V mom
                    CVx = dPsidy+eps*lambda*mu^2*dTdx;
                    if CVx>=0
                        dV_dx = (3*V(2,i)-4*V(2,i-1)+V(2,i-2))*0.5*hi; 
                    else
                        dV_dx = -(3*V(2,i)-4*V(2,i+1)+V(2,i+2))*0.5*hi; 
                    end
                    LapV = (V(1,i)+V(2,i-1)-4*V(2,i)+V(2,i+1)+V(3,i))*hi^2;
                    Svh(2,i) = (CVx*dV_dx + eps*V(j,i)*dPsidy - eps*mu*LapV) + Vres(2,i);
                    
                    % Heat
                    if dPsidy>=0
                        dT_dx = (3*T(2,i)-4*T(2,i-1)+T(2,i-2))*0.5*hi; 
                    else
                        dT_dx = -(3*T(2,i)-4*T(2,i+1)+T(2,i+2))*0.5*hi; 
                    end
                    LapT = (T(1,i)+T(2,i-1)-4*T(2,i)+T(2,i+1)+T(3,i))*hi^2;
                    Sth(2,i) = (dPsidy*dT_dx - eps/Pr*LapT) + Tres(2,i);
                    
                else
                    Svh(2,Nx+1) = (V(2,Nx+1)*dPsidy - (V(1,Nx+1)+V(2,Nx)-4*V(2,Nx+1)+V(2,Nx+2)+V(3,Nx+1))*hi^2) + Vres(2,Nx+1);
                    Sth(2,Nx+1) = (T(1,Nx+1)+T(2,Nx)-4*T(2,Nx+1)+T(2,Nx+2)+T(3,Nx+1)) + Tres(2,Nx+1);
                end
            end
        end
    end
end