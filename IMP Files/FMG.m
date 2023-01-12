%% FMG-FAS Alg. to solve the non-linear system
% Inputs: 
% Psi - Initial Stream Function (SF) guess     W - Initial Vorticity guess 
% V - Initial azimuthal velocity guess         T - Initial Temperature guess
% BC - function that contains the Boundary Condition information
% H - size of fine grid                        Re - Reynolds Number
% Pr - Prandtl number                          lambda - non-dim sensitivity constant
% Outputs:
% Psi - solution field for SF                  W - solution field for Vorticity
% V - solution field for azimuthal velocity    T - solution field for Temperature
% res - final residual of solution
function [Psi,W,V,T,res] = FMG(Psi,W,V,T,BC,H,Re,Pr,lambda)
    %% Pre-Condition and Restrict to Coarse Grid
    tol = 1e-5; % desired accuracy/tolerance 
    a = 0.9; % under-relaxation parameter
    h = H; str1 = fprintf('Restricting to Coarse Grid...\n');
    while h<1
        % Relax
        Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp; St = Sp;
        [Psi,W,V,T] = Relax(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,9,a);
        
        % Restrict
        Psi = restrict(Psi); W = restrict(W); 
        V = restrict(V); T = restrict(T);
        h = 2*h; tol = 4*tol;
        
        % Update BCs
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h; 
        [Pxxt,Pt,VI,TI] = BC(h,Re,Pr,lambda);
        % Psi
        Psi(Ny+1,2:Nx+1) = Pt(:)';
        % W
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % V
        V(Ny+1,2:Nx+1) = VI;
        % T
        T(Ny+1,2:Nx+1) = TI;
        % Extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = 3*Psi(2,2:Nx+1) - 3*Psi(3,2:Nx+1) + Psi(4,2:Nx+1);
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

    %% Solve on Coarse Grid
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
    Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp; St = Sp; pr = 1e12;
    str2 = fprintf('residual = %f, %.2f%% Complete',[NaN,0]);
    while pr>tol
        fprintf(repmat('\b',1,str2))
        str2 = fprintf('residual = %f, %.2f%% Complete',[pr,tol/pr*100]);
        [Psi,W,V,T] = Relax(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,50,0.9);
        [Pres,Wres,Vres,Tres] = residual(Psi,W,V,T,Sp,Sw,Sv,St,h,Re,Pr,lambda); 
        pr = max([abs(Pres);abs(Wres);abs(Vres);abs(Tres)],[],'all'); 
    end
    fprintf(repmat('\b',1,str2))

    %% FMG Alg.
    while h>H
        % Interpolate 
        Psi = interp(Psi); W = interp(W); V = interp(V); T = interp(T); 
        h = h/2; tol = tol/4;
        % Update BCs
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
        fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
        hi = 1/h; [Pxxt,Pt,VI,TI] = BC(h,Re,Pr,lambda);
        % Psi
        Psi(Ny+1,2:Nx+1) = Pt(:)';
        % W
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % V
        V(Ny+1,2:Nx+1) = VI; V(2:Ny+1,2) = 1 - 0.5/Re^2*(0:h:ceil(sqrt(Re))).^2; 
        % T
        T(Ny+1,2:Nx+1) = TI;
        % Extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = 3*Psi(2,2:Nx+1) - 3*Psi(3,2:Nx+1) + Psi(4,2:Nx+1);
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
        % Post-Relax
        Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp; St = Sp;
        [Psi,W,V,T] = Relax(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,6,a);
        % FAS Alg.
        [Psi,W,V,T] = FAS(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,a,tol);
    end
    
    %% Determine Residual
    fprintf(repmat('\b',1,str1));
    [Pres,Wres,Vres,Tres] = residual(Psi,W,V,T,Sp,Sw,Sv,St,h,Re,Pr,lambda); 
    res = max([abs(Pres);abs(Wres);abs(Vres);abs(Tres)],[],'all');
end