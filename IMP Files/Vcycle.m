%% Recursive V-Cycles to improve approximation of solutions
function [Psih,Wh,Vh,Th] = Vcycle(Psih,Wh,Vh,Th,Sph,Swh,Svh,Sth,BC,h,Re,Pr,lambda,a)
    %% Condition current grid
    % Relax
    [Psih,Wh,Vh,Th] = Relax(Psih,Wh,Vh,Th,Sph,Swh,Svh,Sth,BC,h,Re,Pr,lambda,9,a);

    % Compute Residual Errors
    [Presh,Wresh,Vresh,Tresh] = residual(Psih,Wh,Vh,Th,Sph,Swh,Svh,Sth,h,Re,Pr,lambda);

    % Restrict to Coarse grid
    Psi2h = restrict(Psih); W2h = restrict(Wh); 
    V2h = restrict(Vh); T2h = restrict(Th);
    Pres2h = restrict(Presh); Wres2h = restrict(Wresh); 
    Vres2h = restrict(Vresh); Tres2h = restrict(Tresh);
    
    % Update BCs
    N2h = size(Psi2h); N2x = N2h(2)-2; N2y = N2h(1)-2; hii = 1/(2*h); 
    [Pxxt2h,~] = BC(2*h,Re,Pr,lambda);
    W2h(N2y+1,2:N2x+1) = 2*(Psi2h(N2y+1,2:N2x+1)-Psi2h(N2y,2:N2x+1))*hii^2 - Pxxt2h(:)';
    W2h(2:N2y+1,2) = -2*Psi2h(2:N2y+1,3)*hii^2;
    W2h(3:N2y,N2x+1) = (2*Psi2h(3:N2y,N2x+1)-Psi2h(2:N2y-1,N2x)-Psi2h(4:N2y+1,N2x))*hii^2;
    % Extrapolated points
    Psi2h(2:N2y+1,1) = Psi2h(2:N2y+1,3);
    Psi2h(2:N2y+1,N2x+2) = Psi2h(2:N2y+1,N2x);
    Psi2h(1,2:N2x+1) = 3*Psi2h(2,2:N2x+1) - 3*Psi2h(3,2:N2x+1) + Psi2h(4,2:N2x+1);
    Psi2h(N2y+2,2:N2x+1) = Psi2h(N2y,2:N2x+1);
    W2h(2:N2y+1,1) = 3*W2h(2:N2y+1,2) - 3*W2h(2:N2y+1,3) + W2h(2:N2y+1,4);
    W2h(2:N2y+1,N2x+2) = 3*W2h(2:N2y+1,N2x+1) - 3*W2h(2:N2y+1,N2x) + W2h(2:N2y+1,N2x-1);
    W2h(1,2:N2x+1) = -3*W2h(3,2:N2x+1) + W2h(4,2:N2x+1);
    W2h(N2y+2,2:N2x+1) = 3*W2h(N2y+1,2:N2x+1) - 3*W2h(N2y,2:N2x+1) + W2h(N2y-1,2:N2x+1);
    V2h(2:N2y+1,1) = 3*V2h(2:N2y+1,2) - 3*V2h(2:N2y+1,3) + V2h(2:N2y+1,4);
    V2h(2:N2y+1,N2x+2) = V2h(2:N2y+1,N2x);
    V2h(1,2:N2x+1) = V2h(3,2:N2x+1);
    V2h(N2y+2,2:N2x+1) = 3*V2h(N2y+1,2:N2x+1) - 3*V2h(N2y,2:N2x+1) + V2h(N2y-1,2:N2x+1);
    T2h(2:N2y+1,1) = 3 - 3*T2h(2:N2y+1,3) + T2h(2:N2y+1,4);
    T2h(2:N2y+1,N2x+2) = T2h(2:N2y+1,N2x);
    T2h(1,2:N2x+1) = T2h(3,2:N2x+1);
    T2h(N2y+2,2:N2x+1) = 3*T2h(N2y+1,2:N2x+1) - 3*T2h(N2y,2:N2x+1) + T2h(N2y-1,2:N2x+1);
    
    % form RHS
    [Sp2h,Sw2h,Sv2h,St2h] = RHS(Psi2h,W2h,V2h,T2h,Pres2h,Wres2h,Vres2h,Tres2h,2*h,Re,Pr,lambda);

    % stop recursion at smallest grid size, otherwise continue recursion
    if h>0.4
        [Ap2h,Aw2h,Av2h,At2h] = Relax(Psi2h,W2h,V2h,T2h,Sp2h,Sw2h,Sv2h,St2h,BC,2*h,Re,Pr,lambda,50,0.9);
    else 
        [Ap2h,Aw2h,Av2h,At2h] = Vcycle(Psi2h,W2h,V2h,T2h,Sp2h,Sw2h,Sv2h,St2h,BC,2*h,Re,Pr,lambda,a);
    end

    %% Interpolate upto fine grid
    % determine error
    Ep2h = Ap2h - Psi2h; Ew2h = Aw2h - W2h;
    Ev2h = Av2h - V2h; Et2h = At2h - T2h;
    % interpolate error
    Eph = interp(Ep2h); Ewh = interp(Ew2h);
    Evh = interp(Ev2h); Eth = interp(Et2h);
    % correct
    Psih = Psih + Eph; Wh = Wh + Ewh;
    Vh = Vh + Evh; Th = Th + Eth;
    
    % Update BCs
    N = size(Psih); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h; 
    [Pxxth,~] = BC(h,Re,Pr,lambda);
    Wh(Ny+1,2:Nx+1) = 2*(Psih(Ny+1,2:Nx+1)-Psih(Ny,2:Nx+1))*hi^2 - Pxxth(:)';
    Wh(2:Ny+1,2) = -2*Psih(2:Ny+1,3)*hi^2;
    Wh(3:Ny,Nx+1) = (2*Psih(3:Ny,Nx+1)-Psih(2:Ny-1,Nx)-Psih(4:Ny+1,Nx))*hi^2;
    % Extrapolated points
    Psih(2:Ny+1,1) = Psih(2:Ny+1,3);
    Psih(2:Ny+1,Nx+2) = Psih(2:Ny+1,Nx);
    Psih(1,2:Nx+1) = 3*Psih(2,2:Nx+1) - 3*Psih(3,2:Nx+1) + Psih(4,2:Nx+1);
    Psih(Ny+2,2:Nx+1) = Psih(Ny,2:Nx+1);
    Wh(2:Ny+1,1) = 3*Wh(2:Ny+1,2) - 3*Wh(2:Ny+1,3) + Wh(2:Ny+1,4);
    Wh(2:Ny+1,Nx+2) = 3*Wh(2:Ny+1,Nx+1) - 3*Wh(2:Ny+1,Nx) + Wh(2:Ny+1,Nx-1);
    Wh(1,2:Nx+1) = -3*Wh(3,2:Nx+1) + Wh(4,2:Nx+1);
    Wh(Ny+2,2:Nx+1) = 3*Wh(Ny+1,2:Nx+1) - 3*Wh(Ny,2:Nx+1) + Wh(Ny-1,2:Nx+1);
    Vh(2:Ny+1,1) = 3*Vh(2:Ny+1,2) - 3*Vh(2:Ny+1,3) + Vh(2:Ny+1,4);
    Vh(2:Ny+1,Nx+2) = Vh(2:Ny+1,Nx);
    Vh(1,2:Nx+1) = Vh(3,2:Nx+1);
    Vh(Ny+2,2:Nx+1) = 3*Vh(Ny+1,2:Nx+1) - 3*Vh(Ny,2:Nx+1) + Vh(Ny-1,2:Nx+1);
    Th(2:Ny+1,1) = 3 - 3*Th(2:Ny+1,3) + Th(2:Ny+1,4);
    Th(2:Ny+1,Nx+2) = Th(2:Ny+1,Nx);
    Th(1,2:Nx+1) = Th(3,2:Nx+1);
    Th(Ny+2,2:Nx+1) = 3*Th(Ny+1,2:Nx+1) - 3*Th(Ny,2:Nx+1) + Th(Ny-1,2:Nx+1);

    % Post-Smoothing
    [Psih,Wh,Vh,Th] = Relax(Psih,Wh,Vh,Th,Sph,Swh,Svh,Sth,BC,h,Re,Pr,lambda,6,a);
end