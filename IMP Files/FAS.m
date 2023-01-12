%% FAS Alg.
% Inputs:
% Psi - Approximation for SF                 W - Approximation for Vort.
% V - Approximation for azi. vel.            T - Approximation for Temp.        
% Sp - RHS of SF eqn                         Sw - RHS of Vort. eqn              
% Sv - RHS of V eqn                          St - RHS of T eqn
% BC - Bound. Cond. function                 h - size of grid
% Re - Reynolds Number                       Pr - Prandtl number
% lambda - sensitivity constant              a - under-relaxation parameter 
% tol - tolerance
% Outputs:
% Psi - Solution for SF at grid size h       W - Solution for Vort. at grid size h
% V - Solution for azi. vel. at grid size h  T - Solution for Temp. at grid size h
function [Psi,W,V,T] = FAS(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,a,tol)
    %% First Residual
    [Pres,Wres,Vres,Tres] = residual(Psi,W,V,T,Sp,Sw,Sv,St,h,Re,Pr,lambda); 
    r = max([abs(Pres);abs(Wres);abs(Vres);abs(Tres)],[],'all');
    str = fprintf('residual = %f, %.2f%% Complete',[NaN,0]);
    
    %% Iterate until desired accuarcy
    while r>tol
        fprintf(repmat('\b',1,str))
        str = fprintf('residual = %f, %.2f%% Complete',[r,tol/r*100]);
        
        % V-Cycle
        [Psi,W,V,T] = Vcycle(Psi,W,V,T,Sp,Sw,Sv,St,BC,h,Re,Pr,lambda,a);

        % Residual 
        [Pres,Wres,Vres,Tres] = residual(Psi,W,V,T,Sp,Sw,Sv,St,h,Re,Pr,lambda); 
        r = max([abs(Pres);abs(Wres);abs(Vres);abs(Tres)],[],'all');
    end
    fprintf(repmat('\b',1,str))
end