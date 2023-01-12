%% Solves the Impinging region equations with a temperature dependant viscosity 
% of the form mu(T)=1/(1+lambda*T)
% Inputs:
% Re - Reynolds Number (sqrt) 
% Pr - Prandtl number
% lambda - non-dim sensitivity constant

function IMPRegion(Re,Pr,lambda)
    %% Initialise
    str1 = fprintf('Initialising...\n');
    % step size
    h = 2^-4;
    % set eta
    eta = 0:h:30; Neta = length(eta); 
    % set beta
    betamax = ceil(sqrt(Re)); beta = 0:h:betamax; Nbeta = length(beta);  

    % assign variable storage
    Psi = zeros(Nbeta+2,Neta+2); Omega = zeros(Nbeta+2,Neta+2); 
    V = zeros(Nbeta+2,Neta+2); T = zeros(Nbeta+2,Neta+2);

    % Obtain and update BC's
    hi = 1/h; [Pxxt,Pt,VI,TI] = BC(h,Re,Pr,lambda);
    % Psi
    Psi(Nbeta+1,2:Neta+1) = Pt(:)'; 
    % W
    Omega(Nbeta+1,2:Neta+1) = 2*(Psi(Nbeta+1,2:Neta+1)-Psi(Nbeta,2:Neta+1))*hi^2 - Pxxt(:)';
    Omega(2:Nbeta+1,2) = -2*Psi(2:Nbeta+1,3)*hi^2;
    Omega(3:Nbeta,Neta+1) = (2*Psi(3:Nbeta,Neta+1)-Psi(2:Nbeta-1,Neta)-Psi(4:Nbeta+1,Neta))*hi^2;
    % V
    V(Nbeta+1,2:Neta+1) = VI; V(2:Nbeta+1,2) = 1 - 0.5/Re^2*beta.^2; 
    % T
    T(Nbeta+1,2:Neta+1) = TI; T(2:Nbeta+1,2) = 1;
    % Extrapolated points
    Psi(2:Nbeta+1,1) = Psi(2:Nbeta+1,3);
    Psi(2:Nbeta+1,Neta+2) = Psi(2:Nbeta+1,Neta);
    Psi(1,2:Neta+1) = 3*Psi(2,2:Neta+1) - 3*Psi(3,2:Neta+1) + Psi(4,2:Neta+1);
    Psi(Nbeta+2,2:Neta+1) = Psi(Nbeta,2:Neta+1);
    Omega(2:Nbeta+1,1) = 3*Omega(2:Nbeta+1,2) - 3*Omega(2:Nbeta+1,3) + Omega(2:Nbeta+1,4);
    Omega(2:Nbeta+1,Neta+2) = 3*Omega(2:Nbeta+1,Neta+1) - 3*Omega(2:Nbeta+1,Neta) + Omega(2:Nbeta+1,Neta-1);
    Omega(1,2:Neta+1) = -3*Omega(3,2:Neta+1) + Omega(4,2:Neta+1);
    Omega(Nbeta+2,2:Neta+1) = 3*Omega(Nbeta+1,2:Neta+1) - 3*Omega(Nbeta,2:Neta+1) + Omega(Nbeta-1,2:Neta+1);
    V(2:Nbeta+1,1) = 3*V(2:Nbeta+1,2) - 3*V(2:Nbeta+1,3) + V(2:Nbeta+1,4);
    V(2:Nbeta+1,Neta+2) = V(2:Nbeta+1,Neta);
    V(1,2:Neta+1) = V(3,2:Neta+1);
    V(Nbeta+2,2:Neta+1) = 3*V(Nbeta+1,2:Neta+1) - 3*V(Nbeta,2:Neta+1) + V(Nbeta-1,2:Neta+1);
    T(2:Nbeta+1,1) = 3 - 3*T(2:Nbeta+1,3) + T(2:Nbeta+1,4);
    T(2:Nbeta+1,Neta+2) = T(2:Nbeta+1,Neta);
    T(1,2:Neta+1) = T(3,2:Neta+1);
    T(Nbeta+2,2:Neta+1) = 3*T(Nbeta+1,2:Neta+1) - 3*T(Nbeta,2:Neta+1) + T(Nbeta-1,2:Neta+1);

    %% FMG-FAS Alg.
    fprintf(repmat('\b',1,str1)); str1 = fprintf('FMG-FAS Algorithm...\n');
    [Psi,Omega,V,T,~] = FMG(Psi,Omega,V,T,@BC,h,Re,Pr,lambda);
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solution Converged.\n');

    %% POST PROCESS
    str2 = fprintf('Post-Processing...\n');
    % load Boundary Layer flow 
    filename = ['../Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat'];
    load(filename,'VelBL','eta','theta'); UB = VelBL{1}; Eta = 0:h:30;  
    % Determine U inlet
    Tm = -betamax/Re + pi/2; [~,i] = min(abs(theta-Tm));
    UI = spline(eta,-UB(i,:),Eta); % U in
    % Determine W via finite differences of SF
    W = zeros(Nbeta,Neta); 
    W(2:Nbeta-1,:) = (Psi(4:Nbeta+1,2:Neta+1)-Psi(2:Nbeta-1,2:Neta+1))*0.5*hi;
    W(1,:) = -(3*Psi(2,2:Neta+1)-4*Psi(3,2:Neta+1)+Psi(4,2:Neta+1))*0.5*hi;
    W(Nbeta,:) = 0; W(:,1) = 0;  
    % Determine U via finite differences of SF
    U = zeros(Nbeta,Neta);
    U(:,2:Neta-1) = -(Psi(2:Nbeta+1,4:Neta+1)-Psi(2:Nbeta+1,2:Neta-1))*0.5*hi;
    U(:,Neta) = 0; U(:,1) = 0; U(Nbeta,:) = UI; U(1,:) = 0; 
    
    % Solve Poisson Equation for Pressure
    PB = VelBL{4}; PI = Re^(-1/2)*spline(eta,PB(i,:),Eta); % P in
    P = Pressure(U,V(2:end-1,2:end-1),W,T(2:end-1,2:end-1),PI,h,Re,lambda); 

    % save data to file
    VelIMP{1} = U; VelIMP{2} = V(2:end-1,2:end-1); VelIMP{3} = W; VelIMP{4} = T(2:end-1,2:end-1);
    VelIMP{5} = Psi(2:end-1,2:end-1); VelIMP{6} = Omega(2:end-1,2:end-1); VelIMP{7} = P;
    filename = ['../Flows/IMP/lambda=',num2str(lambda),'/IMP_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    if exist(['../Flows/IMP/lambda=',num2str(lambda)],'dir')==0; mkdir('../Flows/IMP',['lambda=',num2str(lambda)]); end
    eta = Eta; save(filename, 'VelIMP', 'eta', 'beta')
    fprintf(repmat('\b',1,str2)); str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end
%% BOUNDARY CONDITIONS
function [Pxx,Psit,VI,TI] = BC(h,Re,Pr,lambda)
    hi = 1/h;
    % load Boundary Layer flow 
    filename = ['../Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat'];
    load(filename,'VelBL','eta','theta'); 
    UB = VelBL{1}; VB = VelBL{2}; TB = VelBL{5}; Eta = 0:h:30;
    betamax = ceil(sqrt(Re)); %betamax = ceil(Re^(4/7)); 
    % Determine betamax in theta 
    Tm = -betamax/Re + pi/2; [~,i] = min(abs(theta-Tm));
    % Determine U inlet
    UI = spline(eta,-UB(i,:),Eta);
    % Determine V inlet 
    VI = spline(eta,VB(i,:),Eta);
    % Determine TI
    TI = spline(eta,TB(i,:),Eta); % T in
    % Vorticity BC's
    Pxx = zeros(size(Eta)); Pxx(2:end-1) = (UI(3:end)-UI(1:end-2))*0.5*hi;
    Pxx(1) = -(3*UI(1)-4*UI(2)+UI(3))*0.5*hi; Pxx(end) = (3*UI(end)-4*UI(end-1)+UI(end-2))*0.5*hi; Pxx=-Pxx;
    % SF BC's
    Psit = -cumtrapz(Eta,UI); % Psi Top
end