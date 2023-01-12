%% Determine Pressure distribution
% via the Pressure Laplace equation
function P = Pressure(U,V,W,T,PI,h,Re,lambda)
    %% Initialise
    [Nbeta,Neta] = size(U); hi = 1/h; eps = 1/Re; mu = 1./(1+lambda*T);
    
    % determine derivatives
    dUdeta = zeros(Nbeta,Neta); dUdbeta = dUdeta; 
    dWdeta = dUdeta; dWdbeta = dUdeta; dVdeta = dUdeta;
    dTdeta = dUdeta; dTdbeta = dUdeta;
    d2Udeta2 = dUdeta; d2Udbeta2 = dUdeta;
    d2Wdeta2 = dUdeta; d2Wdbeta2 = dUdeta;
    d2Tdeta2 = dUdeta; d2Tdbeta2 = dUdeta; D2T = dUdeta;
    
    % U
    dUdeta(:,1) = -(3*U(:,1)-4*U(:,2)+U(:,3))*0.5*hi;
    dUdeta(:,2:Neta-1) = (U(:,3:Neta)-U(:,1:Neta-2))*0.5*hi; 
    dUdeta(:,Neta) = (3*U(:,Neta)-4*U(:,Neta-1)+U(:,Neta-2))*0.5*hi;
    
    dUdbeta(1,:) = -(3*U(1,:)-4*U(2,:)+U(3,:))*0.5*hi;
    dUdbeta(2:Nbeta-1,:) = (U(3:Nbeta,:)-U(1:Nbeta-2,:))*0.5*hi;
    dUdbeta(Nbeta,:) = (3*U(Nbeta,:)-4*U(Nbeta-1,:)+U(Nbeta-2,:))*0.5*hi;
    
    d2Udeta2(:,1) = (2*U(:,1)-5*U(:,2)+4*U(:,3)-U(:,4))*hi^2;
    d2Udeta2(:,2:Neta-1) = (U(:,1:Neta-2)-2*U(:,2:Neta-1)+U(:,3:Neta))*hi^2;
    d2Udeta2(:,Neta) = (2*U(:,Neta)-5*U(:,Neta-1)+4*U(:,Neta-2)-U(:,Neta-3))*hi^2;
    
    d2Udbeta2(1,:) = (2*U(1,:)-5*U(2,:)+4*U(3,:)-U(4,:))*hi^2;
    d2Udbeta2(2:Nbeta-1,:) = (U(1:Nbeta-2,:)-2*U(2:Nbeta-1,:)+U(3:Nbeta,:))*hi^2;
    d2Udbeta2(Nbeta,:) = (2*U(Nbeta,:)-5*U(Nbeta-1,:)+4*U(Nbeta-2,:)-U(Nbeta-3,:))*hi^2;
    
    % W
    dWdeta(:,1) = -(3*W(:,1)-4*W(:,2)+W(:,3))*0.5*hi;
    dWdeta(:,2:Neta-1) = (W(:,3:Neta)-W(:,1:Neta-2))*0.5*hi; 
    dWdeta(:,Neta) = (3*W(:,Neta)-4*W(:,Neta-1)+W(:,Neta-2))*0.5*hi;
    
    dWdbeta(1,:) = -(3*W(1,:)-4*W(2,:)+W(3,:))*0.5*hi;
    dWdbeta(2:Nbeta-1,:) = (W(3:Nbeta,:)-W(1:Nbeta-2,:))*0.5*hi;
    dWdbeta(Nbeta,:) = (3*W(Nbeta,:)-4*W(Nbeta-1,:)+W(Nbeta-2,:))*0.5*hi;
    
    d2Wdeta2(:,1) = (2*W(:,1)-5*W(:,2)+4*W(:,3)-W(:,4))*hi^2;
    d2Wdeta2(:,2:Neta-1) = (W(:,1:Neta-2)-2*W(:,2:Neta-1)+W(:,3:Neta))*hi^2;
    d2Wdeta2(:,Neta) = (2*W(:,Neta)-5*W(:,Neta-1)+4*W(:,Neta-2)-W(:,Neta-3))*hi^2;
    
    d2Wdbeta2(1,:) = (2*W(1,:)-5*W(2,:)+4*W(3,:)-W(4,:))*hi^2;
    d2Wdbeta2(2:Nbeta-1,:) = (W(1:Nbeta-2,:)-2*W(2:Nbeta-1,:)+W(3:Nbeta,:))*hi^2;
    d2Wdbeta2(Nbeta,:) = (2*W(Nbeta,:)-5*W(Nbeta-1,:)+4*W(Nbeta-2,:)-W(Nbeta-3,:))*hi^2;
    
    % V
    dVdeta(:,1) = -(3*V(:,1)-4*V(:,2)+V(:,3))*0.5*hi;
    dVdeta(:,2:Neta-1) = (V(:,3:Neta)-V(:,1:Neta-2))*0.5*hi; 
    dVdeta(:,Neta) = (3*V(:,Neta)-4*V(:,Neta-1)+V(:,Neta-2))*0.5*hi;
    
    % T
    dTdeta(:,1) = -(3*T(:,1)-4*T(:,2)+T(:,3))*0.5*hi;
    dTdeta(:,2:Neta-1) = (T(:,3:Neta)-T(:,1:Neta-2))*0.5*hi; 
    dTdeta(:,Neta) = (3*T(:,Neta)-4*T(:,Neta-1)+T(:,Neta-2))*0.5*hi;
    
    dTdbeta(1,:) = -(3*T(1,:)-4*T(2,:)+T(3,:))*0.5*hi;
    dTdbeta(2:Nbeta-1,:) = (T(3:Nbeta,:)-T(1:Nbeta-2,:))*0.5*hi;
    dTdbeta(Nbeta,:) = (3*T(Nbeta,:)-4*T(Nbeta-1,:)+T(Nbeta-2,:))*0.5*hi;
    
    d2Tdeta2(:,1) = (2*T(:,1)-5*T(:,2)+4*T(:,3)-T(:,4))*hi^2;
    d2Tdeta2(:,2:Neta-1) = (T(:,1:Neta-2)-2*T(:,2:Neta-1)+T(:,3:Neta))*hi^2;
    d2Tdeta2(:,Neta) = (2*T(:,Neta)-5*T(:,Neta-1)+4*T(:,Neta-2)-T(:,Neta-3))*hi^2;
    
    d2Tdbeta2(1,:) = (2*T(1,:)-5*T(2,:)+4*T(3,:)-T(4,:))*hi^2;
    d2Tdbeta2(2:Nbeta-1,:) = (T(1:Nbeta-2,:)-2*T(2:Nbeta-1,:)+T(3:Nbeta,:))*hi^2;
    d2Tdbeta2(Nbeta,:) = (2*T(Nbeta,:)-5*T(Nbeta-1,:)+4*T(Nbeta-2,:)-T(Nbeta-3,:))*hi^2;
    
    D2T(2:Nbeta-1,2:Neta-1) = (T(1:Nbeta-2,1:Neta-2)-T(1:Nbeta-2,3:Neta) - ...
                               T(3:Nbeta,1:Neta-2)+T(3:Nbeta,3:Neta))*0.25*hi^2;
    D2T(:,1) = -(3*dTdbeta(:,1)-4*dTdbeta(:,2)+dTdbeta(:,3))*0.5*hi;
    D2T(:,Neta) = (3*dTdbeta(:,Neta)-4*dTdbeta(:,Neta-1)+dTdbeta(:,Neta-2))*0.5*hi;

    %% Determine f
    LapU = d2Udeta2 + d2Udbeta2; LapW = d2Wdeta2 + d2Wdbeta2;
    dmudeta = -lambda*mu.^2.*dTdeta; dmudbeta = -lambda*mu.^2.*dTdbeta;
    d2mudeta2 = 2*lambda^2*mu.^3.*dTdeta.^2 - lambda*mu.^2.*d2Tdeta2;
    d2mudbeta2 = 2*lambda^2*mu.^3.*dTdbeta.^2 - lambda*mu.^2.*d2Tdbeta2;
    D2mu = 2*lambda^2*mu.^3.*dTdeta.*dTdbeta - lambda*mu.^2.*D2T;
    visc = dmudeta.*LapW + dmudbeta.*LapU + ...
           d2mudeta2.*dWdeta + d2mudbeta2.*dUdbeta + D2mu.*(dUdeta+dWdbeta);
    F = 2*sqrt(Re)*(dWdeta.*dUdbeta - dUdeta.*dWdbeta + eps*V.*dVdeta + eps*visc)*h^2; 
    F(:,1) = F(1:end,1) + 2*h/sqrt(Re); F(Nbeta,:) = PI;

    %% Laplace operator 
    e = ones(Neta,1); 
    D = spdiags([e,-4*e,e],-1:1,Neta,Neta); D(1,2) = 2; D(Neta,Neta-1) = 2;
    I = speye(Neta);
    Z = sparse(Neta,Neta); 
    A = [D, 2*I, repmat(Z,[1,Nbeta-2])];
    for j = 2:Nbeta-1
        A = [A; repmat(Z,[1,j-2]), I, D, I, repmat(Z,[1,Nbeta-3-j+2])];
    end
    A = [A; repmat(Z,[1,Nbeta-1]), I];

    % reorganise B into a vector
    f = zeros(Neta*Nbeta,1);
    for j=1:Nbeta
       f(1+(j-1)*Neta:j*Neta) = F(j,:); 
    end

    %% Solve linear system
    Pv = A\f;

    % reorganise solution into matrix form
    P = zeros(Nbeta,Neta);
    for j=1:Nbeta
        P(j,:) = Pv(1+(j-1)*Neta:j*Neta); 
    end
end