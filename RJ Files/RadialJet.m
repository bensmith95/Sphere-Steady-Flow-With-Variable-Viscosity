%% Solves the boundary layer equations with a temperature dependant viscosity 
% of the form mu(T)=1/(1+lambda*T)
% Inputs:
% Re - square root of Reynolds number
% Pr - Prandtl number
% lambda - non-dim sensitivity constant

function RadialJet(Re,Pr,lambda)
    %% Initialise
    % Load Inlet
    filename = ['../Flows/IMP/lambda=',num2str(lambda),'/IMP_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    load(filename,'VelIMP', 'beta', 'eta'); 
    UIMP = VelIMP{1}; VIMP = VelIMP{2}; WIMP = VelIMP{3}; TIMP = VelIMP{4}; Eta=eta; Beta=beta;
    ii = max(-Re*UIMP,[],1); ii = find(ii>3); ii = ii(end); % Determine etamax in eta
    jj = length(find(WIMP(:,ii)>=0)); % Determine betamax in beta
    
    % Space marching parameters
    beta = -flip(beta(1:jj)); Nbeta = length(beta);
    filename = ['../Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat']; load(filename,'eta');
    [~,i] = min(abs(eta-Eta(end))); % Determine Etamax in eta
    r = 1+[Eta(ii-1:end) eta(i+1:end)]/Re; dr = r(end)-r(end-1);
    r = [r logspace(log10(r(end)+dr),log10(20),200)]; Nr = length(r);  

    % Initialise flow field
    Urj = zeros(Nbeta,Nr); Vrj = Urj; Wrj = Urj; Trj = Urj;

    % Inlet conditions
    if abs(Beta(end)+beta(1))~=0
        % cubic "extrapolation"
        dVdr = 0.5*(3*VIMP(jj,ii)-4*VIMP(jj,ii-1)+VIMP(jj,ii-2))/(Eta(ii)-Eta(ii-1));
        dTdr = 0.5*(3*TIMP(jj,ii)-4*TIMP(jj,ii-1)+TIMP(jj,ii-2))/(Eta(ii)-Eta(ii-1));
        a = Eta(ii); b = Eta(end); 
        A = [VIMP(jj,ii),TIMP(jj,ii)]; B = 1e-8; dA = [dVdr,dTdr];
        P = (A.*log(B./A) + (a-b)*dA)./(A*(a-b)^2);
        Q = (-2*a*A.*log(B./A) + (b^2-a^2)*dA)./(A*(a-b)^2);
        R = (A.*(a^2*log(B)+b*log(A)*(b-2*a)) + a*b*(a-b)*dA)./(A*(a-b)^2);
        VBC = exp(P(1)*Eta(ii:end).^2 + Q(1)*Eta(ii:end) + R(1));
        TBC = exp(P(2)*Eta(ii:end).^2 + Q(2)*Eta(ii:end) + R(2));
    else
        VBC = VIMP(end,ii:end); TBC = TIMP(end,ii:end);
    end
    Uin = -flip(UIMP(1:jj,ii-1:ii))*Re; Vin = flip(VIMP(1:jj,ii-1:ii));
    Win = flip(WIMP(1:jj,ii-1:ii)); Tin = flip(TIMP(1:jj,ii-1:ii));
    Urj(:,1:2) = Uin; Wrj(:,1:2) = Win; 
    Vrj(:,1:2) = Vin; Vrj(1,2:length(VBC)+1) = VBC;
    Trj(:,1:2) = Tin; Trj(1,2:length(TBC)+1) = TBC;

    %% NEWTON-PICARD MULTIVARIATE SCHEME
    str1 = fprintf('Solving r =  %.3f\n',r(2));
    for i=3:Nr
        fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving r =  %.3f\n',r(i));

        % Picard Method
        [Urj,Vrj,Wrj,Trj] = picard(Urj,Vrj,Wrj,Trj,Pr,lambda,r,beta,i,0.5);

        % Newton Method
        [Urj,Vrj,Wrj,Trj] = newton(Urj,Vrj,Wrj,Trj,Pr,lambda,r,beta,i,0.95);

    end
    Urj = Urj(:,2:end); Vrj = Vrj(:,2:end); 
    Trj = Trj(:,2:end); Wrj = Wrj(:,2:end); r = r(2:end);
    fprintf(repmat('\b',1,str1)); 
    str1 = fprintf('All radial distances solved for.\n'); 

    % save file
    VelRJ{1} = Urj; VelRJ{2} = Vrj; VelRJ{3} = Wrj; VelRJ{4} = Trj;
    filename = ['../Flows/RJ/lambda=',num2str(lambda),'/RJ_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    if exist(['../Flows/RJ/lambda=',num2str(lambda)],'dir')==0; mkdir('../Flows/RJ',['lambda=',num2str(lambda)]); end
    save(filename,'VelRJ','r','beta');
    str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end