%% Obtains and assembles the full flow 
% Inputs:
% Re - Renolds number (sqrt)
% Pr - Prandtl number
% lambda - non-dim sensitivity constant
function BaseFlow(Re,Pr,lambda)
    fprintf('\nObtaining Base Flow:\n')
    %% Obtain Boundary Layer region 
    % Check if data aleady exists
    BLfile = ['Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat'];
    if exist(BLfile,'file')
        str1 = fprintf('1) Boundary Layer Region File Found.\n'); pause(1)
    else
        % Solve 
        str1 = fprintf('1) Solving Boundary Layer Region...\n');
        cd 'BL Files'; BLRegion(Pr,lambda); cd ..
    end
    fprintf(repmat('\b',1,str1)); fprintf('1) Obtained Boundary Layer Region.\n')

    %% Obtain Impinging region
    % Check if data aleady exists
    IMPfile = ['Flows/IMP/lambda=',num2str(lambda),'/IMP_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    if exist(IMPfile,'file')
        str2 = fprintf('2) Impinging Region File Found.\n'); pause(1)
    else
        % Solve
        str2 = fprintf('2) Solving Impinging Region...\n');
        cd 'IMP Files'; IMPRegion(Re,Pr,lambda); cd ..
    end
    fprintf(repmat('\b',1,str2)); fprintf('2) Obtained Impinging Region.\n')

    %% Obtain Radial Jet region
    % Check if data aleady exists
    RJfile = ['Flows/RJ/lambda=',num2str(lambda),'/RJ_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    if exist(RJfile,'file')
        str3 = fprintf('3) Radial Jet Region File Found.\n'); pause(1)
    else
        % Solve
        str3 = fprintf('3) Solving Radial Jet Region...\n');
        cd 'RJ Files'; RadialJet(Re,Pr,lambda); cd ..
    end
    fprintf(repmat('\b',1,str3)); fprintf('3) Obtained Radial Jet Region.\n')
    
    %% Assemble Full Base Flow
    % Check if data aleady exists
    FULLfile = ['Flows/FULL/lambda=',num2str(lambda),'/FULL_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
    if exist(FULLfile,'file')
        str4 = fprintf('4) Full Base Flow File Found.\n'); pause(1)
    else
        str4 = fprintf('4) Building Full Base Flow...\n');
        % load Boundary Layer region
        load(BLfile); UBL = VelBL{1}; VBL = VelBL{2}; 
        WBL = VelBL{3}; TBL = VelBL{5}; BLeta = eta;
        % load Impinging region
        load(IMPfile); UIMP = -flip(VelIMP{1}); VIMP = flip(VelIMP{2}); 
        WIMP = flip(VelIMP{3}); TIMP = flip(VelIMP{4});
        IMPeta = eta; IMPbeta = beta;
        % load Radial Jet region
        beta = []; load(RJfile); URJ = VelRJ{1}; VRJ = VelRJ{2}; 
        WRJ = VelRJ{3}; TRJ = VelRJ{4};
        
        % Sigmoid function
        re = 2*(1+1/sqrt(Re))-r(1); [~,re] = min(abs(r-re));
        ii = max(Re*UIMP,[],1); ii = find(ii>3); ii = ii(end); % start of RJ
        jj = length(beta); % truncate to remove O(eps) reverse flow 
        sigmoid = (2*r(1:re).^3 - 3*(r(1)+r(re))*r(1:re).^2 + 6*r(1)*r(re)*r(1:re) + r(1)^2*(r(1)-3*r(re)))/(r(1)-r(re))^3;
        for k=1:length(sigmoid)
            %Wsig(:,k) = (1-sigmoid(k))*WIMP(end-jj+1:end,ii+k-1) + sigmoid(k)*WRJ(:,k);
            Usig(:,k) = (1-sigmoid(k))*UIMP(end-jj+1:end,ii+k-1) + sigmoid(k)*URJ(:,k)/Re;
            %Vsig(:,k) = (1-sigmoid(k))*VIMP(end-jj+1:end,ii+k-1) + sigmoid(k)*VRJ(:,k);
            %Tsig(:,k) = (1-sigmoid(k))*TIMP(end-jj+1:end,ii+k-1) + sigmoid(k)*TRJ(:,k);
        end
        %WRJ = [WIMP(end-jj+1:end,1:ii-1) Wsig WRJ(:,re+1:end)];
        WRJ = [WIMP(end-jj+1:end,1:ii-1) WRJ];
        URJ = [UIMP(end-jj+1:end,1:ii-1) Usig URJ(:,re+1:end)/Re];
        %VRJ = [VIMP(end-jj+1:end,1:ii-1) Vsig VRJ(:,re+1:end)];
        VRJ = [VIMP(end-jj+1:end,1:ii-1) VRJ];
        %TRJ = [TIMP(end-jj+1:end,1:ii-1) Tsig TRJ(:,re+1:end)];
        TRJ = [TIMP(end-jj+1:end,1:ii-1) TRJ];
        r = [IMPeta(1:ii)/Re+1 r(2:end)];

        % beta->Inf extrapolation
        if abs(-IMPbeta(end)-beta(1))~=0
            % mesh grid
            x = repmat(IMPeta,length(IMPbeta),1); y = -repmat(flip(IMPbeta(:)),1,length(IMPeta));
            % add IMP region
            Sv = zeros(size(x)); Sv(:,1:ii) = [VIMP(1:end-jj,1:ii); VRJ(:,1:ii)];
            Sw = zeros(size(x)); Sw(:,1:ii) = [WIMP(1:end-jj,1:ii); WRJ(:,1:ii)];
            St = zeros(size(x)); St(:,1:ii) = [TIMP(1:end-jj,1:ii); TRJ(:,1:ii)];
            X1 = x(:,1:ii); X1 = X1(:); Y1 = y(:,1:ii); Y1 = Y1(:); 
            Sv1 = Sv(:,1:ii); Sv1 = Sv1(:); Sw1 = Sw(:,1:ii); Sw1 = Sw1(:); St1 = St(:,1:ii); St1 = St1(:);
            % add replacement RJ in IMP region 
            Sv(end-jj+1:end,ii+1:end) = VRJ(:,ii+1:length(IMPeta)); Sv(1,ii+1:end) = VIMP(1,ii+1:end);
            Sw(end-jj+1:end,ii+1:end) = WRJ(:,ii+1:length(IMPeta)); Sw(1,ii+1:end) = WIMP(1,ii+1:end);
            St(end-jj+1:end,ii+1:end) = TRJ(:,ii+1:length(IMPeta)); St(1,ii+1:end) = TIMP(1,ii+1:end);
            X2 = [x(end-jj+1:end,ii+1:end); x(1,ii+1:end)]; X2 = X2(:);   
            Y2 = [y(end-jj+1:end,ii+1:end); y(1,ii+1:end)]; Y2 = Y2(:); 
            Sv2 = [Sv(end-jj+1:end,ii+1:end); Sv(1,ii+1:end)]; Sv2 = Sv2(:);
            Sw2 = [Sw(end-jj+1:end,ii+1:end); Sw(1,ii+1:end)]; Sw2 = Sw2(:);
            St2 = [St(end-jj+1:end,ii+1:end); St(1,ii+1:end)]; St2 = St2(:);
            X3 = x(2:end-jj,end); Y3 = y(2:end-jj,end); 
            Sv3 = Sv(2:end-jj,end); Sw3 = Sw(2:end-jj,end); St3 = St(2:end-jj,end);
            %  extrapolation using 'natural' Delaunay triangulation(?)
            I = scatteredInterpolant([X1;X2;X3],[Y1;Y2;Y3],[Sv1;Sv2;Sv3],'natural'); 
            Sv(2:end-jj,ii+1:end-1) = I(x(2:end-jj,ii+1:end-1),y(2:end-jj,ii+1:end-1));
            I = scatteredInterpolant([X1;X2;X3],[Y1;Y2;Y3],[Sw1;Sw2;Sw3],'natural'); 
            Sw(2:end-jj,ii+1:end-1) = I(x(2:end-jj,ii+1:end-1),y(2:end-jj,ii+1:end-1));
            I = scatteredInterpolant([X1;X2;X3],[Y1;Y2;Y3],[St1;St2;St3],'natural'); 
            St(2:end-jj,ii+1:end-1) = I(x(2:end-jj,ii+1:end-1),y(2:end-jj,ii+1:end-1));
            % replace O(eps) flow with extrapolated
            Z = zeros(length(IMPbeta)-jj,length(r)-length(IMPeta));
            WRJ = [Sw [Z; WRJ(:,length(IMPeta)+1:end)]];
            VRJ = [Sv [Z; VRJ(:,length(IMPeta)+1:end)]];
            TRJ = [St [Z; TRJ(:,length(IMPeta)+1:end)]];

            % Quadratic "extrapolation" for U
            Su  = zeros(length(IMPbeta(jj:end)),length(r(ii+1:end)));
            for n=ii+1:length(r)
               dUdb = -0.5*(3*URJ(1,n)-4*URJ(2,n)+URJ(3,n))/(IMPbeta(jj-1)-IMPbeta(jj));
               % Quadratic
               Su(:,n-ii) = ( 0.5*dUdb*IMPbeta(jj:end).^2 - IMPbeta(end)*dUdb*IMPbeta(jj:end) + ...
                              URJ(1,n)*(IMPbeta(jj)-IMPbeta(end)) - 0.5*IMPbeta(jj)*dUdb*(IMPbeta(jj)-2*IMPbeta(end)) )/(IMPbeta(jj)-IMPbeta(end));

            end
            URJ = [UIMP(1:end-jj+1,1:ii) flip(Su); URJ(2:end,:)];
        end

        % interpolate to same r
        [~,idx] = min(abs((r-1)*Re-30)); % determine etamax in r
        UBL = spline(BLeta,UBL,(r(1:idx)-1)*Re); VBL = spline(BLeta,VBL,(r(1:idx)-1)*Re);
        WBL = spline(BLeta,WBL,(r(1:idx)-1)*Re); TBL = spline(BLeta,TBL,(r(1:idx)-1)*Re);
        % determine beta->Inf in theta
        IMPbeta = -flip(IMPbeta);
        thetaIN = -max(abs(IMPbeta))/Re + pi/2; [~,j] = min(abs(theta-thetaIN)); 
        UBL1 = (spline(theta(j:end),UBL(j:end,ii:end)',IMPbeta/Re+pi/2))';
        WBL1 = (spline(theta(j:end),WBL(j:end,:)',IMPbeta/Re+pi/2))';
        theta = [theta(1:j-1) IMPbeta/Re+pi/2]; 
        
        % BL Sigmoid function
        Wsig = []; Usig = [];
        WRJ(1,1:idx) = WBL(j,:)/Re; m = 1+1/(IMPbeta(2)-IMPbeta(1));
        sigmoid = ( 2*IMPbeta(1:m).^3 - 3*(IMPbeta(1)+IMPbeta(m))*IMPbeta(1:m).^2 + ...
                    6*IMPbeta(1)*IMPbeta(m)*IMPbeta(1:m) + IMPbeta(1)^2*(IMPbeta(1)-3*IMPbeta(m)) )/(IMPbeta(1)-IMPbeta(m))^3;
        for k=1:length(sigmoid)
            Wsig(k,:) = (1-sigmoid(k))*WBL1(k,:)/Re + sigmoid(k)*WRJ(k,1:idx); %j+k-1
            Usig(k,:) = (1-sigmoid(k))*UBL1(k,:) + sigmoid(k)*URJ(k,ii:idx);
        end
        WRJ(1:m,1:idx) = Wsig; URJ(1:m,ii:idx) = Usig;
        
        % attach RJ region to BL region
        UFULL = zeros(length(theta),length(r)); VFULL = UFULL; WFULL = UFULL; TFULL = UFULL;
        UFULL(1:j-1,1:idx) = UBL(1:j-1,:); UFULL(j:end,:) = URJ; 
        VFULL(1:j-1,1:idx) = VBL(1:j-1,:); VFULL(j:end,:) = VRJ;
        WFULL(1:j-1,1:idx) = WBL(1:j-1,:)/Re; WFULL(j:end,:) = WRJ;
        TFULL(1:j-1,1:idx) = TBL(1:j-1,:); TFULL(j:end,:) = TRJ;
        % save data
        VelFULL{1} = UFULL; VelFULL{2} = VFULL; VelFULL{3} = WFULL; VelFULL{4} = TFULL;
        if exist(['Flows/FULL/lambda=',num2str(lambda)],'dir')==0; mkdir('Flows/FULL',['lambda=',num2str(lambda)]); end
        save(FULLfile,'VelFULL','r','theta');
        fprintf(repmat('\b',1,str4)); str4 = fprintf('Flow saved in %s\n', FULLfile); pause(1) 
    end
    fprintf(repmat('\b',1,str4)); fprintf('4) Obtained Full Base Flow.\n\n')
end