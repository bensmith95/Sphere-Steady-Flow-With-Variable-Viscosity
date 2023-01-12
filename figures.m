%% PLOT FIGURES

clear all; clc; close all

% Reynolds number
Re = sqrt(1e4);

% Prandtl Number
Pr = 0.72;

% Non-dim sensitivity contant
lambda = 0.5;

% display info
set(0,'units','pixels'); disp = get(0,'ScreenSize');
Lx = disp(3); Ly = disp(4); 

%% Plot Boundary Layer Region
% Note - not to scale

% load file
BLfile = ['Flows/BL/BL_Pr=',num2str(Pr),'_lambda=',num2str(lambda),'.mat']; load(BLfile); 
UBL = VelBL{1}; VBL = VelBL{2}; WBL = VelBL{3}; 
PBL = VelBL{4}; TBL = VelBL{5}; dPBL = VelBL{14};
BLeta = eta; BLtheta = theta;

% initialise co-ords
rBL = 1+BLeta/50;
XBL = repmat(sin(BLtheta)',1,length(rBL)); YBL = repmat(cos(BLtheta)',1,length(rBL));
for i=1:length(rBL); XBL(:,i) = rBL(i)*XBL(:,i); YBL(:,i) = rBL(i)*YBL(:,i); end

figure(1); t1 = tiledlayout(2,3); 
TT = ['Heated Boundary Layer Flow: $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
title(t1,TT,'interpreter','latex');

% plot U
h1(1) = nexttile(t1); hold on; set(gca,'Color','none'); x = 0:1e-3:1; area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,UBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]);
set(get(cb,'label'),'string','(i) $U$','interpreter','latex');

% plot V
h1(2) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,VBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
set(get(cb,'label'),'string','(ii) $V$','interpreter','latex');

% plot P
h1(3) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,PBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([-3 0]);
set(get(cb,'label'),'string','$\bar{P}$','interpreter','latex');

% plot W
h1(4) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,WBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([-1 2.2]); 
set(get(cb,'label'),'string','(iii) $\overline{W}$','interpreter','latex');

% plot T
h1(5) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,TBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]);
set(get(cb,'label'),'string','(iv) $T$','interpreter','latex');

% plot dP/dtheta
h1(6) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,dPBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); 
set(get(cb,'label'),'string','$\partial_\theta\bar{P}$','interpreter','latex');

% settings
set(h1,'Colormap', jet)
t1.TileSpacing = 'tight'; t1.Padding = 'compact'; 
set(gcf, 'Position',  [0.075*Lx, 0.075*Ly, 0.825*Lx, 0.825*Ly])

% labelling
grey = [0.4 0.4 0.4]; ticks = [1 1+7.5/50 1+15/50 1+22.5/50 1+30/50]; 
for i=1:6
    nexttile(i); 
    % axis ticks
    set(gca,'Xtick',[]); 
    yticks(ticks); yticklabels({'0','7.5','15','22.5','30'}); 
    yh = ylabel('$\eta$','interpreter','latex'); 
    yh.Position(1)=-0.15; yh.Position(2)=1+13/50; yh.Rotation=0;
    % grid
    plot([0 0], [0 ticks(end)], 'color', grey);
    plot([0 ticks(end)], [0 0], 'color', grey);
    for j=ticks
        x = 0:1e-3:j; plot(x,sqrt(j^2-x.^2),'color', grey)
    end
    for k=15:15:75
        plot([cos(k*pi/180) ticks(end)*cos(k*pi/180)], [sin(k*pi/180) ticks(end)*sin(k*pi/180)], 'color', grey)
        text(0.96*sin(k*pi/180), 0.96*cos(k*pi/180), [num2str(k),'^o'],'HorizontalAlignment','right')
    end
end

%% Plot Impinging Region

% Load flow
IMPfile = ['Flows/IMP/lambda=',num2str(lambda),'/IMP_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
load(IMPfile); UIMP = VelIMP{1}; VIMP = VelIMP{2}; WIMP = VelIMP{3}; TIMP = VelIMP{4}; 
Psi = VelIMP{5}; Omega = VelIMP{6}; PIMP = VelIMP{7}; IMPeta = eta; IMPbeta = beta;

% Plot U-W Velocity Magnitude & vector field
figure(); TT = ['Impinging Region Velocity Field: $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
fn = pcolor(IMPeta,beta,sqrt(WIMP.^2+UIMP.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
fc = colorbar('location','southoutside'); caxis([0,0.35]); colormap('jet');
set(get(fc,'label'),'string','$\sqrt{W^2+U^2}$','interpreter','latex');
hold on; ii = 16; jj = 16; 
quiver(IMPeta(1:ii:end),beta(1:jj:end),WIMP(1:jj:end,1:ii:end),UIMP(1:jj:end,1:ii:end),'color','k'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]); 

% Plot Stream Function
figure(); TT = ['$\psi(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
contourf(IMPeta,beta,Psi,15,'LineColor','none' ); 
colorbar('location','southoutside'); caxis([0,0.8]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot Vorticity
figure(); TT = ['$\omega(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
contourf(IMPeta,beta,Omega,15,'LineColor','none'); 
colorbar('location','southoutside'); caxis([-0.15,0.1]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot V velocity component 
figure(); TT = ['$V_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
fn = pcolor(IMPeta,beta,VIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot W velocity component 
figure(); TT = ['$W_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
Wp = WIMP; Wp(Wp<0)=NaN; % remove negative/reverse flow
fn = pcolor(IMPeta,beta,Wp); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.35]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot U velocity component 
figure(); TT = ['$U_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
fn = pcolor(IMPeta,beta,-UIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.15]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot T 
figure(); TT = ['$T_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
fn = pcolor(IMPeta,beta,TIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot Pressure 
figure(); TT = ['$\bar{P}_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];
fn = pcolor(eta,beta,PIMP/sqrt(Re)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([-0.015,0.01]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

%% Plot radial jet figures

% load radial jet
RJfile = ['Flows/RJ/lambda=',num2str(lambda),'/RJ_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
load(RJfile); Urj = VelRJ{1}; Vrj = VelRJ{2}; Wrj = VelRJ{3}; Trj = VelRJ{4}; r1 = r;

% Plot W component
figure(); TT = ['$W_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)]; 
fn = pcolor(r1,beta,Wrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,0.35]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 10]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot U component
figure(); TT = ['$\overline{U}_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)]; 
fn = pcolor(r1,beta,Urj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,3]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 2]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot V component
figure(); TT = ['$V_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)]; 
fn = pcolor(r1,beta,Vrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,1]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 2]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

% Plot T component
figure(); TT = ['$T_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)]; 
fn = pcolor(r1,beta,Trj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,1]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 2]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

%% plot r-> InF Solutions
figure(); t = tiledlayout(2,2); 
TT = ['$r\rightarrow\infty$ solutions at $\sqrt{R_e}=$',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)];

% Calculate relevant constants
M = trapz(beta,r1(end)^2*Wrj(:,end).^2); N = trapz(beta,r1(end)^3*Wrj(:,end).*Vrj(:,end));
K = sqrt(2)/(M*3)^(1/3); sigma = N/M; 
C = (gamma(0.5+Pr)*(0.5+Pr)/(sqrt(pi)*gamma(1+Pr)))*sqrt(2)/K*trapz(beta,Trj(:,end).*sech(beta/(sqrt(2)*K))'.^2);

% Plot U(r->inf,beta)
nexttile(t); plot(beta,Urj(:,end),'-b'); hold on;
Uinf = -sqrt(2)/(K*r1(end))*tanh(beta/(sqrt(2)*K));
plot(beta(1:4:end),Uinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','southwest'); 
title('(a) $\overline{U}_\infty$','interpreter','latex')
xlabel('\beta'); 

% Plot V(r->inf,beta)
nexttile(t); plot(beta,Vrj(:,end),'-b'); hold on;
Vinf = sigma/(K^2*r1(end)^2)*sech(beta/(sqrt(2)*K)).^2;
plot(beta(1:4:end),Vinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(b) $V_\infty$','interpreter','latex')
xlabel('\beta');

% Plot W(r->inf,beta)
nexttile(t); plot(beta,Wrj(:,end),'-b'); hold on;
Winf = 1/(K^2*r1(end))*sech(beta/(sqrt(2)*K)).^2;
plot(beta(1:4:end),Winf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(c) $W_\infty$','interpreter','latex')
xlabel('\beta');

% Plot T(r->inf,beta)
nexttile(t); plot(beta,Trj(:,end),'-b'); hold on;
Tinf = C*sech(beta/(sqrt(2)*K)).^(2*Pr);
plot(beta(1:4:end),Tinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(d) $T_\infty$','interpreter','latex')
xlabel('\beta');
 
t.TileSpacing = 'compact'; t.Padding = 'tight';
title(t,TT,'interpreter','latex');
set(gcf,'Position',  [0.25*Lx, 0.15*Ly, 0.5*Lx, 0.7*Ly])

%% Plot full flow figures

% load full flow
FULLfile = ['Flows/FULL/lambda=',num2str(lambda),'/FULL_Pr=',num2str(Pr),'_Re=',num2str(Re),'_lambda=',num2str(lambda),'.mat'];
load(FULLfile); U = VelFULL{1}; V = VelFULL{2}; W = VelFULL{3}; T = VelFULL{4};

% transform to cartesian coords
X = repmat(sin(theta)',1,length(r)); Y = repmat(cos(theta)',1,length(r));
for i=1:length(r); X(:,i) = r(i)*X(:,i); Y(:,i) = r(i)*Y(:,i); end

% Plot velocity components
figure(); t1 = tiledlayout(2,2); 
title(t1,['Velocity Components at $\sqrt{R_e}$=',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)],'interpreter','latex');

% plot U
h(1) = nexttile(t1);
fn = pcolor(X,Y,U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]); 
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(i) $U$','interpreter','latex');

% plot V
h(2) = nexttile(t1); 
fn = pcolor(X,Y,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(ii) $V$','interpreter','latex');

% plot W
h(3) = nexttile(t1); 
fn = pcolor(X,Y,W); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 0.35]) 
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(iii) $W$','interpreter','latex');

% plot T
h(4) = nexttile(t1); 
fn = pcolor(X,Y,T); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(iv) $T$','interpreter','latex');

% settings
set(h,'Colormap', jet)
t1.TileSpacing = 'compact'; t1.Padding = 'compact'; 
set(gcf,'Position',  [0.2*Lx, 0.1*Ly, 0.6*Lx, 0.8*Ly])

%% Plot Velocity magnitudes and vector field

figure(); t2 = tiledlayout(1,2); 
title(t2,['Flow at $\sqrt{R_e}=$',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)],'interpreter','latex')

% plot planar magnitude
h(1) = nexttile(t2);
fn = pcolor(X,Y,sqrt(U.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
xlim([0 2]); ylim([0 1.25]);
set(get(cb,'label'),'string','$\sqrt{U^2+W^2}$','interpreter','latex');

% plot velcoity magnitude and vector field
h(2) = nexttile(t2); 
fn = pcolor(X,Y,sqrt(U.^2+V.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.2]); ylim([0 1.05]);
set(get(cb,'label'),'string','$||${\boldmath$u$}$||$','interpreter','latex');
hold on;
% convert to cartesian coords
u = W.*sin(theta') + U.*cos(theta'); v = W.*cos(theta') - U.*sin(theta');
quiver(X(1:40:end,1:40:end),Y(1:40:end,1:40:end),u(1:40:end,1:40:end),v(1:40:end,1:40:end),'color',[0.7 0.7 0.7],'autoscalefactor',0.1)

% settings
set(h,'Colormap', jet)
t2.TileSpacing = 'compact'; t2.Padding = 'compact'; 
set(gcf,'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])

%% Plot Velocities near Equator

figure(); t3 = tiledlayout(2,2); 
title(t3,['Equatorial Velocity Components at $\sqrt{R_e}$=',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)],'interpreter','latex');

% plot U
h(1) = nexttile(t3);
fn = pcolor(X,Y,U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(i) $U$','interpreter','latex');

% plot V
h(2) = nexttile(t3); 
fn = pcolor(X,Y,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(ii) $V$','interpreter','latex');

% plot W
h(3) = nexttile(t3); 
fn = pcolor(X,Y,W); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 0.35]) 
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(iii) $W$','interpreter','latex');

% plot T
h(4) = nexttile(t3); 
fn = pcolor(X,Y,T); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(iv) $T$','interpreter','latex');

% settings
set(h,'Colormap', jet)
t3.TileSpacing = 'compact'; t3.Padding = 'compact'; 
set(gcf,'Position',  [0.2*Lx, 0.1*Ly, 0.6*Lx, 0.8*Ly])

%% Plot Velocity magnitudes and vector field near Equator 

figure(); t4 = tiledlayout(1,2); 
title(t4,['Equatorial Flow at $\sqrt{R_e}=$',num2str(Re),', $P_r=$ ',num2str(Pr),', $\lambda=$ ',num2str(lambda)],'interpreter','latex')

% plot planar magnitude
h(1) = nexttile(t4);
fn = pcolor(X,Y,sqrt(U.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','$\sqrt{U^2+W^2}$','interpreter','latex');

% plot velcoity magnitude and vector field
h(2) = nexttile(t4); 
fn = pcolor(X,Y,sqrt(U.^2+V.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','$||${\boldmath$u$}$||$','interpreter','latex');
hold on;
quiver(X(1:20:end,1:20:end),Y(1:20:end,1:20:end),u(1:20:end,1:20:end),v(1:20:end,1:20:end),'color','k','autoscalefactor',0.05)

% settings
set(h,'Colormap', jet)
t4.TileSpacing = 'compact'; t4.Padding = 'compact'; 
set(gcf, 'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])
