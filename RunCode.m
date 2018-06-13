% simple 2D finite volume code
clearvars; close all; clc;

% Mesh Information
% domain from from x=a to x=b and y=c to y=d
a = 0;
b = 54;
Nx = 64; % 256; % Number of x grid cells
dx = (b-a)/Nx; % x grid cell width
c = 0;
d = 1;
Ny = 32; % 128; % Number of y grid cells
dy = (d-c)/Ny; % y grid cell width
[xgrid,ygrid]=meshgrid((a+(dx/2)):dx:(b-(dx/2)),(c+dy/2):dy:(d-(dy/2)));
N = Nx*Ny; % Total number of entries

% Time step information
T = 3; % Final time
dt = 0.25*min(dx,dy).^2;
Nt = T/dt; % Number of timesteps
tgrid=linspace(0,T,Nt);

% % Parameters % %
loadparameters;
Qup = 10^(-7);
accumulation = params.Pe*ones(Nx,Ny+1);
downstreamvel = (params.uvel./params.accumulation)*ones(Nx+1,Ny);
Tsurf = -ones(Nx,1);
Tbed = zeros(Nx,1);
etamat = ones(Nx,Ny);
shearheat = 2*(params.A.^(-1/params.n)).*...
    ((polyval(params.edotxy,xgrid')).^((params.n+1)/params.n))*(params.H^2)./(params.DeltaT*params.K);

% Diffusion matrix
lambda.matrix = ones(Nx,Ny);
Tran.xplus = -(2*dy/dx)*([lambda.matrix(2:Nx,:);Inf*ones(1,Ny)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.xminus = -(2*dy/dx)*([Inf*ones(1,Ny);lambda.matrix(1:(Nx-1),:)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.yplus = -(2*dx/dy)*([lambda.matrix(:,2:Ny),Inf*ones(Nx,1)].^(-1)+lambda.matrix.^(-1)).^(-1);
Tran.yminus = -(2*dx/dy)*([Inf*ones(Nx,1),lambda.matrix(:,1:(Ny-1))].^(-1)+lambda.matrix.^(-1)).^(-1);

h.matrix = zeros(size(xgrid')); basalN = 0; drainQ = zeros(size(xgrid(1,:)))'; SavePressure = ones(Nx,Ny)';
for n = 1:Nt
    % enthalpy
    Temp.matrix = min(h.matrix,0); Temp.real = params.Tmelt + params.DeltaT*Temp.matrix;
    phi.matrix = max(h.matrix,0); phi.real = params.epsilon*100*phi.matrix;
    % vertical advection
    adfv = (accumulation(:,1:Ny).*h.matrix-accumulation(:,2:(Ny+1)).*[h.matrix(:,2:Ny),h.matrix(:,Ny)])*dx;
    % horizontal advection
    adfh = -(downstreamvel(2:(Nx+1),:).*h.matrix-downstreamvel(1:Nx,:).*[h.matrix(1,:); h.matrix(1:(Nx-1),:)])*dy;
    % diffusion
    uxp = Tran.xplus.*([Temp.matrix(2:Nx,:); Temp.matrix(Nx,:)]-Temp.matrix);
    uxm = Tran.xminus.*(Temp.matrix-[Temp.matrix(1,:); Temp.matrix(1:(Nx-1),:)]);
    uyp = Tran.yplus.*([Temp.matrix(:,2:Ny),Tsurf]-Temp.matrix);
    uym = Tran.yminus.*(Temp.matrix-[Tbed,Temp.matrix(:,1:(Ny-1))]);
    % Effective pressure computation
    Iphi = reshape((phi.matrix~=0)',N,1);
    if (n*dt)<=(1/2)
        basalN = ((1/2)*(xgrid(1,:).^2));
        [EffP,drainQ] = Nfun2d(dx,dy,Nx,Ny,Iphi,zeros(Nx,Ny),phi.matrix,params.alpha,etamat,params.delta,params.g,params.kappa,basalN');
    elseif ~mod(n,50)
        basalN = subglacialhydrology(xgrid(1,:),abs(drainQ)*params.qscale,100*1000,Qup)./params.N0; basalN(isnan(basalN))=0;
        [EffP,drainQ] = Nfun2d(dx,dy,Nx,Ny,Iphi,zeros(Nx,Ny),phi.matrix,params.alpha,etamat,params.delta,params.g,params.kappa,basalN');
        % warning('off')
    else
        [EffP,drainQ] = Nfun2d(dx,dy,Nx,Ny,Iphi,zeros(Nx,Ny),phi.matrix,params.alpha,etamat,params.delta,params.g,params.kappa,basalN');
    end
    % compute enthalpy
    h.matrix = h.matrix-(dt/(dx*dy))*(uxp-uxm+uyp-uym)+(dt/(dx*dy))*(adfh+adfv)+...
        dt*(shearheat - (EffP'.*phi.matrix./etamat));
    if ~mod(n,4000)
        % difference
        pressure_error = sum(sum(sqrt((SavePressure-EffP).^2)./mean(mean(EffP))));
        SavePressure = EffP;
        % figure
        subplot(2,1,1)
        surf(xgrid,ygrid,(EffP)*params.N0/1000,'facecolor','interp','edgecolor','none','facelighting','flat')
        title(['t = ',num2str(n*dt)],'interpreter','latex','fontsize',20)
        view(2); grid off; hc = colorbar; ylabel(hc,'$N$ (kPa)','interpreter','latex','fontsize',20)
        axis([0 54 0 1]); set(gca,'fontsize',18);
        ylabel('height, $z$ (km)','interpreter','latex','fontsize',20)
        subplot(2,1,2)
        plot(xgrid(1,:),params.qscale*10*1000*abs(drainQ),'k','linewidth',2)
        title(['error = ',num2str(pressure_error)],'interpreter','latex','fontsize',20)
        hc = colorbar; set(hc,'visible','off')
        xlabel('downstream distance, $x$ (km)','interpreter','latex','fontsize',20)
        ylabel('flux, Q (m$^3$/s)','interpreter','latex','fontsize',20)
        set(gca,'fontsize',18)
        axis([0 54 0 10^(-5)])
        drawnow;
        saveas(gcf,['ds_effp',num2str(n),'.jpg'])
    end
end
EffPbottom = subglacialhydrology(xgrid(1,:),abs(drainQ)*params.qscale,100*1000,Qup); % EffPbottom(isnan(EffPbottom))=0;
% % % % % % % % % % % % % % % % % % %
save('hydrologydataset')
% % % % % % % % % % % % % % % % % % %
figure(4)
surf(xgrid,ygrid,EffP*params.N0/1000,'facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$N$ (kPa)','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'xlim',[0 55]); caxis([0,15])
xlabel('$x$ (km)','interpreter','latex','fontsize',20)
ylabel('$z$ (km)','interpreter','latex','fontsize',20)
saveas(4,'downstream_effectivepressure.jpg')
% % % % % % % % % % % % % % % % % % %
figure(5)
surf(xgrid,ygrid,(Temp.matrix')*params.DeltaT,'facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$T$ ($^{\circ}$C)','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'xlim',[0 55])
xlabel('$x$ (km)','interpreter','latex','fontsize',20)
ylabel('$z$ (km)','interpreter','latex','fontsize',20)
% % plot height of temperate zone % %
zTmat = (Temp.matrix') == 0;
for i = 1:size(zTmat,2)
    if ~isempty(find(zTmat(:,i),1))
        temperateheight(i) = ygrid(find(zTmat(:,i),1,'last'),i);
    end
end
hold on;
plot3(xgrid(1,:),temperateheight,100*ones(size(temperateheight)),'k','linewidth',4)
saveas(5,'downstream_temperature.jpg')
% % % % % % % % % % % % % % % % % % %
figure(6)
surf(xgrid,ygrid,phi.real','facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$\phi$ (\%)','interpreter','latex','fontsize',20)
set(gca,'fontsize',18,'xlim',[0 55])
xlabel('$x$ (km)','interpreter','latex','fontsize',20)
ylabel('$z$ (km)','interpreter','latex','fontsize',20)
saveas(6,'downstream_porosity.jpg')