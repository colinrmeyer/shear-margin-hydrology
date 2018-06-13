% load DAMTP runs and make new plots
clearvars; close all; clc;
load('damtphydrologydataset.mat');
% figure(3);
% EffPbottom(EffPbottom==0)=NaN;
% plot(xgrid(1,:),EffPbottom./1000)
% return
figure(4)
EffP(EffP==0)=NaN;
surf(xgrid,ygrid,EffP*params.N0/1000,'facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$N$ (kPa)','interpreter','latex','fontsize',20); caxis([0,25])
axis([0 55 0 1])
% % plot height of temperate zone % %
zphimat = phi.matrix' ~= 0; yup = ygrid(:,1);
temperateheight = zeros(size(zphimat,2),1);
for i = 1:size(zphimat,2)
    if ~isempty(yup(zphimat(:,i)))
        temperateheight(i) = max(yup(zphimat(:,i)));
    end
end
hold on;
plot3(xgrid(1,:),temperateheight,100*ones(size(xgrid(1,:))),'k','linewidth',4)
set(gca,'fontsize',18)
xlabel('downstream distance, $x$ (km)','interpreter','latex','fontsize',20)
ylabel('height, $z$ (km)','interpreter','latex','fontsize',20)
saveas(4,'downstream_effectivepressure.jpg')

figure(5)
surf(xgrid,ygrid,(Temp.matrix')*params.DeltaT,'facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$T$ ($^{\circ}$C)','interpreter','latex','fontsize',20); %caxis([-25,0])
set(gca,'fontsize',18)
axis([0 55 0 1])
xlabel('downstream distance, $x$ (km)','interpreter','latex','fontsize',20)
ylabel('height, $z$ (km)','interpreter','latex','fontsize',20)
% % plot height of temperate zone % %
hold on;
plot3(xgrid(1,:),temperateheight,100*ones(size(xgrid(1,:))),'k','linewidth',4)
% % Theoretical height of the temperate zone
% eta_diff = max(1 - sqrt( (2*params.K*params.DeltaT) ./ (polyval(params.pW,xgrid(1,:)).*(params.H^2)) ),0);
% Pe = abs(params.Pe); Q = ( (params.K*params.DeltaT) ./ (polyval(params.pW,xgrid(1,:)).*(params.H^2)) ).^(-1);
% Gamma = Pe./Q;
% a = Pe.*Gamma+1;
% eta_acc = max((1 - ((lambertw(-exp(-a))+a)./Pe)),0);
% plot3(xgrid(1,:),eta_diff,100*ones(size(temperateheight)),'y','linewidth',2)
% plot3(xgrid(1,:),eta_acc,100*ones(size(temperateheight)),'g','linewidth',2)
% % % % % % % % % % % % % % % % % % %
saveas(5,'downstream_temperature.jpg')

figure(6)
surf(xgrid,ygrid,(phi.matrix')*params.epsilon*100,'facecolor','interp','edgecolor','none','facelighting','flat')
colormap jet;
view(2); grid off; yc = colorbar; ylabel(yc,'$\phi$ (\%)','interpreter','latex','fontsize',20); caxis([0,5])
hold on;
plot3(xgrid(1,:),temperateheight,100*ones(size(xgrid(1,:))),'k','linewidth',4)
set(gca,'fontsize',18)
axis([0 55 0 1])
xlabel('downstream distance, $x$ (km)','interpreter','latex','fontsize',20)
ylabel('height, $z$ (km)','interpreter','latex','fontsize',20)
saveas(6,'downstream_porosity.jpg')