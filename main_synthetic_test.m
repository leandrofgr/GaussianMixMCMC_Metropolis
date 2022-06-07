
addpath('.\functions\')
load('.\data\data_synth_3layers_oil_water.mat')
load('.\data\cmaps.mat')

% Number of iterations
n_it = 5000;

% Trace to be inverted from the 2D model
trace = 25;

%% Input data
real_seismic = real_seismic_aki; 

SNR = 10;
real_seismic(:,:,1) = real_seismic(:,:,1) + sqrt(mean(var(real_seismic(:,:,1)))/SNR)*noise_mean0_std1(:,:,1);
real_seismic(:,:,2) = real_seismic(:,:,2) + sqrt(mean(var(real_seismic(:,:,2)))/SNR)*noise_mean0_std1(:,:,2);
real_seismic(:,:,3) = real_seismic(:,:,3) + sqrt(mean(var(real_seismic(:,:,3)))/SNR)*noise_mean0_std1(:,:,3);
real_seismic(:,:,4) = real_seismic(:,:,4) + sqrt(mean(var(real_seismic(:,:,4)))/SNR)*noise_mean0_std1(:,:,4);


real_seismic1d(:,:) = real_seismic(:,trace,:);
real_log_vp = log(real_vp(1:end-1,trace));
real_log_vs = log(real_vs(1:end-1,trace));
real_log_rho = log(real_rho(1:end-1,trace));
real_facies_well = real_facies(1:end-1,trace);

I = size(real_log_vp,1);
prob_map = ones(I,1,length(PRIOR_elasticLog))/length(PRIOR_elasticLog);



%% Input parameters
SNR = SNR;
SNR_par = SNR*[1 1 1 1]

% Transition matrix:
P = [0.90    0.055    0.035;
    0.20    0.8    0;
    0.2    0.2    0.60];

PRIOR_ = PRIOR_elasticLog;


[ INVERSION ] = GaussianMixMCMC_metropolis(real_seismic1d, theta, SNR_par, wavelet, PRIOR_, n_it, prob_map, P);


%%   DISPLAY/SHOW RESULTS
INVERSION_ =  INVERSION ;

time_well = [2000:4:2000+(I-1)*4]';
time_well_fine = [2000:2:2000+(I-1)*4+2]';
time = time_well;



figure
ax1 = subplot(1,6,1)
pcolor([1 2],[time_well time_well],[real_facies_well real_facies_well])
shading flat
ylim([time(1) time(end)])
set(gca,'Ydir','reverse')
colormap(ax1,cmap_3facies)
title('Reference facies')
xticks([])
ylabel('Time (ms)','FontSize',12)


ax1 = subplot(1,6,2)
imagesc(INVERSION_.FACIES.samples)
colormap(ax1,cmap_3facies)
title('MChain')

ax1 = subplot(1,6,3)
pcolor([1 2],[time_well time_well],[INVERSION_.FACIES.likely INVERSION_.FACIES.likely])
title('Estimated Facies')
shading flat
set(gca,'Ydir','reverse')
xticks([])
colormap(ax1,cmap_3facies)
c = colorbar;
c.Label.String = 'Shale                       Brine                       Oil';
set(c,'YTick',[])
c.Label.FontSize = 12;

ax1 = subplot(1,6,4)
plot((real_log_vp),time)
pcolor(exp(INVERSION_.VP.axis)',[time_well],[INVERSION_.VP.probability])
%colormap(ax1,cmap_prop)
colormap(ax1,mycmap_prob)
shading flat
hold all
plot(exp(real_log_vp),time,'k','linewidth',2)
plot(exp(INVERSION_.VP.map),time,'b','linewidth',2)
plot(exp(INVERSION_.VP.mean),time,'c','linewidth',2)
plot( exp(PRIOR_(2).MU(1))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(1))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(2).MU(1) + 2*sqrt(PRIOR_(2).C(1,1)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(2).MU(1) - 2*sqrt(PRIOR_(2).C(1,1)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(1) + 2*sqrt(PRIOR_(1).C(1,1)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(1).MU(1) - 2*sqrt(PRIOR_(1).C(1,1)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)

set(gca,'Ydir','reverse')
ylim([time(1) time(end)])
grid
xlabel('P-impedance (m/s g/cm^3)','FontSize',12)
yticks([])

ax1 = subplot(1,6,5)
plot((real_log_vs),time)
pcolor(exp(INVERSION_.VS.axis)',[time_well],[INVERSION_.VS.probability])
shading flat
hold all
plot(exp(real_log_vs),time,'k','linewidth',2)
plot(exp(INVERSION_.VS.map),time,'b','linewidth',2)
plot(exp(INVERSION_.VS.mean),time,'c','linewidth',2)
plot( exp(PRIOR_(2).MU(2))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(2))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(2).MU(2) + 2*sqrt(PRIOR_(2).C(2,2)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(2).MU(2) - 2*sqrt(PRIOR_(2).C(2,2)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(2) + 2*sqrt(PRIOR_(1).C(2,2)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(1).MU(2) - 2*sqrt(PRIOR_(1).C(2,2)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
set(gca,'Ydir','reverse')
ylim([time(1) time(end)])
grid
xlabel('S-impedance (m/s g/cm^3)','FontSize',12)
%colormap(ax1,cmap_prop)
colormap(ax1,mycmap_prob)

ax1 = subplot(1,6,6)
%pcolor(exp(axis_rhob_fine)',[time_well_fine],[probability_rhob_fine])
pcolor(exp(INVERSION_.RHOB.axis)',[time_well],[INVERSION_.RHOB.probability])
shading flat
hold all
plot(exp(real_log_rho),time,'k','linewidth',2)
plot(exp(INVERSION_.RHOB.map),time,'b','linewidth',2)
plot(exp(INVERSION_.RHOB.mean),time,'c','linewidth',2)
plot( exp(PRIOR_(1).MU(3))*ones(size(time)),time,'--')
plot( exp(PRIOR_(2).MU(3))*ones(size(time)),time,'--')
set(gca,'Ydir','reverse')
xlabel('Density (g/cm^3)','FontSize',12)
ylim([time(1) time(end)])
grid
colormap(ax1,mycmap_prob)
c = colorbar;
c.Label.String = 'Probability';
c.Label.FontSize = 12;
set(c,'YTick',[])
plot( exp(PRIOR_(2).MU(3))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(3))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(2).MU(3) + 2*sqrt(PRIOR_(2).C(3,3)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(2).MU(3) - 2*sqrt(PRIOR_(2).C(3,3)))*ones(size(time)),time,'--','color',[1 1 0.15],'linewidth',2)
plot( exp(PRIOR_(1).MU(3) + 2*sqrt(PRIOR_(1).C(3,3)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)
plot( exp(PRIOR_(1).MU(3) - 2*sqrt(PRIOR_(1).C(3,3)))*ones(size(time)),time,'--','color',[0.35 0.8 0.35],'linewidth',2)

legend('Posterior ','Reference model','MAP','Mean','Shale prior','Sand prior')


%%%%%%%%%%%%%%%%
figure
Id = I-1;

n_synth = 500;
n_it = size(INVERSION_.SYNTHETIC(1:Id ,:),2);
iterations = round(  rand(n_synth,1)*(n_it-0.8*n_it) + 0.8*n_it );

sy_mean = mean(INVERSION_.SYNTHETIC(1:Id ,iterations),2);
ax1 = subplot(1,3,1)
plot(real_seismic1d(:,1),time,'k','linewidth',2)
hold all
plot(sy_mean,time(1:end-1),'c','linewidth',2)
grid
%plot(real_seismic(:,1) - INVERSION_.SYNTHETIC(1:I,iterations),repmat(time,1,n_synth),'r','linewidth',1)
plot(INVERSION_.SYNTHETIC(1:Id,iterations),repmat(time(1:end-1),1,n_synth),'r','linewidth',1)
plot(sy_mean,time(1:end-1),'c','linewidth',2)
plot(real_seismic1d(:,1),time,'k','linewidth',2)
ylim([time(1) time(end)])
%xlim([ -0.1 0.1 ])
set(gca,'Ydir','reverse')
title('Angle-stack 8^o')
ylabel('Time (ms)','FontSize',12)

sy_mean = mean(INVERSION_.SYNTHETIC(Id +1:2*Id ,iterations),2);
ax1 = subplot(1,3,2)
plot(real_seismic1d(:,2),time,'k','linewidth',2)
hold all
plot(sy_mean,time(1:end-1),'c','linewidth',2)
grid
%plot(real_seismic(:,2) - INVERSION_.SYNTHETIC(I+1:2*I,iterations),repmat(time,1,n_synth),'r','linewidth',1)
plot(INVERSION_.SYNTHETIC(Id +1:2*Id ,iterations),repmat(time(1:end-1),1,n_synth),'r','linewidth',1)
plot(sy_mean,time(1:end-1),'c','linewidth',2)
plot(real_seismic1d(:,2),time,'k','linewidth',2)
ylim([time(1) time(end)])
%xlim([ -0.1 0.1 ])
set(gca,'Ydir','reverse')
title('Angle-stack 18^o')
%yticks([])

sy_mean = mean(INVERSION_.SYNTHETIC(2*Id +1:3*Id ,iterations),2);
ax1 = subplot(1,3,3)
plot(real_seismic1d(:,3),time,'k','linewidth',2)
hold all
plot(sy_mean,time(1:end-1),'c','linewidth',2)
grid
%plot(real_seismic(:,3) - INVERSION_.SYNTHETIC(2*I+1:3*I,iterations),repmat(time,1,n_synth),'r','linewidth',1)
plot(INVERSION_.SYNTHETIC(2*Id+1:3*Id,iterations),repmat(time(1:end-1),1,n_synth),'r','linewidth',1)
plot(sy_mean,time(1:end-1),'c','linewidth',2)
plot(real_seismic1d(:,3),time,'k','linewidth',2)
ylim([time(1) time(end)])
%xlim([ -0.1 0.1 ])
set(gca,'Ydir','reverse')
title('Angle-stack 28^o')
legend('Reference seismic', 'Synthetics mean' ,'Synthetic seismic')
%yticks([])


figure
semilogx(INVERSION_.log_likelihood,'LineWidth',2)
grid
xlabel('MCMC Steps/Iteration')
ylabel('Log Likelihood ')
