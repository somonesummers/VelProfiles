% Starting sandbox for making synthetic data set of known slope and speeds.
clear
%% Input parameters
depth = 1e3; % Depth of Domain[m]
velMax = 50/(3600*24*365); % [m/s] Max Speed
slopeMax = 1e-1; % [ ] Max slope, assumed to be at bed 
dr = .3; %[m] wavelength of system
n = 100; %[ ] number of modeling depth bins
dz = depth/(n - 1); %[m] size of model depth bins
dt = 24*3600/2; %[s] sample period
nt = 365*1; %must be odd
timeMax = dt*(nt-1);
%% Make synthetic data
% linear layers of known slope. n^4 velocity profile
z = 0:dz:depth;
s = slopeMax/1000:slopeMax/(n):slopeMax; %Linear change in slope
% s = slopeMax:-slopeMax/(n):(slopeMax/1000); %Linear change in slope
% s = slopeMax*(sin(0:2*pi/(n-1):2*pi)) + slopeMax/100*randn(size(z)); %Sine changes in slope
s_clean = s;
s = s + slopeMax/100*randn(size(z));
v_clean = velMax*(z/depth).^(4);
v = v_clean + velMax/100*randn(size(z));
t = (0:dt:timeMax)';

phi = rand(size(v));
dPhi = 2*pi * s .* v ./ dr .* t + 2*pi*phi;
x_clean = cos(dPhi) + 1i*sin(dPhi);
x = cos(dPhi) + 1i*sin(dPhi) + .2*(randn(size(dPhi)) + 1i*randn(size(dPhi)));

%% Solve for velocity
[v_star, m2_v] = fitVelocity(x,z,t,s,velMax,dr);
[v_star2, m2_v2] = fitVelocity_2(x,z,t,s,velMax,dr,m2_v(1,:));

%% Solve for slope
[s_star, m2_s] = fitSlope(x,z,t,slopeMax,v,dr); 

%% Solve for slopeVelocity
[sv_star, m3_s, F1] = fitSV(x,z,t,slopeMax,velMax,dr); 
[sv_star_2, m4_s, F2] = fitSV_2(x,z,t,slopeMax,velMax,dr,movmean(m3_s(1,:),floor(n/5))); 
%% Plot it UP

figure(2)
clf
nfigs  = 6;
di = floor((n)/nfigs);
for i = 1:nfigs
    subplot(nfigs,2,2*i-1)
    plot(t,real(x(:,i*di)),'.','color',rgb('gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m3_s(1,i*di) + 2*pi*m3_s(2,i*di)),'-','color',rgb('dark rose'),'lineWidth',2) 
%     if(i == 1)
%         title('first fits')
%      else
%         title("Layer " +z(i*di))
%     end
end
for i = 1:nfigs
    subplot(nfigs,2,2*i)
    plot(t,real(x(:,i*di)),'.','color',rgb('gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m4_s(1,i*di) + 2*pi*m4_s(2,i*di)),'-','color',rgb('blue'),'lineWidth',2) 
%     if(i == 1)
%        title('second fits')
%     else 
%        title("Layer " +z(i*di))
%     end
end
setFontSize(16)
% savePng('fitn')
% 
%% Plot slopes
figure(3)
clf
plot(abs(s_clean),z,'--','color',rgb('black'),'lineWidth',3)
xlabel('Slope of layer')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
setFontSize(16)
% savePng('figs/slope1')
hold on
plot(abs(s),z,'.','color',rgb('gray'),'markersize',15)
% savePng('figs/slope2')
plot(s_star,z,'.','color',rgb('lilac'),'markersize',15)
plot(movmean(s_star,floor(n/10)),z,'--','color',rgb('lilac'),'lineWidth',3)
% savePng('figs/slope3')

%% Plot Velocities
figure(4)
clf
plot(v,z,'.','color',rgb('gray'),'markerSize',15)
xlabel('Velocity')
xlim([0 velMax])
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
setFontSize(16)
% savePng('figs/vel1')
hold on
plot(v_clean,z,'--','color',rgb('black'),'linewidth',3);
% savePng('figs/vel2')
plot(v_star,z,'.','color',rgb('burnt orange'),'markerSize',15)
plot(movmean(v_star,floor(n/10)),z,'--','color',rgb('burnt orange'),'lineWidth',3)
% savePng('figs/vel3')



%% Plot SV product
figure(5)
clf
plot(abs(s_clean.*v_clean),z,'--','color',rgb('black'),'lineWidth',5)
xlabel('Slope velocity product')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
xlim ([0 max(abs(s_clean.*v_clean))])
setFontSize(16)
% savePng('figs/sv1')
hold on
plot(abs(s.*v),z,'.','color',rgb('gray'),'MarkerSize',25)
% savePng('figs/sv2')
plot(sv_star,z,'.','color',rgb('light rose'),'MarkerSize',25)
plot(movmean(sv_star,floor(n/10)),z,'--','color',rgb('dark rose'),'lineWidth',5)
% savePng('figs/sv3')
plot(sv_star_2,z,'.','color',rgb('baby blue'),'MarkerSize',25)
plot(movmean(sv_star_2,floor(n/10)),z,'--','color',rgb('blue'),'lineWidth',5)
% savePng('figs/sv4')

figure(6)
clf
plot(F1,z,'.','color',rgb('black'),'MarkerSize',25)
xlabel('misfit')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
setFontSize(16)
% savePng('figs/sv1')
hold on
plot(F2,z,'.','color',rgb('gray'),'MarkerSize',25)
