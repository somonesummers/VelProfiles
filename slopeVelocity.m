% Starting sandbox for making synthetic data set of known slope and speeds.
clear
%% Input parameters
depth = 1e3; % Depth of Domain[m]
velMax = 50/(3600*24*365); % [m/s] Max Speed
slopeMax = 1e-1; % [ ] Max slope, assumed to be at bed 
dr = .3; %[m] Rough Range resolution (from radar processsing)
n = 100; %[ ] number of modeling depth bins
dz = depth/(n - 1); %[m] size of model depth bins
dt = 24*3600/2; %[s] sample freq
nt = 365*3; %must be odd
timeMax = dt*(nt-1);
%% Make synthetic data
% linear layers of known slope. n^4 velocity profile
z = 0:dz:depth;
% s = 0:slopeMax/(n-1):slopeMax; %Linear change in slope
% s = -1*(0:slopeMax/(n-1):slopeMax); %Linear change in slope
s = slopeMax*(sin(0:2*pi/(n-1):2*pi)) + slopeMax/100*randn(size(z)); %Sine changes in slope
v = velMax*(z/depth).^(4) + velMax/100*randn(size(z));
t = (0:dt:timeMax)';

phi = rand(size(v));
dPhi = 2*pi * s .* v ./ dr .* t + 2*pi*phi;
x_clean = cos(dPhi) + 1i*sin(dPhi);
x = cos(dPhi) + 1i*sin(dPhi) + .1*(randn(size(dPhi)) + 1i*randn(size(dPhi)));

%% Solve for velocity
[v_star, m2_v] = fitVelocity(x,z,t,s,velMax,dr); 

%% Solve for slope
[s_star, m2_s] = fitSlope(x,z,t,slopeMax,v,dr); 

%% Solve for slopeVelocity
[sv_star, m3_s] = fitSV(x,z,t,slopeMax,velMax,dr); 

%% Plot it UP
figure(1)
clf
plot(v*pi*1e7,z,'--','color',rgb('black'),'linewidth',2,'HandleVisibility','off')
hold on
% plot(v_star2*pi*1e7, z,'*-','color',rgb('rose'),'linewidth',2)
% plot(v_star *pi*1e7, z,'p-','color',rgb('dark lime'),'linewidth',2)
plot(v_star*pi*1e7, z,'s-','color',rgb('coral'),'linewidth',2)
xlabel('Velocity')
ylabel('Depth')
set(gca, 'YDir','reverse')
xlim([-velMax*.1*pi*1e7 velMax*1.1*pi*1e7])
% legend('Linear','FFT','Constrained Direct Fit') 

figure(2)
clf
nfigs  = 6;
di = floor((n)/nfigs);
for i = 1:nfigs
    subplot(nfigs,2,2*i-1)
     if(i == 1)
        title('velocity fits')
    end
    plot(t,real(x(:,i*di)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m2_v(1,i*di) + 2*pi*m2_v(2,i*di)),'-','color',rgb('coral'),'lineWidth',2) 
    title("Layer " +z(i*di))
end
for i = 1:nfigs
    subplot(nfigs,2,2*i)
    if(i == 1)
        title('slope fits')
    end
    plot(t,real(x(:,i*di)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m2_s(1,i*di) + 2*pi*m2_s(2,i*di)),'-','color',rgb('lilac'),'lineWidth',2) 
    title("Layer " +z(i*di))
end

%% Plot slopes
xx = velMax * t;
figure(3)
clf
plot(abs(s),z,'--','color',rgb('black'))
hold on
plot(s_star,z,'s-','color',rgb('lilac'),'lineWidth',2)
xlabel('Slope of layer')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')

%% Plot SV product
xx = velMax * t;
figure(4)
clf
plot(abs(s.*v),z,'--','color',rgb('black'))
hold on
plot(sv_star,z,'s-','color',rgb('dark rose'),'lineWidth',2)
xlabel('Slope velocity product')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')