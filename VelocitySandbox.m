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

if( n < 25) %Plot wiggles for small n only
    figure(1)
    clf
    wiggle(real(x'),'hk')
end

%% Process Data
tic
X = fft(x,[],1);
f = (1/dt)/nt*(-(nt-1)/2:(nt-1)/2);
X = circshift(X,(nt+1)/2);
% figure(2)
% clf
% plot(f,abs(X(:,1)/nt),'*-')
% hold on
% plot(f,abs(X(:,11)/nt),'*-')
% plot(f,abs(X(:,21)/nt),'*-')
% legend
% xlim([0 1e-6])


%% Extract freq
[M,I] = max(X,[],1);
f_max = f(I);
v_star = f_max*dr./s;
disp("FFT takes " + toc + " seconds") 
%% Fit with line
tic
v_star2 = zeros(size(v_star));
m1 = zeros(2,size(v_star,2));
for i = 1:11 %size(x,2)
    yy = real(x(:,i));
    xx = t;
    G = [xx ones(size(xx))];
    m = G\yy;
    y_star = G*m;
    v_star2(i) = abs(m(1))*dr./s(i);
    m1(:,i) = m;
end
disp("Linear Fit takes " + toc + " seconds") 
%% Fit with sine
tic
[v_star3, m2] = fitVelocity(x,z,t,s,velMax,dr); 
disp("Direct Fit takes " + toc + " seconds") 

%% Plot vels
figure(1)
clf
plot(v*pi*1e7,z,'--','color',rgb('black'),'linewidth',2,'HandleVisibility','off')
hold on
% plot(v_star2*pi*1e7, z,'*-','color',rgb('rose'),'linewidth',2)
% plot(v_star *pi*1e7, z,'p-','color',rgb('dark lime'),'linewidth',2)
plot(v_star3*pi*1e7, z,'s-','color',rgb('blue'),'linewidth',2)
xlabel('Velocity')
ylabel('Depth')
xlim([-velMax*.1*pi*1e7 velMax*1.1*pi*1e7])
% legend('Linear','FFT','Constrained Direct Fit') 

figure(2)
clf
nfigs  = 6;
di = floor((n)/nfigs);
for i = 1:nfigs
    subplot(nfigs,2,i)
    plot(t,real(x(:,i*di)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,t*m1(1,i*di) + m1(2,i*di),'-','color',rgb('rose'),'lineWidth',2)
    plot(t,sin(2*pi*t*f_max(i*di) + 1*pi*m2(2,i*di)),'color',rgb('dark lime'),'lineWidth',2) %borrow phase from fit as FFT doesn't save that
    plot(t,sin(2*pi*t*m2(1,i*di) + 2*pi*m2(2,i*di)),'-','color',rgb('blue'),'lineWidth',2)
    legend('Linear','FFT','Constrained Direct Fit') 
    title("Layer " +z(i*di))
end
% for i = 1:nfigs
%     subplot(nfigs,2,i)
%     plot(t,real(x(:,i*di)),'.','color',rgb('light gray'),'HandleVisibility','off')
%     hold on
%     plot(t,real(x_clean(:,i*di)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
%     plot(t,t*m1(1,i*di) + m1(2,i*di),'-','color',rgb('rose'),'lineWidth',2)
%     plot(t,sin(2*pi*t*f_max(i*di) + 1*pi*m2(2,i*di)),'color',rgb('dark lime'),'lineWidth',2) %borrow phase from fit as FFT doesn't save that
%     plot(t,sin(2*pi*t*m2(1,i*di) + 2*pi*m2(2,i*di)),'-','color',rgb('blue'),'lineWidth',2)
%     legend('Linear','FFT','Constrained Direct Fit') 
%     title("Layer " +z(i*di))
% end

%% Plot slopes
xx = velMax * t;
figure(3)
clf
plot(s,z,'color',rgb('black'))
