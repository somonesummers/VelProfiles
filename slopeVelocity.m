% Main file to run for IGARSS Paper
% Creates a synthetic data set of layer motion over time based on vertical
% velocity and layer motion due to horizontal motion of sloped layers. 

clear
rng(125) %125 for fig in paper

%% Input parameters
depth = 1e3; % Depth of Domain[m]
v_zb = 1.5/(3600*24*365); % [m/s] vertical speed at bed
velMax = 50/(3600*24*365); % [m/s] Max Speed
slopeMax = 1e-1; % [ ] Max slope, assumed to be at bed 
dr = .3; %[m] wavelength of system
n = 500; %[ ] number of modeling depth bins
dz = depth/(n - 1); %[m] size of model depth bins
dt = 24*3600/2; %[s] sample period
nt = 365*1; %must be odd
timeMax = dt*(nt-1);
%% Make synthetic data
% linear layers of known slope. n^4 velocity profile
z = 0:dz:depth;
v_z = z/depth*v_zb;
% s = .05*ones(size(z));
% s = slopeMax/100:slopeMax/(n):slopeMax; %Linear change in slope
% s = slopeMax:-slopeMax/(n):(slopeMax/1000); %Linear change in slope
% s = slopeMax*(sin(0:2*pi/(n-1):2*pi)) + slopeMax/100*randn(size(z)); %Sine changes in slope
s = movmean(cumsum(randn(size(z))*10e-2/sqrt(n)),floor(n/100)); %Random walk of slope;
s_clean = s;
% s = s + slopeMax/10*randn(size(z));
v_clean = velMax*(z/depth).^(4);
v = v_clean;% + velMax/10*randn(size(z));
t = (0:dt:timeMax)';

phi = rand(size(v));
dPhi = 2*pi * (s .* v  + v_z)./ dr .* t + 2*pi*phi;
x_clean = cos(dPhi) + 1i*sin(dPhi);
x = cos(dPhi) + 1i*sin(dPhi) + .5*(randn(size(dPhi)) + 1i*randn(size(dPhi)));

%% Solve for velocity
% [v_star, m2_v] = fitVelocity(x,z,t,s,velMax + (v_zb),dr);
% [v_star2, m2_v2] = fitVelocity_2(x,z,t,s,velMax,dr,m2_v(1,:));
% 
%% Solve for slope
% [s_star, m2_s] = fitSlope(x,z,t,slopeMax,v,dr); 

%% Solve for slopeVelocity
[sv_star, m3_s, F1] = fitSV(x,z,t,slopeMax,velMax,dr,v_zb); 

for i = 1:21
    [sv_star_2, m4_s, F2] = fitSV_2(x,z,t,slopeMax,velMax,dr,movmean(m3_s(1,:),floor(n/5))); 
    disp("Loop " + i + ": res " + abs((sum(F2) - sum(F1))/sum(F1)));
    if(abs((sum(abs(F2-F1)))/sum(F2)) < 1e-2 || i == 10)
        disp("Broke on loop " + i + " with res " + abs((sum(F2) - sum(F1))/sum(F1)));
        break;
    end
    F1 = F2;
    m3_s = m4_s;
end

%% Plot it UP
figure(1)
clf
plot(s,z)
xlabel('Slope')
ylabel('Depth')
set(gca, 'YDir','reverse')

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
    ylabel('Re(x)');
   if(i == nfigs)
        xlabel('Time [s]')  
    end
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
    ylabel('Re(x)');
    if(i == nfigs)
        xlabel('Time [s]')  
    end
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
% figure(3)
% clf
% plot(abs(s_clean),z,'--','color',rgb('black'),'lineWidth',3)
% xlabel('Slope of layer')
% ylabel('Depth')
% set(gca, 'YDir','reverse')
% setFontSize(16)
% % savePng('figs/slope1')
% hold on
% plot(abs(s),z,'.','color',rgb('gray'),'markersize',15)
% % savePng('figs/slope2')
% plot(s_star,z,'.','color',rgb('lilac'),'markersize',15)
% plot(movmean(s_star,floor(n/10)),z,'--','color',rgb('lilac'),'lineWidth',3)
% % savePng('figs/slope3')

%% Plot Velocities

% G = [z' z.^4'];
% m_v = G\(v_star');
% 
% m_sv = G\((sv_star./s)');
% 
% 
% figure(4)
% clf
% plot(v,z,'.','color',rgb('gray'),'markerSize',15)
% xlabel('Velocity')
% xlim([0 velMax])
% ylabel('Depth')
% set(gca, 'YDir','reverse')
% setFontSize(16)
% % savePng('figs/vel1')
% hold on
% plot(v_clean,z,'--','color',rgb('black'),'linewidth',4);
% % savePng('figs/vel2')
% % plot(v_star,z,'.','color',rgb('soft green'),'markerSize',15)
% % plot(movmean(v_star,floor(n/10)),z,'--','color',rgb('dark green'),'lineWidth',4)
% 
% plot(sv_star_2./s,z,'.','color',rgb('baby blue'),'MarkerSize',10)
% % plot(movmean(sv_star_2,floor(n/10))./s,z,'--','color',rgb('blue'),'lineWidth',4)
% 
% % plot(G*m_v,z,'--','color',rgb('green'),'lineWidth',4)
% % plot(m_v(1)*G(:,1),z,':','color',rgb('lime'),'lineWidth',4)
% % plot(m_v(2)*G(:,2),z,':','color',rgb('dark lime'),'lineWidth',4)
% plot(G*m_sv,z,'--','color',rgb('turquoise'),'lineWidth',4)
% plot(m_sv(1)*G(:,1),z,':','color',rgb('light aqua'),'lineWidth',4)
% plot(m_sv(2)*G(:,2),z,':','color',rgb('aqua'),'lineWidth',4)
% legend
% savePng('figs/vel3')



%% Plot SV product

%Fit full data set
slopeSmooth = n/10;
s_smooth = movmean(s,slopeSmooth);
G    = [z' (z.^4 .* (s_smooth))'];
% m_sv = G\((sv_star_2)');

fit = @(b,zz,ss)  abs(b(1)*zz + b(2) * zz.^4 .* ss);    % Function to fit
fcn = @(b) sum((fit(b,z,s_smooth) - sv_star_2).^2);
OPTIONS = optimset('Display','none','TolX',1e-12);
m_fit = fminsearchbnd(fcn, [1e-8 -1e-8], [],[], OPTIONS);
if(m_fit(1) < 0)
    m_fit = m_fit*(-1);
end

% Echo Free processing
ef_i = floor(n*4/6);
G_ef = G(1:ef_i,:);

fit_ef = @(b,zz,ss)  abs(b(1)*zz + b(2) * zz.^4 .* ss);    % Function to fit
fcn_ef = @(b) sum((fit_ef(b,z(1:ef_i),s_smooth(1:ef_i)) ...
            - sv_star_2(1:ef_i)).^2);
OPTIONS = optimset('Display','none','TolX',1e-12);
m_fit_ef = fminsearchbnd(fcn_ef, [1e-8 -1e-8], [],[], OPTIONS);
if(m_fit_ef(1) < 0)
    m_fit_ef = m_fit_ef*(-1);
end
%% Outputs

%% make m/yr
% m_fit = m_fit * 3.154e7;
% m_fit_ef = m_fit_ef * 3.154e7;


disp("Full profile")
disp("Vz err: (|"  + v_zb + " - " + m_fit(1)*depth + "|) => " + abs(v_zb -  m_fit(1)*depth)/v_zb*100 + "%");
disp("Vz normalized resid: " + norm(v_z - m_fit(1)*z)/n + " => " + norm(v_z - m_fit(1)*z)/n*3.154e7 + "m/yr")
disp("Vx err: (|"  + velMax + " - " + m_fit(2)*depth.^4 + "|) => " + abs(velMax -  m_fit(2)*depth.^4)/velMax*100 + "%");
disp("Vx normalized resid: " + norm(v_clean - m_fit(2).*z.^4)/n + " => " + norm(v_clean - m_fit(2).*z.^4)/n*3.154e7 + "m/yr" )
disp("Echo Free Zone exluded")
disp("Vz err: (|"  + v_zb + " - " + m_fit_ef(1)*depth + "|) => " + abs(v_zb -  m_fit_ef(1)*depth)/v_zb*100 + "%");
disp("Vx err: (|"  + velMax + " - " + m_fit_ef(2)*depth.^4 + "|) => " + abs(velMax -  m_fit_ef(2)*depth.^4)/velMax*100 + "%");
disp("Vx normalized resid: " + norm(v_clean - m_fit_ef(2)*z.^4)/n )

resid_vxs = norm(v_clean - m_fit(2).*z.^4)/n;
resid_vz = norm(v_z - m_fit(1)*z)/n;
%%
figure(4)
clf
subplot(131)
plot(v_clean,z,'color',rgb('red'))
hold on
plot(m_fit(2)*z.^4,z,'color',rgb('blue'),'lineWidth',5)
plot(m_fit_ef(2)*z.^4,z,'color',rgb('light blue'),'lineWidth',5)
xlabel('Velocity [m/s]')
ylabel('Depth')
set(gca, 'YDir','reverse')
legend('Vel_x','Vel_x fit','Vel_x fit Echo Free')
subplot(132)
plot(s,z,'color',rgb('red'))
hold on
plot(s_smooth,z,'color',rgb('blue'),'lineWidth',5)
xlabel('Slope []')
ylabel('Depth')
set(gca, 'YDir','reverse')
legend('Random Slope','Smoothed Slope')
subplot(133)
plot(v.*s,z,'color',rgb('red'))
hold on
plot(m_fit(2)*z.^4.*s_smooth,z,'color',rgb('blue'),'lineWidth',5)
xlabel('Vel*S [m/s]')
ylabel('Depth')
set(gca, 'YDir','reverse')
legend('V_x*s','V_x*s fit')
%%
figure(5)
clf
plot(sv_star,z,'o','color',rgb('light rose'),'MarkerSize',10,'DisplayName', 'Initial Bin Fits')
% plot(movmean(sv_star,floor(n/10)),z,'--','color',rgb('dark rose'),'lineWidth',5,...
%     'DisplayName', 'Data ')
xlabel('Velocity [m/s]')
ylabel('Depth')
set(gca, 'YDir','reverse')
legend
% xlim ([0 max(abs(s_clean.*v_clean+v_z))])
setFontSize(16)
% savePng('figs/sv1')
hold on
% savePng('figs/sv2')
% savePng('figs/sv3')
plot(sv_star_2,z,'.','color',rgb('baby blue'),'MarkerSize',25,'DisplayName', 'Final Bin Fits')
% plot(movmean(sv_star_2,floor(n/10)),z,'--','color',rgb('blue'),'lineWidth',5)

plot(abs(s_clean.*v_clean+v_z),z,'-','color',rgb('dark gray'),'lineWidth',5,...
    'DisplayName', 'Total Signal (V_r)')


plot(v_z,z,'-','color',rgb('light gray'),'lineWidth',4,'DisplayName', 'V_z')
plot(v_clean.*s,z,'-','color',rgb('gray'),'lineWidth',4, 'DisplayName', 'V_x * S')

% plot(G*m_sv,z,'-.','color',rgb('turquoise'),'lineWidth',4,'DisplayName', 'Total Fit')
% plot(m_sv(1)*G(:,1),z,':','color',rgb('dark mint'),'lineWidth',4,'DisplayName', 'Fit V_z')
% plot(m_sv(2)*G(:,2),z,'--','color',rgb('mint'),'lineWidth',4,'DisplayName', 'Fit V_x * S')

plot(fit(m_fit,z,s_smooth),z,'-.','color',rgb('red'),'lineWidth',4,'DisplayName', 'Total Fit (V_r)')
plot(m_fit(1)*G(:,1),z,':','color',rgb('light red'),'lineWidth',4,'DisplayName', 'Fit V_z')
plot(m_fit(2)*G(:,2),z,'--','color',rgb('scarlet'),'lineWidth',4,'DisplayName', 'Fit V_x * S')

% plot(fit_ef(m_fit_ef,z,s_smooth),z,'-.','color',rgb('turquoise'),'lineWidth',4,'DisplayName', 'Total Fit Echo Free')
% plot(m_fit_ef(1)*G_ef(:,1),z(1:ef_i),':','color',rgb('dark mint'),'lineWidth',4,'DisplayName', 'Fit V_z Echo Free')
% plot(m_fit_ef(2)*G_ef(:,2),z(1:ef_i),'--','color',rgb('mint'),'lineWidth',4,'DisplayName', 'Fit V_x * S Echo Free')

% savePng('figs/svFits')
% saveVect('figs/svFits')

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
