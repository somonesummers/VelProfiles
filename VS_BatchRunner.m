clc; clear;
nl = 500;
r_vz = zeros(1,nl);
r_vx = zeros(1,nl);
mm = zeros(2,nl);
ss = zeros(500,nl);
parfor j = 1:nl
%     clearvars -except j nl r_vz r_vxs
    
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

    %% Solve for slopeVelocity
    [sv_star, m3_s, F1] = fitSV(x,z,t,slopeMax,velMax,dr,v_zb); 

    for i = 1:21
        [sv_star_2, m4_s, F2] = fitSV_2(x,z,t,slopeMax,velMax,dr,movmean(m3_s(1,:),floor(n/5))); 
%         disp("Loop " + j + "." + i + ": res " + abs((sum(F2) - sum(F1))/sum(F1)));
        if(abs((sum(abs(F2-F1)))/sum(F2)) < 1e-2 || i == 21)
            disp("Broke on loop " + i + " with res " + abs((sum(F2) - sum(F1))/sum(F1)));
            break;
        end
        F1 = F2;
        m3_s = m4_s;
    end
    
    %% Fit full data set
    slopeSmooth = n/10;
    s_smooth = movmean(s,slopeSmooth);
    G    = [z' (z.^4 .* (s_smooth))'];

    fit = @(b,zz,ss)  abs(b(1)*zz + b(2) * zz.^4 .* ss);    % Function to fit
    fcn = @(b) sum((fit(b,z,s_smooth) - sv_star_2).^2);
    OPTIONS = optimset('Display','none','TolX',1e-12);
    m_fit = fminsearchbnd(fcn, [1e-8 -1e-8], [],[], OPTIONS);
    if(m_fit(1) < 0)
        m_fit = m_fit*(-1);
    end

    %% Echo Free processing
%     ef_i = floor(n*4/6);
%     G_ef = G(1:ef_i,:);
% 
%     fit_ef = @(b,zz,ss)  abs(b(1)*zz + b(2) * zz.^4 .* ss);    % Function to fit
%     fcn_ef = @(b) sum((fit_ef(b,z(1:ef_i),s_smooth(1:ef_i)) ...
%                 - sv_star_2(1:ef_i)).^2);
%     OPTIONS = optimset('Display','none','TolX',1e-12);
%     m_fit_ef = fminsearchbnd(fcn_ef, [1e-8 -1e-8], [],[], OPTIONS);
%     if(m_fit_ef(1) < 0)
%         m_fit_ef = m_fit_ef*(-1);
%     end
%% Outputs

r_vx(j)  = norm(v_clean - m_fit(2).*z.^4)/n;
r_vz(j)  = norm(v_z - m_fit(1)*z)/n;
mm(:,j)  = m_fit;
ss(:,j)   = s;
end
save dataTrials2

%% 
clear
load dataTrials500_0.5.mat
r_vx_2 = r_vx(r_vx ~=0);
r_vz_2 = r_vz(r_vz ~=0);
figure
subplot(211)

histogram(r_vx_2,nl/5)
set(gca, 'YScale', 'log')
hold on
scatter(mean(r_vx_2),1,'b*')
title("Full")
scatter(std(r_vx_2),1,'r*')
subplot(212)
histogram(r_vz_2,nl/5)
hold on
scatter(mean(r_vz_2),1,'b*')
scatter(std(r_vz_2),1,'r*')

r_vx_2_ef = r_vx_ef(r_vx_ef ~=0);
r_vz_2_ef = r_vz_ef(r_vz_ef ~=0);
figure
subplot(211)
histogram(r_vx_2_ef,nl/5)
set(gca, 'YScale', 'log')
hold on
scatter(mean(r_vx_2_ef),1,'b*')
scatter(std(r_vx_2_ef),1,'r*')
title("ef")
subplot(212)
histogram(r_vz_2_ef,nl/5)
hold on
scatter(mean(r_vz_2_ef),1,'b*')
scatter(std(r_vz_2_ef),1,'r*')

disp("Full: Vz: " + mean(r_vz_2) + " Vx: " + mean(r_vx_2));
disp("Ef: Vz: " + mean(r_vz_2_ef) + " Vx: " + mean(r_vx_2_ef));