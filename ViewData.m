clear
load('ImageP2.mat')
load('ImageP2.mat','timeInDays')
layers = size(RawImage,1);
[m, bed] = max(mean(abs(RawImage(2*1500:end,:)),2));
bed = 2*1500+bed;
% figure(1)
% clf
% plot(log(abs(RawImage(:,50))))

in = 1:floor(bed/15):bed;
buff = 3;
dwnSample = 1;
x = RawImage(1:dwnSample:bed,:)./mean(sqrt(2)*abs(RawImage(1:dwnSample:bed,:)),2);
x_clean = movmean(x,buff,2);
z = Rcoarse(1:dwnSample:bed);
n = size(z,2);
t = timeInDays*24*3600;
slopeMax = .01; %Max Slope
velMax = 50/(3600*24*365); % [m/s] Max Speed
dr = 284e-3; %wavelength in medium

[sv_star  , m3_s,F1] = fitSV(  x',z,t,slopeMax,velMax,dr);
[sv_star_2, m4_s,F2] = fitSV_2(x',z,t,slopeMax,velMax,dr,movmean(m3_s(1,:),floor(bed/50))); 
[sv_star_3, m5_s,F3] = fitSV_2(x',z,t,slopeMax,velMax,dr,movmean(m4_s(1,:),floor(bed/50))); 

sv_star   = sv_star   * 3.154e7;
sv_star_2 = sv_star_2 * 3.154e7;
sv_star_3 = sv_star_3 * 3.154e7;


%% Plotting
figure(1)
clf
plot(movmean(real(RawImage(in(1),:)),buff)/mean(sqrt(2)*abs(RawImage(in(1),:))))
hold on
for i = 2:length(in)
    plot(movmean(real(RawImage(in(i),:)),buff)/mean(sqrt(2)*abs(RawImage(in(i),:))))
end

figure(2)
clf
nfigs  = 12;
di = floor((bed)/nfigs/dwnSample);
for i = 1:nfigs
    subplot(nfigs,1,i)
    plot(t,real(x(i*di,:)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(i*di,:)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m3_s(1,i*di) + 2*pi*m3_s(2,i*di)),'-','color',rgb('coral'),'lineWidth',2) 
    plot(t,sin(2*pi*t*m4_s(1,i*di) + 2*pi*m4_s(2,i*di)),'-','color',rgb('light lime'),'lineWidth',2) 
    plot(t,sin(2*pi*t*m5_s(1,i*di) + 2*pi*m5_s(2,i*di)),'-','color',rgb('sky blue'),'lineWidth',2) 
    title("Layer at " +(z(i*di)))
%     axis off
    hold off
end

figure(3)
clf
nfigs  = 12;
zi = floor(4000/dwnSample);
di = 1;
for i = 1:nfigs
    subplot(nfigs,1,i)
    plot(t,real(x(zi+(i-1)*di,:)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(zi+(i-1)*di,:)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m3_s(1,zi+(i-1)*di) + 2*pi*m3_s(2,zi+(i-1)*di)),'-','color',rgb('coral'),'lineWidth',2) 
    plot(t,sin(2*pi*t*m4_s(1,zi+(i-1)*di) + 2*pi*m4_s(2,zi+(i-1)*di)),'-','color',rgb('light lime'),'lineWidth',2) 
    plot(t,sin(2*pi*t*m5_s(1,zi+(i-1)*di) + 2*pi*m5_s(2,zi+(i-1)*di)),'-','color',rgb('sky blue'),'lineWidth',2) 
    title("Layer " +z(zi+(i-1)*di))
    hold off
end

%%
figure(4)
clf
plot(sv_star,z,'.','color',rgb('light rose'),'MarkerSize',10,'HandleVisibility','off')
setFontSize(16)
xlim([0 2])
xlabel('Slope velocity product [rad*m/yr]')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
hold on
% savePng('Willans1')
plot(movmean(sv_star,floor(n/20)),z,'--','color',rgb('dark rose'),'lineWidth',4)
legend('1st Fit')
% savePng('Willans2')
plot(sv_star_2,z,'.','color',rgb('light lime'),'MarkerSize',10,'HandleVisibility','off')
plot(movmean(sv_star_2,floor(n/20)),z,'--','color',rgb('dark lime'),'lineWidth',4)
legend('1st Fit','2nd Fit')
% savePng('Willans3')
plot(sv_star_3,z,'.','color',rgb('baby blue'),'MarkerSize',10,'HandleVisibility','off')
plot(movmean(sv_star_3,floor(n/20)),z,'--','color',rgb('blue'),'lineWidth',4)
legend('1st Fit','2nd Fit','3rd fit')
% savePng('Willans4')
hold off

figure(5)
clf
plot(F1,z,'.','color',rgb('light rose'),'MarkerSize',10,'HandleVisibility','off')
setFontSize(16)
xlabel('Misfit')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
hold on
% savePng('Willans1')
plot(movmean(F1,floor(n/20)),z,'--','color',rgb('dark rose'),'lineWidth',4)
legend('1st Fit')
% savePng('Willans2')
plot(F2,z,'.','color',rgb('light lime'),'MarkerSize',10,'HandleVisibility','off')
plot(movmean(F2,floor(n/20)),z,'--','color',rgb('dark lime'),'lineWidth',4)
legend('1st Fit','2nd Fit')
% savePng('Willans3')
plot(F3,z,'.','color',rgb('baby blue'),'MarkerSize',10,'HandleVisibility','off')
plot(movmean(F3,floor(n/20)),z,'--','color',rgb('blue'),'lineWidth',4)
legend('1st Fit','2nd Fit','3rd fit')
% savePng('Willans4')
hold off

figure(6)
clf
plot(sv_star(F1<movmean(F1,floor(n/20))),z(F1<movmean(F1,floor(n/20))),'.','color',rgb('light rose'),'MarkerSize',10,'HandleVisibility','on')
setFontSize(16)
xlim([0 2])
xlabel('Slope velocity product [rad*m/yr]')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')
hold on
% savePng('Willans1')
legend('1st Fit')
% savePng('Willans2')
plot(sv_star_2(F2<movmean(F2,floor(n/20))),z(F2<movmean(F2,floor(n/20))),'.','color',rgb('light lime'),'MarkerSize',10,'HandleVisibility','on')
legend('1st Fit','2nd Fit')
% savePng('Willans3')
plot(sv_star_3(F3<movmean(F3,floor(n/20))),z(F3<movmean(F3,floor(n/20))),'.','color',rgb('baby blue'),'MarkerSize',10,'HandleVisibility','on')
ln = length(sv_star_3(F3<movmean(F3,floor(n/20))));
plot(movmean(sv_star_3(F3<movmean(F3,floor(n/20))),floor(ln/20)),z(F3<movmean(F3,floor(n/20))),'--',...
    'color',rgb('blue'),'lineWidth',4,'HandleVisibility','off')
legend('1st Fit','2nd Fit','3rd fit')
title('Moving average with worst 50% fits excluded')
savePng('Willans7')
hold off