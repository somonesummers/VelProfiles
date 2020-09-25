clear
load('Image2.mat')
load('ImageP2.mat','timeInDays')
layers = size(RawImage,1);
[m, bed] = max(mean(abs(RawImage(2*1500:end,:)),2));
bed = 2*1500+bed;
figure
plot(log(abs(RawImage(:,50))))

in = 1:floor(bed/15):bed;
buff = 3;

x = RawImage(1:bed,:)./mean(sqrt(2)*abs(RawImage(1:bed,:)),2);
x_clean = movmean(x,buff,2);
z = Rcoarse(1:bed);
t = timeInDays*24*3600;
slopeMax = .01; %Max Slope
velMax = 50/(3600*24*365); % [m/s] Max Speed
dr = Rcoarse(2)-Rcoarse(1); %TODO, correct this value

[sv_star, m3_s] = fitSV(x',z,t,slopeMax,velMax,dr);
%%
[sv_star_2, m4_s] = fitSV_2(x',z,t,slopeMax,velMax,dr,movmean(m3_s(2,:),250)); 
% Plotting
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
di = floor((bed)/nfigs);
for i = 1:nfigs
    subplot(nfigs,1,i)
     if(i == 1)
        title('velocity fits')
    end
    plot(t,real(x(i*di,:)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(t,real(x_clean(i*di,:)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
    plot(t,sin(2*pi*t*m3_s(1,i*di) + 2*pi*m3_s(2,i*di)),'-','color',rgb('coral'),'lineWidth',2) 
    plot(t,sin(2*pi*t*m4_s(1,i*di) + 2*pi*m4_s(2,i*di)),'-','color',rgb('light lime'),'lineWidth',2) 
    title("Layer " +(i*di))
end

xx = velMax * t;
figure(4)
clf
plot(sv_star,z,'s','color',rgb('light rose'),'lineWidth',2)
hold on
plot(movmean(sv_star,50),z,'--','color',rgb('dark rose'),'lineWidth',2)
plot(sv_star_2,z,'s','color',rgb('light lime'),'lineWidth',2)
plot(movmean(sv_star_2,50),z,'--','color',rgb('dark lime'),'lineWidth',2)

xlabel('Slope velocity product')
ylabel('Depth')
set(gca, 'YDir','reverse')
set(gca, 'YDir','reverse')