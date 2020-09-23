clear
load('ImageP2.mat')
layers = size(RawImage,1);
[m, bed] = max(mean(abs(RawImage(3500:end,:)),2));
bed = 3500+bed;
figure
plot(log(abs(RawImage(:,50))))

in = 1:floor(bed/15):bed;
buff = 3;

x = RawImage./mean(sqrt(2)*abs(RawImage),2);
x_clean = movmean(x,buff,2);

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
    plot(real(x(i*di,:)),'.','color',rgb('light gray'),'HandleVisibility','off')
    hold on
    plot(real(x_clean(i*di,:)),'--','color',rgb('black'),'lineWidth',2,'HandleVisibility','off')
%     plot(t,sin(2*pi*t*m2_v(1,i*di) + 2*pi*m2_v(2,i*di)),'-','color',rgb('coral'),'lineWidth',2) 
    title("Layer " +(i*di))
end