% Fit rough sin data
clear
x = [0:1e4:5e7]';
f = 1.6516e-08;
phi = rand;
y_prime = sin(2*pi*x*f+2*pi*phi);
y = y_prime + .01*randn(size(x));
% 
% 
% 
yu = max(y);
yl = min(y);
yr = (yu-yl)/2;                               % Range of ?y?
yz = movmean(y-yu+(yr),floor(length(y)/50));
zx = x(yz .* circshift(yz,[1 0]) <= 0);     % Find zero-crossings
if length(zx) > 2
    per = 2*mean(diff(zx));                     % Estimate period
else
    per = 2*max(x)-min(x);
end
ym = mean(y);                               % Estimate offset
fit = @(b,x)  1.*(sin(2*pi*x.*b(1) + 2*pi*b(2))) + 0;    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);
OPTIONS = optimset('Display','none');% Least-Squares cost function
s = fminsearch(fcn, [1/per;  0],OPTIONS);                       % Minimise Least-Squares
xp = linspace(min(x),max(x));
disp([[5;phi] s])
% y = Gm

m = [0 0];
for i = 1:10
    m_old = m;
    G = [x-m(1)^2*(x).^3/6+m(1)^4*x.^5/120 ones(size(x))];
    m = G\y;
    y_star = G*m;
    if (norm(m-m_old) < 1e-4)
        break
    end
end


figure(1)
clf
plot(x,y,'.')
hold on
plot(x,yz,'-.')
plot(x,y_prime,'--k')
% plot(x,y_star,'r')
plot(xp, fit(s,xp),'b')
