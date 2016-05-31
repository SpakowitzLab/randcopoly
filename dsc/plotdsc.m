% %%
clear;
close all
t = csvread('test.csv');

% find numerical second derivative
x = t(:,2);
y = t(:,3);
nt = length(t);

dy = zeros(nt-1,1);
d2y = zeros(nt-2,1);

% calculate first derivatives (forward differences)
for ii = 2:length(x)
    h = x(ii)-x(ii-1);
    dy(ii-1) = (y(ii)-y(ii-1))/h;
end

% calculate second derivatives (central differences)
for ii = 2:length(x)-1
    h1 = x(ii)-x(ii-1);
    h2 = x(ii+1)-x(ii);
    d2y(ii-1) = 2* ( (y(ii+1)-y(ii))/h2 - (y(ii)-y(ii-1))/h1 ) / (h1+h2);
end

% filter second derivatives
% d2y(abs(dy)>1)=0;
xp = x(2:nt);
xp2 = x(2:nt-1);

figure;hold
for ii = 1:10:nt
    col = (ii-1)/(nt-1);
    color = [col 0 1-col];
    plot(t(ii,2),t(ii,3),'.','color',color);
end

figure;hold
for ii = 1:10:nt-1
    col = (ii-1)/(nt-1);
    color = [col 0 1-col];
    plot(xp(ii),dy(ii),'.','color',color);
end
ylim([-.02,.03])
plot([-22.92,-22.92],[-1,1],'k-','linewidth',2)

figure;hold
for ii = 1:10:nt-2
    col = (ii-1)/(nt-1);
    color = [col 0 1-col];
    plot(xp2(ii),d2y(ii),'.','color',color);
end
ylim([-.2,.3])

% figure;plot(xp,dy);
% figure;plot(xp2,abs(d2y));