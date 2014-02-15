close all;

dx = 0.1;
dt = 1e-4;
T = 1000;
r = 1/10;

x = 0:dx:1;

c = -(x-0.5).^2 + 0.25;

c2 = zeros(1, length(c));

h1 = plot(x,c,'b');
ylim([0,0.25]);
hold all;

for t=1:T
    for i=2:length(c)-1
        c2(i) = c(i) + r*(c(i-1) - 2*c(i) + c(i+1));
    end
    c = c2;
    pause(1/60);
    delete(h1);
    h1 = plot(x,c,'b');
    hold all;
end