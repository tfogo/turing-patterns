close all;

dx = 0.1;
dt = 1e-4;
T = 1000;

x = 0:dx:1;
x = [(0-dx),x,(1+dx)];

a = -0.5:dx:0.5;
c = -a.^2 + 0.25;
c = [0,c,0];

c2 = zeros(1, length(c));

h1 = plot(x,c,'b');
ylim([0,0.25]);
xlim([0,1]);
hold all;

for t=1:T
    for i=3:length(c)-2
        c2(i) = c(i) + (1/10)*(c(i-1) - 2*c(i) + c(i+1));
    end
    c = c2;
    pause(1/60);
    delete(h1);
    h1 = plot(x,c,'b');
    hold all;
end


