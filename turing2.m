%dirichlet

close all;

dx = 0.1;
dy = 0.1;
dt = 1e-2;
T = 10000;
D = 0.02;

x = 0:dx:1;
y = 0:dy:1;

[x,y] = meshgrid(x,y);

c = -x.^2 -y.^2 + 2;

h1 = surf(x,y,c);
caxis([0,2]);
xlim([0,1]);
ylim([0,1]);
% view(2);
hold all;

c2 = zeros(length(c)-2);

k1 = D*dt/dx^2;
k2 = D*dt/dy^2;

for t=1:T
    for i=2:length(c)-1
        for j = 2:length(c)-1
            c2(i-1,j-1) = k1*(c(i-1,j) - 2*c(i,j) + c(i+1,j)) + ...
                k2*(c(i,j-1) - 2*c(i,j) + c(i,j+1)) + c(i,j);
        end
    end
    c(2:length(c)-1,2:length(c)-1) = c2;
    pause(1/60);
    delete(h1);
    h1 = surf(x,y,c);
    %view(2);   
    hold all;
end