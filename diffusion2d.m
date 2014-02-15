close all;

dx = 0.1;
dy = 0.1;
dt = 1e-4;
T = 10000;
D = 0.02;

x = (0-dx):dx:(1+dx);
y = (0-dy):dy:(1+dy);

[x,y] = meshgrid(x,y);

c = -x.^2 -y.^2 + 2;

h1 = surf(x,y,c);
caxis([0,2]);
xlim([0,1]);
ylim([0,1]);
% view(2);
hold all;

c2 = zeros(length(c));

k1 = D*dt/dx^2;
k2 = D*dt/dy^2;

for t=1:T
    for i=3:length(c)-2
        for j = 3:length(c)-2
            c2(i,j) = k1*(c(i-1,j) - 2*c(i,j) + c(i+1,j)) + ...
                k2*(c(i,j-1) - 2*c(i,j) + c(i,j+1)) + c(i,j);
        end
    end
    c = c2;
    pause(1/60);
    delete(h1);
    h1 = surf(x,y,c);
%     view(2);   
    hold all;
end