function turingpdc

close all;

dx = 0.02;
dy = 0.02;
dt = 1e-4;
T = 10000;
D = 0.2;

x = 0:dx:1;
y = 0:dy:1;

[x,y] = meshgrid(x,y);

c = zeros(length(x));
c(6:8,:) = 10;
c(:,6:8) = 10;

function p = pindex(i)
if i == length(c)
    p = i;
elseif i == 0
    p = length(c);
elseif i == length(c)+1
    p = 1;
else 
    p = i;
end
end

h1 = surf(x,y,c,'edgecolor','none');
colormap gray;
caxis([0,5]);
xlim([0,1]);
ylim([0,1]);
%view(2);
hold all;

c2 = zeros(length(c));

k1 = D*dt/dx^2;
k2 = D*dt/dy^2;

for t=1:T
    for i=1:length(c)
        for j = 1:length(c)
            if i == 1 
                c2(i,pindex(j)) = k1*(c(2,j) - 4*dy*c(1,j) - 2*c(i,j) + ...
                    c(i+1,j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                    c(i,pindex(j+1))) + c(i,j);
            elseif i == length(c)
                c2(i,pindex(j)) = k1*(c(i-1,j) - 2*c(i,j) + ...
                    c(i-1,j) - 4*dy*c(i,j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                    c(i,pindex(j+1))) + c(i,j);
            else
                c2(i,pindex(j)) = k1*(c(i-1,j) - 2*c(i,j) + ...
                    c(i+1,j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                    c(i,pindex(j+1))) + c(i,j);
            end
        end
    end
    c = c2;
    pause(1/60);
    delete(h1);
    h1 = surf(x,y,c,'edgecolor','none');
    %view(2);
    hold all;
end
end