function turingpdc2

close all;

dx = 0.02;
dy = 0.02;
dt = 1e-4;
T = 10000;
D = 0.122;
delta = 2;

alpha = 0.398;
beta = -0.4;
gamma = -alpha;
r1 = 3.5;
r2 = 0;

x = 0:dx:1;
y = 0:dy:1;

[x,y] = meshgrid(x,y);

activatorConc = zeros(length(x));
for k=1:5
    activatorConc(k*10-1:k*10,:) = 10;
end   


inhibitorConc = zeros(length(x));
for k=1:5
    inhibitorConc(:,k*10-1:k*10) = 10;
end  

function p = pindex(i)
    if i == length(activatorConc)
        p = i;
    elseif i == 0
        p = length(activatorConc);
    elseif i == length(activatorConc)+1
        p = 1;
    else 
        p = i;
    end
end

function c2 = diffusionStep(c,D)
    c2 = zeros(length(c));
    
    k1 = D*dt/dx^2;
    k2 = D*dt/dy^2;
    
    for i=1:length(c)
        for j = 1:length(c)
            if i == 1 
                c2(i,pindex(j)) = k1*(c(2,j) - 2*c(i,j) + ...
                    c(i+1,j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                    c(i,pindex(j+1))) + c(i,j);
            elseif i == length(c)
                c2(i,pindex(j)) = k1*(c(i-1,j) - 2*c(i,j) + ...
                    c(i-1,j)) + k2*(c(i,pindex(j-1))...
                    - 2*c(i,j) + c(i,pindex(j+1))) + c(i,j);
            else
                c2(i,pindex(j)) = k1*(c(i-1,j) - 2*c(i,j) + ...
                    c(i+1,j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                    c(i,pindex(j+1))) + c(i,j);
            end
        end
    end     
end

function c2 = reactionStep1(u,v)
    c2 = zeros(length(u));
    
    for i=1:length(u)
        for j = 1:length(u)
            c2(i,j) = dt*(alpha*u(i,j)*(1 - r1*v(i,j)^2) + ...
                v(i,j)*(1 - r2*u(i,j)));
        end
    end
end

function c2 = reactionStep2(u,v)
    c2 = zeros(length(u));
    
    for i=1:length(u)
        for j = 1:length(u)
            c2(i,j) = dt*(beta*v(i,j)*(1 + (alpha*r1/beta)*u(i,j)*v(i,j))...
                + u(i,j)*(gamma + r2*v(i,j)));
        end
    end
end

subplot(1,2,1),
h1 = surf(x,y,activatorConc,'edgecolor','none');
set(gcf,'renderer','painters');
% colormap gray;
caxis([0,5]);
xlim([0,1]);
ylim([0,1]);
view(2);
title('Activator Concentration');
hold all;

subplot(1,2,2),
h2 = surf(x,y,inhibitorConc,'edgecolor','none');
set(gcf,'renderer','painters');
% colormap gray;
caxis([0,5]);
xlim([0,1]);
ylim([0,1]);
view(2);
title('Inhibitor Concentration');
hold all;




for t=1:T
    
    a1 = activatorConc;
    i1 = inhibitorConc;
    
    activatorConc = diffusionStep(activatorConc, D*delta) + reactionStep1(a1,i1);
    
    inhibitorConc = diffusionStep(inhibitorConc, delta) + ...
        reactionStep2(a1,i1);
    
    pause(1/60);
    
    subplot(1,2,1, 'replace'),
    h1 = surf(x,y,activatorConc,'edgecolor','none');
    set(gcf,'renderer','painters');
    view(2);
    title('Activator Concentration');
    
    subplot(1,2,2,'replace'),
    h2 = surf(x,y,inhibitorConc,'edgecolor','none');
    set(gcf,'renderer','painters');
    view(2);
    title('Inhibitor Concentration');
    
    hold all;
end
end