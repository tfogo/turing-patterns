function turingAllPdc

close all;


dx = 0.04;
dy = 0.04;
dt = 1e-4;
T = 10000;
D = 0.516;
delta = 2;

alpha = .899;
beta = -.91;
gamma = -alpha;
r1 = .35;
r2 = 0;

x = 0:dx:10;
y = 0:dy:10;

[x,y] = meshgrid(x,y);

activatorConc = rand(length(x));

      
inhibitorConc = rand(length(x));

    
 
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
            c2(pindex(i),pindex(j)) = k1*(c(pindex(i-1),j) - 2*c(i,j) + ...
                c(pindex(i+1),j)) + k2*(c(i,pindex(j-1)) - 2*c(i,j) + ...
                c(i,pindex(j+1))) + c(i,j);
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
%caxis([0,5]);
xlim([0,1]);
ylim([0,1]);
view(2);
title('Activator Concentration');
hold all;

subplot(1,2,2),
h2 = surf(x,y,inhibitorConc,'edgecolor','none');
set(gcf,'renderer','painters');
% colormap gray;
% caxis([0,5]);
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
    
    pause(1/100);
    
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