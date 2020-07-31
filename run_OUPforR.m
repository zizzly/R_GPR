

close all; clear all; 

description='OUP'; 

x0=0; 
c=1; 
g=0.5; 
dt=0.1; 
T=1000; 
Nx=1000; 

[x] = OrsnUhl1Dw(x0,c,g,T,dt);

xfit=linspace(min(x),max(x),Nx)'; dx=xfit(2)-xfit(1);

t=0:dt:T; t=t(:); 

save(strcat('oupForR',description)); 
 
writematrix(x,strcat('oupX',description,'.txt'))
writematrix(xfit,strcat('oupXfit',description,'.txt'))
figure; plot(t,x); xlabel('t'); ylabel('x') 