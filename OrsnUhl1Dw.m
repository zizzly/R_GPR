%ER, June 21, 2019

function [xSDE] = OrsnUhl1Dw(x0,c,g,T,dt)

% solves  dx= -cxdt + gdW and plots
% xSDE - solve linear SDE using euler method/numerically
% xINT - solve using integrated/analytic solution
% xclean - noise free solution 
% T - max time
% dt - time interval 

N = round(T/dt);
t= [0:dt:T]'; %time vector

xSDE=zeros(N,1); %initialize

w=randn(N,1)*sqrt(dt); %noise term for each sample, gaussian with Var=dt
W=zeros(N,1);       % initialize - sum for integral of analytic solution    
xSDE(1)=x0;         % set initial condition 

%iterate over each time step 
for i=1:N
    xSDE(i+1)= xSDE(i) + -c*xSDE(i)*dt + g*w(i);  %Euler method
    W(i+1)=W(i)+  exp(c*t(i))*w(i);         %accumulated noise for analytic
end

%xInt=x0*exp(-c*t) + g*exp(-c*t).*W; %analytic solution 

%xClean=x0*exp(-c*t);                %no noise solution 

end 


