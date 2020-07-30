Upoint=0.5; 

Dx=x(2:end)-x(1:end-1); x=x(1:end-1); 

%% find point
[~,bu]=min(abs(Ubin-Upoint));
bu=8; 
%% extract data 
indU=binU==bu;
DxU=Dx(indU); xU=x(indU);
D1xU=(1/dt)*DxU;          % calculate drifts *** where x binninng would be added  
D2xU=(1/dt)*DxU.*DxU;     % calculate diffusions 
D2xUlog=log(D2xU);


[~,indS]=min(abs(xfit-min(xU)));[~,indE]=min(abs(xfit-max(xU))); %% find parts of xfit vector corresponding to xU range
    dx=xfit(2)-xfit(1);

xfitu=xfitsave{bu}; nu=length(xfitu);
d1xfit=DriftFit(1:nu,bu); 
    %% find zeros
     [x0,stb]=ZeroGPR(xfitu,d1xfit,err,count_max,Nx); x0=sort(x0);  %find zero's
     n0=length(x0);y0=zeros(n0,1); 
     


%% refit as models aren't saved --- solutions wont exactly match, only affects uncertainty
          itts=1000;     
        %itts=round(length(xU)/20);
         
     sigN=std(xU); sigF=KDrift(1,bu); sigL=KDrift(2,bu); 
     [x0Dist,x0min,x0max,nNaN,Lstar] = x0uncertaintyR(xU,xfitu,d1xfit,x0,sigN,sigL,sigF,itts,err,count_max,prc,reg); 
   
     
     Nx1=length(xfitu);
   D1xnew1 = d1xfit + Lstar*randn(Nx1,1);  %new realizationn 
   D1xnew2 = d1xfit + Lstar*randn(Nx1,1);  %new realizationn 

figure 
hold on; grid on 
plot(xU,D1xU,'.')
plot(xfitu,d1xfit)
plot(xfitu,D1xnew1)
plot(xfitu,D1xnew2);
z0=zeros(length(x0),1);
scatter(x0,z0,[],nNaN,'filled')
errorbar(x0,z0,x0min-x0,x0max-x0,'horizontal')
legend('data','fit -mean', 'realization 1','realization 2','zero','errors')

   
   