close all; clear all; 

description='SecondAttempt'; 
%% load 
xsave=readmatrix(strcat('vwXBin',description,'.txt')); xsave=xsave'; 
xfit=readmatrix(strcat('vwXfit',description,'.txt')); dx=xfit(2)-xfit(1); 
U=readmatrix(strcat('vwU',description,'.txt'));
Ubin=readmatrix(strcat('vwU',description,'.txt')); 

%% parameters 
Nu=30; Nx=1000;
dt=1;
err=1e-5; count_max=10; %for zero finding iteration code
prc=97.5; itts=1000; reg=1e-12;

%% initialize 
Uplot=1; x0Plot=1; 
D1xfit=nan(Nx,Nu); D2xfit=nan(Nx,Nu);
 x0MIN=1; x0MAX=1;
 x0Dist=cell(Nu,1); 
 nNaN=1;
 STB=1; 

 for bu=1:Nu

     xU=xsave(:,bu); 
     bad=isnan(xU); 
     xU(bad)=[]; 
      DxU=xU(2:end)-xU(1:end-1); xU=xU(1:end-1);
     p=99; 
     [~,indS]=min(abs(xfit-prctile(xU,100-p)));[~,indE]=min(abs(xfit-prctile(xU,p))); %% find parts of xfit vector corresponding to xU range
        
      %xfitu=xfit(indS:indE); 
      m1=length(xfitu); m2=length(indS:indE); 
     if m1 ~=m2
         indE=indE-1;
     end
 
     %DxU=xU(2:end)-xU(1:end-1); %pull corresponding Dx, xU points at that Uvalue
       
     D1xU=(1/dt)*DxU;          % calculate drifts *** where x binninng would be added  
     D2xU=(1/dt)*DxU.*DxU;     % calculate diffusions 
     D2xUlog=log(D2xU);
   
     model= fitrgp(xU(:), D1xU,'BasisFunction','none','ActiveSetMethod','likelihood','PredictMethod','sr','Holdout',0.2);   %create GPR fit for drift
    
     model2=fitrgp(xU(:),D2xUlog,'BasisFunction','constant','ActiveSetMethod','likelihood','PredictMethod','sr','Holdout',0.2);      
     model=model.Trained{1};
     model2=model2.Trained{1};
     
     fprintf('have fit models for bin %f / %f \n',bu,Nu)
     D1xfit(indS:indE,bu)=predict(model,xfitu); %calculate drift predictions at xfit - defined in parameters section
     sfit=predict(model2,xfitu);   %calculate diffusion predictions
     D2xfit(indS:indE,bu)=exp(sfit);
     %% find zeros
     [x0,stb]=ZeroGPR(xfitu,D1xfit(indS:indE,bu),err,count_max,Nx); x0=sort(x0);  %find zero's
     n0=length(x0);
     
      fprintf('have zeros for bin %f / %f',bu,Nu)
     
     %% uncertainty estimates 
     
        %itts=round(length(xU)/20);
        
     [x0Dist{bu},x0min,x0max,nnn] = x0uncertainty3(xfitu,D1xfit(indS:indE,bu),model,x0,itts,err,count_max,prc);
      
    fprintf('have uncertainty for bin %f / %f \n',bu,Nu)
     
    nNaN=[nnn;nNaN];
     %Test=[Test;test];
     %Ptest=[Ptest;ptest];
     x0MIN=[x0MIN;x0min];
     x0MAX=[x0MAX;x0max];
     u=Ubin(bu)*ones(n0,1);
     Uplot=[Uplot;u];
     x0Plot=[x0Plot;x0];
     
     STB=[STB;stb]; 
     end

   Uplot=Uplot(2:end); x0Plot=x0Plot(2:end);x0MIN=x0MIN(2:end); x0MAX=x0MAX(2:end);%Test=Test(2:end); %Ptest=Ptest(2:end);
 nNaN=nNaN(2:end); STB=STB(2:end);
 Eqms=[Uplot(:),x0Plot(:),x0MIN,x0MAX,nNaN(:),STB(:)];
 Eqms=sortrows(Eqms,1);
 
istable=find(Eqms(:,6)<0);
 iun=find(Eqms(:,6)>=0);
 
 EqmS=Eqms(istable,:);
 EqmU=Eqms(iun,:);


 %% noiseless model run
 
 cd = 1.3e-3;                % drag coeff
lambda = 4e-4;              % normalized surface feedback parameter from vdW
Q = 1.5e-5;                 % normalized net radiation parameter from vdW


 
 ds=0.001;Smax=5;  dx=0.0001; dU=0.001;  itts=50; x01=max(x0Plot); 
 
[Umodel,xmodel,~]= vw_pseudoArclength(Q,cd,lambda,x01,dx,dU,ds,Smax,itts);
save('MatLabEquivOfRrun')  

figure 
grid on 
hold on
plot(Umodel,xmodel,'r.')
errorbar(Eqms(:,1),Eqms(:,2),Eqms(:,3)-Eqms(:,2),Eqms(:,4)-Eqms(:,2),'ko')
scatter(EqmS(:,1),EqmS(:,2),[],EqmS(:,5),'filled')
scatter(EqmU(:,1),EqmU(:,2),[],EqmU(:,5),'filled','h')
legend('Model','extracted','location','best')
xlabel('U -scaled'); ylabel('x0')
caxis([0,1]);
colormap(flipud(gray))
title('Matlab run') 
axis([0,2.5,-0.01,0.08])
h=colorbar; ylabel(h,'fraction of zeros obtained')
str1=strcat('eqmsVW',description,'.png');


figure 
imagesc(Ubin,xfit,D1xfit)
colormap jet
caxis([-0.015,0.015])
axis([0,2.5,-0.01,0.08])
xlabel('U - scaled'); ylabel('x')
set(gca,'YDir','normal')
colorbar
title('Drift - Matlab')


figure 
imagesc(Ubin,xfit,sqrt(D2xfit))
colormap jet
set(gca,'YDir','normal')
%caxis([0,5e-3])
colorbar
axis([0,2.5,-0.01,0.08])
xlabel('U - scaled'); ylabel('x')
title('Diffusion - Matlab')
set(gca,'YDir','normal')
str3=strcat('driftVW',description,'.png');


