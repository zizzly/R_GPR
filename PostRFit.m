close all; clear all

description='SecondAttempt';


err=1e-5; count_max=10; %for zero finding iteration code
prc=97.5; itts=1000; reg=1e-12;


load(strcat('vwRunForR',description))
 
 
DiffFit=readmatrix('DiffFit.csv');

DriftFit=readmatrix('DriftFit.csv');

BasisDiff=readmatrix('BasisDiff.csv'); 

KDiff=readmatrix('KDiff.csv'); 

KDrift=readmatrix('KDrift.csv'); 

DiffFit=DiffFit(:, 2:30);  DriftFit=DriftFit(:, 2:30);  BasisDiff=BasisDiff(:, 2:30); 
KDiff=KDiff(:, 2:30); KDrif=KDrift(2:end); 

nNaN=1; x0MIN=1; x0MAX=1; Uplot=1; x0Plot=1; STB=1; K_R=nan(Nu,3);
for k=1:Nu-1
     xU=xsave{k}; xU=xU(:); 
      p=99; 
     [~,indS]=min(abs(xfit-prctile(xU,100-p)));[~,indE]=min(abs(xfit-prctile(xU,p))); %% find parts of xfit vector corresponding to xU range
     
     d1xfit=DriftFit(indS:indE,k);
    
       
      xfitu=xfit(indS:indE); 
   
    sigN=std(xU); sigF=KDrift(1,k); sigL=KDrift(2,k); 
   
    K_R(k,1)=sigL;
    K_R(k,2)=sigF;
    K_R(k,3)=sigN; 
    
    %bad=isnan(D1xfit); D1xfit(bad)=[]; xfitu(bad)=[]; 
    
    
    [x0,stb]=ZeroGPR(xfitu,d1xfit,err,count_max,Nx); x0=sort(x0);  %find zero's
     n0=length(x0);
     
      fprintf('have zeros for bin %f / %f',k,Nu)
     
     %% uncertainty estimates 

     [x0Dist,x0min,x0max,nnn,Leff] = x0uncertaintyR(xU,xfitu,d1xfit,x0,sigN,sigL,sigF,itts,err,count_max,prc,reg); 
   
    fprintf('have uncertainty for bin %f / %f \n',k,Nu)
     
    nNaN=[nNaN;nnn];

     x0MIN=[x0MIN;x0min];
     x0MAX=[x0MAX;x0max];
     u=Ubin(k)*ones(n0,1);
     Uplot=[Uplot;u];
     x0Plot=[x0Plot;x0];
     STB=[STB;stb]; 
end
     
nNaN=nNaN(2:end); x0MIN=x0MIN(2:end); x0MAX=x0MAX(2:end); Uplot=Uplot(2:end); x0Plot=x0Plot(2:end); STB=STB(2:end); 
   
Eqms=[Uplot(:),x0Plot(:),x0MIN,x0MAX,nNaN(:),STB(:)];
 Eqms=sortrows(Eqms,1);
 
istable=find(Eqms(:,6)<0);
 iun=find(Eqms(:,6)>=0);
 
 EqmS=Eqms(istable,:);
 EqmU=Eqms(iun,:);
 
 
ds=0.001;Smax=5;  dx=0.0001; dU=0.001;  itts=50; x01=max(x0Plot); 
 
[Umodel,xmodel,~]= vw_pseudoArclength(Q,cd,lambda,x01,dx,dU,ds,Smax,itts);
 
figure 
grid on 
hold on
plot(Umodel,xmodel,'k.')
errorbar(Eqms(:,1),Eqms(:,2),Eqms(:,3)-Eqms(:,2),Eqms(:,4)-Eqms(:,2),'ko')
scatter(EqmS(:,1),EqmS(:,2),[],EqmS(:,5),'filled','s')
scatter(EqmU(:,1),EqmU(:,2),180,EqmU(:,5),'filled','p')%legend('Model','extracted','location','best')
xlabel('U -scaled'); ylabel('x0')
caxis([0,1]);
colormap(flipud(winter))
axis([0,2.5,-0.01,0.08])
h=colorbar; ylabel(h,'fraction of zeros obtained')


 
 

%drift=figure('visible','off'); 
figure
imagesc(Ubin(1:end-1),xfit,DriftFit)
colormap jet
caxis([-0.015,0.015])
xlabel('U - scaled'); ylabel('x')
set(gca,'YDir','normal')
colorbar
title('Drift')
%saveas(drift,'driftAdam.png');

%diff=figure('visible','off');
figure
imagesc(Ubin(1:end-1),xfit,sqrt(DiffFit))
colormap jet
set(gca,'YDir','normal')
caxis([0,5e-3])
colorbar
xlabel('U - scaled'); ylabel('x')
title('Diffusion')
set(gca,'YDir','normal')
%saveas(diff,'diffAdam.png');

load('MatLabEquivOfRrun')  

figure
hold on 
grid on
plot(Ubin,K_R(:,1),'o')
plot(Ubin,KDriftSave(:,1),'o')
xlabel('U')
ylabel('sig L')
legend('R','Matlab')


figure
hold on 
grid on
plot(Ubin,K_R(:,2),'o')
plot(Ubin,KDriftSave(:,2),'o')
xlabel('U')
ylabel('sig F')
legend('R','Matlab')


figure
hold on 
grid on
plot(Ubin,K_R(:,3),'o')
plot(Ubin,KDriftSave(:,3),'o')
xlabel('U')
ylabel('sig N')
legend('R','Matlab')



