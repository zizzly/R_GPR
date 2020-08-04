close all; clear all

description='FirstAttempt';


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

nNaN=1; x0MIN=1; x0MAX=1; Uplot=1; x0Plot=1; STB=1; 
for k=1:Nu-1
    
    D1xfit=DriftFit(:,k);
    xfitu=xfitsave{k}; nu=length(xfitu); D1xfit=D1xfit(1:nu); %puts nans at bottom
    xU=xsave{k}; xU=xU(:); 
    sigN=std(xU); sigF=KDrift(1,k); sigL=KDrift(2,k); 
   
    
    %bad=isnan(D1xfit); D1xfit(bad)=[]; xfitu(bad)=[]; 
    
    
    [x0,stb]=ZeroGPR(xfitu,D1xfit,err,count_max,Nx); x0=sort(x0);  %find zero's
     n0=length(x0);
     
      fprintf('have zeros for bin %f / %f',k,Nu)
     
     %% uncertainty estimates 

     [x0Dist,x0min,x0max,nnn,Leff] = x0uncertaintyR(xU,xfitu,D1xfit,x0,sigN,sigL,sigF,itts,err,count_max,prc,reg); 
   
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

