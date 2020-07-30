close all; clear all; 

description='FirstAttempt'; 

cd = 1.3e-3;                % drag coeff
lambda = 4e-4;              % normalized surface feedback parameter from vdW
Q = 1.5e-5;                 % normalized net radiation parameter from vdW

dt=30; sigma = .0003; N=1e5; nstep=10; Umn=1; Usig=0.7; 

[U,x] = vdw_SDEmodel(cd, Q, lambda, sigma,dt,N,nstep,Umn,Usig); 
U=U{1}(:); x=x{1}(:); 

Nu=30; Nx=1000; 

xfit=linspace(min(x),max(x),Nx)'; dx=xfit(2)-xfit(1);
[~,Ubin,binU] = histcounts(U,Nu); %bin over X
 
%convert Ubin edges to centers - histcounts outputs bin edges
 Ubin=Ubin(:); Bu=length(Ubin); du=(Ubin(2:Bu)-Ubin(1:Bu-1));Ubin=Ubin(1:Bu-1) + 0.5*du; Bu=length(Ubin); 
 
%% intialize 
 xsave=cell(1,Nu);
 xfitsave=cell(Nu,1); 
%% calculate GPR model fit for each U, output D1x, D2x, x0's for each  
 for bu=1:Bu
     %% extract uconst points
     indU=binU==bu;          %index of points in data that have Ubin(bu)
   xsave{bu}=x(indU); %pull corresponding Dx, xU points at that Uvalue
   [~,indS]=min(abs(xfit-min(xsave{bu})));[~,indE]=min(abs(xfit-max(xsave{bu}))); %% find parts of xfit vector corresponding to xU range 
   xfitsave{bu}=[xfit(indS):dx:xfit(indE)]';
   
   
 end 

 Usave{1}=U;
 Usave{2}=Ubin; 
 
 

save(strcat('vwRunForR',description)); 
 
writecell(xsave',strcat('vwXBin',description,'.txt'))
writecell(xfitsave,strcat('vwXfitBin',description,'.txt'))

 