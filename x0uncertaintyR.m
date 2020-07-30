function [x0Dist,x0min,x0max,nNaN,Leff] = x0uncertaintyR(xU,xfit,yfit,x0,sigN,sigL,sigF,itts,err,count_max,prc,reg)


%% extract model parameters 
badn=1e20;
bad=NaN;
Nx=length(xfit);

if length(xU)>1000
    
xsparse=downsample(xU,100);
else 
    xsparse=xU;
end
N=length(xsparse); 



n0=length(x0); 

%% find covariance 


    K_star=sigF^2*exp( (-0.5*(xfit(:)' - xfit(:)).^2)/(sigL^2));
    K=sigF^2*exp(-0.5*(xsparse(:)' - xsparse(:)).^2/(sigL^2)) + sigN^2*eye(N);
    K_x_star=sigF^2*exp(-0.5*(xfit(:)' - xsparse(:)).^2/(sigL^2));
    
    Kinv=inv(K); 
    
   
        
    Keff= K_star; 
    
 
    
    Leff=chol(Keff+ reg*eye(Nx))'; 
    

%% iterate over new realizations 
x0New=1; %initialize 

 for k= 1:itts 
   D1xnew = yfit + Leff*randn(Nx,1); 
   x0new=ZeroGPR(xfit,D1xnew,err,count_max,Nx);  %zero of new realization 
   
   if isempty(x0new) ==1  %if can't find a zero 
       x0new= bad;
   end
   
   x0New=[x0New;x0new]; %stack zeros found 
  
 end

 x0New=x0New(2:end);
 
 
 %% assign to each x0
 
nNaN=(itts-length(find(x0New==bad)))/itts; 
nNaN=nNaN*ones(n0,1);
 

 if length(x0New)>=1
 
 ind=knnsearch(x0,x0New);
 x0Dist=cell(n0,1);
 
 
 x0min=zeros(n0,1); x0max=zeros(n0,1);

 for n=1:n0
     x0Dist{n}=x0New(ind==n);
     itn=length(x0Dist{n});
     nNaN(n)=(itn)/length(x0New);
     x0min(n)=prctile(x0Dist{n},100-prc); x0max(n)=prctile(x0Dist{n},prc); 
 end
 else   %found no zero's for other itterations, so uncertainty is ~inf
     
     x0min=zeros(n0,1); x0max=zeros(n0,1);
     for n=1:n0
     x0Dist{n}=NaN;
     nNaN(n)=0;
     x0min(n)=-badn;
     x0max(n)=badn; 
     end
 
 end
end

 %PLim=[pmin(:), pmax(:)];
 %legend(Legend)

