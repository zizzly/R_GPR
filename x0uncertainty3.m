function [x0Dist,x0min,x0max,nNaN,Lstar] = x0uncertainty3(xfit,yfit,model,x0,itts,err,count_max,prc)

% PLim gives percentile data goes out to 

%% extract model parameters 
bad=NaN;
Nx=length(xfit);
kparam=model.KernelInformation.KernelParameters; 

sigL=kparam(1); sigF=kparam(2);
%sigN=model.Sigma;


n0=length(x0); 

    
    
%% find covariance 
    Kstar=sigF^2*exp( (-0.5*(xfit(:)' - xfit(:)).^2)/(sigL^2));
    
    eigKs=eig(Kstar); neig=eigKs(eigKs<0); reg=10*max(abs(neig));
    if isempty(reg)==1
        reg=1e-12;
    end
    
    Kstar_reg=Kstar + reg*eye(Nx); 
    Lstar=chol(Kstar_reg)'; 

%% iterate over new realizations 
x0New=1; %initialize 

 for k= 1:itts 
   D1xnew = yfit + Lstar*randn(Nx,1);  %new realizationn 
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


