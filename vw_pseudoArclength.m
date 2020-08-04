

function [U,x,turn]= vw_pseudoArclength(Q,cd,lambda,x01,dx,dU,ds,Smax,itts)

s=0:ds:Smax; Ns=length(s); 

turn=zeros(Ns,1);   Gam_u=zeros(Ns,1); U=zeros(Ns,1); x=zeros(Ns,1); 

G= @(x,U) vwDE(x, U, Q, cd, lambda);

U(1)=0;   %initial guess

x(1)= fzero(@(x)G(x,U(1)),x01); %get true zero for initial guess 



for k=2:Ns
    
    dGdx= ( G(x(k-1)+dx,U(k-1)) - G(x(k-1),U(k-1)))/dx;
    dGdU= ( G(x(k-1),U(k-1)+dU) - G(x(k-1),U(k-1)))/dU; 
    
    %% set direction s steps forward in
    if k==2  
        if dGdx<0
            duds_init=-1;
            Gam_u(1)= duds_init; 
        else
            duds_init=1;
            Gam_u(1)= duds_init;
        end
    end
    %% 
    
    norm=sqrt(1/(dGdx^2 + dGdU^2));
    Gam_u(k)= duds_init*dGdx*norm;
    Gam_x= -duds_init*dGdU*norm; 
    
    if Gam_u(k-1)/Gam_u(k)<0
        turn(k)=1; 
    end
    
    

    
    xnew= x(k-1) + Gam_x*ds; 
    U(k)=U(k-1) + Gam_u(k)*ds; 
    
    dG=@(x)( G(x,U(k)) - G(x,U(k-1)) )/Gam_u(k)*ds;
    %% newton method - itterate to find zero 
    for j=1:itts
        dGnew=dG(xnew);
        if dGnew >0
        xnew = xnew - G(xnew,U(k))/dGnew; 
        end 
    end
    
    x(k)=xnew ; 
end

end 
    
    
        