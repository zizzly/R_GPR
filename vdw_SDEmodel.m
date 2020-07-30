function [U,x] = vdw_SDEmodel(cd, Q, lambda, sigma,dt,N,nstep,Umn,Usig)
% uses normalized equations defined by adam 
%Q      = net radiation 
%lambda = earth coupling parameter in model 
%cd     = neutral drag coefficient  1.3e-3;
%N      = mumber of saved values 1e5
%nstep  = number of steps per saved value 10
%dt     = timestep 30
%Umn    =  mean U  , 1
%Usig   = standard deviation of U,  0.7
% sigma = additive noise .0003;  


randn('state',sum(100*clock));
t = (0:N-1)*dt*nstep;

tau_u = .1*max(t);          % autocorrelation e-folding time of U

T = max(t);

for jj=1:2                  %jj=1 red noise U run ; jj=2 monotonically increasing U run
    if (jj==2)
        U_sort = sort(U{1});
    end
    x{jj} = zeros(size(t));     % x is normalized inversion strength
    U_save = zeros(size(t));
    x_old = Q/lambda;
    t_ctr = 0;
    u_old = 0;
    v_old = 0;
    for j=1:N
        for k=1:nstep
            if (jj==1)
                u_new = u_old+dt*(-u_old/tau_u)+sqrt(2*dt/tau_u)*randn(1,1);
                v_new = v_old+dt*(-v_old/tau_u)+sqrt(2*dt/tau_u)*randn(1,1);
                u_old = u_new;
                v_old = v_new;
                U_fine = Usig*sqrt((Umn+u_new)^2+v_new^2);
                %U = (((j-1)*nstep+k)*dt/T)*2.5;
            else
                U_fine = U_sort(j);
            end
            Rb = x_old/U_fine^2;
            dxdt = Q-lambda*x_old-cd*U_fine*fm(Rb)*x_old;
            x_new = x_old+dt*dxdt+sqrt(dt)*sigma*randn(1,1);
            x_old = x_new;
            t_ctr = t_ctr+dt;
        end
        x{jj}(j) = x_new;
        U_save(j) = U_fine;
    end
    
    U{jj} = U_save;
end


