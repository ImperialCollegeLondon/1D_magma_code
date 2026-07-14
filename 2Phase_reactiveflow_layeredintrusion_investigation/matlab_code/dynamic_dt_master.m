function [dt_new,iter] = dynamic_dt_master(dt_way,iter,dt_old,Min_dt,Max_dt,N_dt,Courant, dz, u_all,N,inc_N)


% dt_way=0 means decrease timestep
% dt_way=1 means increase timestep

%dt_steps=linspace(Min_dt,Max_dt,N_dt);

%% decrease timestep
if dt_way==0
    dt_steps=logspace(log10(Min_dt),log10(Max_dt),N_dt);
    
    for i=N_dt:-1:1
        if dt_old>dt_steps(i)
            dt_new = dt_steps(i);
            %disp('Timestep decreased')
            break
        end
        if i==1
            dt_new=Min_dt;
        end
    end

    if dt_new == Min_dt && dt_old~= Min_dt
        disp('You have reached the minimum dt')
       
    end

    if dt_old ~= Min_dt
        iter=0;
    else
        dt_new=dt_old;
    end

%% increase timestep       
elseif dt_way==1

    dt_steps=logspace(log10(Min_dt),log10(Max_dt),N_dt*inc_N);
    
    for i=1:1:N_dt*inc_N
        if dt_old<dt_steps(i)
            dt_new = dt_steps(i);
            %disp('Timestep increased')
            break
        end
        if i==N_dt*inc_N
            dt_new=Max_dt;
        end
    end

    if dt_old == Max_dt
        dt_new = Max_dt;
    end

    if dt_new == Max_dt && dt_old ~= Max_dt
        disp('You have reached the maximum dt')
    end

    
   
end


%% Check dt is smaller than the CFL condition, if not set as CFL condition.
if max(abs(u_all))~=0
    courant_dt = min(abs([Courant*dz'./u_all(1:N); Courant*dz'./u_all(N+2:2*N+1)]));

    dt_new = min(dt_new,courant_dt);
    
end







end