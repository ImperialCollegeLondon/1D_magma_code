function phi_star=CV_saturation_solve_explicit_master(u_all,phi_old,phi,dz,dt,use_new)
% For solving contribution of compaction, divergence part explicitly given
% use_new defined from 0 to 1, which is from [phi_old to phi]
% H.H 5-14-2020
    N=length(dz);
    Nu=length(u_all)/2;
    use_new=max(0,use_new);
    use_new=min(1,use_new);
    use_old=1-use_new;
    phi_use=phi*use_new+phi_old*use_old;
    flux=zeros(N,1);
    for i=1:N-1
        if u_all(i+1)>0
            flux(i)  = flux(i)  -u_all(i+1)*phi_use(i)/dz(i);
            flux(i+1)= flux(i+1)+u_all(i+1)*phi_use(i)/dz(i+1);
        else
            flux(i)  = flux(i)  +u_all(i+1)*phi_use(i+1)/dz(i);
            flux(i+1)= flux(i+1)-u_all(i+1)*phi_use(i+1)/dz(i+1);
        end        
    end
%     if u_all(1)>0
        flux(1)=flux(1)+u_all(1)*phi_use(1)/dz(1);
%     end
%     if u_all(N+1)<0
        flux(N)=flux(N)-u_all(N+1)*phi_use(N)/dz(N);
%     else
%         flux(N)=flux(N)-u_all(N+1)*phi_use(N)/dz(N);
%     end
    phi_star=phi_old(1:N)+flux*dt;
end