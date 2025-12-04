function Temp_R = CV_reactiveflow(dt,dz,C_all,phi, u_all,phi_old,u_all_old,Weight_new)

N = length(dz);
Cs_dx = zeros(N,1);
Cl_dx = zeros(N,1);
temp_cl = zeros(N,1);
temp_cs = zeros(N,1);

for i=1:N-1
    Cs_dx(i) = (C_all(N+i+1)-C_all(N+i))/dz(i);
    Cl_dx(i) = (C_all(i+1)-C_all(i))/dz(i);
end
Cs_dx(N) = C_all(2*N);
Cl_dx(N) = C_all(N);

uf=u_all(1:N+1);
um=u_all(N+2:2*N+2);
uf_old=u_all_old(1:N+1);
um_old=u_all_old(N+2:2*N+2);

for i=1:N
   % if uf(i+1)>0
    temp_cl(i) = (uf(i)*phi(i))*Weight_new + (uf_old(i)*phi_old(i))*(1-Weight_new);
    %else
    %    temp_cl(i) = uf(i+1)*phi(i+1);
    %end
    %if um(i+1)>0
    temp_cs(i) = (um(i)*(1-phi(i)))*Weight_new + (um_old(i)*(1-phi_old(i)))*(1-Weight_new);
    %else
     %   temp_cs(i) = um(i+1)*(1-phi(i+1));
    %end
end

%temp_cs(N) = um(N+1)*(1-phi(N));
%temp_cl(N) = uf(N+1)*phi(N);

Temp_R = (temp_cs.*Cs_dx + temp_cl.*Cl_dx)*dt;

%Temp_R = ((phi(N+1:2*N).*u_all(N+1:2*N)).*Cs_dx + (phi(1:N).*u_all(1:N)).*Cl_dx)*dt;

end