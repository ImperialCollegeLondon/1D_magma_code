% HH 2/2/2020
% Solving the enthalpy equation based on the modified fomulation
function [H,H2]=CV_enthalpy_solve_source_master(H_old,phi,u_all,u_bar,dz,dt,kt,cp,Lf,BC_type,BC_value,source)
% A control volumn enthalpy transport solver 
% Boundary type: 1, fixed temperature (not enthalpy)  2, smooth 3.no
% enthalpy flux 4.fixed enthalpy 5. no sensitive(temperature) heat flux
% 6. smooth temperature
% 
% H.H 2-22-2020
    N=length(dz);
    uf=u_all(1:N+1);
    um=u_all(N+2:2*N+2);   

    kkk=H_old;
  
    Source1=zeros(N,1); % The right hand side of enthalpy transport equation, advection of latent heat
    for i=1:N-1
        if uf(i+1)>0            
            Flux=uf(i+1)*phi(i)^2; 
        else
            Flux=uf(i+1)*phi(i+1)^2;            
        end
        if um(i+1)>0
            Flux=Flux+um(i+1)*(1-phi(i)^2);
        else
            Flux=Flux+um(i+1)*(1-phi(i+1)^2);
        end
        Source1(i)=Source1(i)    -Flux/dz(i);
        Source1(i+1)=Source1(i+1)+Flux/dz(i+1); 
    end
    
    
    dz=dz';
    kk=kt;%(kt(1:N-1).*dz(1:N-1)+kt(2:N).*dz(2:N))./(dz(1:N-1)+dz(2:N));
    %index1=find(phi(1:N)==0);
    %kk(index1)=kt(index1);
    %kk(index1(2:end)-1)=kt(index1(2:end)-1);

    coe=2;
    phi_old=phi;
    Source2=zeros(N,1); % The right hand side of enthalpy transport equation, diffusion of the difference between the enthalpy and sensitive heat 
    Source2(2:N-1)=coe*(-kk(1:N-2).*(-phi_old(1:N-2)+phi_old(2:N-1))./(dz(1:N-2)+dz(2:N-1))+kk(2:N-1).*(-phi_old(2:N-1)+phi_old(3:N))./(dz(3:N)+dz(2:N-1)))./dz(2:N-1);
    Source2(1)=coe*kk(1)*(phi_old(2)-phi_old(1))./(dz(1)+dz(2))/dz(1);
    Source2(N)=-coe*kk(N-1)*(phi_old(N)-phi_old(N-1))./(dz(N)+dz(N-1))/dz(N);
    dz=dz';
%     Source2(N)=2*Source2(N-1)-Source2(N-2);
%     Source2(1)=2*Source2(2)-Source2(3);   
%     Source2=zeros(N,1); % The right hand side of enthalpy transport equation, diffusion of the difference between the enthalpy and sensitive heat 
%     for i=1:N-1
%         Flux=2*(kt(i+1)*dz(i+1)+kt(i)*dz(i))*(phi(i+1)-phi(i))/(dz(i)+dz(i+1))^2;  
%         Source2(i)=Source2(i)    +Flux/dz(i);
%         Source2(i+1)=Source2(i+1)-Flux/dz(i+1);         
%     end
switch BC_type(1)  %Default is for no flux boundary 
    case 1
%     Source2(1)=2*Source2(2)-Source2(3);   %for Dirichlet 
%     Source1(1)=2*Source1(2)-Source1(3);
%     if uf(1)>0
%         Source1(1)=Source1(1)+uf(1)*(2*phi(1)-phi(2))^2/dz(1);   %phi at bc =2 phi(1)-phi(2)
%     end
%     if um(1)>0
%         Source1(1)=Source1(1)+um(1)*(1-(2*phi(1)-phi(2))^2)/dz(1);
%     end
end

switch BC_type(2)
    case 1
%     Source2(N)=2*Source2(N-1)-Source2(N-2);
%     Source1(N)=2*Source1(N-1)-Source1(N-2);
%     if uf(N+1)<0
%         Source1(N)=Source1(N)+uf(N+1)*(2*phi(N)-phi(N-1))^2/dz(N);
%     end
%     if um(N+1)<0
%         Source1(N)=Source1(1)+um(N)*(1-(2*phi(N)-phi(N-1))^2)/dz(N);
%     end
end

    Source1=Source1*Lf;
    Source2=Source2*Lf./cp;
    
    bc_type_loc=BC_type;
    bc_value_loc=BC_value;
    switch BC_type(1)
        case 1
            bc_value_loc(1)=BC_value(1)*cp+Lf*phi(1);
        case 4
            bc_type_loc(1)=1;
        case 5
            bc_type_loc(1)=4;  % Newmann type for the transport solver
            bc_value_loc(1)=Lf*(phi(2)-phi(1))/(dz(1)+dz(2))*2;
        case 6
            bc_type_loc(1)=2;  % smooth
            bc_value_loc(1)=(2*phi(2)-phi(1)-phi(3))*Lf;
        case 7
            bc_type_loc(1)=4;  % Newmann type for the transport solver
            bc_value_loc(1)=BC_value(1);
    end
    
    switch BC_type(2)
        case 1
            bc_value_loc(2)=BC_value(2)*cp+Lf*phi(N);
        case 4
            bc_type_loc(2)=1;
        case 5
            bc_type_loc(2)=4;  % Newmann type for the transport solver
            bc_value_loc(2)=Lf*(phi(N)-phi(N-1))/(dz(N)+dz(N-1))*2;
        case 6
            bc_type_loc(2)=2;  % smooth
            bc_value_loc(2)=(2*phi(N-1)-phi(N)-phi(N-2))*Lf;
        case 7
            bc_type_loc(2)=4;  % Newmann type for the transport solver
            bc_value_loc(2)=BC_value(2);
    end
    
%         H=CV_transport(dz,H_old,    u_bar ,kt/cp ,dt,Source1+Source2,bc_type_loc,bc_value_loc);
    One_field=ones(N,1);
%     H2=0;
%     H=CV_transport_scaled(dz,H_old,One_field,u_bar,kt/cp,dt,Source1+Source2,bc_type_loc,bc_value_loc);
    H=CV_transport_scaled_v2_master(dz,H_old,One_field,One_field,u_bar,kt/cp,dt,Source1+Source2-source,bc_type_loc,bc_value_loc);  %V2 more accurate!!

    aaa=1;
end