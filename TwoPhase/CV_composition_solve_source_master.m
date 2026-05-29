% HH 2/2/2020
% solving the component transportation
function Cphi=CV_composition_solve_source_master(C_all,phi,phi_old,Cphi_old,u_all,dz,dt,kc,BC_type,BC_value,Method,source)
% A control volumn composition solver 
% solves phi*C instead of C!!
% Currently only consider no flux boundary condition for composition. 
% BC_type, 1.fixed compoenent (C*phi) 2. smooth 3,no_flux 4,fixed component density (C) 
% Method: 1, Conservative 2, non-conservative
% H.H 2-22-2020
%If code or comments are written by Catherine Booth, it is annotated CAB

%CAB - sorting out variables
    N=length(dz);
    uf=u_all(1:N+1);
    um=u_all(N+2:2*N+2);   
    
%     Cphi_all=zeros(2*N,1);
%     Cphi_old=C_all.*phi;

%     Sourcef=zeros(N,1); %Source term to correct the diffusion term
%     Sourcem=zeros(N,1); 

% CAB - should this be an input?
    kc_m2f=1e-5;    %ratio between solid and fluid chemical diffusivity
    
%     for i=1:N-1     
%         Flux=2*(kc(i+1)*dz(i)+kc(i)*dz(i+1))*(phi(i+1)-phi(i))/(dz(i)+dz(i+1))^2;         
%         Sourcef(i)=Sourcef(i)    +Flux*C_all(i+1)/dz(i);
%         Sourcef(i+1)=Sourcef(i+1)-Flux*C_all(i+1)/dz(i+1);            
% 
%         Sourcem(i)=Sourcem(i)    -Flux*C_all(N+i+1)/dz(i)*kc_m2f;
%         Sourcem(i+1)=Sourcem(i+1)+Flux*C_all(N+i+1)/dz(i+1)*kc_m2f;
%     end
%     
    
    Cphi_all=zeros(N,4);
% CAB - boundary conditions
    BC_type_loc=BC_type;
    BC_value_loc1=BC_value(1,:);
    BC_value_loc2=BC_value(2,:);
    switch BC_type(1) %Lower boundary
        case 4
            BC_type_loc(1)=1;
            BC_value_loc1(1)=BC_value(1,1)*phi(1); %(2*phi(1)-phi(2));
            BC_value_loc2(1)=BC_value(2,1)*phi(N+1);% (2*phi(N+1)-phi(N+2));
        case 2
            BC_type_loc(1)=1;
            BC_value_loc1(1)=2*Cphi_old(1)-Cphi_old(2); 
            BC_value_loc2(1)=2*Cphi_old(N+1)-Cphi_old(N+2);
    end
    
    switch BC_type(2) %Upper boundary
        case 4
            BC_type_loc(2)=1;
            BC_value_loc1(2)=BC_value(1,2)*phi(N);%(2*phi(N)-phi(N-1));
            BC_value_loc2(2)=BC_value(2,2)*phi(2*N);%(2*phi(2*N)-phi(2*N-1));
        case 2
            BC_type_loc(2)=1;
            BC_value_loc1(2)=2*Cphi_old(N)-Cphi_old(N-1); 
            BC_value_loc2(2)=2*Cphi_old(2*N)-Cphi_old(2*N-1);            
    end
    % CAB - source 1,2
    One_field=ones(N,1);
    Cphi_all(:,1)=CV_transport_scaled_v2_master(dz,Cphi_old(1:N),One_field, One_field,   uf,kc,       dt, -source(:,1),BC_type_loc,BC_value_loc1);
    Cphi_all(:,2)=CV_transport_scaled_v2_master(dz,Cphi_old(N+1:2*N),One_field,One_field,um,kc*kc_m2f,dt, -source(:,2),BC_type_loc,BC_value_loc2);
    % CAB - source 3,4
    if Method==2
        Cphi_all(:,3)=CV_transport_scaled_v2_master(dz,phi_old(1:N)-Cphi_old(1:N),One_field, One_field,       uf,kc,       dt,-source(:,3),BC_type_loc,BC_value_loc1);
        Cphi_all(:,4)=CV_transport_scaled_v2_master(dz,phi_old(N+1:2*N)-Cphi_old(N+1:2*N),One_field,One_field,um,kc*kc_m2f,dt,-source(:,4),BC_type_loc,BC_value_loc2);
    end

switch Method
    case 1
    SUM=ones(N,1);
    case 2
    SUM=sum(Cphi_all(:,1:4),2);
end     
%     Index=SUM<0.18;
%     SUM(Index)=1;
    Cphi=zeros(2*N,1);
    Cphi(1:N)=Cphi_all(:,1)./SUM;
    Cphi(N+1:2*N)=Cphi_all(:,2)./SUM;
    
%     Cphi=max(0,Cphi);
%     Cphi=min(1,Cphi); 
    aaa=1;
end