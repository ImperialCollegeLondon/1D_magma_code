function H=Threephase_enthalpy_solver3(Mass_data, Mass_data_old,T_old, H_old, u_all,cp, k ,Lf, dz, cellz, nodez, dt, source, rho_BC, BC_flux, N_component)
% A three phase enthalpy solver extended from the two phase version
% H.H 5-12-2022
% cp should be [cpl, cps, cpg]
% k should be [kl,ks,kg]
% Mass_data: [m1, n1, v1, m2, n2,v2, v3]
Boundary=0;  %0: fixed enthalpy, 1: conservative

N=length(dz);
uf=u_all(1:N+1);
um=u_all(N+2:2*N+2);
if length(u_all)>2*N+2
    ug=u_all(2*N+3:3*N+3);
else
    ug=zeros(N+1,1);
end
m1=Mass_data(:,1);
n1=Mass_data(:,2);
v1=Mass_data(:,3);
m2=Mass_data(:,4);
n2=Mass_data(:,5);
v2=Mass_data(:,6);
v3=Mass_data(:,7);

mm1=interp1(cellz,m1,nodez,'linear','extrap');
nn1=interp1(cellz,n1,nodez,'linear','extrap');
vv1=interp1(cellz,v1,nodez,'linear','extrap');
mm2=interp1(cellz,m2,nodez,'linear','extrap');
nn2=interp1(cellz,n2,nodez,'linear','extrap');
vv2=interp1(cellz,v2,nodez,'linear','extrap');
vv3=interp1(cellz,v3,nodez,'linear','extrap');

if N_component==5
    l1=Mass_data(:,7);
    l2=Mass_data(:,8);
    l3=Mass_data(:,9);
    ll1=interp1(cellz,l1,nodez,'linear','extrap');
    ll2=interp1(cellz,l2,nodez,'linear','extrap');
    ll3=interp1(cellz,l3,nodez,'linear','extrap');
end



Mc=(m1+n1+v1)*cp(1)+(m2+n2+v2)*cp(2)+v3*cp(3);
MMc=interp1(cellz,Mc,nodez,'linear','extrap');
if N_component==3
    u_bar=(mm1+nn1+vv1).*cp(1,1).*uf+(mm2+nn2+vv2).*cp(2,2).*um+vv3.*cp(3,3).*ug;
    k_bar=(m1+n1+v1)*k(1)+(m2+n2+v2)*k(2)+v3*k(3);
else
    u_bar=(mm1+nn1+vv1+ll1).*cp(1,1).*uf+(mm2+nn2+vv2+ll2).*cp(2,2).*um+(vv3+ll3).*cp(3,3).*ug;
    k_bar=(m1+n1+v1+l1)*k(1)+(m2+n2+v2+l2)*k(2)+(v3+l3)*k(3);
end



u_right=u_bar./MMc-uf;
Latent=m1*Lf(1)+n1*Lf(2)+v1*Lf(2);

  
Source1=zeros(N,1); % The right hand side of enthalpy transport equation, advection of latent heat
for i=1:N-1
    if u_right(i+1)>0
        Flux=u_right(i+1)*Latent(i);
    else
        Flux=u_right(i+1)*Latent(i+1);
    end
    Source1(i)=Source1(i)    -Flux/dz(i);
    Source1(i+1)=Source1(i+1)+Flux/dz(i+1);
end
    
    

% kk=k_bar; %(k_bar(1:N-1).*dz(1:N-1)'+k_bar(2:N).*dz(2:N)')./(dz(1:N-1)'+dz(2:N)');
phi_old=Latent./Mc;

% coe=2;
% Source2=zeros(N,1); % The right hand side of enthalpy transport equation, diffusion of the difference between the enthalpy and sensitive heat
% Source2(2:N-1)=coe*(-kk(1:N-2).*(-phi_old(1:N-2)+phi_old(2:N-1))./(dz(1:N-2)+dz(2:N-1))+kk(2:N-1).*(-phi_old(2:N-1)+phi_old(3:N))./(dz(3:N)+dz(2:N-1)))./dz(2:N-1);
% Source2(1)=coe*kk(1)*(phi_old(2)-phi_old(1))./(dz(1)+dz(2))/dz(1);
% Source2(N)=-coe*kk(N-1)*(phi_old(N)-phi_old(N-1))./(dz(N)+dz(N-1))/dz(N);

Source3=zeros(N,1); 
for i=2:N
    Flux=2*(phi_old(i)-phi_old(i-1))/(dz(i)+dz(i-1))*k_bar(i);
    Source3(i)=Source3(i)-Flux/dz(i);
    Source3(i-1)=Source3(i-1)+Flux/dz(i-1);
end
% Source3(1)=coe*kk(1)*(phi_old(2)-phi_old(1))./(dz(1)+dz(2))/dz(1);
% Source3(N)=-coe*kk(N-1)*(phi_old(N)-phi_old(N-1))./(dz(N)+dz(N-1))/dz(N);

dz=dz';


% if N_component==3
%     Wl=Mass_data_old(:,1)+Mass_data_old(:,2)+Mass_data_old(:,3);
%     Ws=Mass_data_old(:,4)+Mass_data_old(:,5)+Mass_data_old(:,6);
%     Wg=Mass_data_old(:,7);
% else
%     Wl=sum(Mass_data_old(:,[1:3 8]),2);
%     Ws=sum(Mass_data_old(:,[4:6 9]),2);
%     Wg=sum(Mass_data_old(:,[7 10]), 2);
% end

% Hl_old=Wl.*T_old*cp(1)+    Mass_data_old(:,1)*Lf(1)+Mass_data_old(:,2)*Lf(2)+Mass_data_old(:,3)*Lf(2);  %!!!
% Hs_old=Ws.* T_old*cp(2);
% Hg_old=Wg.* T_old*cp(3);
% H_old=Hl_old+Hs_old+Hg_old;
% H=CV_transport_scaled_v2(dz,H_old,One_field,One_field,u_bar./MMc,k_bar,dt,Source1+Source2-source,[3,3],[0,0]);  %V2 more accurate!!

H=CV_transport_varified(dz,H_old,Mc, u_bar/MMc ,k_bar,dt,-(Source1+Source3-source), rho_BC*T_old(1)*cp(3), ug, BC_flux,Boundary);

%% corrections
% flux_b=(T_old(2)-T_old(1))/(cellz(2)-cellz(1))*k(1)*sum(Mass_data_old(1,:))*dt;
% flux_t=(T_old(N-1)-T_old(N))/(cellz(N)-cellz(N-1))*k(1)*sum(Mass_data_old(N,:))*dt;   %% fix k
% 
% Increase=-flux_b-flux_t;
% H(2:end-1)=H(2:end-1)+(Increase-sum(H.*dz')+sum(H_old.*dz'))/sum(dz(2:end-1));
end

function results=CV_transport_varified(dz,value_old,W,u,k,dt,Source, extra_H_BC, ug,BC_flux,Boundary)
    % this function solves a varified transport equation with saturation
    % d(val*phi)/dt+d(val*phi*u)/dz=d(W*k d (val/W)dz)dz
    % k and cp should be difined on node instead of on cell
    % in a control volumn manner
    % HH 9-10-2022
    % u: speed, k: diffusivity
    % BC_type, containing two values,   1: fixed,  2: smooth 3: no flux
    N=length(dz);
    value_old=reshape(value_old,[N,1]);
    indexi=zeros(1,(N-2)*3);
    indexj=zeros(1,(N-2)*3);
    value=zeros(1,(N-2)*3);
    
    coe=1;  % theta=1: implicit; theta=0.5: Crank-Nicolson 
%   WK are k_bar values on nodes, starts from the second node, ends on the last but one
    WK=(k(1:N-1).*dz(2:N)+k(2:N).*dz(1:N-1))./(dz(1:N-1)+dz(2:N));
    for i=0:N-3
        indexi(i*3+1)=i+2; indexj(i*3+1)=i+1; value(i*3+1)=      -coe*2*WK(i+1)/W(i+1)/(dz(i+1)+dz(i+2))/dz(i+2);
        indexi(i*3+2)=i+2; indexj(i*3+2)=i+2; value(i*3+2)=1/dt +coe*2*WK(i+1)/W(i+2)/(dz(i+1)+dz(i+2))/dz(i+2)+coe*2*WK(i+2)/W(i+2)/(dz(i+2)+dz(i+3))/dz(i+2);
        indexi(i*3+3)=i+2; indexj(i*3+3)=i+3; value(i*3+3)=      -coe*2*WK(i+2)/W(i+3)/(dz(i+2)+dz(i+3))/dz(i+2);

        if u(i+2)>0
            value(i*3+1)=value(i*3+1)-coe*u(i+2)/dz(i+2);   %CHECK
        else
            value(i*3+2)=value(i*3+2)-coe*u(i+2)/dz(i+2);
        end
      if u(i+3)>0
          value(i*3+2)=value(i*3+2)+coe*u(i+3)/dz(i+2);
      else
          value(i*3+3)=value(i*3+3)+coe*u(i+3)/dz(i+2);
      end      
    end   
  %Lower boundary
  if Boundary==1
      indexi=[indexi, 1,1]; indexj=[indexj,1,2];
      value=[value, 0, 0];
      if u(2)>0
          value(end-1)=1/dt+coe*u(2)/dz(1)+coe*2*WK(1)/W(1)/(dz(1)+dz(2))/dz(1);
          value(end)=-coe*2*WK(1)/W(2)/(dz(1)+dz(2))/dz(1);
      else
          value(end-1)=1/dt+coe*2*WK(1)/W(1)/(dz(1)+dz(2))/dz(1);
          value(end)=coe*u(2)/dz(1)-coe*2*WK(1)/W(2)/(dz(1)+dz(2))/dz(1);
      end

      %Upper boundary
      indexi=[indexi, N,N]; indexj=[indexj,N-1,N];
      value=[value, 0, 0];
      if u(N)>0
          value(end-1)=-coe*u(N)/dz(N)-coe*2*WK(N-1)/W(N-1)/(dz(N)+dz(N-1))/dz(N);
          value(end)=1/dt+coe*2*WK(N-1)/W(N)/(dz(N)+dz(N-1))/dz(N);
      else
          value(end-1)=-coe*2*WK(N-1)/W(N-1)/(dz(N)+dz(N-1))/dz(N);
          value(end)=1/dt-coe*u(N)/dz(N)+coe*2*WK(N-1)/W(N)/(dz(N)+dz(N-1))/dz(N);
      end
  else
      indexi=[indexi, 1]; indexj=[indexj,1]; value=[value, 1];
      indexi=[indexi, N]; indexj=[indexj,N]; value=[value, 1];
  end

  Q=sparse(indexi,indexj,value,N,N);

  % NEEDS TO BE UPDATED FOR coe not 1
%   temp=zeros(N,1);

  rhs=value_old/dt+Source;
  if ug(1)>0 && BC_flux==1
     rhs(1)=rhs(1)+extra_H_BC*ug(1);
  end  
  if Boundary==0
    rhs(1)=value_old(1);
    rhs(N)=value_old(N);
  end

  results=Q\rhs;
  %     aaa=1;
end