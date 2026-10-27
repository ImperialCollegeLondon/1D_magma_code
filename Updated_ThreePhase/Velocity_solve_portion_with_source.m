function [U_new, Portion]=Velocity_solve_portion_with_source(phi,u_all_old, Mass_data, ZerosU_values, C_values,mu_f,mu_m,rhof,rhom,cellz, nodez, dz,g,Source)
% A velocity solver carried on nodes, suitable for CV version of codes. 
    N=size(Mass_data,1)+1;
    % mu_m=reshape(mu_m,size(phi));
    % mu_f=reshape(mu_f,size(phi));
    % 
    cap=1e-12;
    % eps=1e-12;


    m1=Mass_data(:,1);
    n1=Mass_data(:,2);
    v1=Mass_data(:,3);
    % m2=Mass_data(:,4);
    % n2=Mass_data(:,5);
    % v2=Mass_data(:,6);
    % v3=Mass_data(:,7);

    % 
    Ml=interp1(cellz,m1+n1+v1,nodez,'linear','extrap');
    % Ms=interp1(cellz,m2+n2+v2,nodez,'linear','extrap');
    % Mg=interp1(cellz,v3,nodez,'linear','extrap');
    % M_all=Mg+Ml+Ms;
    % 
    % 
    % solid_fraction_cap=Ms./M_all;

    ZerosU=[1,N];
    % ZerosU=union(ZerosU, find(solid_fraction_cap>1-cap));
    ZerosU=union(ZerosU,find(phi<cap));
    Nonzero=setxor(1:N,ZerosU);

    if isempty(Nonzero)
        U_new=zeros(N*2,1);
        Portion=zeros(N,3);
        return
    end 
    
    if length(Nonzero)<=2
        U_new=zeros(N*2,1);
        Portion=zeros(N,3);
        return
    end

    indexi=zeros(length(Nonzero)*3,1);
    indexj=zeros(length(Nonzero)*3,1);
    value=zeros(length(Nonzero)*3,1);    
    
    
    APHI=(mu_f).*phi;
    for i=1:length(Nonzero)   
%           Flux1=  APHI(Nonzero(i)-1)/dz(Nonzero(i)-1)*2/(dz(Nonzero(i)-1)+dz(Nonzero(i)));
%           Flux2=  APHI(Nonzero(i))/dz(Nonzero(i))*2/(dz(Nonzero(i)-1)+dz(Nonzero(i)));          
%           indexi(i*3-2)=Nonzero(i); indexj(i*3-2)=Nonzero(i)-1; value(i*3-2)= Flux1;
%           indexi(i*3-1)=Nonzero(i); indexj(i*3-1)=Nonzero(i);   value(i*3-1)=-(Flux1+Flux2);
%           indexi(i*3)  =Nonzero(i); indexj(i*3)=Nonzero(i)+1;   value(i*3)  = Flux2;
          Flux1=(APHI(Nonzero(i)+1)+APHI(Nonzero(i)))/dz(Nonzero(i))  /(dz(Nonzero(i)-1)+dz(Nonzero(i)));
          Flux2=(APHI(Nonzero(i)-1)+APHI(Nonzero(i)))/dz(Nonzero(i)-1)/(dz(Nonzero(i)-1)+dz(Nonzero(i)));
          indexi(i*3-2)=Nonzero(i); indexj(i*3-2)=Nonzero(i)-1; value(i*3-2)= Flux2;
          indexi(i*3-1)=Nonzero(i); indexj(i*3-1)=Nonzero(i);   value(i*3-1)=-(Flux1+Flux2);
          indexi(i*3)  =Nonzero(i); indexj(i*3)=Nonzero(i)+1;   value(i*3)  = Flux1;
    end
    Qf=sparse(indexi,indexj,value,N,N);
    
%     Coef=log(0.5)*(1-Dynamic_matrix_vis_coefficient)/(-Dynamic_matrix_vis_coefficient);
%     mu_m=mu_m.*exp(-Coef*phi(1:N)./phi(N+1:2*N))+0.1;
    
    %     APHI=(mu_m*4/3+zeta_m(power,zeta0,phi0,phi(1:N),Precision)).*phi(N+1:2*N);
%     APHI=zeta_m(power,mu_m,phi(1:N),Precision).*phi(N+1:2*N);
    APHI=mu_m.*(1-phi);
    for i=1:length(Nonzero)

          Flux1=(APHI(Nonzero(i)+1)+APHI(Nonzero(i)))/dz(Nonzero(i))  /(dz(Nonzero(i)-1)+dz(Nonzero(i)));
          Flux2=(APHI(Nonzero(i)-1)+APHI(Nonzero(i)))/dz(Nonzero(i)-1)/(dz(Nonzero(i)-1)+dz(Nonzero(i)));
          indexi(i*3-2)=Nonzero(i); indexj(i*3-2)=Nonzero(i)-1; value(i*3-2)= Flux2;
          indexi(i*3-1)=Nonzero(i); indexj(i*3-1)=Nonzero(i);   value(i*3-1)=-(Flux1+Flux2);
          indexi(i*3)  =Nonzero(i); indexj(i*3)=Nonzero(i)+1;   value(i*3)  = Flux1;
    end
    Qm=sparse(indexi,indexj,value,N,N);
    C=sparse(Nonzero,Nonzero,C_values(Nonzero),N,N);
    Mphi=sparse(Nonzero,Nonzero,phi(Nonzero),N,N);
    One_m_Mphi=sparse(Nonzero,Nonzero,1-phi(Nonzero),N,N);

%This part uses the bulk velocity instead of the original continuity
%equation to obtain a more stable solution. This means U_source is the
%integral of the velocity source instead of the source itself

%     double_Nonzero=zeros(1,length(Nonzero)*2);
%     double_Nonzero(1:2:end-1)=Nonzero;
%     double_Nonzero(2:2:end)=Nonzero;
%     
%     double_Nonzero2=zeros(1,length(Nonzero)*2);
%     double_Nonzero2(1:2:end-1)=Nonzero-1;
%     double_Nonzero2(2:2:end)=Nonzero+1;
%     
%     dz2=dz(Nonzero)+dz(Nonzero-1);
%     double_dz2=ones(1,length(Nonzero)*2);
%     double_dz2(1:2:end-1)=-dz2;
%     double_dz2(2:2:end)=dz2;
% 
%     Mdphi=sparse(double_Nonzero,double_Nonzero2,phi(double_Nonzero2)./double_dz2',N,N);
%     One_m_Mdphi=sparse(double_Nonzero,double_Nonzero2,phi(double_Nonzero2+N)./double_dz2',N,N);



    LHS=[One_m_Mphi*Qf-C -Mphi*Qm+C; Mphi One_m_Mphi];
    rhs=[(phi.*(1-phi).*(rhof-rhom))*g; Source];
    %% Imposing boundary counditions v=0 and positions where phi=0 or 1
    for i=1:length(ZerosU)
        if phi(ZerosU(i))>0
            rhs(ZerosU(i))=ZerosU_values(i);
            rhs(ZerosU(i)+N)=0;
        else
            rhs(ZerosU(i))=0;
            rhs(ZerosU(i)+N)=ZerosU_values(i);
        end
    end
    rhs(N+1)=0;
    rhs(2*N)=0;
    LHS=LHS+sparse([ZerosU N+ZerosU],[ZerosU N+ZerosU],ones(length(ZerosU)*2,1),2*N,2*N);
    
    
    U_new=LHS\rhs;
    
    Portion=[One_m_Mphi*Qf*U_new(1:N), -Mphi*Qm*U_new(N+1:2*N), C*(U_new(1:N)-U_new(N+1:2*N))];
    
    U_new(find(Ml==0)+N)=0;
%     u_all=gmres(LHS,rhs);
%     [u_all,flag]=bicgstab(LHS,rhs);
%     if flag ~=0
%         disp(flag)
%     end
    aaa=1;
end