function results=CV_transport_scaled_one_sided(dz,value_old,phi,phi_old,u,k,dt,Source,BC_type,BC_value)
    % this function solves transport equation with one-sided diffusion
    % d(val*phi)/dt+d(val*phi*u)/dz=dk*phi*dval/dx^2+Source
    % in a control volumn manner
    % HH 2-20-2020
    % u: speed, k: diffusivity
    % BC_type, containing two values,   1: fixed,  2: smooth 3: no flux
    N=length(dz);
    value_old=reshape(value_old,[N,1]);
    phi=reshape(phi,[N,1]);
    indexi=zeros(1,(N-2)*3);
    indexj=zeros(1,(N-2)*3);
    value=zeros(1,(N-2)*3);
    
    coe=1;  % theta=1: implicit; theta=0.5: Crank-Nicolson 
%     dz=dz';
    kk=(k(1:N-1).*dz(2:N)+k(2:N).*dz(1:N-1))./(dz(1:N-1)+dz(2:N));
    for i=0:N-3   
        % if value_old(i+1)>value_old(i+2)
            Flux1=-2*kk(i+1)/(dz(i+1)+dz(i+2))/dz(i+2);
        % else
        %     Flux1=0;
        % end
        % 
        % if value_old(i+2)>value_old(i+3)
            Flux2=-2*kk(i+2)/(dz(i+2)+dz(i+3))/dz(i+2);
        % else
        %     Flux2=0;
        % end
      indexi(i*3+1)=i+2; indexj(i*3+1)=i+1; value(i*3+1)=      Flux1*phi(i+1); %-coe*2*kk(i+1)*phi(i+1)/(dz(i+1)+dz(i+2))/dz(i+2);
      indexi(i*3+2)=i+2; indexj(i*3+2)=i+2; value(i*3+2)=phi(i+2)/dt -Flux1*phi(i+2)-Flux2*phi(i+2); %+coe*2*kk(i+1)*phi(i+2)/(dz(i+1)+dz(i+2))/dz(i+2)+coe*2*kk(i+2)*phi(i+2)/(dz(i+2)+dz(i+3))/dz(i+2);
      indexi(i*3+3)=i+2; indexj(i*3+3)=i+3; value(i*3+3)=      Flux2*phi(i+3);
      
      if u(i+2)>0
          value(i*3+1)=value(i*3+1)-coe*u(i+2)*phi(i+1)/dz(i+2);
      else
          value(i*3+2)=value(i*3+2)-coe*u(i+2)*phi(i+2)/dz(i+2);
      end
      if u(i+3)>0
          value(i*3+2)=value(i*3+2)+coe*u(i+3)*phi(i+2)/dz(i+2);
      else
          value(i*3+3)=value(i*3+3)+coe*u(i+3)*phi(i+3)/dz(i+2);
      end      
    end   
    switch BC_type(1) %Lower boundary
        case 1
            if u(2)>0
                indexi=[indexi, 1,1]; indexj=[indexj,1,2]; value=[value, phi(1)/dt+coe*u(2)*phi(1)/dz(1)+coe*2*kk(1)*phi(1)/(dz(1)+dz(2))/dz(1)+coe*k(1)*phi(1)/dz(1)^2,  -coe*2*kk(1)*phi(2)/(dz(1)+dz(2))/dz(1)];   %Lower no flux
            else
                indexi=[indexi, 1,1]; indexj=[indexj,1,2]; value=[value,            phi(1)/dt+coe*2*kk(1)*phi(1)/(dz(1)+dz(2))/dz(1)+coe*k(1)*phi(1)/dz(1)^2,coe*u(2)*phi(2)/dz(1)-coe*2*kk(1)*phi(2)/(dz(1)+dz(2))/dz(1)];   %Lower no flux  
            end    
            if u(1)<0
                indexi=[indexi, 1]; indexj=[indexj,1]; value=[value,-coe*u(1)*phi(1)/dz(1)];
            end
        case 2
            indexi=[indexi, 1,1,1]; indexj=[indexj,1,2,3]; value=[value, 1,-2,1];   %Lower smooth
        case 3
            if u(2)>0
                indexi=[indexi, 1,1]; indexj=[indexj,1,2]; value=[value, phi(1)/dt+coe*u(2)*phi(1)/dz(1)+coe*2*kk(1)*phi(1)/(dz(1)+dz(2))/dz(1),  -coe*2*kk(1)*phi(2)/(dz(1)+dz(2))/dz(1)];   %Lower no flux
            else
                indexi=[indexi, 1,1]; indexj=[indexj,1,2]; value=[value, phi(1)/dt+coe*2*kk(1)*phi(1)/(dz(1)+dz(2))/dz(1),   coe*u(2)*phi(2)/dz(1)-coe*2*kk(1)*phi(2)/(dz(1)+dz(2))/dz(1)];   %Lower no flux  
            end
            value(end-1)=value(end-1)-coe*u(1)*phi(1)/dz(2);
    end
    
    switch BC_type(2) %Upper boundary
        case 1
            if u(N)>0
                indexi=[indexi, N,N]; indexj=[indexj,N-1,N]; value=[value, -coe*u(N)*phi(N-1)/dz(N)-coe*2*kk(N-1)*phi(N-1)/(dz(N)+dz(N-1))/dz(N),           phi(N)/dt+coe*2*kk(N-1)*phi(N)/(dz(N)+dz(N-1))/dz(N)+coe*k(N)*phi(N)/dz(N)^2];   %Upper no flux
            else
                indexi=[indexi, N,N]; indexj=[indexj,N-1,N]; value=[value,            -coe*2*kk(N-1)*phi(N-1)/(dz(N)+dz(N-1))/dz(N),phi(N)/dt-coe*u(N)*phi(N)/dz(N)+coe*2*kk(N-1)*phi(N)/(dz(N)+dz(N-1))/dz(N)+coe*k(N)*phi(N)/dz(N)^2];   %Upper no flux
            end
            if u(N+1)>0
                indexi=[indexi, N]; indexj=[indexj,N]; value=[value,coe*u(N+1)*phi(N)/dz(N)];
            end
        case 2
            indexi=[indexi, N,N,N]; indexj=[indexj,N,N-1,N-2]; value=[value, 1,-2,1];   %Upper smooth
        case 3
     
            if u(N)>0
                indexi=[indexi, N,N]; indexj=[indexj,N-1,N]; value=[value, -coe*u(N)*phi(N-1)/dz(N)-coe*2*kk(N-1)*phi(N-1)/(dz(N)+dz(N-1))/dz(N),           phi(N)/dt+coe*2*kk(N-1)*phi(N)/(dz(N)+dz(N-1))/dz(N)];   %Upper no flux
            else
                indexi=[indexi, N,N]; indexj=[indexj,N-1,N]; value=[value,            -coe*2*kk(N-1)*phi(N-1)/(dz(N)+dz(N-1))/dz(N),phi(N)/dt-coe*u(N)*phi(N)/dz(N)+coe*2*kk(N-1)*phi(N)/(dz(N)+dz(N-1))/dz(N)];   %Upper no flux
            end
            value(end)=value(end)+coe*u(N+1)*phi(N)/dz(N);
    end
    
    Q=sparse(indexi,indexj,value,N,N);
    
    temp=zeros(N,1);  
    if coe~=1
        temp(2:N-1)=(1-coe)*(2*kk(1:N-2).*(phi_old(1:N-2).*value_old(1:N-2)-phi_old(2:N-1).*value_old(2:N-1))./(dz(1:N-2)+dz(2:N-1))-2*kk(2:N-1).*(phi_old(2:N-1).*value_old(2:N-1)-phi_old(3:N).*value_old(3:N))./(dz(3:N)+dz(2:N-1)))./dz(2:N-1);
        temp(1)=(1-coe)*2*kk(1)*(phi_old(2).*value_old(2)-phi_old(1).*value_old(1))./(dz(1)+dz(2))/dz(1);
        temp(N)=-(1-coe)*2*kk(N-1)*(phi_old(N).*value_old(N)-phi_old(N-1).*value_old(N-1))./(dz(N)+dz(N-1))/dz(N);
        
        for i=2:N-1
            if u(i)>0
                temp(i)=temp(i)+(1-coe)*u(i)*phi_old(i-1)*value_old(i-1)/dz(i);
            else
                temp(i)=temp(i)+(1-coe)*u(i)*phi_old(i)*value_old(i)/dz(i);
            end
            if u(i+1)>0
                temp(i)=temp(i)-(1-coe)*u(i+1)*phi_old(i)*value_old(i)/dz(i);
            else
                temp(i)=temp(i)-(1-coe)*u(i+1)*phi_old(i+1)*value_old(i+1)/dz(i);
            end
        end
        if u(2)>0
            temp(1)=temp(1)-(1-coe)*u(2)*phi(1)*value_old(1)/dz(1);
        else
            temp(1)=temp(1)-(1-coe)*u(2)*phi(2)*value_old(2)/dz(1);
        end
            temp(1)=temp(1)+(1-coe)*u(1)*phi(1)*value_old(1)/dz(1);


        if u(N)>0
            temp(N)=temp(N)+(1-coe)*u(N)*phi(N-1)*value_old(N-1)/dz(N);
        else
            temp(N)=temp(N)+(1-coe)*u(N)*phi(N)*value_old(N)/dz(N);
        end
        temp(N)=temp(N)-(1-coe)*u(N+1)*phi(N)*value_old(N)/dz(N);
    end
    
    rhs=value_old.*phi_old/dt+Source+temp;
    switch BC_type(1) %Lower boundary
        case 1
            rhs(1)= rhs(1)+k(1)*phi(1)*BC_value(1)/dz(1)^2-(1-coe)*k(1)*phi(1)*value_old(1)/dz(1)^2;         %Lower fix
            if u(1)>0
                rhs(1)=rhs(1)+u(1)*phi(1)*BC_value(1)/dz(1);   %value(i*3+1)=value(i*3+1)-u(i+2)*phi(i+1)/dz(i+2);
            end
        case 2
            rhs(1)= 0;   %Lower smooth
        case 3
%             rhs(1)=rhs(1)-(1-coe)*k(1)*phi(1)*value_old(1)/dz(1)^2;
    end

    switch BC_type(2) %Upper boundary
        case 1
            rhs(N)= rhs(N)+k(N)*phi(N)*BC_value(2)/dz(N)^2-(1-coe)*k(N)*phi(N)*value_old(N)/dz(N)^2;         %Lower fix
            if u(N+1)<0
                rhs(N)=rhs(N)-u(N+1)*phi(N)*BC_value(2)/dz(N);
            end
        case 2
            rhs(N)= 0;   %Lower smooth
        case 3
%             rhs(N)=rhs(N)-(1-coe)*k(N)*phi(N)*value_old(N)/dz(N)^2;
            
    end
    
    results=Q\rhs;
%     aaa=1;
    end
    