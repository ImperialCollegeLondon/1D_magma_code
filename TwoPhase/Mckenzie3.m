% clear;

g=9.81;
drho=500;
mu_f=1e3;
mu_m_mix=1e14;
d=2.75e-3;      %grain size
n=3;
b=125;

phi0=0.5;

K0=phi0^n*d^2/b;
um_0=K0/mu_f*drho*g*(1-phi0); %
sigma=(K0*mu_m_mix/mu_f)^0.5;   %Compaction length
tao0=sigma/um_0;    %time scale

Year=365.25*24*3600;


% To_show=(10:10:200)*Year/tao0;
To_show=(10:10:180)*Year/tao0;
% To_show=[0.5,1,2,4,10];

N_to_show=length(To_show);
pos=1;
Z_len=500;

N=500+1;
theta=1;
fixed_dt=1;
dt=To_show(1)/10;%0.005;
End_time=To_show(end)+0.1;


Precision=1e-6;  %Saturation calculating precision
z=linspace(0,Z_len,N);
phi=ones(N,1);
% phi(round(N/5):round(2*N/5))=1.5;
um=zeros(N,1);

indexi=zeros((N-2)*3+3,1);
indexj=zeros((N-2)*3+3,1);
value=zeros((N-2)*3+3,1); 

To_save=1;
if To_save==1 && fixed_dt==1
   phi_data=zeros(N_to_show,N); 
   um_data=zeros(N_to_show,N); 
end
Time=0;
temp1=zeros(N,1);
counter=0;
%%

figure(13)
clf;
subplot(1,2,1)
hold on;
Handel2=plot(phi(1:N),z);
H_T=text(0.2,Z_len*0.8, "T=0",'fontsize',12);
H_dt=text(0.2,Z_len*0.7, "dT=0",'fontsize',12);
xlim([0,1.1])

subplot(1,2,2)
hold on;
Handel3=plot(um(1:N),z);
xlim([-1,0])

Zeroelements=0;
Zeroelements_old=0;
Top=0;
while Time<End_time 
    um_old=um;
    phi_old=phi;
    improve=1;
    iter=0;
    dz=z(2)-z(1);
    
    while improve>Precision && iter<50
        indexi=zeros((N-2-Zeroelements_old)*3+3,1);
        indexj=zeros((N-2-Zeroelements_old)*3+3,1);
        value=zeros((N-2-Zeroelements_old)*3+3,1); 
        
        theta2=1;
        
        indexi(1)=1; indexj(1)=1; value(1)=1;  %boundary condition at the bottom
        for i=0:N-3-Zeroelements_old 
              indexi(i*3+2)=i+2; indexj(i*3+2)=i+1; value(i*3+2)=theta2;
              indexi(i*3+3)=i+2; indexj(i*3+3)=i+2; value(i*3+3)=-2*theta2;
              indexi(i*3+4)=i+2; indexj(i*3+4)=i+3; value(i*3+4)=theta2;
        end
%         for i=0:N-3-Zeroelements_old 
%               indexi(i*3+2)=i+2; indexj(i*3+2)=i+1; value(i*3+2)=theta2*((1-phi0*phi(i+1))+(1-phi0*phi(i+2)))/(1-phi0)/2;
%               indexi(i*3+3)=i+2; indexj(i*3+3)=i+2; value(i*3+3)=-theta2*((1-phi0*phi(i+1))+2*(1-phi0*phi(i+2))+(1-phi0*phi(i+3)))/(1-phi0)/2;
%               indexi(i*3+4)=i+2; indexj(i*3+4)=i+3; value(i*3+4)=theta2*((1-phi0*phi(i+3))+(1-phi0*phi(i+2)))/(1-phi0)/2;
%         end        
        indexi(end-1)=N-Zeroelements_old; indexj(end-1)=N-Zeroelements_old;   value(end-1)=1;   %boundary condition at the top
        indexi(end)=N-Zeroelements_old; indexj(end)=N-Zeroelements_old-1;     value(end)=-1;   %boundary condition at the top
%         for i=1:Zeroelements_old
%            indexi(end-i+1)=N-i+1; indexj(end-i+1)=N-i+1; value(end-i+1)=1; 
%         end 
%         temp=zeros(N-Zeroelements_old,1);
%         temp(2:N-1-Zeroelements_old)=(um_old(3:N-Zeroelements_old)-2*um_old(2:N-1-Zeroelements_old)+um_old(1:N-2-Zeroelements_old))*(1-theta)/dz^2-(1-theta)./phi(2:N-1-Zeroelements_old).^n.*um_old(2:N-1-Zeroelements_old);    
        
        
        Q=sparse(indexi,indexj,value,N-Zeroelements_old,N-Zeroelements_old);
        K=sparse(2:N-Zeroelements_old-1,2:N-Zeroelements_old-1,theta2./phi(2:N-Zeroelements_old-1).^n,N-Zeroelements_old,N-Zeroelements_old);
        LHS=Q/dz^2-K;
        rhs=((1-phi0*phi)/(1-phi0));
        rhs=rhs(1:end-Zeroelements_old);
%         rhs=rhs-temp;
        % Boundary condition
        rhs(end)=0;
        rhs(1)=0;
        um(1:end-Zeroelements_old)=LHS\rhs;  
        um(end-Zeroelements_old+1:end)=0;
        %Set dt
        if ~fixed_dt
            dt=min(min(abs(dz./um)),0.05);
            if pos<=N_to_show
                if Time+dt>To_show(pos)
                    dt=To_show(pos)-Time;
                end
            end
        end


        %%
        %Solving phi
        indexi2=zeros((N-2-Zeroelements)*3+3,1);
        indexj2=zeros((N-2-Zeroelements)*3+3,1);
        value2= zeros((N-2-Zeroelements)*3+3,1); 
        
        indexi2(1)=1; indexj2(1)=1; value2(1)=-um(1)/dz*theta+1/dt;  
        indexi2(2)=1; indexj2(2)=2; value2(2)=um(2)/dz*theta;        
        for i=0:N-3-Zeroelements   
              indexi2(i*3+3)=i+2; indexj2(i*3+3)=i+1; value2(i*3+3)=-um(i+1)/dz/2*theta;
              indexi2(i*3+4)=i+2; indexj2(i*3+4)=i+2; value2(i*3+4)=1/dt;
              indexi2(i*3+5)=i+2; indexj2(i*3+5)=i+3; value2(i*3+5)=um(i+3)/dz/2*theta;
        end
        indexi2((N-3-Zeroelements)*3+6)=N-Zeroelements; indexj2((N-3-Zeroelements)*3+6)=N-Zeroelements;   value2((N-3-Zeroelements)*3+6)=1;   %boundary condition at the top
        indexi2((N-3-Zeroelements)*3+7)=N-Zeroelements; indexj2((N-3-Zeroelements)*3+7)=N-Zeroelements-1;   value2((N-3-Zeroelements)*3+7)=-1;   %boundary condition at the top
%         for i=1:Zeroelements
%            indexi2(end-i+1)=N-i+1; indexj2(end-i+1)=N-i+1; value2(end-i+1)=1; 
%         end   
        
        LHS=sparse(indexi2,indexj2,value2,N-Zeroelements,N-Zeroelements);
        temp1=zeros(N,1);
        if theta<=1
            temp1(2:N-1-Zeroelements)=-(um_old(3:N-Zeroelements).*phi_old(3:N-Zeroelements)-um_old(1:N-2-Zeroelements).*phi_old(1:N-2-Zeroelements)-...
                (um_old(3:N-Zeroelements)-um_old(1:N-2-Zeroelements))/phi0)/2/dz*(1-theta)+(um(3:N-Zeroelements)-um(1:N-2-Zeroelements))/phi0/2/dz*theta;
            temp1(1)=-(um_old(2).*phi_old(2)-um_old(1).*phi_old(1)-(um_old(2)-um_old(1))/phi0)/dz*(1-theta)+(um(2)-um(1))/phi0/dz*theta; 
        end
%         temp1(Zeroelements+2:N-1)=1/phi0*(um(Zeroelements+3:N)-um(Zeroelements+1:N-2))/2/dz;
        rhs=phi_old/dt+temp1;
%         rhs(N-Zeroelements+1:N)=1;
        rhs(N-Zeroelements)=0; %boundary condition at the top: phi=1 or smooth
        phi_old_nonlinear=phi;
        phi(1:end-Zeroelements)=LHS\rhs(1:end-Zeroelements);
        improve=max(abs(phi-phi_old_nonlinear));
        iter=iter+1;           
    end

    Time=Time+dt;
    disp([counter, iter])
    counter=counter+1;
    set(Handel2, 'xdata',phi(1:end-Zeroelements),'ydata',z(1:end-Zeroelements))
    set(Handel3, 'xdata',um(1:end-Zeroelements),'ydata',z(1:end-Zeroelements))
    eval("set(H_T, 'string', string('T=')+num2str(Time))")
    eval("set(H_dt, 'string', string('dT=')+num2str(dt))")
    if pos<=N_to_show
        if abs((Time-To_show(pos)))<1e-4
            subplot(1,2,1)
            plot(phi(1:N-Zeroelements-2),z(1:N-Zeroelements-2),'linewidth',1,'color','g');
            subplot(1,2,2)
            plot(um(1:N-Zeroelements),z(1:N-Zeroelements),'linewidth',1,'color','g');
            if To_save==1
                phi_data(pos,:)=phi;
                um_data(pos,:)=um;
            end
            pos=pos+1;
        end
    end

    Top=Top-um(end-Zeroelements)*dt;
%     Moving upper boundary.....this procedure is soo fucking ugly!
    Zeroelements=ceil((Top+1e-7)/dz);
    Zeroelements_old=Zeroelements;
    drawnow;

end

if To_save==1
   writematrix(phi_data','MacKenzie.txt')
end
