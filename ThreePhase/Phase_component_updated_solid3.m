function [phi,S,T,rho,iter,T_S_region,Ts,Ts_local,Mass_data,S_cap,Tl_local,Case, Saturation, Partition_CL_CU]=Phase_component_updated_solid3(Input_quantities,P, T_guess, T_S_region, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant, cp,   rho_constant, Lf, Partition_coef, Ts, Tl, min_por, max_melt, dz, Data_point,Data_y,Data_point2,Data_y2,PD_range, tol, Conservation_type, simplified_TS, force_solid)
% solve the saturations and component concentrations in a three phase system based on a eutectic melting, partitian coefficient and melt saturation relationship
% HH 5/11/2022
% Sys_constant contains coefficient [a1,b1,c1; a2,b2,c2; aa,bb,cc; cap_A, cap_B, 0]; 
% Jacobians and Rhs:  Jacobian matrces and corresponding Rhs for situations at 1. between solidus and liquidus 2. at solidus; eacth contain situation, s0l1, s0l0, s1l1, s1l0
% cp is [cpl, cps, cpg]
% rho is [rhol, rhos, rhog];
% rho_constant is [rhol_min, rhol_max; rhos_min, rhos_max]
% Lf is [Lf_most, Lf_least]
% K is the partition coefficient
% min_por: minimal porosity of the solid
% New version, decouple the desity, porosity and saturation with the rest
% v2: guarantees that volatile phase disappear when the melt phase: option 2
% v3: adding H2O dependant solidus and liquidus
% v4: Ts based on H2O of the 'first melt'
% v5: Adding minimal solid porosity
% disapperar
% option 3: give a solid saturation
% updated 24-05-2023: A new version of the code, making it more compact and easier to extend for extra components 
% From this update, the chemical_solver_initiator must be called first before this function can be used


% Currently the latent heat from the water is ignored!!
if length(Input_quantities)==4
    Add_CLCU=0;
else
    Add_CLCU=1;
    CL=Input_quantities(5);
    CU=Input_quantities(6);
end

H=Input_quantities(1);
M=Input_quantities(2);
N=Input_quantities(3);
V=Input_quantities(4);

K=Partition_coef(1);  %m/s water
D2=Partition_coef(2); %m/s Cl
D3=Partition_coef(3); %m/s Cu

% P is in kbar 
P2=P*1e3; % Pressure in bar instead of kbar
P3=P*1e2; % P in Mpa 
PT=[-4e-4,-5e-4,-6e-4,-13e-4, -15.5e-4,-17e-4,-16e-4,-5e-4, 0, 26e-4,5e-3, 5e-3];  % coefficient from (ref) for the temperature dependency of water saturation 
PTx=[0,   0.12    0.2  0.3    0.5       1       2       3   4  5,    11, 20];
ind=find(P<=PTx,1,'first');
if isempty(ind)
    ind=length(11);
end

dSdT=(PT(ind)*(P-PTx(ind-1))+PT(ind-1)*(PTx(ind)-P))/(PTx(ind)-PTx(ind-1));

% [~,~, Solid_saturation]=Cal_liquidus(19,P2,Data_point2,Data_y2);
% Solid_saturation=0.05;%(Solid_saturation+2)/100;
% a1=Sys_constant(1,1);
% b1=Sys_constant(1,2);
% c1=Sys_constant(1,3);
% 
% a31=Sys_constant(3,1);
% b31=Sys_constant(3,2);
% c31=Sys_constant(3,3);
% 
% a32=Sys_constant(4,1);
% b32=Sys_constant(4,2);
% c32=Sys_constant(4,3);
% 
aa=Sys_constant(2,:);
bb=Sys_constant(3,:);
cc=Sys_constant(4,:);
cap_A=Sys_constant(5,1);
cap_B=Sys_constant(5,2);

rholmin=rho_constant(1,1);
rholmax=rho_constant(1,2);
rhosmin=rho_constant(2,1);
rhosmax=rho_constant(2,2);

rho_mean=Sys_constant(6,3);


D0P=0.005654/(1+exp((172.8-P3)/40.42));
Coea=0.05749*100;
S1=PD_range(2);
S0=PD_range(1);
Coeb=0.1099*S1;
Coec=0.1099*S0;

%% Calculate the solidus and liquidus based on current H2O and update the liquidus coefficient
Ts0=Ts;
Tl0=Tl;

Ts=Cal_solidus((V)/(M+N+V)*100/K,P2,Data_point,Data_y);
Tl=Cal_liquidus((V)/(M+N+V)*100/K,P2,Data_point2,Data_y2);

% Ts=740.98;
% Tl=1020.1;

Sys_constant(1,1)=Sys_constant(1,1)*(Tl-Ts)/(Tl0-Ts0);
Sys_constant(1,2)=Sys_constant(1,2)*(Tl-Ts)/(Tl0-Ts0);
Sys_constant(1,3)=Tl;%(Sys_constant(1,3)-Ts0)*(Tl-Ts)/(Tl0-Ts0)+Ts;

%%
%Main constraints:
% 1. Conservation of mass M: (m1+m2-M)*1000=0
% 2. Conservation of mass N: n1+n2-N=0
% 3. Conservation of mass V: v1+v2+v3-V=0
% 4. Conservation of enthalpy H: [(m1+n1+v1)*cp1+(m2+n2+v2)*cp2+v3*cp3]*T+(m1*Lf1+n1*Lf2)-H=0  
% 5. Major component distribution from the melting data: 
%    (1) above liquidus: m2=0; n2=0;
%    (2) between solidus and liquidus:%%    a1*(m1/(m1+n1))^2+b1*(m1/(m1+n1))+c1-T=0 ;  m2=0;
%    (3) at solidus: n1=0; T=Ts
%    (4) below solidus: m1=0; n1=0; v1=0;
% 6. For solid solution, solid data satisfy  (alpha*(1-x).^n+(1-alpha)*(1-x.^(1/n)))*(Tl-Ts)+Ts-T=0
% 7. Water partitioning between volatile and melt phases: % v1/(m1+n1+v1)-2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5=0 or v3=0 
% 8. Water partitioning between solid and melt phases:  v2/(m2+n2+v2)-Solid_cap=0   or v2*(m1+n1+v1)-D1*v1*(m2+n2+v2)=0;
% 9. Solid water saturation: cap_A*m2+cap_B*n2-Solid_cap*(m2+n2)=0;
    

% x=[Cl1,Cs1,Cl2,Cs2,phi,S,T,rho(1),rho(2), rho(3)]';


Saturation=0;
T=T_guess;

cb=M/(M+N);
Tl_local=Sys_constant(1,1)*cb^2+Sys_constant(1,2)*cb+Sys_constant(1,3);
if simplified_TS==1
    Ts_local=Sys_constant(6,2)*(1-cb).^Sys_constant(6,1)*(Tl-Ts)+Ts;
else
    Ts_local=(Sys_constant(6,2)*(1-cb).^Sys_constant(6,1)+(1-Sys_constant(6,2))*(1-cb.^(1/Sys_constant(6,1))))*(Tl-Ts)+Ts;
end


Saturation0=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100; %Pressure dependant Saturation
iter=0;


Steps=0.9;
Steps0=0.15;
Max_iter=300; %Max iteration for the Newton's method

Found=0;
Case=0;
% Next_test=T_S_region(1);
while (Found~=1)
    for i=1:6    % 1. below solidus 2. above liquidus 3. between and s0l0 4. s1l0 5. s0l1 6. s1l1
        switch i
            case 1
                sub_iter=0;
                m1=0;
                n1=0;
                m2=M;
                n2=N;
                %             Saturation=Saturation0+(T-800)*dSdT/100;
                S_cap=(M*cap_A+N*cap_B)/(M+N)+0.04;
                if V/(M+N+V)<=S_cap
                    v1=0;
                    v2=V;
                    v3=0;
                else
                    v1=0;
                    v2=S_cap*(M+N)/(1-S_cap);
                    v3=V-v2;
                end
                T=(H-Lf(1)*m1-Lf(2)*n1-Lf(2)*v1)/((m1+n1+v1)*cp(1)+(m2+n2+v2)*cp(2)+v3*cp(3));
                if T<Ts_local+3 || force_solid==1
                    T_S_region(1)=1;
                    Found=1;
                    Case=1;
                    break
                end
            case 2
                sub_iter=0;
                m2=0;
                n2=0;
                m1=M;
                n1=N;
                
                error=10;
                while(error>1)
                Saturation=Saturation0+(T-800)*dSdT/100;
                    if V/(m1+n1+V)<Saturation
                        v1=V;
                        v2=0;
                        v3=0;
                    else
                        v1=(m1+n1)*Saturation/(1-Saturation);
                        v2=0;
                        v3=V-v1;
                    end
                    T_new=(H-Lf(1)*m1-Lf(2)*n1-Lf(2)*v1)/((m1+n1+v1)*cp(1)+(m2+n2+v2)*cp(2)+v3*cp(3));
                    error=abs(T-T_new);
                    T=T_new;
                end
                if T>Tl_local
                    S_cap=0;
                    T_S_region(1)=2;
                    Found=1;
                    Case=2;
                    break
                end
            case 3 %s0l0
                Input=[M/2,N/2,V/2,M/2,N/2,V/2, 0, Tl_local*0.7+Ts_local*0.3, cap_A, M,N,V,H,cp,Lf,P3,cap_A,cap_B, Sys_constant(1,:), K, Sys_constant(6,2), Sys_constant(6,1),  Tl, Ts, dSdT];
                [Output,sub_iter]=subsolve(Constant_index{1}, Input, Jacobians{1}, Rhs{1}, tol, Steps, Max_iter, 1);                

                if sub_iter<Max_iter && any(Output~=0)
                else
                    [Output,sub_iter]=subsolve(Constant_index{1}, Input, Jacobians{1}, Rhs{1}, tol, Steps/2, Max_iter, 1);
                end

                Saturation=Saturation0+(Output(8)-800)*dSdT/100;

                if Output(3)/sum(Output(1:3))<=Saturation*1.01 %l0
                    if Output(6)/sum(Output(4:6))<=Output(9)*1.01 %s0
                        T_S_region(1)=3;
                        Found=1;
                        iter=iter+sub_iter;
                        m1=Output(1); n1=Output(2); v1=Output(3); m2=Output(4); n2=Output(5); v2=Output(6); v3=Output(7); T=Output(8); S_cap=Output(9);
                        Case=3;
                        break
                    end
                end
                
            case 5 %s1l0
                Input=[M/2,N/2,V/2,M/2,N/2,V/2, 0, Tl_local*0.8+Ts_local*0.2, cap_A, M,N,V,H,cp,Lf,P3,cap_A,cap_B, Sys_constant(1,:), K, Sys_constant(6,2), Sys_constant(6,1),  Tl, Ts, dSdT];
                [Output,sub_iter]=subsolve(Constant_index{2}, Input, Jacobians{2}, Rhs{2}, tol, Steps, Max_iter, 1);

                if sub_iter<Max_iter && any(Output~=0)
                else
                    [Output,sub_iter]=subsolve(Constant_index{2}, Input, Jacobians{2}, Rhs{2}, tol, Steps/2, Max_iter, 1);
                end

                Saturation=Saturation0+(Output(8)-800)*dSdT/100;              
                if Output(3)/sum(Output(1:3))<=Saturation*1.01 && sub_iter<Max_iter
                    T_S_region(1)=4;
                    Found=1;
                    iter=iter+sub_iter;
                    m1=Output(1); n1=Output(2); v1=Output(3); m2=Output(4); n2=Output(5); v2=Output(6); v3=Output(7); T=Output(8); S_cap=Output(9);
                    Case=5;
                    break
                end


            case 4 %s0l1
                Input=[M/2,N/2,V/3,M/2,N/2,V/3,V/3,Tl_local*1.+Ts_local*0., cap_A, M,N,V,H,cp,Lf,P3,cap_A,cap_B,Sys_constant(1,:), K, Sys_constant(6,2), Sys_constant(6,1), Tl, Ts, dSdT];
                [Output,sub_iter]=subsolve(Constant_index{3}, Input, Jacobians{3}, Rhs{3}, tol, Steps*0.8, Max_iter, 1);
                
                if sub_iter<Max_iter && any(Output~=0)
                else
                    [Output,sub_iter]=subsolve(Constant_index{3}, Input, Jacobians{3}, Rhs{3}, tol, Steps0, Max_iter, 0);
                end

                if Output(6)/sum(Output(4:6))<=Output(9)*1.1 && sub_iter<Max_iter
                    T_S_region(1)=5;
                    Found=1;
                    iter=iter+sub_iter;
                    m1=Output(1); n1=Output(2); v1=Output(3); m2=Output(4); n2=Output(5); v2=Output(6); v3=Output(7); T=Output(8); S_cap=Output(9);
                    Case=4;
                    break
                end

            case 6 %s1l1
                Input=[M/2,N/2,V/3,M/2,N/2,V/3,V/3,Tl_local*0.2+Ts_local*0.8, cap_A, M,N,V,H,cp,Lf,P3,cap_A,cap_B,Sys_constant(1,:), K, Sys_constant(6,2), Sys_constant(6,1), Tl, Ts, dSdT];
                [Output,sub_iter]=subsolve(Constant_index{4}, Input, Jacobians{4}, Rhs{4}, tol, Steps0, Max_iter, 1);
                if sub_iter<Max_iter && any(Output~=0)
                    Saturation=Saturation0+(Output(8)-800)*dSdT/100;
                    if Output(3)/sum(Output(1:3))>=Saturation*0.9 && Output(6)/sum(Output(4:6))>=Output(9)*0.9
                        T_S_region(1)=6;
                        Found=1;
                        iter=iter+sub_iter;
                        m1=Output(1); n1=Output(2); v1=Output(3); m2=Output(4); n2=Output(5); v2=Output(6); v3=Output(7); T=Output(8); S_cap=Output(9);
                        Case=6;
                        break
                    end
                else
                    Input=[M/2,N/2,V/3,M/2,N/2,V/3,V/3,Tl_local*1.+Ts_local*0., cap_A, M,N,V,H,cp,Lf,P3,cap_A,cap_B,Sys_constant(1,:), K, Sys_constant(6,2), Sys_constant(6,1), Tl, Ts, dSdT];
                    [Output,sub_iter]=subsolve(Constant_index{4}, Input, Jacobians{4}, Rhs{4}, tol, Steps0/2, 500, 1);
                    Found=1;
                    m1=Output(1); n1=Output(2); v1=Output(3); m2=Output(4); n2=Output(5); v2=Output(6); v3=Output(7); T=Output(8); S_cap=Output(9);
                    Case=6;
                    break
                end
        end
        iter=iter+sub_iter;
    end
end

% update the density to its reference value
% if m1+n1>0
%     Cl1=m1/(m1+n1);    
% else
%     Cl1=1;
% end
% 
% if m2+n2>0
%     Cs1=m2/(m2+n2);
% else
%     Cs1=0;
% end
m1=max(m1,0); n1=max(n1,0); v1=max(0,v1); m2=max(m2,0); n2=max(n2,0); v2=max(0,v2); v3=max(v3,0);

% update the density to its real value for MC scheme or update the mass to have
% the original density for VC scheme
if Conservation_type==1
    rho(1)=(rholmin*m1+rholmax*n1)/max(m1+n1,1e-3); %not considering the influence of water
    rho(2)=(rhosmin*m2+rhosmax*n2)/max(m2+n2,1e-3);
    rho(3)=(aa(1)* max(T,300).^bb(1)+aa(2)*max(P2,400).^bb(2)+aa(3)*max(T,300).^bb(3).*max(P2,400).^cc(1))*1000;    %pressure in bar (NOT in Mpa,  HUBER2017 paper is incorrect!)
    V10=(m1+n1+v1)/max(rho(1), 1);
    V20=(m2+n2+v2)/max(rho(2),1);
    V30=v3/max(rho(3),1);

    V_sum=V10+V20+V30;
    V1=V10/V_sum;
    V2=V20/V_sum;
    V3=V30/V_sum;

    if V2>1-min_por %If the porosity is less than min_por, force the solid phase to have min_por porosity. The extra space is occupied by the volatile phase.
        rho(2)=rho(2)*V2/(1-min_por);
        rho(3)=rho(3)*V3/(V3+V2-(1-min_por));
        V3=V3+V2-(1-min_por);
        V2=1-min_por;
    end

    if V1>max_melt
        dm=m1*(1-max_melt/V1);
        dn=n1*(1-max_melt/V1);
        m1=m1-dm;
        n1=n1-dn;
        m2=m2+dm;
        n2=n2+dn;
        V2=V2+(V1-max_melt);
        V1=max_melt;
    end

    phi=(V1+V3);

    if V1+V3>0
        S=V3/(V1+V3);
    else
        S=0;
    end
    if V1>0
        rho(1)=(m1+n1+v1)/V1;
    end
    if V2>0
        rho(2)=(m2+n2+v2)/V2;
    end
    if V3>0
        rho(3)=        v3/V3;
    end
else
    rho(1:3)=rho_mean;
    phi=(m1+n1+v1+v3)/(M+N+V);
    S=v3/max(m1+n1+v1+v3,1e-5);
end

if Add_CLCU==1
    Max_iter=200;
    sub_step=0.5;
    if v3>0
        if m1+n1>0
            if m2+n2>0

                funcL=Extra_func{1};
                funcR=Extra_func{2};

                found=0;

                while found==0 && sub_step>0.3
                    l1=CL/3;            l2=CL/3;            l3=CL/3;
                    u1=CU/3;            u2=CU/3;            u3=CU/3;
                    Dcl=1;
                    Dcu=1;
                    error=10;
                    sub_iter=0;
                    while (error>tol && sub_iter <Max_iter)
                        LHS=funcL(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                        RHS=funcR(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                        dx=LHS\RHS;
                        l1=max(l1-dx(1)*sub_step,0); l2=max(l2-dx(2)*sub_step,0); l3=max(l3-dx(3)*sub_step,0);
                        u1=max(u1-dx(4)*sub_step,0); u2=max(u2-dx(5)*sub_step,0); u3=max(u3-dx(6)*sub_step,0);
                        Dcl=max(Dcl-dx(7)*sub_step,0); Dcu=max(Dcu-dx(8)*sub_step,0);
                        error=max(abs(dx(:)));
                        sub_iter=sub_iter+1;
                    end
                    iter=iter+sub_iter;
                    if [l1 l2 l3 u1 u2 u3 Dcl Dcu]>0
                        found=1;
                    else
                        sub_step=sub_step/2;
                        Max_iter=Max_iter*2;
                    end
                end
                if sub_iter>=Max_iter/2-1
                    l2=CL/(1+D2);
                    l1=CL*D2/(1+D2);
                    l3=0;
                    u1=CU/(1+D3);
                    u2=CU*D3/(1+D3);
                    u3=0;
                end
            else % melt +volatile
                funcL=Extra_func{5};
                funcR=Extra_func{6};
                error=1;
                sub_iter=0;
                l1=CL/2;            l2=0;            l3=CL/2;
                u1=CU/2;            u2=0;            u3=CU/2;
                Dcl=3;
                Dcu=10;
                while (error>tol && sub_iter <Max_iter)
                    LHS=funcL(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                    RHS=funcR(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                    dx=LHS\RHS;
                    l1=max(l1-dx(1)*sub_step,0);  l3=max(l3-dx(2)*sub_step,0);
                    u1=max(u1-dx(3)*sub_step,0);  u3=max(u3-dx(4)*sub_step,0);
                    Dcl=max(Dcl-dx(5)*sub_step,0); Dcu=max(Dcu-dx(6)*sub_step,0);
                    error=max(abs(dx(:)));
                    sub_iter=sub_iter+1;
                end
                iter=iter+sub_iter;
            end
        else %solid+volatile
            funcL=Extra_func{3};
            funcR=Extra_func{4};
            error=1;
            sub_iter=0;
            l1=0;            l2=CL/2;            l3=CL/2;
            u1=0;            u2=CU/2;            u3=CU/2;
            Dcl=0; Dcu=0;
            while (error>tol && sub_iter <Max_iter)
                LHS=funcL(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                RHS=funcR(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                dx=LHS\RHS;
                l2=max(l2-dx(1),0);  l3=max(l3-dx(2),0);
                u2=max(u2-dx(3),0);  u3=max(u3-dx(4),0);
                error=max(abs(dx(:)));
                sub_iter=sub_iter+1;
            end
            Dcl=nan;
            Dcu=nan;
            iter=iter+sub_iter;
        end
    else
        if m1+n1>0 % solid+melt
            if m2+n2>0
                funcL=Extra_func{7};
                funcR=Extra_func{8};
                error=1;
                sub_iter=0;
                l1=CL/2;            l2=CL/2;            l3=0;
                u1=CU/2;            u2=CU/2;            u3=0;
                Dcl=0; Dcu=0;
                while (error>tol && sub_iter <Max_iter)
                    LHS=funcL(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                    RHS=funcR(m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3,Dcl,Dcu,CL,CU,D2,D3,Coea,Coeb,Coec,D0P);
                    dx=LHS\RHS;
                    l1=max(l1-dx(1)*sub_step,0);  l2=max(l2-dx(2)*sub_step,0);
                    u1=max(u1-dx(3)*sub_step,0);  u2=max(u2-dx(4)*sub_step,0);
                    error=max(abs(dx(:)));
                    sub_iter=sub_iter+1;
                end
                Dcl=nan;
                Dcu=nan;
                iter=iter+sub_iter;
            else
                l1=CL;            l2=0;            l3=0;
                u1=CU;            u2=0;            u3=0;
                Dcl=0; Dcu=0;
            end
        else
            l2=CL; l1=0; l3=0;
            u2=CU; u1=0; u3=0;
            Dcl=nan;
            Dcu=nan;
        end
    end
end



% Mass conservation correction
err=M-m1-m2;
fraction=m1/max(m1+m2,1e-5);
m1=m1+err*fraction;
m2=m2+err*(1-fraction);

err=N-n1-n2;
fraction=n1/max(n1+n2,1e-5);
n1=n1+err*fraction;
n2=n2+err*(1-fraction);

err=V-v1-v2-v3;
fraction=v1/max(v1+v2+v3,1e-5);
fraction2=v2/max(v1+v2+v3,1e-5);
v1=v1+err*fraction;
v2=v2+err*fraction2;
v3=v3+err*(1-fraction-fraction2);

if Add_CLCU==1
    err=CL-l1-l2-l3;
    fraction=l1/max(l1+l2+l3,1e-5);
    fraction2=l2/max(l1+l2+l3,1e-5);
    l1=l1+err*fraction;
    l2=l2+err*fraction2;
    l3=l3+err*(1-fraction-fraction2);

    err=CU-u1-u2-u3;
    fraction=u1/max(u1+u2+u3,1e-5);
    fraction2=u2/max(u1+u2+u3,1e-5);
    u1=u1+err*fraction;
    u2=u2+err*fraction2;
    u3=u3+err*(1-fraction-fraction2);

    if v3+l3+u3<0.1
        v2=v2+v3;
        v3=0;
        l2=l2+l3;
        l3=0;
        u2=u2+u3;
        u3=0;
    end
    Mass_data=max([m1,n1,v1,m2,n2,v2,v3,l1,l2,l3,u1,u2,u3],0);
    Partition_CL_CU=[Dcl, Dcu];
else
    if v3<0.1
        v2=v2+v3;
        v3=0;
    end
    Mass_data=[m1,n1,v1,m2,n2,v2,v3];    
    Partition_CL_CU=[0,0];
end

T=(H-Lf(1)*m1-Lf(2)*n1-Lf(2)*v1)/((m1+n1+v1)*cp(1)+(m2+n2+v2)*cp(2)+v3*cp(3));  %also need to correction for T! ffs....
S_cap=(m2*cap_A+n2*cap_B)/max(m2+n2,1e-2);



end

function [Output,iter]=subsolve(Index_const, in, Jacobian, Rhs, tol, step_size,Max_iter, Force_positive)
    % A simple Newton's method
    index=1:9;
    index(Index_const)=[];
    error=1;
    iter=0;
    while (error>tol && iter <Max_iter)
        LHS=Jacobian(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15),in(16),in(17),in(18),in(19),in(20),in(21),in(22),in(23),in(24),in(25),in(26),in(27),in(28),in(29),in(30));%feval(Jacobian, Input); 
        RHS=Rhs(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15),in(16),in(17),in(18),in(19),in(20),in(21),in(22),in(23),in(24),in(25),in(26),in(27),in(28),in(29),in(30));%feval(Rhs, Input);
        dx=LHS\RHS;
        % if ~isreal(dx)
        %     error('Model definition failure')
        % end
        in(index)=in(index)'-dx*step_size;
        if Force_positive==1
            in(1:9)=max(in(1:9),0);
        end
%         if Dynamic_steps==2
%             in(1:3)=max(0,in(1:3));
%         end
        error=max(abs(dx(:)));
        iter=iter+1;
    end
    Output=in(1:9);
end

