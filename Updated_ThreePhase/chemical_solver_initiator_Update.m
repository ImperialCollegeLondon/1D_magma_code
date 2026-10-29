function [Jacobians,Rhs, Constant_index, Extra_func]=chemical_solver_initiator(is_solid_solution,N_component,simplify_Ts, options)
% Calculate the Jacobian matrix entries which are used for the chemical models based on number of components and given conditions  
% this part should be called first before the main program runs and not in the main non-linear iteration, since it
% relys on Symbolic Math toolbox which is expensive to run 
% Haiayng Hu 24/5/2023
% is_solid_solution: 0 for eutectic 1 for solid solution
% N_component: 3 for A, B, H2O;  4 for Cl as extra
% Cu is treated as trace element
% simplify_Ts: if true, ignore the 1/n term in the solidus formulation

if options==1
    simplified_5p=1;
    add_CLCUS=0;
else
    simplified_5p=0;
    add_CLCUS=1;
end

if N_component~=3 && N_component~=5
    error('Wrong system component number, currently can only be 3 or 5.')
end

% Variables:
% 1. m1, component A in the melt (mass)
% 2. n1, component B in the melt 
% 3. v1, water in the melt
% 4. m2, component A in the solid
% 5. n2, component B in the solid
% 6. v2, water in the solid
% 7. v3, water in volatile phase
% 8. T,  temperature
% 9. S_cap, maximum solid saturation 

% For 5 components, extra variables:  
%10. l1, Cl in the melt
%11. l2, Cl in the solid
%12. l3, Cl in the volatile
%13.Dcl, f/m partitioning coefficient

% Constraints:
% 1. Conservation of mass M: m1+m2-M=0
% 2. Conservation of mass N: n1+n2-N=0
% 3. Conservation of mass V: v1+v2+v3-V=0
%10. Conservation of mass L: l1+l2+l3-CL=0;
% 4. Conservation of enthalpy H: [(m1+n1+v1)*cp1+(m2+n2+v2)*cp2+v3*cp3]*T+(m1*Lf1+n1*Lf2)-H=0  
% 4b:
% ((m1+n1+v1+l1)*cp1+(m2+n2+v2+l2)*cp2+(v3+l3)*cp3)*T+(m1*Lf1+n1*Lf2)-H; 
% 5. Major component distribution from the melting data: 
%    (1) above liquidus: m2=0; n2=0;
%    (2) between solidus and liquidus:%%    a1*(m1/(m1+n1))^2+b1*(m1/(m1+n1))+c1-T=0 ;  m2=0;
%    (3) at solidus: n1=0; T=Ts
%    (4) below solidus: m1=0; n1=0; v1=0;
% 6. For solid solution, solid data satisfy  (a2*(1-x).^n+(1-a2)*(1-x.^(1/n)))*(Tl-Ts)+Ts-T=0
% 6(b) for simplified solution, (1-x).^n*(Tl-Ts)+Ts-T=0
% 7. Water partitioning between volatile and melt phases: % v1/(m1+n1+v1)-(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-(T-800)*dSdT/100; or v3=0 
% 8. Water partitioning between solid and melt phases:  v2/(m2+n2+v2)-S_cap=0   or v2*(m1+n1+v1)-D1*v1*(m2+n2+v2)=0;
% 9. Solid water saturation: cap_A*m2+cap_B*n2-(Solid_cap-cap_por)*(m2+n2)=0; it has 4% fractured water content
%  
%
% 5-component constraints:
%11. Cl partition coefficient (Tattitch 2021): Coea*l3*(m1+n1+v1+l1)+(Coeb*m1+Clc*n1)*(v3+l3)+(Coec2-log(Dcl))*(v3+l3)*(m1+n1+v1+l1)=0
%12. Cl s/m partitioning: D2*l1*(m2+n2+v2+l2)-l2*(m1+n1+v1+l1)=0;
%13. Cl f/m partitioning: l3*(m1+n1+v1+l1)-Dcl*l1*(v3+l3)
%variables
if N_component==3
    syms m1 n1 v1 m2 n2 v2 v3 T S_cap 
    Variables=[m1 n1 v1 m2 n2 v2 v3 T S_cap];
    %constants and coefficients
    syms M N V H cp1 cp2 cp3 Lf1 Lf2 P3 cap_A cap_B a1 b1 c1 D1 a2 n Tl Ts dSdT cap_por
    Constants=[M N V H cp1 cp2 cp3 Lf1 Lf2 P3 cap_A cap_B a1 b1 c1 D1 a2 n Tl Ts dSdT cap_por];
    %constraints 1_3
    eqn1=(m1+m2-M)/1000;
    Jac1=jacobian(eqn1, Variables);
    % eqn1=matlabFunction(eqn1,'Vars', [Variables Constants]);

    eqn2=(n1+n2-N)/1000;
    Jac2=jacobian(eqn2, Variables);
    % eqn2=matlabFunction(eqn2,'Vars', [Variables Constants]);

    eqn3=(v1+v2+v3-V)/1000;
    Jac3=jacobian(eqn3, Variables);
    % eqn3=matlabFunction(eqn3,'Vars', [Variables Constants]);

    %constraint 4
    eqn4=(((m1+n1+v1)*cp1+(m2+n2+v2)*cp2+v3*cp3)*T+(m1*Lf1+n1*Lf2+v1*Lf2)-H)/1e7;
    Jac4 = jacobian(eqn4, Variables);
    % eqn4=matlabFunction(eqn4,'Vars', [Variables Constants]);
    %constraint 5,6

    eqn5=a1*m1^2+b1*m1*(m1+n1)+(c1-T)*(m1+n1)^2;
    Jac5= jacobian(eqn5, Variables);
    % eqn5=matlabFunction(eqn5,'Vars', [Variables Constants]);

    if simplify_Ts==1
        eqn6=((1-m2/(m2+n2)).^n)*(Tl-Ts)+Ts-T;
    else
        eqn6=(a2*(1-m2/(m2+n2)).^n+(1-a2)*(1-(m2/(m2+n2))^(1/n)))*(Tl-Ts)+Ts-T;
    end
    Jac6= jacobian(eqn6, Variables);

    % adjust ment
    if simplify_Ts==0
        Jac6(4)=-(Tl - Ts)*((max(m2/(m2 + n2),0)^(1/n - 1)*(m2/(m2 + n2)^2 - 1/(m2 + n2))*(a2 - 1))/n - a2*n*(m2/(m2 + n2)^2 - 1/(m2 + n2))*(1 - m2/(m2 + n2))^(n - 1));
        Jac6(5)=(Tl - Ts)*((a2*m2*n*(1 - m2/(m2 + n2))^(n - 1))/(m2 + n2)^2 - (m2*max(m2/(m2 + n2),0)^(1/n - 1)*(a2 - 1))/(n*(m2 + n2)^2));
        eqn6=(a2*(1-m2/(m2+n2)).^n+(1-a2)*(1-norm((m2/(m2+n2))^(1/n))))*(Tl-Ts)+Ts-T;
    end


    %constraint 7
    eqn7=v1/(m1+n1+v1)-(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-(T-800)*dSdT/100;
    Jac7=jacobian(eqn7, Variables);
    % eqn7=matlabFunction(eqn7,'Vars', [Variables Constants]);

    %constraint 8
    eqn8a=v2/(m2+n2+v2)-S_cap;
    Jac8a=jacobian(eqn8a, Variables);
    % eqn8a=matlabFunction(eqn8a,'Vars', [Variables Constants]);

    eqn8b=v2*(m1+n1+v1)-D1*v1*(m2+n2+v2);
    Jac8b=jacobian(eqn8b, Variables);
    % eqn8b=matlabFunction(eqn8b,'Vars', [Variables Constants]);

    %constraint 9
    eqn9=cap_A*m2+cap_B*n2-(S_cap-cap_por)*(m2+n2);
    Jac9=jacobian(eqn9, Variables);
    % eqn9=matlabFunction(eqn9,'Vars', [Variables Constants]);

elseif N_component==5
    syms m1 n1 v1 m2 n2 v2 v3 T S_cap l1 l2 l3 Dcl
    Variables=[m1 n1 v1 m2 n2 v2 v3 T S_cap l1 l2 l3 Dcl ];
    %constants and coefficients
    syms M N V CL H cp1 cp2 cp3 Lf1 Lf2 P3 cap_A cap_B a1 b1 c1 D1 a2 n Tl Ts dSdT Coea Coeb Coec Coec2 D2 
    Constants=[M N V CL H cp1 cp2 cp3 Lf1 Lf2 P3 cap_A cap_B a1 b1 c1 D1 a2 n Tl Ts dSdT Coea Coeb Coec Coec2 D2];

    %constraints 1_4
    eqn1=(m1+m2-M)/1000;
    Jac1=jacobian(eqn1, Variables);
    % eqn1=matlabFunction(eqn1,'Vars', [Variables Constants]);

    eqn2=(n1+n2-N)/1000;
    Jac2=jacobian(eqn2, Variables);
    % eqn2=matlabFunction(eqn2,'Vars', [Variables Constants]);

    eqn3=(v1+v2+v3-V)/100;
    Jac3=jacobian(eqn3, Variables);
    % eqn3=matlabFunction(eqn3,'Vars', [Variables Constants]);

    eqn10=(l1+l2+l3-CL)/100;
    Jac10=jacobian(eqn10, Variables);

    %constraint 4
    eqn4=(((m1+n1+v1+l1)*cp1+(m2+n2+v2+l2)*cp2+(v3+l3)*cp3)*T+(m1*Lf1+n1*Lf2)-H)/1e9;
    Jac4 = jacobian(eqn4, Variables);
    % eqn4=matlabFunction(eqn4,'Vars', [Variables Constants]);
    %constraint 5,6

    eqn5=a1*m1^2+b1*m1*(m1+n1)+(c1-T)*(m1+n1)^2;
    Jac5= jacobian(eqn5, Variables);
    % eqn5=matlabFunction(eqn5,'Vars', [Variables Constants]);

    if simplify_Ts==1
        eqn6=((1-m2/(m2+n2)).^n)*(Tl-Ts)+Ts-T;
    else
        eqn6=(a2*(1-m2/(m2+n2)).^n+(1-a2)*(1-(m2/(m2+n2))^(1/n)))*(Tl-Ts)+Ts-T;
    end
    Jac6= jacobian(eqn6, Variables);

    % adjust ment
    if simplify_Ts==0
        Jac6(4)=-(Tl - Ts)*((max(m2/(m2 + n2),0)^(1/n - 1)*(m2/(m2 + n2)^2 - 1/(m2 + n2))*(a2 - 1))/n - a2*n*(m2/(m2 + n2)^2 - 1/(m2 + n2))*(1 - m2/(m2 + n2))^(n - 1));
        Jac6(5)=(Tl - Ts)*((a2*m2*n*(1 - m2/(m2 + n2))^(n - 1))/(m2 + n2)^2 - (m2*max(m2/(m2 + n2),0)^(1/n - 1)*(a2 - 1))/(n*(m2 + n2)^2));
        eqn6=(a2*(1-m2/(m2+n2)).^n+(1-a2)*(1-norm((m2/(m2+n2))^(1/n))))*(Tl-Ts)+Ts-T;
    end

    %constraint 7
    eqn7=v1/(m1+n1+v1+l1)-(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-(T-800)*dSdT/100;
    Jac7=jacobian(eqn7, Variables);
    % eqn7=matlabFunction(eqn7,'Vars', [Variables Constants]);

    %constraint 8
    eqn8a=v2/(m2+n2+v2+l2)-S_cap;
    Jac8a=jacobian(eqn8a, Variables);
    % eqn8a=matlabFunction(eqn8a,'Vars', [Variables Constants]);

    eqn8b=v2*(m1+n1+v1+l1)-D1*v1*(m2+n2+v2+l2);
    Jac8b=jacobian(eqn8b, Variables);
    % eqn8b=matlabFunction(eqn8b,'Vars', [Variables Constants]);

    %constraint 9
    eqn9=cap_A*m2+cap_B*n2-S_cap*(m2+n2);
    Jac9=jacobian(eqn9, Variables);
    % eqn9=matlabFunction(eqn9,'Vars', [Variables Constants]);

    %constraint 11
%11. Cl s/m partitioning: Kcl*l1*(m2+n2+v2+l2)-l2*(m1+n1+v1+l1)=0;
%     eqn11=Coea*l3*(m1+n1+v1+l1)+(Coeb*m1+Coec*n1)*(v3+l3)+(Coec2-log(Dcl))*(v3+l3)*(m1+n1+v1+l1);
    eqn11=exp(Coea*l3/(v3+l3)+(Coeb*m1+Coec*n1)/(m1+n1+v1+l1))*Coec2-Dcl;
    Jac11=jacobian(eqn11, Variables);

    %constraint 12  
    eqn12=D2*l1*(m2+n2+v2+l2)-l2*(m1+n1+v1+l1);
    Jac12=jacobian(eqn12, Variables);

    %constraint 13
    eqn13=(l3/(v3+l3)*(m1+n1+v1+l1)-Dcl*l1)/(M+N+V+CL);
    Jac13=jacobian(eqn13, Variables);
end

%% Between solidus and liquidus, 
if is_solid_solution==0  %Eutectic
    % solid unsaturated, melt saturated 
    % m2=0; m1=M;
    Jacob_SL_s0_l1=[Jac2; Jac3; Jac4; Jac5; Jac7; Jac8b; Jac9]; %order 7
    Jacob_SL_s0_l1(:,[1,4])=[];
    Jacob_SL_s0_l1 = matlabFunction(Jacob_SL_s0_l1, 'Vars', [Variables Constants]);
    rhs_SL_s0_l1=[eqn2; eqn3;eqn4;eqn5;eqn7;eqn8b;eqn9];
    rhs_SL_s0_l1=matlabFunction(rhs_SL_s0_l1, 'Vars', [Variables Constants]);
    
    
    % solid unsaturated, melt unsaturated
    % m2=0; m1=M; v3=0;
    Jacob_SL_s0_l0=[Jac2; Jac3; Jac4; Jac5; Jac8b; Jac9]; %order 6
    Jacob_SL_s0_l0(:,[1,4,7])=[];
    Jacob_SL_s0_l0 = matlabFunction(Jacob_SL_s0_l0, 'Vars', [Variables Constants]);
    rhs_SL_s0_l0=[eqn2; eqn3;eqn4;eqn5;eqn8b;eqn9];
    rhs_SL_s0_l0=matlabFunction(rhs_SL_s0_l0, 'Vars', [Variables Constants]);
    
    % solid saturated, melt saturated 
    % m2=0; m1=M;
    Jacob_SL_s1_l1=[Jac2; Jac3; Jac4; Jac5; Jac7; Jac8a; Jac9]; %order 7
    Jacob_SL_s1_l1(:,[1,4])=[];
    Jacob_SL_s1_l1 = matlabFunction(Jacob_SL_s1_l1, 'Vars', [Variables Constants]);
    rhs_SL_s1_l1=[eqn2; eqn3;eqn4;eqn5;eqn7;eqn8a;eqn9];
    rhs_SL_s1_l1=matlabFunction(rhs_SL_s1_l1, 'Vars', [Variables Constants]);
    
    % solid saturated, melt unsaturated
    % m2=0; m1=M; v3=0;
    Jacob_SL_s1_l0=[Jac2; Jac3; Jac4; Jac5; Jac8a; Jac9]; %order 6
    Jacob_SL_s1_l0(:,[1,4,7])=[];
    Jacob_SL_s1_l0 = matlabFunction(Jacob_SL_s1_l0, 'Vars', [Variables Constants]);
    rhs_SL_s1_l0=[eqn2; eqn3;eqn4;eqn5;eqn8a;eqn9];
    rhs_SL_s1_l0=matlabFunction(rhs_SL_s1_l0, 'Vars', [Variables Constants]);
else %Solid solution
    if N_component==3
        % solid unsaturated, melt saturated
        Jacob_SL_s0_l1=[Jac1; Jac2; Jac3; Jac4; Jac5; Jac6; Jac7; Jac8b; Jac9]; %order 9
        %     Jacob_SL_s0_l1(:,[])=[];
        Jacob_SL_s0_l1 = matlabFunction(Jacob_SL_s0_l1, 'Vars', [Variables Constants]);
        rhs_SL_s0_l1=[eqn1; eqn2; eqn3;eqn4;eqn5; eqn6; eqn7;eqn8b;eqn9];
        rhs_SL_s0_l1=matlabFunction(rhs_SL_s0_l1, 'Vars', [Variables Constants]);

        % solid unsaturated, melt unsaturated
        % v3=0;
        Jacob_SL_s0_l0=[Jac1; Jac2; Jac3; Jac4; Jac5; Jac6; Jac8b; Jac9]; %order 8
        Jacob_SL_s0_l0(:,7)=[];
        Jacob_SL_s0_l0 = matlabFunction(Jacob_SL_s0_l0, 'Vars', [Variables Constants]);
        rhs_SL_s0_l0=[eqn1; eqn2; eqn3;eqn4;eqn5; eqn6; eqn8b;eqn9];
        rhs_SL_s0_l0=matlabFunction(rhs_SL_s0_l0, 'Vars', [Variables Constants]);

        % solid saturated, melt saturated
        Jacob_SL_s1_l1=[Jac1; Jac2; Jac3; Jac4; Jac5; Jac6; Jac7; Jac8a; Jac9]; %order 9
        %     Jacob_SL_s1_l1(:,[])=[];
        Jacob_SL_s1_l1 = matlabFunction(Jacob_SL_s1_l1, 'Vars', [Variables Constants]);
        rhs_SL_s1_l1=[eqn1; eqn2; eqn3;eqn4;eqn5;eqn6;eqn7;eqn8a;eqn9];
        rhs_SL_s1_l1=matlabFunction(rhs_SL_s1_l1, 'Vars', [Variables Constants]);

        % solid saturated, melt unsaturated
        % v3=0;
        Jacob_SL_s1_l0=[Jac1; Jac2; Jac3; Jac4; Jac5; Jac6; Jac8a; Jac9]; %order 8
        Jacob_SL_s1_l0(:,7)=[];
        Jacob_SL_s1_l0 = matlabFunction(Jacob_SL_s1_l0, 'Vars', [Variables Constants]);
        rhs_SL_s1_l0=[eqn1;eqn2; eqn3;eqn4;eqn5;eqn6;eqn8a;eqn9];
        rhs_SL_s1_l0=matlabFunction(rhs_SL_s1_l0, 'Vars', [Variables Constants]);
    elseif N_component==5
        % solid unsaturated, melt saturated
        Jacob_SL_s0_l1=[Jac1; Jac2; Jac3; Jac10; Jac4; Jac5; Jac6; Jac7; Jac8b; Jac9; Jac11; Jac12; Jac13]; %order 13
        %     Jacob_SL_s0_l1(:,[])=[];
        Jacob_SL_s0_l1 = matlabFunction(Jacob_SL_s0_l1, 'Vars', [Variables Constants]);
        rhs_SL_s0_l1=[eqn1; eqn2; eqn3; eqn10; eqn4;eqn5; eqn6; eqn7;eqn8b;eqn9;eqn11;eqn12;eqn13];
        rhs_SL_s0_l1=matlabFunction(rhs_SL_s0_l1, 'Vars', [Variables Constants]);

        % solid unsaturated, melt unsaturated
        % v3=0;
        Jacob_SL_s0_l0=[Jac1; Jac2; Jac3; Jac4; Jac10; Jac5; Jac6; Jac8b; Jac9;  Jac12]; %order 10
        Jacob_SL_s0_l0(:,[7 12 13])=[];
        Jacob_SL_s0_l0 = matlabFunction(Jacob_SL_s0_l0, 'Vars', [Variables Constants]);
        rhs_SL_s0_l0=[eqn1; eqn2; eqn3;eqn4;eqn10; eqn5; eqn6; eqn8b;eqn9; eqn12];
        rhs_SL_s0_l0=matlabFunction(rhs_SL_s0_l0, 'Vars', [Variables Constants]);

        % solid saturated, melt saturated
        Jacob_SL_s1_l1=[Jac1; Jac2; Jac3; Jac4; Jac10; Jac5; Jac6; Jac7; Jac8a; Jac9; Jac11;Jac12;Jac13]; %order 13
        %     Jacob_SL_s1_l1(:,[])=[];
        Jacob_SL_s1_l1 = matlabFunction(Jacob_SL_s1_l1, 'Vars', [Variables Constants]);
        rhs_SL_s1_l1=[eqn1; eqn2; eqn3;eqn4;eqn10;eqn5;eqn6;eqn7;eqn8a;eqn9;eqn11;eqn12;eqn13];
        rhs_SL_s1_l1=matlabFunction(rhs_SL_s1_l1, 'Vars', [Variables Constants]);

        % solid saturated, melt unsaturated
        % v3=0;
        Jacob_SL_s1_l0=[Jac1; Jac2; Jac3; Jac4; Jac10; Jac5; Jac6; Jac8a; Jac9; Jac12]; %order 10
        Jacob_SL_s1_l0(:,[7 12 13])=[];
        Jacob_SL_s1_l0 = matlabFunction(Jacob_SL_s1_l0, 'Vars', [Variables Constants]);
        rhs_SL_s1_l0=[eqn1;eqn2; eqn3;eqn4;eqn10;eqn5;eqn6;eqn8a;eqn9;eqn12];
        rhs_SL_s1_l0=matlabFunction(rhs_SL_s1_l0, 'Vars', [Variables Constants]);
    end
end
%% At solidus
% solid unsaturated, melt saturated 
%n1=0; n2=N; T=Ts
if is_solid_solution==0
    Jacob_S_s0_l1=[Jac1; Jac3; Jac4; Jac7; Jac8b; Jac9]; %order 6
    Jacob_S_s0_l1(:,[2,5,8])=[];
    Jacob_S_s0_l1=matlabFunction(Jacob_S_s0_l1, 'Vars', [Variables Constants]);
    rhs_S_s0_l1=[eqn1; eqn3;eqn4;eqn7;eqn8b;eqn9];
    rhs_S_s0_l1=matlabFunction(rhs_S_s0_l1, 'Vars', [Variables Constants]);
    
    
    % solid unsaturated, melt unsaturated
    % m2=0; m1=M; v3=0; T=Ts;
    Jacob_S_s0_l0=[Jac1; Jac3; Jac4; Jac8b; Jac9]; %order 5
    Jacob_S_s0_l0(:,[2,5,7,8])=[];
    Jacob_S_s0_l0=matlabFunction(Jacob_S_s0_l0, 'Vars', [Variables Constants]);
    rhs_S_s0_l0=[eqn1; eqn3;eqn4;eqn8b;eqn9];
    rhs_S_s0_l0=matlabFunction(rhs_S_s0_l0, 'Vars', [Variables Constants]);
    
    % solid saturated, melt saturated 
    %n1=0; n2=N; T=Ts
    Jacob_S_s1_l1=[Jac1; Jac3; Jac4; Jac7; Jac8a; Jac9]; %order 6
    Jacob_S_s1_l1(:,[2,5,8])=[];
    Jacob_S_s1_l1=matlabFunction(Jacob_S_s1_l1, 'Vars', [Variables Constants]);
    rhs_S_s1_l1=[eqn1; eqn3;eqn4;eqn7;eqn8a;eqn9];
    rhs_S_s1_l1=matlabFunction(rhs_S_s1_l1, 'Vars', [Variables Constants]);
    
    % solid saturated, melt unsaturated
    % m2=0; m1=M; v3=0; T=Ts;
    Jacob_S_s1_l0=[Jac1; Jac3; Jac4; Jac8a; Jac9]; %order 5
    Jacob_S_s1_l0(:,[2,5,7,8])=[];
    Jacob_S_s1_l0=matlabFunction(Jacob_S_s1_l0, 'Vars', [Variables Constants]);
    rhs_S_s1_l0=[eqn1; eqn3;eqn4;eqn8a;eqn9];
    rhs_S_s1_l0=matlabFunction(rhs_S_s1_l0, 'Vars', [Variables Constants]);
end

%%
if is_solid_solution==0
    Jacobians={Jacob_SL_s0_l0, Jacob_SL_s1_l0, Jacob_SL_s0_l1, Jacob_SL_s1_l1, Jacob_S_s0_l0, Jacob_S_s1_l0, Jacob_S_s0_l1, Jacob_S_s1_l1};
    Rhs={rhs_SL_s0_l0, rhs_SL_s1_l0, rhs_SL_s0_l1, rhs_SL_s1_l1, rhs_S_s0_l0,rhs_S_s1_l0,rhs_S_s0_l1,rhs_S_s1_l1 };
    Constant_index={[1,4,7],[1,4,7],[1,4], [1,4], [2,5,7,8], [2,5,7,8],[2,5,8],[2,5,8] };
else    
    Jacobians={Jacob_SL_s0_l0, Jacob_SL_s1_l0, Jacob_SL_s0_l1, Jacob_SL_s1_l1};
    Rhs={rhs_SL_s0_l0, rhs_SL_s1_l0, rhs_SL_s0_l1, rhs_SL_s1_l1};
    if N_component==3
        Constant_index={[7],[7],[], []};
    elseif N_component==5
        Constant_index={[7 12 13],[7 12 13],[], []};
    end
end

% m1=in(1); n1=in(2); v1=in(3); m2=in(4); n2=in(5); v2=in(6); v3=in(7); T=in(8); S_cap=in(9); l1=in(10); l2=in(11); l3=in(12); Dcl=in(13); P3=in(24); dSdT=in(35);

Extra_func=cell(1,4);
if N_component==5
    syms D0P S_cap D2
    %5p below solidus, normal version
    eqn=[l2*(v3+l3)*20-l3*(M+N+v2+l2)*D2;  l2+l3-CL; v2+v3-V; v2/(M+N+v2+l2)-S_cap];
    Jac=jacobian(eqn, [v2 v3 l2 l3]);
    Extra_func{1}=matlabFunction(Jac,'Vars', [v2 v3 l2 l3 M N V CL Coea Coeb S_cap D0P D2]);
    Extra_func{2}=matlabFunction(eqn,'Vars', [v2 v3 l2 l3 M N V CL Coea Coeb S_cap D0P D2]);

    %5p above liquidus
    eqn=[l1*(v3+l3)*exp(Coea*l3/(v3+l3)+(Coeb*M+Coec*N)/(M+N+v1+l1))*Coec2-l3*(M+N+v1+l1);  l1+l3-CL; v1+v3-V; v1/(M+N+v1+l1)-(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-(T-800)*dSdT/100; (((M+N+v1+l1)*cp1+(v3+l3)*cp3)*T+(M*Lf1+N*Lf2)-H)/1e9];
    Jac=jacobian(eqn, [v1 v3 l1 l3 T]);
    Extra_func{3}=matlabFunction(Jac, 'Vars', [v1 v3 l1 l3 T M N V CL H Coea Coeb Coec Coec2 P3 dSdT D0P cp1 cp2 cp3 Lf1 Lf2]);
    Extra_func{4}=matlabFunction(eqn,'Vars',  [v1 v3 l1 l3 T M N V CL H Coea Coeb Coec Coec2 P3 dSdT D0P cp1 cp2 cp3 Lf1 Lf2]);
end

if N_component==3 && simplified_5p==1
    % simplified version, which treats Cl and Cu as decoupled from the rest three components    
    syms l1 l2 l3 u1 u2 u3 CL CU Dcl Dcu D2 D3 Coea Coeb Coec D0P
    % all three phases exist    
    eqn=[l1+l2+l3-CL; u1+u2+u3-CU; l1/(m1+n1+v1+l1+u1)*D2-l2/(m2+n2+v2+l2+u2); u1/(m1+n1+v1+l1+u1)*D3-u2/(m2+n2+v2+l2+u2); l3/(v3+l3+u3)-l1/(m1+n1+v1+l1+u1)*Dcl; u3/(v3+l3+u3)-u1/(m1+n1+v1+l1+u1)*Dcu;   exp(Coea*l3/(v3+l3+u3)+(Coeb*m1+Coec*n1)/(m1+n1+v1+l1+u1))*D0P-Dcl; 8e4*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))^2*(1-((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55)))^14+380*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))+0.8-Dcu];
    Jac=jacobian(eqn,[l1 l2 l3 u1 u2 u3 Dcl Dcu]);
    Extra_func{1}=matlabFunction(Jac,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    Extra_func{2}=matlabFunction(eqn,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    % Solid and volatile
    eqn=[l2+l3-CL; u2+u3-CU; l2/(m2+n2+v2+l2+u2)*20-l3/(v3+l3+u3)*D2; u2/(m2+n2+v2+l2+u2)*100-u3/(v3+l3+u3)*D3];
    Jac=jacobian(eqn,[l2 l3 u2 u3]);
    Extra_func{3}=matlabFunction(Jac,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    Extra_func{4}=matlabFunction(eqn,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    % melt and volatile
    eqn=[l1+l3-CL; u1+u3-CU; l1/(m1+n1+v1+l1+u1)*Dcl-l3/(v3+l3+u3); u1/(m1+n1+v1+l1+u1)*Dcu-u3/(v3+l3+u3); exp(Coea*l3/(v3+l3+u3)+(Coeb*m1+Coec*n1)/(m1+n1+v1+l1+u1))*D0P-Dcl;8e4*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))^2*(1-((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55)))^14+380*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))+0.8-Dcu];
    Jac=jacobian(eqn,[l1 l3 u1 u3 Dcl Dcu]);
    Extra_func{5}=matlabFunction(Jac,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    Extra_func{6}=matlabFunction(eqn,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    % solid and melt
    eqn=[l1+l2-CL; u1+u2-CU; l2/(m2+n2+v2+l2+u2)-l1/(m1+n1+v1+l1+u1)*D2; u2/(m2+n2+v2+l2+u2)-u1/(m1+n1+v1+l1+u1)*D3];
    Jac=jacobian(eqn,[l1 l2 u1 u2]);
    Extra_func{7}=matlabFunction(Jac,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
    Extra_func{8}=matlabFunction(eqn,'Vars', [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 Dcl Dcu CL CU D2 D3 Coea Coeb Coec D0P]);
end

if N_component==3 && add_CLCUS==1
    % simplified version, which treats Cl, Cu and S decoupled from the rest three components  
    % For Cl
    syms l1 l2 l3 u1 u2 u3 u4 s10 s11 s2 s3 CL CU Sul Dcl Dcu D2 D3 D32 D4 Coea Coeb Coec D0P S6_ratio SCSS SCAS
    % all three phases exist    
    eqn=[l1+l2+l3-CL; l1/(m1+n1+v1+l1)*D2-l2/(m2+n2+v2+l2);  l3/(v3+l3)-l1/(m1+n1+v1+l1)*Dcl; exp(Coea*l3/(v3+l3)+(Coeb*m1+Coec*n1)/(m1+n1+v1+l1))*D0P-Dcl];
    
    sys_variables=[m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3 u4 s10 s11 s2 s3 Dcl Dcu CL CU Sul D2 D3 D32 D4 Coea Coeb Coec D0P S6_ratio SCSS SCAS];

    Jac=jacobian(eqn,[l1 l2 l3 Dcl]);
    Extra_func{1}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{2}=matlabFunction(eqn,'Vars', sys_variables);
    % Solid and volatile
    eqn=[l2+l3-CL; l2/(m2+n2+v2+l2)*20-l3/(v3+l3)*D2];
    Jac=jacobian(eqn,[l2 l3]);
    Extra_func{3}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{4}=matlabFunction(eqn,'Vars', sys_variables);
    % melt and volatile
    eqn=[l1+l3-CL; l1/(m1+n1+v1+l1)*Dcl-l3/(v3+l3); exp(Coea*l3/(v3+l3)+(Coeb*m1+Coec*n1)/(m1+n1+v1+l1))*D0P-Dcl];
    Jac=jacobian(eqn,[l1 l3 Dcl]);
    Extra_func{5}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{6}=matlabFunction(eqn,'Vars', sys_variables);
    % solid and melt
    eqn=[l1+l2-CL; l2/(m2+n2+v2+l2)-l1/(m1+n1+v1+l1)*D2];
    Jac=jacobian(eqn,[l1 l2]);
    Extra_func{7}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{8}=matlabFunction(eqn,'Vars', sys_variables);

    % For Cu
    % all three phases exist. D3 is now solfide phase/melt partitioning instead of solid/melt partition
    % sulfide phase exist, s2>0
    eqn=[u1+u2+u3+u4-CU; u1/(m1+n1+v1+l1+u1)*D3-u2/(s2/0.3495); u1/(m1+n1+v1+l1+u1)*D32-u4/(m2+n2+v2+l2+u4) ;u3/(v3+l3+u3)-u1/(m1+n1+v1+l1+u1)*Dcu;   8e4*(((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))^2*(1-((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55)))^14*(1+(s3/34.0809)/(l3/58.44+v3/18.02+u3/63.55+s3/34.0809)))+380*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))+0.8-Dcu];
    Jac=jacobian(eqn,[u1 u2 u3 u4 Dcu]);
    Extra_func{9}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{10}=matlabFunction(eqn,'Vars', sys_variables);
    % sulfide phase not exit, s2=0
    eqn=[u1+u2+u3+u4-CU; u1/(m1+n1+v1+l1+u1)*D32-u4/(m2+n2+v2+l2+u4) ;u3/(v3+l3+u3)-u1/(m1+n1+v1+l1+u1)*Dcu;   8e4*(((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))^2*(1-((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55)))^14*(1+(s3/34.0809)/(l3/58.44+v3/18.02+u3/63.55+s3/34.0809)))+380*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))+0.8-Dcu];
    Jac=jacobian(eqn,[u1 u3 u4 Dcu]);
    Extra_func{101}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{102}=matlabFunction(eqn,'Vars', sys_variables);

    % Solid and volatile
    eqn=[u2+u3-CU; u2/(m2+n2+v2+l2+u2)*100-u3/(v3+l3+u3)*D3];
    Jac=jacobian(eqn,[u1 u2 u3]);
    Extra_func{11}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{12}=matlabFunction(eqn,'Vars', sys_variables);
    % melt and volatile
    eqn=[u1+u3-CU; u1/(m1+n1+v1+l1+u1)*Dcu-u3/(v3+l3+u3); 8e4*(((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))^2*(1-((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55)))^14*(1+(s3/34.0809)/(l3/58.44+v3/18.02+u3/63.55+s3/34.0809)))+380*((l3/58.44)/(l3/58.44+v3/18.02+u3/63.55))+0.8-Dcu];
    Jac=jacobian(eqn,[u1 u2 u3 Dcu]);
    Extra_func{13}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{14}=matlabFunction(eqn,'Vars', sys_variables);
    % solid and melt
    eqn=[u1+u2-CU; u2/(m2+n2+v2+l2+u2)-u1/(m1+n1+v1+l1+u1)*D3];
    Jac=jacobian(eqn,[u1 u2]);
    Extra_func{15}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{16}=matlabFunction(eqn,'Vars', sys_variables);

    %For Sulfur S
    % all three phases exist. 
    % D3 is now solfide phase/melt partitioning instead of solid/melt partition 
    % sulfide phase exist
    eqn=[s10+s11+s2+s3-Sul; s10+s11-SCSS*(m1+n1+v1+l1); s11/(s11+s10)-S6_ratio; s3/(v3+l3+s3)-(s10+s11)/(m1+n1+v1+l1+s10+s11)*D4];
    Jac=jacobian(eqn,[s10 s11 s2 s3]);
    Extra_func{17}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{18}=matlabFunction(eqn,'Vars', sys_variables);
    % no sulfide or sulfate phase, s2=0
    eqn=[s10+s11+s3-Sul;  s11/(s11+s10)-S6_ratio; s3/(v3+l3+s3)-(s10+s11)/(m1+n1+v1+l1+s10+s11)*D4];
    Jac=jacobian(eqn,[s10 s11 s3]);
    Extra_func{19}=matlabFunction(Jac,'Vars', sys_variables);
    Extra_func{20}=matlabFunction(eqn,'Vars', sys_variables);
end

end