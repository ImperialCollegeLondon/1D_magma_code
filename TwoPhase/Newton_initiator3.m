%% A general generator for Jacobians used in the 1D magma model
%% HH 2026-03-13
function [Jac_mom, rhs_mom, Jac_con, rhs_con, Jac_ct, rhs_ct, Jac_ent, rhs_ent, Jac_solidus, rhs_solidus, Jac_liquidus, rhs_liquidus, Jac_ct2,rhs_ct2, Jac_ssat, rhs_ssat, Jac_lsat, rhs_lsat]=Newton_initiator...
    (Has_volatile)
% system valriables contain:
% solid and melt velocities: us, uf
% conservative quantities: melt fraction (phi),  temperature (T), solid major composition (Cs), melt major composition (Cl)
% if volatile is included, it contains volatile fractoin (S), solid volatile volatile composition (Cs2), melt volatile composition (Cl2).
% all parameters (e.g. viscosities, densities) should be functions of phi, T, Cs, Cl, S, Cs2, Cl2 and corresponding linear look-up table be established beofore assembly

% cell values at cell i-1, i, i+1 are labled as 0, 1, 2, respectively
% nodal values (only velocity) at nodes i-1, i, i+1 are labled as 0,1, respectively
% node 1 and 2 are nodes of cell 1; node 0 and 1 are nodes of cell 0.


syms phi_0 phi_1 phi_2 T_0 T_1 T_2 cs_0 cs_1 cs_2 cl_0 cl_1 cl_2
% volatile related variables
syms S_0 S_1 S_2 cs2_0 cs2_1 cs2_2 cl2_0 cl2_1 cl2_2

syms umi_0 umi_1 umi_2  
ufi = sym('ufi');
syms ufi2

% cell lengths
syms dzi_0 dzi_1 dzi_2

% if Has_volatile==1
%     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2; S_0; S_1; S_2; cs2_0; cs2_1; cs2_2; cl2_0; cl2_1; cl2_2];
% else
%     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
% end

%% Parameters linear coefficients
% linear, bylinear, trilinear function coefficients
% linear: y=ax_1+bx
% bilinear: y= ax_1+bx_2+cx_1*x_2+d
% trilinear: y= ax_1+ bx_2 +c x_3 + dx_1*x_2 +ex_1*x_3+ fx_2*x_3+gx_1*x_2*x_3+h
% syms a0 b0 c0 d0 e0 f0 g0 h0
% syms a1 b1 c1 d1 e1 c1 g1 h1


%% Momentum equation
% parameters for momentum equation
% Solid viscosity is a function of phi
syms a0 b0 a1 b1
syms a2 b2 c2 d2 e2 f2 g2 h2
mus_0=a0*phi_0+b0;
mus_1=a1*phi_1+b1;

% melt viscosity is a function of Cl and Cl2 and used only in coupling term,
% so we directly model coupling term C
% coupling coefficient is a function of muf and phi
if Has_volatile==2
    Coupling=a2*cl_1+b2*cl2_1+c2*phi_1+d2*cl_1*cl2_1+e2*cl_1*phi_1+f2*cl2_1*phi_1+g2*cl_1*cl2_1*phi_1+h2;
else
    Coupling=a2*cl_1+b2*phi_1+c2*cl_1*phi_1+d2;
end


% densities are function of compositions
syms a3 b3 a4 b4 a5 b5 c5 d5 a6 b6 c6 d6
rhos_0=a3*cs_0+b3;
rhos_1=a4*cs_1+b4;
if Has_volatile==2
    rhol_0=a5*cl_0+b5*cl2_0+c5*cl_0*cl2_0+d5;
    rhol_1=a6*cl_1+b6*cl2_1+c6*cl_1*cl2_1+d6;   
else
    rhol_0=a5*cl_0+b5;
    rhol_1=a6*cl_1+b6;  
end

% density difference at node i
syms g mum_0
drhog=(rhos_0+rhos_1-rhol_0-rhol_1)/2*g;

eps=1e-12;
momentum=(-((umi_2-umi_1)/dzi_1*(1-phi_1)*mus_1-(umi_1-umi_0)/dzi_0*(1-phi_0)*mus_0)*2/(dzi_0+dzi_1) *2*phi_0*phi_1/(phi_0+phi_1+eps)...
    +phi_0*phi_1/(phi_0+phi_1+eps)*(2-phi_1-phi_0)*drhog-Coupling*(ufi-umi_1))/mum_0;




if Has_volatile==2
    Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
else
    Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1];
end
parameters=[a0; b0; a1; b1; a2; b2; c2; d2; e2; f2; g2; h2; a3; b3; a4; b4; a5; b5; c5; d5; a6; b6; c6; d6;  dzi_0; dzi_1; g; mum_0];  
Jac_mom=jacobian(momentum, Variables);
Jac_mom=simplify(Jac_mom);
Jac_mom= matlabFunction(Jac_mom, 'Vars', {[Variables; parameters]});
rhs_mom= matlabFunction(momentum, 'Vars', {[Variables; parameters]});


%% Continuity equation
continuity=umi_1*(1-phi_1/2-phi_0/2)+ufi*(phi_1/2+phi_0/2);
Variables=[umi_1; ufi; phi_0; phi_1];
Jac_con=jacobian(continuity, Variables);
Jac_con=simplify(Jac_con);
Jac_con= matlabFunction(Jac_con, 'Vars', {Variables});
rhs_con= matlabFunction(continuity, 'Vars', {Variables});


%% Major component transport equation
syms dt
% old bulk composition
syms OLD_com

Variables=[umi_1; umi_2; ufi; ufi2; phi_0; phi_1; phi_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
trans_comp=(phi_1*cl_1+(1-phi_1)*cs_1)-OLD_com-(ufi*(phi_0*cl_0+phi_1*cl_1)/2-ufi2*(phi_1*cl_1+phi_2*cl_2)/2+umi_1*((1-phi_0)/2*cs_0+(1-phi_1)/2*cs_1)-umi_2*((1-phi_1)/2*cs_1+(1-phi_2)/2*cs_2))/dzi_1*dt;
Jac_ct=jacobian(trans_comp, Variables);
Jac_ct= matlabFunction(Jac_ct, 'Vars', {[Variables; dt; dzi_1; OLD_com]});
rhs_ct= matlabFunction(trans_comp, 'Vars', {[Variables; dt; dzi_1; OLD_com]});


%% Enthalpy transport equation
syms cp Lf kt
% old enthalpy
syms OLD_ent
trans_enthalpy=((cp*T_1+Lf*phi_1)-OLD_ent-Lf*((phi_0+phi_1)/2*ufi-(phi_1+phi_2)/2*ufi2)/dzi_1*dt-kt*2/dzi_1*((T_2-T_1)/(dzi_2+dzi_1)-(T_1-T_0)/(dzi_1+dzi_0))*dt)/20/Lf;
Variables=[ufi; ufi2; phi_0; phi_1; phi_2; T_0; T_1; T_2];
Jac_ent=jacobian(trans_enthalpy, Variables);
Jac_ent=simplify(Jac_ent);
Jac_ent=matlabFunction(Jac_ent,'Vars',{[Variables; dzi_0; dzi_1; dzi_2; dt; cp; Lf; kt; OLD_ent]});
rhs_ent= matlabFunction(trans_enthalpy, 'Vars', {[Variables; dzi_0; dzi_1; dzi_2; dt; cp; Lf; kt; OLD_ent]});




%% solidus and liquidus
% solidus and liquidus should be defined in temrs of normalized temperature T'=(T-Ts)/(Tl-Ts) where Tl and Ts is fixed when no volatile is present
% as system parameter which can becomes a variable if volatile component is included in the system, a typical liquidus:
% T=func(cl)*(Tl-Ts)+Ts
% syms Ts Tl
% syms D1 % wather s/m parition coefficient
% if Has_volatile==1
%     % when volatile is present Tl and Ts become functions of bulk water content (thus cb2=cl2*phi+(1-phi)*cs2)
%     % presure dependency is constant for each time step.
%     cb2=cl2_1*phi_1+cs2_1*(1-phi_1)+S_1;
%     Tl=a0*cb2/D1+b0;
%     Ts=a1*cb2/D1+b1;
%     Variables=[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1; ];
% else
%     Tl=a0;
%     Ts=a1;
%     Variables=[phi_1; T_1; cs_1; cl_1];
% end
% T_scaled=(T_1-Ts)/(Tl-Ts);
% 
% syms i2 j2 k2 l2 m2 n2 o2 p2
% % the most general solidus constraint, quadra-linear form G(phi,T_scaled, Cs,Cl)=0
% solidus= a2*phi_1 + b2*T_scaled + c2*cs_1 + d2*cl_1 ...
%     + e2*phi_1.*T_scaled + f2*phi_1.*cs_1 + g2*phi_1.*cl_1 ...
%     + h2*T_scaled.*cs_1 + i2*T_scaled.*cl_1 + j2*cs_1.*cl_1 ...
%     + k2*phi_1.*T_scaled.*cs_1 + l2*phi_1.*T_scaled.*cl_1 ...
%     + m2*phi_1.*cs_1.*cl_1 + n2*T_scaled.*cs_1.*cl_1 ...
%     + o2*phi_1.*T_scaled.*cs_1.*cl_1 + p2;
% 
% 
% 
% Parameters=[a0; b0; a1; b1; ...
%     a2; b2; c2; d2; e2; f2;g2; h2; i2; j2; k2; l2; m2; n2; o2; p2;...
%     D1];
% Jac_solidus=jacobian(solidus/100, Variables);
% Jac_solidus=matlabFunction(Jac_solidus,'Vars',{[Variables; Parameters]});
% rhs_solidus=matlabFunction(solidus,'Vars',{[Variables; Parameters]});
% 
% % In the most general form, liquidus and solidus are of the same format symbolicly 
% Jac_liquidus=Jac_solidus;
% rhs_liquidus=rhs_solidus;
% 
syms A1 n_order D1
eps=1e-12;

if Has_volatile==1
    % when volatile is present Tl and Ts become functions of bulk water content (thus cb2=cl2*phi+(1-phi)*cs2)
    % presure dependency is constant for each time step.
    cb2=cl2_1*phi_1+cs2_1*(1-phi_1)+S_1;
    Ts=a0*cb2/D1/0.13+b0;
    Tl=a1*cb2/D1/0.2+b1;
    Variables=[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1];
else
    Ts=a0;
    Tl=b0;
    Variables=[phi_1; T_1; cs_1; cl_1];
end

C1= Tl;
B1= Ts-A1-Tl;
Parameters=[a0;b0;a1;b1; D1];
Cb=phi_1*cl_1+(1-phi_1)*cs_1;


MINT=((T_1+C1+50)-sqrt((T_1-C1-50)^2+eps))/2;
Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps))/2;
Constrain5=(Cb+SMIN+sqrt((SMIN-Cb)^2+eps))/2;


liquidus=(Constrain5-cl_1)/1e4;
Jac_liquidus=jacobian(liquidus, Variables);
% Jac_liquidus=simplify(Jac_liquidus);
Jac_liquidus=matlabFunction(Jac_liquidus, 'Vars',{[Variables; Parameters; A1]});
rhs_liquidus=matlabFunction(liquidus,'Vars',{[Variables; Parameters; A1]});

Cond4=((Tl-T_1)/(Tl-Ts))^n_order;
SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;


solidus=(Constrain4-cs_1)/1e4;
Jac_solidus=jacobian(solidus, Variables);
% Jac_solidus=simplify(Jac_solidus);
Jac_solidus=matlabFunction(Jac_solidus,'Vars',{[Variables; Parameters; n_order]});
rhs_solidus=matlabFunction(solidus,'Vars',{[Variables; Parameters; n_order]});
%% Volatile (H2O) component transport equation
if Has_volatile==1
    % old volatile bulk composition
    syms OLD_com2
    syms kf
    Variables=[umi_1; umi_2; ufi; ufi2; phi_0; phi_1; phi_2; S_0; S_1; S_2; cs2_0; cs2_1; cs2_2; cl2_0; cl2_1; cl2_2];
    trans_comp2=((phi_1*cl2_1+(1-phi_1)*cs2_1+S_1)-OLD_com2-(ufi*(phi_0*cl2_0+phi_1*cl2_1)/2-ufi2*(phi_1*cl2_1+phi_2*cl2_2)/2+umi_1*((1-phi_0)/2*cs2_0+(1-phi_1)/2*cs2_1)-umi_2*((1-phi_1)/2*cs2_1+(1-phi_2)/2*cs2_2))/dzi_1*dt...
        -kf*2/dzi_1*((S_2-S_1)/(dzi_2+dzi_1)-(S_1-S_0)/(dzi_1+dzi_0))*dt)/20;
    % trans_comp2=((cl2_1+cs2_1+S_1)-OLD_com2-(ufi*(cl2_0+cl2_1)/2-ufi2*(cl2_1+cl2_2)/2+umi_1*(cs2_0+cs2_1)/2-umi_2*(cs2_1+cs2_2)/2)/dzi_1*dt...
    %     -kf*2/dzi_1*((S_2-S_1)/(dzi_2+dzi_1)-(S_1-S_0)/(dzi_1+dzi_0))*dt)/20;
    Jac_ct2=jacobian(trans_comp2, Variables);
    Jac_ct2=simplify(Jac_ct2);
    Jac_ct2= matlabFunction(Jac_ct2, 'Vars',     {[Variables; dzi_0; dzi_1; dzi_2; dt; kf; OLD_com2]});
    rhs_ct2= matlabFunction(trans_comp2, 'Vars', {[Variables; dzi_0; dzi_1; dzi_2; dt; kf; OLD_com2]});

else
    Jac_ct2=[];
    rhs_ct2=[];
end

%% Melt saturation equation
% melt saturation is a function of temperature and presusre. Pressure dependency part is a constant in the constraint
% e.g. Sat=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100+(T-800)*dSdT/100-v2;
eps=1e-12;
% beta=1e4;
if Has_volatile==1
    % cl2_1_cor=(cl2_1+sqrt(cl2_1^2+eps))/2;
    T_capped=(T_1+500+sqrt((T_1-500)^2+eps))/2;
    Sat=a0*T_capped+b0-cl2_1;
    % Sat=sqrt(Sat^2+phi_1^2+eps)-(Sat+phi_1);
    % condition=(sqrt(Sat^2+S_1^2+eps)-(Sat+S_1))/1e4; %condition=0 means either Sat>0 and S=0;  or Sat=0 and S>0    
    % condition=(Sat+S_1-sqrt((Sat-S_1)^2+eps))/1e4;
    condition=(Sat+S_1-(sqrt(Sat-S_1)^2+eps))/1e4;
    % condition=Sat+S_1-sqrt((Sat-S_1)^2+eps);

    Variables=[phi_1; T_1; S_1; cl2_1 ];
    Jac_lsat=jacobian(condition, Variables);
    Jac_lsat= matlabFunction(Jac_lsat, 'Vars', {[Variables; a0; b0]});
    rhs_lsat= matlabFunction(condition, 'Vars', {[Variables; a0; b0]});
else
    Jac_lsat=[];
    rhs_lsat=[];
end

%% Solid saturation equation
% Solid saturation is a function of solid component cs
% either water in solid/melt follow the partition rule or solid is saturated
syms cap_A cap_B
eps=1e-12;
% beta=1e5;
if Has_volatile==1    
    Sat=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1+0.02;
    Par=cl2_1*D1-cs2_1;
    % Par=-1/beta*log(exp(-beta*phi_1^2)+exp(-beta*Par^2)); %log softmin
    % Par=(Par+sqrt(Par^2+eps))/2;
    % Par=sqrt(Par^2+phi_1^2+eps)-(Par+phi_1);
    % condition=(sqrt(Sat^2+Par^2+eps)-(Sat+Par))/1e4; %condition=0 means either Sat>0 and Par=0;  or Sat>0 and Par=0
    condition=(Sat+Par-sqrt((Sat-Par)^2+eps))/1e4;

    Variables=[phi_1; cs_1; cs2_1; cl2_1];
    Jac_ssat=jacobian(condition, Variables);
    Jac_ssat= matlabFunction(Jac_ssat, 'Vars', {[Variables; cap_A; cap_B; D1]});
    rhs_ssat= matlabFunction(condition, 'Vars', {[Variables; cap_A; cap_B; D1]});
else
    Jac_ssat=[];
    rhs_ssat=[];
end
%% Melt and solid saturation, and volatile partitioning constraints
% Constraints:
% 1: S(H(phi)(Satl-Cl2)+(1-H(phi))(Sats-Cs2))=0
% 2: phi*(1-phi)*(Sats-Cs2)(Cs2-D1*Cl2)=0;
% 2a:  phi*(1-phi)*(min(Sats,D1*Cl2)-Cs2)=0;
% syms cap_A cap_B
% 
% tao=1e-6;
% if Has_volatile
%     Variables=[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1];
% 
%     eps=1e-15;
%     Sats=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1/(1-phi_1+eps);
%     T_capped=(T_1+500+sqrt((T_1-500)^2+eps))/2;
%     Satl=a0*T_capped+b0-cl2_1/(phi_1+eps);
% 
%     Par_l=cl2_1/(phi_1+eps)*D1;
%     Par=cs2_1/(1-phi_1+eps)-((cap_A*cs_1+cap_B*(1-cs_1))+Par_l-sqrt(((cap_A*cs_1+cap_B*(1-cs_1))-Par_l+eps)^2))/2;
% 
%     R=(phi_1/(phi_1+eps)*Satl+(1-phi_1/(phi_1+eps))*Sats);
%     constraint=sqrt(S_1^2+R^2+eps)-(S_1+R);
%     Jac_lsat=jacobian(constraint, Variables);
%     Jac_lsat=simplify(Jac_lsat);
%     Jac_lsat= matlabFunction(Jac_lsat, 'Vars', {[Variables; a0; b0; cap_A; cap_B; D1]});
%     rhs_lsat= matlabFunction(constraint, 'Vars', {[Variables; a0; b0; cap_A; cap_B; D1]});
% 
% 
%     constraint=-tao*log(exp(-phi_1^2/tao)+exp(-(1-phi_1)^2/tao)+exp(-Par^2/tao));
%     Jac_ssat=jacobian(constraint, Variables);
%     Jac_ssat=simplify(Jac_ssat);
%     Jac_ssat= matlabFunction(Jac_ssat, 'Vars', {[Variables; a0; b0; cap_A; cap_B; D1]});
%     rhs_ssat= matlabFunction(constraint, 'Vars', {[Variables; a0; b0; cap_A; cap_B; D1]});
% else
%         Jac_lsat=[];
%         rhs_lsat=[];
%         Jac_ssat=[];
%         rhs_ssat=[];
% end
%%


save('Jacobians_for_Newton.mat', 'Jac_mom', 'rhs_mom', 'Jac_con', 'rhs_con', 'Jac_ct', 'rhs_ct', 'Jac_ent', 'rhs_ent', 'Jac_solidus', 'rhs_solidus', 'Jac_liquidus', 'rhs_liquidus', 'Jac_ct2', 'rhs_ct2', 'Jac_ssat','rhs_ssat','Jac_lsat','rhs_lsat');
end