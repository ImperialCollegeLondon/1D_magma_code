%% A general generator for Jacobians used in the 1D magma model
%% HH 2026-03-13
function [Jac_mom, rhs_mom, Jac_con, rhs_con, Jac_ct, rhs_ct, Jac_ent, rhs_ent, Jac_solidus, rhs_solidus, Jac_liquidus, rhs_liquidus, Jac_ct2,rhs_ct2, Jac_ssat, rhs_ssat, Jac_lsat, rhs_lsat, Jac_ct_central , rhs_ct_central, Jac_con_central, rhs_con_central]=Newton_initiator2...
    (Has_volatile, is_eutectic)
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

eps=1e-3;
eps_phi01=3e-2;
momentum=(-((umi_2-umi_1)/dzi_1*(1-phi_1+eps_phi01)*mus_1-(umi_1-umi_0)/dzi_0*(1-phi_0+eps_phi01)*mus_0)*2/(dzi_0+dzi_1) *0.5*(phi_0+phi_1+eps_phi01)...
    +phi_0*phi_1/(phi_0+phi_1+eps)*(2-phi_1-phi_0)*drhog-Coupling*(ufi-umi_1))/mum_0;




% if Has_volatile==2
%     Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
% else
%     Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1];
% end

% parameters=[a0; b0; a1; b1; a2; b2; c2; d2; e2; f2; g2; h2; a3; b3; a4; b4; a5; b5; c5; d5; a6; b6; c6; d6;  dzi_0; dzi_1; g; mum_0];  
% Jac_mom=jacobian(momentum, Variables);
% Jac_mom=simplify(Jac_mom);
% Jac_mom= matlabFunction(Jac_mom, 'Vars', {[Variables; parameters]});
% rhs_mom= matlabFunction(momentum, 'Vars', {[Variables; parameters]});

Jac_mom=jacobian(momentum, [umi_0, umi_1, umi_2, ufi, phi_0, phi_1, cs_0, cs_1, cl_0, cl_1]);
Jac_mom=simplify(Jac_mom);
Jac_mom= matlabFunction(Jac_mom, 'Vars', {umi_0, umi_1, umi_2, ufi, phi_0, phi_1, cs_0, cs_1, cl_0, cl_1, a0, b0, a1, b1, a2, b2, c2, d2, e2, f2, g2, h2, a3, b3, a4, b4, a5, b5, c5, d5, a6, b6, c6, d6,  dzi_0, dzi_1, g, mum_0});
rhs_mom= matlabFunction(momentum, 'Vars', {umi_0, umi_1, umi_2, ufi, phi_0, phi_1, cs_0, cs_1, cl_0, cl_1, a0, b0, a1, b1, a2, b2, c2, d2, e2, f2, g2, h2, a3, b3, a4, b4, a5, b5, c5, d5, a6, b6, c6, d6,  dzi_0, dzi_1, g, mum_0});


%% Continuity equation
syms con_scaling
eps=1e-6;

constraint_mean_por=(phi_1+phi_0)/2;

% constraint_mean_por=(constraint_mean_por+Cap-sqrt((constraint_mean_por-Cap)^2+eps))/2;
% constraint_mean_por=(constraint_mean_por-Cap2-sqrt((constraint_mean_por-Cap2)^2+eps))/2;
% continuity=(umi_1*(1-constraint_mean_por)+ufi*constraint_mean_por)*con_scaling;  %Need to be enlarged
continuity=(umi_1*(1-phi_1)+ufi*phi_0)*con_scaling;  %Need to be enlarged
% Variables=;

continuity=simplify(continuity);
Jac_con=jacobian(continuity, [umi_1, ufi, phi_0, phi_1]);
Jac_con=simplify(Jac_con);

Jac_con= matlabFunction(Jac_con, 'Vars', {umi_1, ufi, phi_0, phi_1, con_scaling});
rhs_con= matlabFunction(continuity, 'Vars', {umi_1, ufi, phi_0, phi_1, con_scaling});

%
continuity=(umi_1*(1-constraint_mean_por)+ufi*constraint_mean_por)*con_scaling;
continuity=simplify(continuity);
Jac_con_central=jacobian(continuity, [umi_1, ufi, phi_0, phi_1]);
Jac_con_central=simplify(Jac_con_central);

Jac_con_central= matlabFunction(Jac_con_central, 'Vars', {umi_1, ufi, phi_0, phi_1, con_scaling});
rhs_con_central= matlabFunction(continuity, 'Vars', {umi_1, ufi, phi_0, phi_1, con_scaling});
%% Major component transport equation
syms dt
% old bulk composition
syms OLD_com k_stable1 k_stable2

cb_0=phi_0*cl_0+(1-phi_0)*cs_0;
cb_1=phi_1*cl_1+(1-phi_1)*cs_1;
cb_2=phi_2*cl_2+(1-phi_2)*cs_2;



%----------------------------------------
% RAW FLUXES (as in your code)
% %----------------------------------------
F_m = ...
    ufi*(phi_0*cl_0 + phi_1*cl_1)/2 + ...
    umi_1*((1-phi_0)/2*cs_0 + (1-phi_1)/2*cs_1);

F_p = ...
    ufi2*(phi_1*cl_1 + phi_2*cl_2)/2 + ...
    umi_2*((1-phi_1)/2*cs_1 + (1-phi_2)/2*cs_2);


% F_m=ufi*phi_0*cl_0+umi_1*(1-phi_1)*cs_1;
% F_p=ufi2*phi_1*cl_1+umi_2*(1-phi_2)*cs_2;

% K_cut=800;
% alpha_p=1./(1+exp(K_cut*(-cb_1+1e-2)));
% alpha_m=1./(1+exp(K_cut*(-cb_0+1e-2)));
%----------------------------------------
% APPLY CONTROL ONLY TO OUTGOING PART
%----------------------------------------
F_m_ctrl =  F_m ;
F_p_ctrl =  F_p ;

%----------------------------------------
% TRANSPORT EQUATION (MODIFIED)
%----------------------------------------
% trans_comp = cb_1 - OLD_com - (F_m_ctrl - F_p_ctrl)/dzi_1 * dt;
trans_comp = cb_1 - OLD_com - (F_m_ctrl - F_p_ctrl)/dzi_1 * dt-2/dzi_1*(k_stable2*(cb_2-cb_1)/(dzi_2+dzi_1)-k_stable1*(cb_1-cb_0)/(dzi_1+dzi_0))*dt;
trans_comp = simplify(trans_comp);


Jac_ct=jacobian(trans_comp, [umi_1, umi_2, ufi, ufi2, phi_0, phi_1, phi_2, cs_0, cs_1, cs_2, cl_0, cl_1, cl_2]);
Jac_ct = simplify(Jac_ct);

Jac_ct_central = matlabFunction(Jac_ct, 'Vars', ...
   {umi_1, umi_2, ufi, ufi2, ...
    phi_0, phi_1, phi_2, ...
    cs_0, cs_1, cs_2, ...
    cl_0, cl_1, cl_2, ... %cb_1    
    dt, dzi_0, dzi_1, dzi_2, OLD_com, k_stable1, k_stable2});

rhs_ct_central = matlabFunction(trans_comp, 'Vars', ...
   {umi_1, umi_2, ufi, ufi2, ...
    phi_0, phi_1, phi_2, ...
    cs_0, cs_1, cs_2, ...
    cl_0, cl_1, cl_2, ... %alpha_m, alpha_p, ...
    dt, dzi_0, dzi_1, dzi_2, OLD_com, k_stable1, k_stable2});


%%%%%%%%%%%%%%
F_m=ufi*phi_0*cl_0+umi_1*(1-phi_1)*cs_1;
F_p=ufi2*phi_1*cl_1+umi_2*(1-phi_2)*cs_2;

%----------------------------------------
% APPLY CONTROL ONLY TO OUTGOING PART
%----------------------------------------
F_m_ctrl =  F_m ;
F_p_ctrl =  F_p ;

%----------------------------------------
% TRANSPORT EQUATION (MODIFIED)
%----------------------------------------
% trans_comp = cb_1 - OLD_com - (F_m_ctrl - F_p_ctrl)/dzi_1 * dt;
trans_comp = cb_1 - OLD_com - (F_m_ctrl - F_p_ctrl)/dzi_1 * dt-2/dzi_1*(k_stable2*(cb_2-cb_1)/(dzi_2+dzi_1)-k_stable1*(cb_1-cb_0)/(dzi_1+dzi_0))*dt;
trans_comp = simplify(trans_comp);


Jac_ct=jacobian(trans_comp, [umi_1, umi_2, ufi, ufi2, phi_0, phi_1, phi_2, cs_0, cs_1, cs_2, cl_0, cl_1, cl_2]);
Jac_ct = simplify(Jac_ct);

Jac_ct= matlabFunction(Jac_ct, 'Vars', ...
   {umi_1, umi_2, ufi, ufi2, ...
    phi_0, phi_1, phi_2, ...
    cs_0, cs_1, cs_2, ...
    cl_0, cl_1, cl_2, ... %cb_1    
    dt, dzi_0, dzi_1, dzi_2, OLD_com, k_stable1, k_stable2});

rhs_ct = matlabFunction(trans_comp, 'Vars', ...
   {umi_1, umi_2, ufi, ufi2, ...
    phi_0, phi_1, phi_2, ...
    cs_0, cs_1, cs_2, ...
    cl_0, cl_1, cl_2, ... %alpha_m, alpha_p, ...
    dt, dzi_0, dzi_1, dzi_2, OLD_com, k_stable1, k_stable2});



%% Enthalpy transport equation
syms cp Lf kt
% old enthalpy
syms OLD_ent H_scaling
trans_enthalpy=((cp*T_1+Lf*phi_1)-OLD_ent-Lf*((phi_0+phi_1)/2*ufi-(phi_1+phi_2)/2*ufi2)/dzi_1*dt-kt*2/dzi_1*((T_2-T_1)/(dzi_2+dzi_1)-(T_1-T_0)/(dzi_1+dzi_0))*dt)/Lf*H_scaling;

% Variables=[ufi; ufi2; phi_0; phi_1; phi_2; T_0; T_1; T_2];
% Jac_ent=jacobian(trans_enthalpy, Variables);
% Jac_ent=simplify(Jac_ent);
% Jac_ent=matlabFunction(Jac_ent,'Vars',{[Variables; dzi_0; dzi_1; dzi_2; dt; cp; Lf; kt; OLD_ent]});
% rhs_ent= matlabFunction(trans_enthalpy, 'Vars', {[Variables; dzi_0; dzi_1; dzi_2; dt; cp; Lf; kt; OLD_ent]});

Jac_ent=jacobian(trans_enthalpy, [ufi, ufi2, phi_0, phi_1, phi_2, T_0, T_1, T_2]);
Jac_ent=simplify(Jac_ent);
Jac_ent=matlabFunction(Jac_ent,'Vars',{ufi, ufi2, phi_0, phi_1, phi_2, T_0, T_1, T_2, dzi_0, dzi_1, dzi_2, dt, cp, Lf, kt, OLD_ent, H_scaling});
rhs_ent= matlabFunction(trans_enthalpy, 'Vars', {ufi, ufi2, phi_0, phi_1, phi_2, T_0, T_1, T_2, dzi_0, dzi_1, dzi_2, dt, cp, Lf, kt, OLD_ent, H_scaling});
%% solidus and liquidus
% solidus and liquidus should be defined in temrs of normalized temperature T'=(T-Ts)/(Tl-Ts) where Tl and Ts is fixed when no volatile is present
% as system parameter which can becomes a variable if volatile component is included in the system, a typical liquidus:
% T=func(cl)*(Tl-Ts)+Ts
syms A1 n_order D1



if Has_volatile==1
    % when volatile is present Tl and Ts become functions of bulk water content (thus cb2=cl2*phi+(1-phi)*cs2)
    % presure dependency is constant for each time step.
    cb2=cl2_1*phi_1+cs2_1*(1-phi_1)+S_1;
    Ts=a0*cb2/D1/0.13+b0;
    Tl=a1*cb2/D1/0.2+b1;
    % Variables=[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1];
else
    Ts=a0;
    Tl=b0;
    % Variables=[phi_1; T_1; cs_1; cl_1];
end

C1= Tl;
B1= Ts-A1-Tl;
% Parameters=[a0;b0;a1;b1; D1];
Cb=phi_1*cl_1+(1-phi_1)*cs_1;



eps=1e-3;
eps2=1e-8;

MINT=((T_1+C1)-sqrt((T_1-C1)^2+eps))/2;
Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps2))/2;
Constrain5=(Cb+SMIN+sqrt((SMIN-Cb)^2+eps2))/2;

% MINT=-1/K*log(exp(-K*T_1/1000)+exp(-K*(C1+50)/1000))*1000;
% Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
% SMIN=-1/K*log(exp(-K*Cond5)+exp(-K*1));
% Constrain5=1/K*log(exp(K*Cb)+exp(K*SMIN));



liquidus=(Constrain5-cl_1)/1e4;

% Jac_liquidus=jacobian(liquidus, Variables);
% % Jac_liquidus=simplify(Jac_liquidus);
% Jac_liquidus=matlabFunction(Jac_liquidus, 'Vars',{[Variables; Parameters; A1]});
% rhs_liquidus=matlabFunction(liquidus,'Vars',{[Variables; Parameters; A1]});

if Has_volatile==1 
    Jac_liquidus=jacobian(liquidus, [phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1]);
    Jac_liquidus=matlabFunction(Jac_liquidus, 'Vars',[phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, a0,b0,a1,b1,D1, A1]);
    rhs_liquidus=matlabFunction(liquidus,'Vars',     [phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, a0,b0,a1,b1,D1, A1]);
else
    Jac_liquidus=jacobian(liquidus, [phi_1, T_1, cs_1, cl_1]);
    Jac_liquidus=matlabFunction(Jac_liquidus, 'Vars',[phi_1, T_1, cs_1, cl_1,       a0,b0,a1,b1,D1, A1]);
    rhs_liquidus=matlabFunction(liquidus,'Vars',     [phi_1, T_1, cs_1, cl_1,       a0,b0,a1,b1,D1, A1]);
end


if is_eutectic==1
    eps2=1e-0;
    solidus=(Cb/2*(1-tanh((T_1-Ts)/eps2))-cs_1)/1e4;
else 
    % K=2e3;
    % cs_min=0;
    eps=1e-8;
    Cond4=((Tl-T_1)/(Tl-Ts))^n_order;
    SMIN1=(Cb+Cond4-sqrt((Cb-Cond4)^2+eps))/2;
    Constrain4=(SMIN1+sqrt(SMIN1^2+eps))/2;

    % SMIN1=-1/K*log(exp(-K*Cb)+exp(-K*Cond4));
    % Constrain4=1/K*log(exp(K*SMIN1)+exp(K*cs_min));
    solidus=(Constrain4-cs_1)/1e4;
end

% Jac_solidus=jacobian(solidus, Variables);
% % Jac_solidus=simplify(Jac_solidus);
% Jac_solidus=matlabFunction(Jac_solidus,'Vars',{[Variables; Parameters; n_order]});
% rhs_solidus=matlabFunction(solidus,'Vars',{[Variables; Parameters; n_order]});

if Has_volatile==1
    Jac_solidus=jacobian(solidus, [phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1]);
    Jac_solidus=matlabFunction(Jac_solidus,'Vars',[phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, a0,b0,a1,b1,D1, n_order]);
    rhs_solidus=matlabFunction(solidus,'Vars',[phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, a0,b0,a1,b1,D1, n_order]);
else
    Jac_solidus=jacobian(solidus, [phi_1, T_1, cs_1, cl_1]);
    Jac_solidus=matlabFunction(Jac_solidus,'Vars',[phi_1, T_1, cs_1, cl_1, a0,b0,a1,b1,D1, n_order]);
    rhs_solidus=matlabFunction(solidus,'Vars',[phi_1, T_1, cs_1, cl_1, a0,b0,a1,b1,D1, n_order]);
end
%% Volatile (H2O) component transport equation
if Has_volatile==1
    % old volatile bulk composition
    syms OLD_com2
    syms kf1 kf2 kf_stable1 kf_stable2
    % Variables=[umi_1; umi_2; ufi; ufi2; phi_0; phi_1; phi_2; S_0; S_1; S_2; cs2_0; cs2_1; cs2_2; cl2_0; cl2_1; cl2_2];
    eps=1e-10;
    
    % Smin=1e-6;
    % dS12=S_2-S_1;
    % eps12=(1e-2*dS12^2+Smin)^2;
    % dS01=S_1-S_0;
    % eps01=(1e-2*dS01^2+Smin)^2;
    
    S12flux=0.5*((S_2-S_1)-sqrt((S_2-S_1)^2+eps));
    S01flux=0.5*((S_1-S_0)-sqrt((S_1-S_0)^2+eps));    
    
    cb2_0=phi_0*cl2_0+(1-phi_0)*cs2_0+S_0;
    cb2_1=phi_1*cl2_1+(1-phi_1)*cs2_1+S_1;
    cb2_2=phi_2*cl2_2+(1-phi_2)*cs2_2+S_2;

    %
    % F_m = ...
    %     ufi*(phi_0*cl2_0 + phi_1*cl2_1)/2 + ...
    %     umi_1*((1-phi_0)/2*cs2_0 + (1-phi_1)/2*cs2_1);
    % 
    % F_p = ...
    %     ufi2*(phi_1*cl2_1 + phi_2*cl2_2)/2 + ...
    %     umi_2*((1-phi_1)/2*cs2_1 + (1-phi_2)/2*cs2_2);
    
    
    F_m=ufi*phi_0*cl2_0+umi_1*(1-phi_1)*cs2_1;
    F_p=ufi2*phi_1*cl2_1+umi_2*(1-phi_2)*cs2_2;
    


    trans_comp2=(cb2_1-OLD_com2- (F_m - F_p)/dzi_1 * dt...
        -2/dzi_1*(kf2*S12flux/(dzi_2+dzi_1)-kf1*S01flux/(dzi_1+dzi_0))*dt-1/dzi_1*(kf_stable2*(cb2_2-cb2_1)/(dzi_2+dzi_1)-kf_stable1*(cb2_1-cb2_0)/(dzi_1+dzi_0))*dt)/1e3; %
    % trans_comp2=((cl2_1+cs2_1+S_1)-OLD_com2-(ufi*(cl2_0+cl2_1)/2-ufi2*(cl2_1+cl2_2)/2+umi_1*(cs2_0+cs2_1)/2-umi_2*(cs2_1+cs2_2)/2)/dzi_1*dt...
    %     -kf*2/dzi_1*((S_2-S_1)/(dzi_2+dzi_1)-(S_1-S_0)/(dzi_1+dzi_0))*dt)/20;

    % Jac_ct2=jacobian(trans_comp2, Variables);
    % Jac_ct2=simplify(Jac_ct2);
    % Jac_ct2= matlabFunction(Jac_ct2, 'Vars',     {[Variables; dzi_0; dzi_1; dzi_2; dt; kf; OLD_com2]});
    % rhs_ct2= matlabFunction(trans_comp2, 'Vars', {[Variables; dzi_0; dzi_1; dzi_2; dt; kf; OLD_com2]});

    trans_comp2=simplify(trans_comp2);
    Jac_ct2=jacobian(trans_comp2, [umi_1, umi_2, ufi, ufi2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2]);
    Jac_ct2=simplify(Jac_ct2);
    Jac_ct2= matlabFunction(Jac_ct2, 'Vars',     {umi_1, umi_2, ufi, ufi2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2, dzi_0, dzi_1, dzi_2, dt, kf1, kf2, OLD_com2, kf_stable1,kf_stable2});
    rhs_ct2= matlabFunction(trans_comp2, 'Vars', {umi_1, umi_2, ufi, ufi2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2, dzi_0, dzi_1, dzi_2, dt, kf1, kf2, OLD_com2, kf_stable1,kf_stable2});
else
    Jac_ct2=[];
    rhs_ct2=[];
end

%% Melt and solid saturation equation
% melt saturation is a function of temperature and presusre. Pressure dependency part is a constant in the constraint
% e.g. Sat=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100+(T-800)*dSdT/100-v2;
syms cap_A cap_B
eps=1e-3;
eps2=1e-8;
% beta=1e4;
if Has_volatile==1
    % cl2_1_cor=(cl2_1+sqrt(cl2_1^2+eps))/2;
    T_capped=(T_1+500+sqrt((T_1-500)^2+eps))/2;
    Satl=a0*T_capped+b0-cl2_1;
    
    % Sats=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1+0.02;
    

    % condition=(Satl+S_1-(sqrt(Satl-S_1)^2+eps))/1e4;
    condition=(Satl+S_1-sqrt((Satl-S_1)^2+eps2))/1e4;

    % Variables=[phi_1; T_1; cs_1; S_1; cs2_1; cl2_1];
    % Jac_lsat=jacobian(condition, Variables);
    % Jac_lsat= matlabFunction(Jac_lsat, 'Vars', {[Variables; a0; b0]}); 
    % rhs_lsat= matlabFunction(condition, 'Vars', {[Variables; a0; b0]});
    Jac_lsat=jacobian(condition, [ T_1,  S_1,  cl2_1]);
    Jac_lsat= matlabFunction(Jac_lsat, 'Vars', {T_1,  S_1,  cl2_1, a0, b0});
    rhs_lsat= matlabFunction(condition, 'Vars',{ T_1,  S_1,  cl2_1, a0, b0});

    % % condition=(Sats+Par-sqrt((Sats-Par)^2+eps))/1e4;
    Par=cl2_1*D1-cs2_1;
    condition=Par/1e4;
    % Variables=[phi_1; T_1; cs_1; S_1; cs2_1; cl2_1];

    % Jac_ssat=jacobian(condition, Variables);
    % Jac_ssat= matlabFunction(Jac_ssat, 'Vars', {[Variables; cap_A; cap_B; D1]});
    % rhs_ssat= matlabFunction(condition, 'Vars', {[Variables; cap_A; cap_B; D1]});
    Jac_ssat=jacobian(condition, [ T_1, cs_1, S_1, cs2_1, cl2_1]);
    Jac_ssat= matlabFunction(Jac_ssat, 'Vars', { T_1, cs_1, S_1, cs2_1, cl2_1, cap_A, cap_B, D1});
    rhs_ssat= matlabFunction(condition, 'Vars', { T_1, cs_1, S_1, cs2_1, cl2_1, cap_A, cap_B, D1});
else
    Jac_lsat=[];
    rhs_lsat=[];
    Jac_ssat=[];
    rhs_ssat=[];
end



save('Jacobians_for_Newton.mat', 'Jac_mom', 'rhs_mom', 'Jac_con', 'rhs_con', 'Jac_ct', 'rhs_ct', 'Jac_ent', 'rhs_ent', 'Jac_solidus', 'rhs_solidus', 'Jac_liquidus', 'rhs_liquidus', 'Jac_ct2', 'rhs_ct2', 'Jac_ssat','rhs_ssat','Jac_lsat','rhs_lsat');
end