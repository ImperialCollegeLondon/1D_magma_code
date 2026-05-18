%cr% HH 2/2/2020
% Input for two phase solver. 
% Settings for a sill intrusion simulation

% If code or comments are written by Catherine Booth, it is annotated CAB

clear all;
fclose('all'); % close all files
%execute='1AA_4M_2_phase_HS_mc_input_update_2.txt';


Inputs= readtable('Input_Files/1AA_2phase_master_input_v6.txt');

Has_volatile=1;


[r,~] = size(Inputs);
names=string(Inputs.Var2);
Number=Inputs.Var3;
for i=1:r
    assignin('base',names(i),Number(i))
end
Precision0=Precision;

if HHJPet==1
    Inputs= readtable('Input_Files/1AA_2phase_HHJPet_master.txt');
    [r,~] = size(Inputs);
    names=string(Inputs.Var2);
    Number=Inputs.Var3;
    for i=1:r
        assignin('base',names(i),Number(i))
    end
    
    is_eutectic=1;
    % set values to zero that are used in sill_intrusion and solid_state
    % but not related to the chosen phase diagram 
    ol=0;
    opx=0;
    cpx=0;
    feld=0;
    ol_mg=0;
    px_mg=0;
    opx_mg=0;
    cpx_mg=0;
    feld_mg=0;
    sol_mn=0;
    ol_mn=0;
    DT=0;
    DC=0;
    lm1=0;
    lc1=0;
    crit_cb4=0;
    crit_cb5=0;
    s5T=0;
    liq_P2_C=0;
    liq_P3_C=0;
    liq_k1=0;
    liq_a1=0;
    liq_b1=0;
    liq_k2=0;
    liq_a2=0;
    liq_b2=0;
    liq_k3=0;
    liq_a3=0;
    liq_b3=0;
    alpha=0;
    n_PD=0;

elseif SSPD==1
    Inputs= readtable('Input_Files/1AA_2phase_SSPD_master.txt');
    is_eutectic=0;
    [r,~] = size(Inputs);
    names=string(Inputs.Var2);
    Number=Inputs.Var3;
    for i=1:r
        assignin('base',names(i),Number(i))
    end
    
    if mod(n_order,2)==0 % make sure n_order is an odd number
        n_order=n_order-1;
    end
    % set values to zero that are used in sill_intrusion and solid_state
    % but not related to the chosen phase diagram 
    ol=0;
    opx=0;
    cpx=0;
    feld=0;
    ol_mg=0;
    px_mg=0;
    opx_mg=0;
    cpx_mg=0;
    feld_mg=0;
    sol_mn=0;
    ol_mn=0;
    DT=0;
    DC=0;
    lm1=0;
    lc1=0;
    crit_cb4=0;
    crit_cb5=0;
    s5T=0;
    liq_P2_C=0;
    liq_P3_C=0;
    liq_k1=0;
    liq_a1=0;
    liq_b1=0;
    liq_k2=0;
    liq_a2=0;
    liq_b2=0;
    liq_k3=0;
    liq_a3=0;
    liq_b3=0;


elseif FourMPD==1
    Inputs= readtable('Input_Files/1AA_2phase_FourMPD_master.txt');
    is_eutectic=2;
    [r,~] = size(Inputs);
    names=string(Inputs.Var2);
    Number=Inputs.Var3;
    for i=1:r
        assignin('base',names(i),Number(i))
    end

    % set values to zero that are used in sill_intrusion and solid_state
    % but not related to the chosen phase diagram 
    A1=0;
    B1=0;
    C1=0;
    alpha=0;
    n_PD=0;

end



if HighTsF==0
    disp('You have chosen NOT to have a higher solidus in the lower crust')
else
    disp('You have chosen to have a higher solidus in the lower crust')
end

if C_type2==1
    disp('You have chosen KC permeability')
else
    disp('You have chosen RG permeability')
end

if perm_min_FLAG==0
    disp('You have chosen to have no minimum effective permeability')
else
    disp('You have chosen to have a minimum effective permeability')
end

if C_type==0
    disp('You have chosen porous only')
else
    disp('You have chosen porous-suspension')
end

if MU_type==0
    disp('You have chosen Jackson et al 2018 mush bulk viscosity forumulation')
else
    disp('You have chosen Costa et al 2019 mush bulk vicosity formulation')
end

if Sill_injection==0
    disp('You have chosen not to intrude sills')
else
    disp('You have chosen to intrude sills')
end

if fixed_dt==0
    disp('Dynamic timesteps: ON')
else
    disp('Dynamic timesteps: OFF')
end

if HHJPet==1
    disp('Using chemical model, solid and melt densities and melt viscosity from Hu et al 2022')
end

if SSPD==1
    disp('Using solid solution chemical model. Solid and melt densities and melt viscosity relationship from Hu et al 2022')
end

if FourMPD==1
    disp('Using layered intrusion chemical model. Solid and melt densities, as well as melt viscosity relationship for layered intrusions.')
end

if ((FourMPD==1)&&(SSPD==1)) || ((FourMPD==1)&& (HHJPet==1)) || ((HHJPet==1)&&(SSPD==1))
    error('More than one chemical model is specified. You should only specify one chemical model')
end

if FourMPD==0 && SSPD==0 && HHJPet==0
    error('You have not specified a chemical model.')
end
    



if SSPD==1 || HHJPet==1
   Range_of_phase_diagram=[min_sio2 max_sio2];
elseif FourMPD ==1
    Range_of_phase_diagram=[min_mgo max_mgo];
end

if SSPD==1 || HHJPet==1
    rho_constant=[rhof_1, rhof_2; rhom_1, rhom_2];
end

Year=3600*24*365.25;

%%
% Thermal data
if HHJPet==1
    TsU=A1+B1+C1; % Solidus is at liquidus of eutectic point
end


kt0=kt0s/DensSill; %3/2850;                %Thermal difusivity of the sill
kt_background=kt0b/DensSill; %3/2850;      %Thermal difusivity of the surrounding 



if C_type2==1
    C_value_A=C_value_A_KC; %125;  % Coefficient for the power low permeability k=1/A*grain_size^(-2)*phi^(-B),    58
    C_value_B=C_value_B_KC;%3;    % Coefficient for the power low permeability,      2.6
else
    C_value_A=grain_size^2/C_value_A_RG_beta;%70;  % beta=70 from Claudio& Alison's talk
    C_value_B=C_value_B_RG;%0;  %not used, but still need to define
end
%C_type=0;  % 1 for porous-suspension type.  0 for porous only 

%plot permeability
phi_test = linspace(0,1,N_vis);


for i=1:N_vis
    C_highmu(i) = c_gen_master(10^mu_f2,grain_size,phi_test(i),C_value_A,C_value_B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF);
    C_lowmu(i) = c_gen_master(10^mu_f1,grain_size,phi_test(i),C_value_A,C_value_B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF);
end



%f4=figure(4);
%semilogy(phi_test,C_highmu,'r')
%hold on
%semilogy(phi_test,C_lowmu,'b')
%xlabel('Melt fration (-)')
%ylabel('Coupling effect between two phases')
%legend({'Max mu f','Min mu f'})
%saveas(f4,'Couplingterm_porous_KC_fig','svg')

%% Phase diagram calcs for when FourMPD is chosen

if FourMPD==1

    % Non dim
    crit_cb1=(crit_cb1_dim - min_mgo)/(max_mgo-min_mgo);
    crit_cb2=(crit_cb2_dim - min_mgo)/(max_mgo-min_mgo);
    crit_cb3=(crit_cb3_dim - min_mgo)/(max_mgo-min_mgo);
    crit_cb4=(crit_cb4_dim - min_mgo)/(max_mgo-min_mgo);
    crit_cb5=(crit_cb5_dim - min_mgo)/(max_mgo-min_mgo);
    liq_P1_C=(liq_P1_C_dim - min_mgo)/(max_mgo-min_mgo);
    liq_P2_C=(liq_P2_C_dim - min_mgo)/(max_mgo-min_mgo);
    liq_P3_C=(liq_P3_C_dim - min_mgo)/(max_mgo-min_mgo);
    liq_P4_C=(liq_P4_C_dim - min_mgo)/(max_mgo-min_mgo);
    mk=mk_dim*(max_mgo-min_mgo);
    ck=ck_dim+mk_dim*min_mgo;
    AC=(AC_dim - min_mgo)/(max_mgo-min_mgo);
    BC=(BC_dim - min_mgo)/(max_mgo-min_mgo);
    CC=(CC_dim - min_mgo)/(max_mgo-min_mgo);
    DC=(DC_dim - min_mgo)/(max_mgo-min_mgo);
    EC=(EC_dim - min_mgo)/(max_mgo-min_mgo);
    FC=(FC_dim - min_mgo)/(max_mgo-min_mgo);
    GC=(GC_dim - min_mgo)/(max_mgo-min_mgo);
    HC=(HC_dim - min_mgo)/(max_mgo-min_mgo);
    IC=(IC_dim - min_mgo)/(max_mgo-min_mgo);
    JC=(JC_dim - min_mgo)/(max_mgo-min_mgo);
    KC=(KC_dim - min_mgo)/(max_mgo-min_mgo);
    LC=(LC_dim - min_mgo)/(max_mgo-min_mgo);
    MC=(MC_dim - min_mgo)/(max_mgo-min_mgo);
    NC=(NC_dim - min_mgo)/(max_mgo-min_mgo);
    OC=(OC_dim - min_mgo)/(max_mgo-min_mgo);
    
    crit_sol_mn_cb = (crit_sol_mn_cb_dim - min_mgo)/(max_mgo-min_mgo);
    
    olpx_C=(olpx_C_dim - min_mgo)/(max_mgo-min_mgo);
    olpxm=olpxm_dim*(max_mgo-min_mgo);
    olpxc=olpxc_dim+olpxm_dim*min_mgo;
    
    copx_C_m = copx_C_m_dim;
    copx_C_c = (copx_C_c_dim - min_mgo + copx_C_m*min_mgo)/(max_mgo-min_mgo);
    
    copx_T_m = copx_T_m_dim*(max_mgo-min_mgo);
    copx_T_c = copx_T_c_dim + copx_T_m_dim*min_mgo;
    
    
    ol_i_m = ol_i_m_dim;
    ol_i_c = (ol_i_m_dim*min_mgo+ol_i_c_dim-min_mgo)/(max_mgo-min_mgo);
    ol_mg_m = ol_mg_m_dim/(max_mgo-min_mgo);
    
    cb_low_high = (cb_low_high_dim - min_mgo)/(max_mgo-min_mgo);
    copx_cb_m_1h = copx_cb_m_l1_h_dim;
    copx_cb_c_1h = (copx_cb_c_l1_h_dim-min_mgo+copx_cb_m_l1_h_dim*min_mgo)/(max_mgo-min_mgo);
    copx_cb_m_1l = copx_cb_m_l1_l_dim;
    copx_cb_c_1l = (copx_cb_c_l1_l_dim - min_mgo + copx_cb_m_l1_l_dim*min_mgo)/(max_mgo-min_mgo);
    
    copx_T_m_1h = copx_T_m_l1_h_dim*(max_mgo-min_mgo);
    copx_T_c_1h = copx_T_c_l1_h_dim + copx_T_m_l1_h_dim*min_mgo; 
    copx_T_m_1l = copx_T_m_l1_l_dim*(max_mgo - min_mgo); 
    copx_T_c_1l = copx_T_c_l1_l_dim + copx_T_m_l1_l_dim*min_mgo;
    
    copx_cb_m_2h = copx_cb_m_l2_h_dim;
    copx_cb_c_2h = (copx_cb_c_l2_h_dim-min_mgo+copx_cb_m_l2_h_dim*min_mgo)/(max_mgo-min_mgo);
    copx_cb_m_2l = copx_cb_m_l2_l_dim;
    copx_cb_c_2l = (copx_cb_c_l2_l_dim - min_mgo + copx_cb_m_l2_l_dim*min_mgo)/(max_mgo-min_mgo);
    
    copx_T_m_2h = copx_T_m_l2_h_dim*(max_mgo-min_mgo);
    copx_T_c_2h = copx_T_c_l2_h_dim + copx_T_m_l2_h_dim*min_mgo; 
    copx_T_m_2l = copx_T_m_l2_l_dim*(max_mgo - min_mgo); 
    copx_T_c_2l = copx_T_c_l2_l_dim + copx_T_m_l2_l_dim*min_mgo;
    
    px_cb_m = px_cb_m_dim/(max_mgo-min_mgo);
    px_cb_c = (px_cb_c_dim-min_mgo)/(max_mgo-min_mgo);
    
    crit_cb4_OG = crit_cb4;
    crit_cb4 = liq_P2_C;
    crit_cb5 = liq_P3_C;
    
    crit_cb4_dim = (crit_cb4)*(max_mgo-min_mgo)+min_mgo;
    crit_cb5_dim = (crit_cb5)*(max_mgo-min_mgo)+min_mgo;
    crit_cb4_OG_dim = (crit_cb4_OG)*(max_mgo-min_mgo)+min_mgo;
    
    
    % Liquid line gradients 
    liq_b1 = ((liq_k1 - liq_P1_T)*liq_P1_C - (liq_k1 - liq_P2_T)*liq_P2_C)/(liq_P2_T - liq_P1_T);
    liq_a1 = (liq_k1 - liq_P1_T)*liq_P1_C - (liq_k1-liq_P1_T)*liq_b1;
    
    liq_b2 = ((liq_k2 - liq_P2_T)*liq_P2_C - (liq_k2 - liq_P3_T)*liq_P3_C)/(liq_P3_T - liq_P2_T);
    liq_a2 = (liq_k2 - liq_P2_T)*liq_P2_C - (liq_k2-liq_P2_T)*liq_b2;
    
    liq_b3 = ((liq_k3 - liq_P3_T)*liq_P3_C - (liq_k3 - liq_P4_T)*liq_P4_C)/(liq_P4_T - liq_P3_T);
    liq_a3 = (liq_k3 - liq_P3_T)*liq_P3_C - (liq_k3-liq_P3_T)*liq_b3;
    
    
    %% Solid line gradients
    %line 1
    lm1 = (AT-BT)/(AC-BC);
    lc1 = BT - lm1*BC;
    
    %line 3
    lm3 = (CT - DT)/(CC - DC);
    lc3 = CT - lm3*CC;
    
    
    %line 8
    lm8 = (JT - KT)/(JC - KC);
    lc8 = liq_k2 +liq_a2/liq_b2; %JT - lm8*JC;
    
    %% Mineral lines
    
    %s1
    s1m = (AC - BC)/(AT - BT);
    s1c = BC - s1m*BT;
    
    %s3
    s3m = (CC - DC)/(CT - DT);
    s3c = CC - s3m*CT;
    
    %s4
    EC = s3m*ET + s3c;
    s4m = (EC - FC)/(ET - FT);
    s4c = EC - s4m*ET;
    
    %s6
    s6m = (GC - HC)/(GT - HT);
    s6c = GC - s6m*GT;
    
    %s7
    s7m = (GC - IC)/(GT - IT);
    s7c = GC - s7m*GT;
    
    %s8
    s8m = 1/lm8;
    s8c = -lc8/lm8;
    
    %s11
    s11m = (LC - MC)/(LT - MT);
    s11c = LC - s11m*LT;
    
    % l13T=(Ts_local+Tl_local)/2;
    l13_start_T = liq_k3 - liq_a3/(crit_cb5 - liq_b3);
    l13_start_cs = s8m*l13_start_T + s8c;
    l13_end_T = OT;
    l13_end_cs = OC;
    
    l13_k = 1125;%sol_k2;%1100;
    l13_b = ((l13_k - l13_start_T)*l13_start_cs - (l13_k-l13_end_T)*l13_end_cs)/(l13_end_T - l13_start_T);
    l13_a = (l13_k - l13_start_T)*l13_start_cs - l13_b*(l13_k - l13_start_T);
    
    % new calcs
    crit_cb3 = (s1m*DT+s1c)+crit_cb2-DC;
    crit_cb3_dim = crit_cb3*(max_mgo-min_mgo)+min_mgo;
    
    %liq_k2 = lc8 - liq_a2/liq_b2;
end
%% Set_up the solid shear and bulk viscosities based on Costa2019
%N_vis=1e6;
phi_range=linspace(0,1,N_vis);
phi_range(phi_range<Precision_BV)=1e-2;
phi_range(phi_range>1-Precision_BV)=1-Precision_BV;


sigma=sigma_max-gamma;
%epsilon=0.5e-4;


%mu_min_FLAG=1;
%mu_min_MF=0.06;
%mu_min_MF_grad = 0.05;
%max_xi_m = 1e20;
if MU_type==0 % MU_type=0 - Nature model, MU_type=1 - Costa 2019 
    mu_m=Ref_Bulk_MN*ones(1,length(phi_range));
    if mu_min_FLAG == 0
        xi_m=mu_m.*(phi_range).^(-0.5);
    else
        for i=1:1:N_vis
            if (phi_range(i))<mu_min_MF_grad
                xi_m(i) = max_xi_m;
            elseif (phi_range(i))<mu_min_MF
                xi_start = mu_m(i).*(mu_min_MF).^(-0.5);
                xi_grad = (max_xi_m - xi_start)/( mu_min_MF_grad - mu_min_MF);
                xi_c = xi_start - xi_grad*mu_min_MF;
                xi_m(i) = xi_grad*(phi_range(i)) + xi_c;
            else
                xi_m(i) = mu_m(i).*(phi_range(i)).^(-0.5);
            end
        end
    end
    mu_all=4/3*mu_m+xi_m;
else

      
    target_fun = @(epsilon) compute_F(0.99, epsilon, phistar, gamma, B_vis, BS_ratio) ...
                      - (ref_shear/ref_shear2) * compute_F(0.5, epsilon, phistar, gamma, B_vis, BS_ratio);

    epsilon0 = 1e-3;   % initial guess
    epsilon = fzero(target_fun, [1e-5 1]);
    F=(1-epsilon)*erf(pi^0.5/2/(1-epsilon)*phi_range/phistar.*(1+(phi_range/phistar).^gamma));

    ref_mu=(1+(phi_range/phistar).^sigma)./(1-F).^(B_vis*phistar);
    ref_xi=ref_mu./max(0.0,phi_range)./(1-phi_range)*BS_ratio;
    ref_all=4/3*ref_mu+ref_xi;
    
    
    mu_all=ref_all*ref_shear/ref_all(end);
    mu_m=ref_mu*ref_shear/ref_all(end);
    xi_m=ref_xi*ref_shear/ref_all(end);
    mu_all=mu_all(end:-1:1);
end




%phi_range=phi_range(end:-1:1);
Font=35;
f3 = figure(3);
clf; set(gcf,'Color','w')
set(gca,'TickDir','out');
semilogy(phi_range,mu_all,'linewidth',2)
hold on
box off
set(gca,'fontsize',Font)
xlabel('Melt fraction (-)', 'FontSize',Font)
ylabel('Solid viscosity (Pa s)', 'FontSize',Font)
line([0.2 0.2], [1 1e20],'linestyle','--', 'color','red','linewidth',1.5)
line([0.6 0.6], [1 1e20],'linestyle','--', 'color','red','linewidth',1.5)
% semilogy(phi_range,mu_m)
% semilogy(phi_range,xi_m, 'x-')

ylim([min(mu_all)/2 max(mu_all)])
% legend({'Sum','Shear','Bulk',})
% saveas(f3,'Shear_bulk_viscosity','svg')

%% Newton parameters initialization
if Use_Newton
    scaling_factor=1;
else
    scaling_factor=sqrt(scaling_factorB^scaling_factorEx);%10^12.5);  %to make the velocity LHS matrix have better conditioning number
end
g=g/scaling_factor;

if Has_volatile==1
    k_stable_value0=1e-10;
    k_stable_value=k_stable_value0;
    H_scaling=1e-2;
else
    k_stable_value0=5e-10;
    k_stable_value=k_stable_value0;
    H_scaling=1e-2;
end
con_scaling=1e5;
%% Defining the sill length and the injection parameters


Injection_depthC=abs(Base_crust-Injection_depth)*1000;

LengthB=abs(Base_crust-Injection_depth)*1000-Sill_length*fine_ratio;   % Below the intrusion depth (metres)
LengthT=abs(Top_crust-Injection_depth)*1000+Sill_length*fine_ratio; %Above the initial intrusion depth (metres)




Show_z=[LengthB LengthB+2*Sill_length*fine_ratio];
min_show_range=Show_z(2)-Show_z(1);



%% define the meshing and mesh adaptivity

part1=linspace(0,LengthB,max(round(min_N/15),5));
part2=linspace(0,2*Sill_length*fine_ratio,min_N);
part3=0:part1(end)-part1(end-1):LengthT-2*Sill_length*fine_ratio;%max(round(min_N/8),5));
nodez=[part1 part1(end)+part2(2:end) part1(end)+part2(end)+part3(2:end)];
N=length(nodez)-1;    


cellz=(nodez(1:end-1)+nodez(2:end))/2;  % cell center points
dz=nodez(2:end)-nodez(1:end-1);
dzF=part2(2)-part2(1);
dzC=part1(2)-part1(1);
%% Sill properties

%Number of nodes per sill
SillNodez=round(Sill_length/(dzF));
%
%
Tot_Sill=Sill_length*Tot_SillN; % Total magma intruded (in metres)

%Intrusion depth - Nodez
depthN = find(nodez<Injection_depthC & nodez>Injection_depthC-dzF);
%Intrusion depth -cellz
depthCN = find(cellz==((nodez(depthN)+nodez(depthN-1))/2));
%Number of nodez/cellz in rand intrusion
rand_injectN=round(Rand_inject/dzF);

%injection parameters
injection_Cb= (Comp_sill-Range_of_phase_diagram(1))/(Range_of_phase_diagram(2)-Range_of_phase_diagram(1)); %bulk composition of injection(s)
injection_phi=MeltFracSill;%injection melt fraction
injection_Cl=injection_Cb./injection_phi; %injection liquid comp
injection_Cs=0.*ones(SillNodez,1); %injection solid comp


% injection T and H needs to be determined after volatile content is determined later



 
if (SSPD==1||HHJPet==1)
    injection_densF =  rhof_1*injection_Cl+(1-injection_Cl)*rhof_2; %injection density - fluid
    injection_densM = rhom_1*injection_Cs+(1-injection_Cs)*rhom_2; %injection density - matrix
elseif FourMPD==1
    if injection_Cl>=crit_mg_si_melt
        injection_densF = rhof_2; %injection density - fluid
    else
        injection_densF = rhof_m*injection_Cl + rhof_c; %injection density - fluid
    end
    injection_densM = rhom_1*injection_Cs+(1-injection_Cs)*rhom_2; %injection density - matrix
end
    
injection_densB = injection_densF.*injection_phi + injection_densM.*(1-injection_phi); %injection density - bulk

%% Setting up time step
%
dt0=dt0Y*Year;
dt_intended=dt0;
Max_dt=Max_dtY*Year;

Min_dt=Min_dtY*Year;

if fixed_dt==0
    dt=Max_dt;
end
End_time=End_timeY*Year;


%% parameters for initial conditions


phi=zeros(N,1);

% Give bulk composition
%Index_Sill=[find(cellz>Sill(1),1,'first') find(cellz<Sill(2),1,'last')];
if HighTsF==1
    DepthTsF=abs(Base_crust-DepthTs)*1000;
    %Soldius depth - Nodez
    depthTsN = find(nodez<DepthTsF & nodez>DepthTsF-dzC);
    %Solidus depth -cellz
    depthTsCN = find(cellz==((nodez(depthTsN)+nodez(depthTsN-1))/2));
else
    DepthTsF=0;
    depthTsN=0;
    depthTsCN=0;
end

Cb_initial=(Cb_initialV-Range_of_phase_diagram(1))/(Range_of_phase_diagram(2)-Range_of_phase_diagram(1)); %crust has 62% SiO2    %0.38;  % Initial sill bulk composition (also used on the back groud rock)
if HighTsF==1
    Cb_initialTs=(Cb_initialVD-Range_of_phase_diagram(1))/(Range_of_phase_diagram(2)-Range_of_phase_diagram(1)); %crust below depth has 50% SiO2
end

Cb=zeros(N,1); %ones(N,1)

if FourMPD==1
    ol=zeros(N,1);
    opx=zeros(N,1);
    cpx=zeros(N,1);
    feld=zeros(N,1);
    ol_mg=zeros(N,1);
    px_mg=zeros(N,1);
    opx_mg=zeros(N,1);
    cpx_mg=zeros(N,1);
    feld_mg=zeros(N,1);
    ol_mn = zeros(N,1);
    sol_mn = zeros(N,1);
end


%% CAB - deciding solidus, liquidus and composition of the crust

    % If a high solidus in the lower crust is wanted
    if HighTsF==1 
        for i=1:1:N

           
        end
        
    %If a high solidus in the lower crust is NOT wanted
    else 
        for i=1:1:N
            if FourMPD==1
                ol(i)=0;
                opx(i)=0;
                cpx(i)=0;
                feld(i)=0;
            end
                
                
                Cb(i)=Cb_initial;
                
                [Ts(i),Tl(i),~,~]=solid_state_Ts_Tl_local_master(HHJPet,SSPD,FourMPD, Cb(i),cp,Lf,DT,DC,lm1,lc1, crit_cb4,crit_cb5,s5T,...
                                                    liq_P2_C, liq_P3_C, liq_k1,liq_a1,liq_b1,liq_k2,liq_a2,liq_b2,liq_k3,liq_a3,liq_b3,A1,B1,C1,alpha, n_PD);
     
        end
    end



%% Setting temperature and corresponding enthalpy

%Geotherm

T(1,1)=T_geo*abs(cellz(1)/1000-Base_crust);

 for i=1:length(cellz)-1
     T(i+1,1)=T(i)-(dz(i)+dz(i+1))/2*T_geo/1000;
 end





Varying_T=0;
  




BC_T_type=[1,1];         %lower and upper boundary condition type for temperature: 1, fixed temperature (not enthalpy)  2, smooth 3. no flux 4.fixed enthalpy
% BC_T=[2*Solidus(0)-1*Liquidus(0.09,A1, B1, C1, A2, B2, C2, ae), 1.5*Solidus(0)-0.5*Liquidus(0.09,A1, B1, C1, A2, B2, C2, ae)]; %lower and upper boundary condition value when fixed. 
BC_T=[T(1), T(end)];

BC_C_type=[3,3];
BC_C_value=[0,0;0,0];

%% Initialize the chemical differentiation number
Lrange1=[6100 6200];  %the range of the studied region 
Lrange2=[6050 6250]; %the range of the calculated region 
CDN_index1=find(cellz>Lrange1(1),1,'first'): find(cellz<Lrange1(2),1,'last');
CDN_index2=find(cellz>Lrange2(1),1,'first'): find(cellz<Lrange2(2),1,'last');
CD0=sum(Cb(CDN_index2).*(cellz(CDN_index2)-Lrange1(1))'.*dz(CDN_index2)');

Cb_all=sum(Cb(CDN_index1).*dz(CDN_index1)');
CD1=0.5*((Lrange1(2)-Lrange1(1))^2-(Lrange1(2)-Lrange1(1)-Cb_all)^2)+sum(Cb(setxor(CDN_index2,CDN_index1)).*(cellz(setxor(CDN_index2,CDN_index1))-Lrange1(1))'.*dz(setxor(CDN_index2,CDN_index1))');
CDN=0;
%% Compute dimensionless compaction length


if HHJPet==1 || SSPD==1
    SiO2_range=Range_of_phase_diagram;
    Ssio2_dim=(Cb_initialV)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1);%si1-si2*tanh(-si3*Cb(1)+si4);
elseif FourMPD==1
    MgO_range=Range_of_phase_diagram;
    SiO2_range=[min_sio2,max_sio2];
    Cb_initialV = Cb_initialV*(MgO_range(2)-MgO_range(1))+MgO_range(1);
    if Cb_initialV>=crit_mg_si_melt
        Ssio2_dim = min_sio2_visc;
    else
        Ssio2_dim = max_sio2_m*Cb_initialV+max_sio2_c;
    end
end

Ssio2_scaled=(Ssio2_dim- SiO2_range(1))/(SiO2_range(2)-SiO2_range(1));

mu_f_dynamics=10.^((mu_f2-mu_f1)*Ssio2_scaled+mu_f1);



% Compaction_length=sqrt(10^mu_m*grain_size^2/C_value_A/10^mu_f2);
% disp(['Compaction Length=',num2str(Compaction_length)])




%% check plots
c_check = linspace(0,1,200);
Ssio2_dim_c=0.*c_check;
rhof_c1=0.*c_check;

% viscosity 
if HHJPet==1 || SSPD==1
    Ssio2_dim_c=(c_check)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1);%si1-si2*tanh(-si3*Cb(1)+si4);
    c_check_dim=Ssio2_dim_c;
elseif FourMPD==1
    c_check_dim = c_check*(MgO_range(2)-MgO_range(1))+MgO_range(1);
    for i=1:200
        if c_check_dim(i)>=crit_mg_si_melt
            Ssio2_dim_c(i) = min_sio2_visc;
        else
            Ssio2_dim_c(i) = max_sio2_m*c_check_dim(i)+max_sio2_c;
        end
    end
end

Ssio2_scaled_c=(Ssio2_dim_c-SiO2_range(1))/(SiO2_range(2)-SiO2_range(1));
mu_f_dynamics_c=10.^((mu_f2-mu_f1)*Ssio2_scaled_c+mu_f1);

% density
if HHJPet==1 || SSPD==1
    rhof_c1=rhof_1*(1-c_check)+(c_check)*rhof_2;
    rhom_c=rhom_1*(1-c_check)+(c_check)*rhom_2;
elseif FourMPD==1
%rhof changes depending on comp
    for i=1:200
        if c_check_dim(i)>=crit_mg_si_melt
            rhof_c1(i) = rhof_2;
        else
            rhof_c1(i) = rhof_m*c_check_dim(i)+rhof_c;
        end
    end
    rhom_c=rhom_1*(1-c_check)+(c_check)*rhom_2;
                  
end
rho_mean=2700;

%f1=figure(1);
%if HHJPet==1 || SSPD==1
%    plot(c_check_dim, mu_f_dynamics_c)
%    xlabel('Melt SiO2 (%)')
%    ylabel('Melt Viscosity (Pas)')
%elseif FourMPD==1
%    subplot(1,3,1)
%    plot(c_check_dim,  Ssio2_dim_c)
%    xlabel('Melt MgO (%)')
%    ylabel('Melt SiO2 (%)')
%    subplot(1,3,2)
%    plot(Ssio2_dim_c, mu_f_dynamics_c)
%    xlabel('Melt SiO2 (%)')
%    ylabel('Melt viscosity (Pas)')
%    subplot(1,3,3)
 %   plot(c_check_dim, mu_f_dynamics_c)
  %  xlabel('Melt MgO (%)')
  %  ylabel('Melt Viscosity (Pas)')
%end
%saveas(f1,'Melt_Viscosity_fig','svg')

%f2=figure(2);
%if HHJPet==1 || SSPD==1
 %   subplot(1,2,1)
  %  plot(c_check_dim, rhof_c1)
   % xlabel('Melt SiO2 (%)')
   % ylabel('Melt density (kgm^{-3}')
   % subplot(1,2,2)
   % plot(c_check_dim, rhom_c)
   
% xlabel('Solid SiO2 (%)')
 %   ylabel('Solid density (kgm^{-3}')
% elseif FourMPD==1
    %subplot(1,2,1)
   % plot(c_check_dim, rhof_c1)
   % xlabel('Melt MgO (%)')
   % ylabel('Melt density (kgm^{-3}')
   % subplot(1,2,2)
   % plot(c_check_dim, rhom_c)
   % xlabel('Solid MgO (%)')
   % ylabel('Solid density (kgm^{-3}')
%end
%saveas(f2,'Density_fig','svg')

%% Data output


Record_data=0; % set to 1 to save the simulation data
Step_counts=0;
Output_Flag=1;
% Record_time=[0:5:240]*Year;
Record_time=[0:Out_time*Year:End_time];%[0:1:50 55:5:230]*Year;  % the time at which the simulation will be saved
% CAB - Plotting of mass error limits
min_cons=0.98;max_cons=1.02;
min_res=0;max_res=1;
% Record_time=[1:5]*Year;


%A compact reusable data ouput
Record_data2=0;
if Record_data2==1
    saveEvery  = 100; %steps contained in each data file 
    if Has_volatile==0
        snapshots = struct('nodez',{},'us', {}, 'ul', {}, 'phi',{}, 'T',{},'Cs', {}, 'Cl',{},                          'Tr',{},'Time',{});
    else
        snapshots = struct('nodez',{},'us', {}, 'ul', {}, 'phi',{}, 'T',{},'Cs', {}, 'Cl',{},'S',{}, 'Cs2',{},'Cl2',{},'Tr',{},'Time',{});
    end
    chunkIdx=0;
    Initial_record_period=20*Year;
    Record_period=Initial_record_period;    
end

%% plotting output
if HHJPet==1 || SSPD==1
    FilenamePlot='1AA_Plot_video_output.txt';
    File_plot_out = fopen(FilenamePlot,'w');
    fprintf(File_plot_out, '%10s %10.3f \n','start_output,',0);
    fprintf(File_plot_out, '%10s %10.3f \n','end_output,',round(End_time/(Out_time*Year)));
    fprintf(File_plot_out, '%10s %10.3f \n','time interval (yr),',Out_time);
    fprintf(File_plot_out, '%10s %10.3f \n','top of plots,',crust_top);
    fprintf(File_plot_out, '%10s %10.3f \n','bottom of plots,',crust_bottom);
    fprintf(File_plot_out, '%10s %10.3f \n','Min SiO2,',min_sio2);
    fprintf(File_plot_out, '%10s %10.3f \n','Max_SiO2,',max_sio2);
    fprintf(File_plot_out, '%10s %10.3f \n','frames per second,',20);
    fprintf(File_plot_out, '%10s %10.3f \n','Contribution outputted,',To_cal_mc);
elseif FourMPD ==1
    FilenamePlot='1AA_Plot_video_output.txt';
    File_plot_out = fopen(FilenamePlot,'w');
    fprintf(File_plot_out, '%10s %10.3f \n','start_output,',0);
    fprintf(File_plot_out, '%10s %10.3f \n','end_output,',round(End_time/(Out_time*Year)));
    fprintf(File_plot_out, '%10s %10.3f \n','time interval (yr),',Out_time);
    fprintf(File_plot_out, '%10s %10.3f \n','top of plots,',crust_top);
    fprintf(File_plot_out, '%10s %10.3f \n','bottom of plots,',crust_bottom);
    fprintf(File_plot_out, '%10s %10.3f \n','Min MgO,',min_mgo);
    fprintf(File_plot_out, '%10s %10.3f \n','Max MgO,',max_mgo);
    fprintf(File_plot_out, '%10s %10.3f \n','Min SiO2,',min_sio2);
    fprintf(File_plot_out, '%10s %10.3f \n','Max_SiO2,',max_sio2);
    fprintf(File_plot_out, '%10s %10.3f \n','frames per second,',20);
    fprintf(File_plot_out, '%10s %10.3f \n','Contribution outputted,',To_cal_mc);
end


%% Calculating the contribution of melt fraction/bulk composition change from heating/compaction/reactive flow.
%To_cal_mc=1;
% Calculating the contribution for the change of melft fraction
if To_cal_mc==1
   Con_T=zeros(N,1);  %contribution from temperature
   Con_H=zeros(N,1);  %contribution from entire enthalpy  
   Con_C=zeros(N,1);  %contribution from composition
   Con_C_dummy=zeros(N,1);
   Con_M=zeros(N,1);  %contribution from compaction

   ACon_T=zeros(N,1); %accumulated contribution from temperature
   ACon_M=zeros(N,1);   
   ACon_C=zeros(N,1); 
   ACon_C_dummy=zeros(N,1);

   CD_M=zeros(N,1); %CD caused by compaction
   CD_R=zeros(N,1);
   CD_R_dummy=zeros(N,1);
end





%% Generate look-up tables for Newton's method

if Use_Newton==1
    disp('Generating all the lookup tables for Newton''s method')
    
    N_mus=N_vis;

    % Solid viscosity are all ready in a lookup table as mu_all, just to create the coefficient
    mus_coef_all = piecewiseFit(mu_all);


    % Coupling term coefficient, as a function of phi, Cl and Cl2 (if there is volatile) 
    N_C=200;
    phi_range=linspace(0,1,N_C);
    Cl_range=linspace(0,1,N_C);
    %calculate the melt viscosity range
    if (HHJPet==1||SSPD==1)
        Ssio2_scaled = Cl_range;
    elseif FourMPD==1
        Ssio2_dim = 0*Cl_range;
        % calculating melt SiO2 from melt MgO
        %dim MgO
        C_all_dim = Cl_range*(MgO_range(2)-MgO_range(1))+MgO_range(1);
        for i=1:N_C
            if C_all_dim(i)>=crit_mg_si_melt
                Ssio2_dim(i) = min_sio2_visc;
            else
                Ssio2_dim(i) = max_sio2_m*C_all_dim(i)+max_sio2_c;
            end
        end
        Ssio2_scaled=(Ssio2_dim- SiO2_range(1))/(SiO2_range(2)-SiO2_range(1));
    else %H2O dependant viscosity from (REFERNCE!) 
        
    end

    if Has_volatile==0
        C_values=zeros(N_C,N_C); %phi, Cl
        muf_range=10.^((mu_f2-mu_f1)*Ssio2_scaled+mu_f1);
        for i=1:N_C
            for j=1:N_C
                C_values(i,j)=c_gen_master(muf_range(j),grain_size,phi_range(i),C_value_A,C_value_B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF);
            end
        end
        C_coef_all= piecewiseFit(C_values);
    else %currently not taking Cl2 into account
        % melt viscosity is a function of both Cl and Cl2
        % vis_b1=  211.2559; vis_b2=-9.8778e+04; vis_b3=  1.4453e+05;
        % vis_c1=0.2519; vis_c2=1.6674e+04; vis_c3=-2.5643e+04;
        % 
        % Cl2_range=linspace(0,1,N_C);
        % SiO2=(Cl_range*74+(1-Cl_range)*47); % SiO2%
        % H2O=Cl2_range*100; % H2O%  !!!
        %     for i=1:N
        %         P3=Pg_real(i)*100;
        %         Saturation=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5);
        %         H2O(i)=H2O(i)/Saturation;
        %     end

        % vis_B=vis_b1*SiO2+vis_b2*H2O+vis_b3*log(1+H2O);
        % vis_C=vis_c1*SiO2+vis_c2*H2O+vis_c3*log(1+H2O);
        % muf_order=-4.55+vis_B./(max(T,500)+273.15-vis_C);

        C_values=zeros(N_C,N_C); %phi, Cl, Cl2
        muf_range=10.^((mu_f2-mu_f1)*Ssio2_scaled+mu_f1);
        for i=1:N_C
            for j=1:N_C
                % for k=1:N_C
                    C_values(i,j)=c_gen_master(muf_range(j),grain_size,phi_range(i),C_value_A,C_value_B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF);
                % end
            end
        end
        C_coef_all= piecewiseFit(C_values);
    end

    
    % mum_0=max(Ref_Bulk_MN,max(max(C_values)));
    mum_0=Ref_Bulk_MN;
    % mum_0=sqrt(Ref_Bulk_MN);
    % Density coefficients
    N_rhos=10;
    N_rhol=100;
    Cl_range=linspace(0,1,N_rhol);
    Cs_range=linspace(0,1,N_rhos);
    if Has_volatile~=2 % currently not using H2O dependant rhof
        if (HHJPet==1 || SSPD==1)
            rhof=rhof_1*(1-Cl_range)+Cl_range*rhof_2;
        elseif FourMPD==1
            C_all_dim = Cl_range*(MgO_range(2)-MgO_range(1))+MgO_range(1);
            %rhof changes depending on comp
            rhof = C_all_dim.*0;
            for i=1:1:length(C_all_dim)
                if C_all_dim(i)>=crit_mg_si_melt
                    rhof(i) = rhof_2;
                else
                    rhof(i) = rhof_m*C_all_dim(i)+rhof_c;
                end
            end    
        end
        rhol_coef_all=piecewiseFit(rhof);
    else
        Cl2_range=linspace(0,0.1,N_rhol);        % 
        Coef_melt_den=[-0.004542610349603, 0.000103097162684, 0.001102980266019, 0.000249259014492];
        % rhof=zeros(N_rhol,N_rhol);
        % for i=1:N_rhol
        %     for j=1:N_rhol
        %         x_other=1-Cl_range(i)-Cl2_range(j);
        % 
        %         VMX=Cl2_range(j).*Cl_range(i).*Coef_melt_den(1)+Cl_range(i).*Coef_melt_den(2)+Cl2_range(j)*Coef_melt_den(3)+Coef_melt_den(4);
        %         V_all=Cl_range(i).*26.86e-6/0.06009+Cl2_range(j).*26.27e-6/0.01802+x_other.*VMX;
        %         rhof(i,j) =1/V_all;
        %     end
        % end
        SiO2=(Cl_range*74+(1-Cl_range)*47)/100; %SiO2%

        [Cl_grid, Cl2_grid] = meshgrid(SiO2, Cl2_range);
        x_other = 1 - Cl_grid - Cl2_grid;
        VMX = Cl2_grid .* Cl_grid * Coef_melt_den(1) + ...
            Cl_grid  * Coef_melt_den(2) + ...
            Cl2_grid * Coef_melt_den(3) + ...
            Coef_melt_den(4);
        V_all = Cl_grid  .* (26.86e-6 / 0.06009) + ...
            Cl2_grid .* (26.27e-6 / 0.01802) + ...
            x_other .* VMX;
        rhof = 1 ./ V_all;
        rhol_coef_all=piecewiseFit(rhof);
    end


    rhom=rhom_1*(1-Cs_range)+Cs_range*rhom_2;
    rhos_coef_all=piecewiseFit(rhom);

    %solidus coefficients
    % if HHJPet == 1
    %     Ts0=A1+B1+C1;
    %     Tl0=C1;
    % end
        
    N_sl=10; %since all varaibles here are normalized, does not need that large entries
    %The range of the scaled temperature should be beyond 0 and 1 for values outside of [Ts Tl] 
    % we force 0 and 1 to be in the range
    k = round((N_sl - 1) / 12);
    h=0.1/k;
    % T_range=-0.0:h:1.0;
    
    T_range=linspace(0,1,N_sl);
    phi_range=linspace(0,1,N_sl);
    Cl_range=linspace(0,1,N_sl);
    Cs_range=linspace(0,1,N_sl);

    % solidus_func=zeros(N_sl, length(T_range), N_sl, N_sl);  %use single precision for memory save
    % liquidus_func=zeros(N_sl, length(T_range), N_sl, N_sl); 
    % %In general solidus is F(phi, T, Cs, Cl)=0 
    % %solidus is cs=softmax(softmin(f,Cb),0);
    % A1_scaled=A1/(Tl0-Ts0);
    % B1_scaled=B1/(Tl0-Ts0);
    % C1_scaled=1;
    % for i=1:N_sl %phi
    %     for j=1:length(T_range) %T
    %         for k=1:N_sl %cs
    %             for p=1:N_sl %cl
    %                 if HHJPet == 1
    %                     cb=Cl_range(p)*phi_range(i)+Cs_range(k)*(1-phi_range(i));                        
    %                     temp=(max(1-T_range(j),0))^n_order;
    %                     temp=min(cb,temp);
    %                     solidus_func(i,j,k,p)=temp-Cs_range(k);
    % 
    %                     temp=max(min(T_range(j),1),0);
    %                     temp=(-B1_scaled-sqrt(B1_scaled^2-4*A1_scaled*(C1_scaled-temp)))/2/A1_scaled;
    %                     temp=max(cb,temp);
    %                     liquidus_func(i,j,k,p)=temp-Cl_range(p);
    %                 elseif SSPD ==1
    % 
    %                 elseif FourMPD==1
    % 
    %                 end
    %             end
    %         end
    %     end
    % end
    % solidus_coef_all=piecewiseFit(solidus_func);
    % liquidus_coef_all=piecewiseFit(liquidus_func);    

    
    %assjust all table numbers to be the cell number to be used later
    N_mus=N_mus-1;
    N_C=N_C-1;
    N_rhos=N_rhos-1;
    N_rhol=N_rhol-1;
    N_sl=N_sl-1;
    N_sl_T=length(T_range)-1;

    K1=0;
end

%% Currently only with Newton's method, volatile and the third phase can be included in the model
Par_v=0.2; % H2O solid/melt partition coefficient
if Has_volatile==1    
    
    %% H2O dependancy for solidus and liquidus - % paper based
    %(0,0)    : (57,326)
    %(1200,0) : (455,326)
    %(0,8000) : (57,70)
    % from 8000 bar to 0 bar, water 13% to 0 %
    %           13  12  11  10    9     8     7    6     5    4    3    2   1     0
    Data_point=[80, 83, 86, 91, 100,  115,  135, 164,  196, 230, 267, 307, 355, 401;... %70
        82, 82, 85,90, 98, 111.5,130.5,158,  189, 220, 255, 293, 343, 388;...   %102
        84,  84, 84,87.5,95.5,108.5,126, 151.5, 180, 211, 243, 280, 329, 375;... %134
        86,  86, 86,   86, 93, 105,  121, 144,   170, 199, 230, 266, 316, 362;... %166
        90.5,  90.5, 90.5,    90.5,90.5,100,  115, 136,   159, 186, 217, 252, 302, 349;... %198
        95,  95, 95,    95,   95, 95,  109, 128.5, 150, 176, 205, 241, 291, 339;... %223
        103,  103, 103,    103,   103,  103,  103, 120,   141, 165, 194, 230, 280, 328;... %248
        113,  113, 113,    113,   113,  113,  113, 113,  133.5,155.5,185,220, 270, 320;... %268
        127,  127, 127,    127,   127,  127,    127,   127,   127,146.5,176, 211.5,261,313;... %285
        141,  141, 141,    141,   141,  141,    141,   141,     141, 141, 170, 206, 256, 309;... %295
        164,  164, 164,    164,   164,  164,    164,   164,     164,   164, 164, 200, 250, 305;... %305
        194,  194, 194,    194,   194,  194,   194,    194,     194,   194, 194, 194, 244, 301;... %316
        240,  240, 240,    240,   240,  240,   240,    240,     240,    240,  240,  240, 240, 298;... %322
        297,  297,  297,   297,   297,  297,    297,   297,     297,   297,   297,   297,   297, 297];   %326
    Data_y=[70, 102,134,166,198,223,248,268,285,295,305,316,322,326];
    Data_y=8000-(Data_y-70)*8000/(326-70);
    Data_point=(Data_point-57)*(1200-600)/(455-57)+600;
    
    N_Ts=100;
    vol_range= linspace(0,13,N_Ts);
    Pressure_range=linspace(0,8000,N_Ts); %Pressure in bar
    F = griddedInterpolant({Data_y(end:-1:1), 0:13}, Data_point(end:-1:1,end:-1:1), 'linear', 'linear');
    [a_grid, b_grid] = ndgrid(Pressure_range, vol_range);
    Ts0_coefficient= piecewiseFit(F(a_grid, b_grid));

    %
    % (0,0) (61,738)
    % (0,1700) (838,738)
    % (0,40000) (61 33)
    %             20   10    5   2    0
    Data_point2=[358, 428,  512, 586, 699;...%120
        321, 397,  477, 547, 653;...%208
        288, 369,  443, 508, 604;...%297
        274, 351,  412, 473, 558;...%390
        340, 340,  391, 447, 517;...%473
        332, 332,  379, 427, 477;...%561
        329, 329,  376, 420, 462;...%596
        374,   374,  374, 414, 438;...%691
        413,   413,    413, 413, 436;...%728
        436,   436,    436, 436, 436];...%738
        Data_y2=[120,208,297,390,473,561,596,691,728,738];
    Data_y2=40000-(Data_y2-120)*40000/(738-120);
    Data_point2=(Data_point2-61)*(1700-800)/(838-61)+800;


    N_Tl=N_Ts;
    vol_range2= linspace(0,20,N_Tl);
    Pressure_range2=linspace(0,40000,N_Tl); %Pressure in bar
    F = griddedInterpolant({Data_y2(end:-1:1), [0 2 5 10 20]}, Data_point2(end:-1:1, end:-1:1), 'linear', 'linear');
    [a_grid, b_grid] = ndgrid(Pressure_range2, vol_range2);
    Tl0_coefficient= piecewiseFit(F(a_grid, b_grid));

    
    N_Ts=N_Ts-1;
    N_Tl=N_Tl-1;

    V_crust=1.2;
    V_sill=8;
    kf=1e-6;

    S_cap=[0.015 0.04];  %Solid water saturation for component A (74%, 0.5-1.5%) and B (47%,  3-5%) 

    

    Pressure=(nodez(end)-cellz)*g*rho_mean/1e5+1; %in bar
    PT=[-4e-4,-5e-4,-6e-4,-13e-4, -15.5e-4,-17e-4,-16e-4,-5e-4, 0, 26e-4,5e-3, 5e-3];  % coefficient from (ref) for the temperature dependency of water saturation (Holtz et al,?) 
    PTx=[0,   0.12    0.2  0.3    0.5       1       2       3   4  5,    11, 20]; %in kbar
    dSdT = interp1(PTx, PT, Pressure/1000, 'linear', 'extrap')';    
    
    P3=Pressure'/10; %Pressure in Mpa
    Lsaturation=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-(max(T,500)-800).*dSdT/100;
   

    S=zeros(N,1);
    Cl2=min(V_crust/Par_v*ones(N,1)/100, Lsaturation);
    Cl2(cellz<9000)=1e-5;
    % Cl2=V_crust/Par_v*ones(N,1)/100;
    % Cs2=V_crust*ones(N,1)/100;
    Cs2=Cl2*Par_v;

    injection_S=0;
    injection_Cl2=V_sill/100.*injection_phi;
    injection_Cs2=V_sill*Par_v/100;
    
    Ssaturation=S_cap(1)*Cs2+S_cap(2)*(1-Cs2);
    Cb2=Cl2.*phi(1:N)+Cs2.*(1-phi(1:N))+S;
    index=find(cellz<3000,1,'last');
    Cb2(1:index)=0;
    % Cb2=Cl2+Cs2+S;

    
    Ts0=zeros(N,1);
    Tl0=zeros(N,1);
    index_v_nts=max(min(floor(Cb2/13*100/Par_v*N_Ts)+1,N_Ts),1);
    index_pressure_nts=max(min(floor(Pressure/8000*N_Ts)+1,N_Ts),1);
    for i=1:N
        coef=Ts0_coefficient(index_pressure_nts(i),index_v_nts(i),:);
        Ts0(i)=coef(1)*Pressure(i)/8000+coef(2)*Cb2(i)/Par_v/13*100+coef(3)*Pressure(i)*Cb2(i)/Par_v/8000/13*100+coef(4);
    end

    index_v_ntl=max(min(floor(Cb2/20*100/Par_v*N_Ts)+1,N_Ts),1);
    index_pressure_ntl=max(min(floor(Pressure/40e3*N_Ts)+1,N_Ts),1);
    for i=1:N
        coef=Tl0_coefficient(index_pressure_ntl(i),index_v_ntl(i),:);
        Tl0(i)=coef(1)*Pressure(i)/40e3+coef(2)*Cb2(i)/Par_v/20*100+coef(3)*Pressure(i)*Cb2(i)/Par_v/400/20+coef(4);
    end
    Ts=Tl0-Cb.^(1/n_order).*(Tl0-Ts0);
    Tl=A1*Cb.^2+(Ts0-Tl0-A1).*Cb+Tl0;
    
    
    Pressure_injection=Injection_depth*g*rho_mean/1e2+1;
    index_v_nts=max(min(floor(V_sill/13/Par_v*N_Ts)+1,N_Ts),1);
    index_pressure_nts=max(min(floor(Pressure_injection/8000*N_Ts)+1,N_Ts),1);
    coef=Ts0_coefficient(index_pressure_nts,index_v_nts,:);
    Ts0_sill=coef(1)*Pressure_injection/8000+coef(2)*V_sill/13+coef(3)*Pressure_injection*V_sill/8000/13+coef(4);
    
    index_v_ntl=max(min(floor(V_sill/20/Par_v*N_Ts)+1,N_Ts),1);
    index_pressure_ntl=max(min(floor(Pressure_injection/40e3*N_Ts)+1,N_Ts),1);
    coef=Tl0_coefficient(index_pressure_ntl,index_v_ntl,:);
    Tl0_sill=coef(1)*Pressure_injection/40e3+coef(2)*V_sill/20+coef(3)*Pressure_injection*V_sill/40e3/20+coef(4);

    Total_cb20=sum(Cb2.*cellz');
else
    Ts0=(A1+B1+C1)*ones(N,1);
    Tl0=C1*ones(N,1);
end

if Has_volatile==1
    if FourMPD==1
        if injection_Cb>=liq_P2_C %injection temperature
            injection_T= liq_k1-(liq_a1/(injection_Cb-liq_b1));
        else
            injection_T= liq_k2 - (liq_a2/(injection_Cb-liq_b1));
        end
    elseif (SSPD==1||HHJPet==1)
        injection_T = A1*injection_Cl.^2 + (Ts0_sill-Tl0_sill-A1)*injection_Cl + Tl0_sill;
    end
else
    if FourMPD==1
        if injection_Cb>=liq_P2_C %injection temperature
            injection_T= liq_k1-(liq_a1/(injection_Cb-liq_b1));
        else
            injection_T= liq_k2 - (liq_a2/(injection_Cb-liq_b1));
        end
    elseif (SSPD==1||HHJPet==1)
        injection_T = A1*injection_Cl.^2 + B1*injection_Cl + C1+50;
    end
end

injection_H= Lf.*injection_phi+cp.*injection_T; %injection enthalpy


disp('All tables generated')


%%
H=Lf*phi+cp*T;
Cl=zeros(N,1);
Cs=zeros(N,1);
TYPE=zeros(N,1);

for i=1:N

    if HHJPet==1
        [~,~,Cl(i), Cs(i),~]=poro_component_solve_JPET_master(Cb(i), H(i),A1, B1, C1, A2, B2, C2, ae, Lf, cp,Precision_PD);
    end

    if SSPD==1
        [~,~,Cl(i), Cs(i),~] = poro_component_solve_solid_master(Cb(i), H(i), A1, Ts0(i)-A1-Tl0(i), Tl0(i), alpha, n_PD ,Lf, cp, Precision_PD,step_size);
    end


    if FourMPD==1
        [~,~,Cl(i), Cs(i),~,CCC(i),DDD(i)]=poro_component_solve_4_components_v2_master(Cb(i), H(i), Lf, cp, min_mgo, max_mgo, crit_T1, crit_T2,...
                                                            crit_cb2, crit_cb3, crit_cb4,crit_cb4_OG, crit_cb5,liq_k1, liq_a1, liq_b1,liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,...
                                                            lm1, lc1, lm3, lc3, lm8, lc8, mk, ck, BT, BC, DT, DC,JC, olpx_C,olpxm, olpxc, s5T, sol_k2,liq_P2_C,s1m,s1c, ...
                                                            lowest_end_T, liq_P3_C,l13_k,l13_a,l13_b,lowest_k,s8m,s8c,crit_T3,liq_P4_C,Precision_PD,step_size);
    end
end




C_all=[Cl;Cs];

%% Initialize the parameters
phi=[phi; 1-phi];

uf=zeros(N+1,1);
um=zeros(N+1,1);
u_all=[uf;um];

C_values=zeros(N,1);  % the coupling term coefficient
kt=kt0*ones(N,1);
kc=kc0*ones(N,1);

% Command window display initialization
prevStr = '';
conservation1=1;
conservation2=1;
%%


if Use_Newton==1

    % Generate the Jacobians for Newton's method or load from the exisitng
    [Jac_mom, rhs_mom, Jac_con, rhs_con, Jac_ct, rhs_ct, Jac_ent, rhs_ent, Jac_solidus, rhs_solidus, Jac_liquidus, rhs_liquidus, Jac_ct2,rhs_ct2, Jac_ssat, rhs_ssat, Jac_lsat, rhs_lsat,  Jac_ct_central , rhs_ct_central]=Newton_initiator2(Has_volatile, is_eutectic);
    
    disp(['All Jacobians generated/loaded']);

    u_all_old=u_all;
    H_old=H;
    T_old=T;
    C_all_old=C_all;
    Cb_old=Cb;

    phi_old=phi;
    if Has_volatile==1
        S_old=S;
        Cs2_old=Cs2;
        Cl2_old=Cl2;
    end
    

    if Has_volatile==1
        nC = size(Tl0_coefficient,3);
        Tl0_coefficient=reshape(Tl0_coefficient, [], nC);

        nC = size(Ts0_coefficient,3);
        Ts0_coefficient=reshape(Ts0_coefficient, [], nC);
        
        % nC = size(rhol_coef_all,3);
        % rhol_coef_all=reshape(rhol_coef_all, [], nC);
    end
    
    nC = size(C_coef_all,3);
    C_coef_all=reshape(C_coef_all, [], nC);
    
    Advance_time=0;
    Time_intended=1*Year;
    Newton_solver3;
end

%% Set up the monitor
%With_monitor=1; 

%Update_frame=1; % How many steps before the monitor  is refreshed

% Chooose what to show on the monitor
% labels for parameters to be plot:
% 1. melt fraction
% 2. liquid velocity
% 3. solid velocity
% 4. enthalpy
% 5. temperature
% 6. liquid composition
% 7. solid composition
% 8. bulk composition
% 9. material state:  0: below solidus, 1: at solidus 2: between solidus and liquidus 3: above liquidus
% 10. force contribution from solid viscous force
% 11. Heat-chemical-compaction contribution, stepwise.
% 12. Heat-chemical-compaction contribution, total (abs).
% 13. Refined recored of bulk composition
% 14. Refined recored of HCC contribution, total.
% 15. Heat-chemical-compaction contribution, total
% 16. Chemical differenciation from compacton and reactive flow, total
% 17. Chemical differenciation from compacton and reactive flow, stepwise
% 18. Chemical differenciation from compacton and reactive flow, total (abs)
% 19. Pressure gradient and the hydrostatic pressure rhofg and rhomg

% 1001. volatiles in phases+melt saturation+ solid saturation
% Plot_Label={'$$\phi$$','uf','um','H','T','Cf','Cm','Cb','St', 'Solid_vis(%)','$$\Delta\phi$$','$$\sum|\Delta\phi|$$'};
% when multiple variable in one figure, temperature will be scaled between solidus and liquidus
% option 11 and 12 must stand alone!

% Plot_configure=...
%     {[1],[2],[4],[5];...
%      [6],[7],[8],[11]};   % a 2*4 plot setting
% Plot_configure=...
%     {[1,5,8],[10],[11]};  


% Plot_configure=  {[1],[5],[6,7,8]};  
if With_monitor==1
    if Has_volatile==1
        Plot_configure=  {[1],[5],[6,7,8],1001};
    else
        if To_cal_mc==1
            Plot_configure=  {[1],[5],[6,7,8],[15],[16]};
        else
            Plot_configure=  {[1],[5],[6,7,8]};
        end
    end    
    Plot_settings_master;
    drawnow
    Update_frequency=0.1; %every X sec
    Start_timer=tic;
    Monitor_frame=0;
end
Adaptive_show_range=1;

%% Mesh adaptivity
Adaptive_mesh=0; %
Adaptive_step_gap=2;
Adaptive_step_time=10*Year;

Last_adapted=0;

Top0=0;
max_adaptive_number=1;
to_adapt=1;

min_change=2e-3;
max_change=4e-2;

min_dx=0.5;   % minimal cell length
max_dx=80; % maximum cell length
min_N=500;  % minimal number of allowed cells
max_N=8000; % maximum number of allowed cells

min_dx0=min_dx;



function ref_all = compute_F(phi_range, epsilon, phistar, gamma, B_vis, BS_ratio)
    sigma=13-gamma;  
    F = (1 - epsilon) * erf(sqrt(pi)/2/(1 - epsilon) ...
        * phi_range/phistar .* (1 + (phi_range/phistar).^gamma));    
    ref_mu=(1+(phi_range/phistar).^sigma)./(1-F).^(B_vis*phistar);
    ref_xi=ref_mu./max(0.0,phi_range)./(1-phi_range)*BS_ratio;
    ref_all=4/3*ref_mu+ref_xi;
end


