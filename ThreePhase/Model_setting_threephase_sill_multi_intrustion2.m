% For solid solution phase diagram
clear;
warning('off')

%A test setting for switching the influence of H2O off. 
% mode 0: all on
% mode 1: only viscosity off
% mode 2: only solidus off
% mode 3: all off
Fix_H2O=0; 

Model_name="D0.2 c1.2 s5 SiO48-60 MDI mu 0.78-17 muf 0 2 4 7 P2F2";




Use_parrallel=1;  %try to use parrallel
if Use_parrallel==1
   Num_core_use=5;  %% set the number of cores to use!
end


rng(13);  %random number seed. For generating consistant result.
%%
PD_range=[47 74]; % the two end members of the phase diagram
N_component=3;
Add_CLCU=1; %add Cl and Cu as trace element

Year=3600*24*365.25;

Conservation_type=2; % 1 for mass conservation, 2 for volume conservation
Transport_method=2;  % 1 transport all component with scaling  2 Transport component A and V only without scaling
Velocity_solver_type=2; % 1: Three phase solver, 2: Two phase solver

Dynamic_scaling=1e7; % To improve the condition number of velocity solver LHS matrix
g=9.81/Dynamic_scaling;

Precision=1e-2;
% Precision2=Precision/1e2;
Max_iter=2;
min_iter=1;

Enhanced_convergence=0;  %


% Coupling term parameters
A_type2=3.5e-3^2/70;  %2.5mm grain size, beta1=70 according to Claudio slides  ( 5)
% grain_size=2e-3;  %    d 1.8e-3 A=300
% C_value_A=200;      % Coefficient for the power low permeability k=1/A*grain_size^(-2)*phi^(-B),    58
% C_value_B=3;        % Coefficient for the power low permeability,      2.6
% C_type=0;  % 1 for porous-suspension type.  0 for porous only 


min_por=0.00; %minimal porosity of the solid
max_melt=1; %maximum allowed melt fraction

Vol_flux=0; %volatile? flux (mass)
S_cap=[0.015 0.04]; %Solid water saturation for component A (74%, 0.5-1.5%) and B (47%,  3-5%) 
%% Viscosities
melt_viscosity_type=1;  %0 for old type, 1 for new
% old model: order [0 5] from least to most evolved melt
% New melt viscosity model, parameters for viscosity of [least/wet, least/dry, most/wet, most/dry] at order 
% melt viscosity follows log10(mul)= A+B/(T-C);
% A=-4.55
if melt_viscosity_type==0
    muf1=0;
    muf2=5;    
else
    %[0 2 3 6]
    % vis_b1= 0.020457395895203e4; vis_b2=-0.660135051778619e4; vis_b3= 1.382974492987449e4;
    % vis_c1=-0.000453961313838e4; vis_c2=0.360872958639910e4; vis_c3=-0.726473385825352e4;
    %[0 2 3 6] NEW!!
    % vis_b1= 241.2067; vis_b2=-1.1208e+05; vis_b3= 1.6578e+05;
    % vis_c1=-4.3207; vis_c2=1.8300e+04; vis_c3=-2.8578e+04;
    %[0 2 4 7] NEW!!
    vis_b1=  211.2559; vis_b2=-9.8778e+04; vis_b3=  1.4453e+05;
    vis_c1=0.2519; vis_c2=1.6674e+04; vis_c3=-2.5643e+04;    
end


%volatile viscosity order
mug=5;
% mum=log(1e13*4/3)/log(10); % Solid viscosity order, remember the 1d/3d conversion

% Set_up the solid shear and bulk viscosities based on Costa2019
N_vis=1e5;
phi_range=linspace(0,1,N_vis);
phi_range(phi_range<1e-2)=1e-2;
phi_range(phi_range>1-2e-2)=1-2e-2;


B_vis=2; % Einstein coefficient

% phistar=0.736;%;
% gamma=3.679;
% sigma=13-gamma;
epsilon=0.5e-4;

phistar= 0.78;
gamma=12.5;
sigma=13-gamma;



%beta
% 18: 5e-5
% 17: 2.2e-4
ref_shear=1e17;
% epsilon=2.2e-4;
% epsilon=1e-2;

BS_ratio=3/5;     %0.25-100


% mu_m=10^14.4*ones(1,N_vis); %(1+(phi_range/phistar).^sigma)./(1-F).^(B_vis*phistar)*ref_shear;
% % mu_m(round(N_vis/5):end)=mu_m(round(N_vis/5):end)*5;
% xi_m=mu_m.*phi_range.^-0.5;    %mu_m./max(0.0,phi_range)./(1-phi_range)*BS_ratio; % xi=mu/phi/(1-phi)*ratio
% mu_all2=4/3*mu_m+xi_m;


F=(1-epsilon)*erf(pi^0.5/2/(1-epsilon)*phi_range/phistar.*(1+(phi_range/phistar).^gamma));

ref_mu=(1+(phi_range/phistar).^sigma)./(1-F).^(B_vis*phistar);
ref_xi=ref_mu./max(0.0,phi_range)./(1-phi_range)*BS_ratio;
ref_all=4/3*ref_mu+ref_xi;


mu_all=ref_all*ref_shear/ref_all(end);
mu_m=ref_mu*ref_shear/ref_all(end);
xi_m=ref_xi*ref_shear/ref_all(end);

% mu_all=ref_all;



mu_all=mu_all(end:-1:1);
phi_range=phi_range(end:-1:1);
%figure(3)
%clf; set(gcf,'color','w');
%set(gca,'TickDir','out');
%semilogy(phi_range,mu_m,'LineWidth',2)
%hold on
%semilogy(phi_range,xi_m,'LineWidth',2)
%semilogy(phi_range(end:-1:1),mu_all,'LineWidth',2)
% semilogy(phi_range(end:-1:1),mu_all2,'LineWidth',2)
% grid on
% ylim([min(mu_all)/2 max(mu_all)])
%xlabel('Melt Fraction (-)',  'FontSize', 13)
%ylabel('Solid viscosity (Pa s)',  'FontSize', 13)
%legend({'Shear','Bulk','Sum','Type2'})
%set(gca,'fontsize',13)
%annotation('textarrow',[0.2 0.13], [0.9 0.85],'String','Ref viscosity','fontsize',13)
%% Density parameters
% Original data: 
%melt
Melt_density_type=2; % type 2 to include the effect of H2O
x1=50.5; y1=2772.04;   x2=70.625; y2=2422.60;
rhof_2=(y2-y1)/(x2-x1)*(47-x1)+y1; %Density fluid, least evolved at 47% 
rhof_1=(y2-y1)/(x2-x1)*(74-x1)+y1; %Density fluid, most evolved AT 74%

% Type2 melt density: Based on (Lesher & Spera,2015)
% VMX=H2O*SiO2*coef1+SiO2*coef2+H2O*coef3+coef4;
% x_other=1-x-H2O
% V=x.*26.86e-6/0.06009+H2O.*26.27e-6/0.01802+x_other.*VMx;
% Den_melt=1/V;
Coef_melt_den=[-0.004542610349603, 0.000103097162684, 0.001102980266019, 0.000249259014492];

%Solid
x1=51.094; y1=2959.43;  x2=70.99; y2=2660.13;
rhom_2=(y2-y1)/(x2-x1)*(47-x1)+y1; %Density solid, least evolved
rhom_1=(y2-y1)/(x2-x1)*(74-x1)+y1; %Density solid, most evolved

rho_constant=[rhof_1, rhof_2; rhom_1, rhom_2];

rho_mean=2800; %the value used when volumn conservation/Bossinessq approximation is used 
%%%%%%  
%% Thermal parameters
Lf=[350e3   600e3];  %Latent heat of component A and B

cp0=[1100 1100 1100];   %Heat capacity of [component A, B, water]
kt0=[2.5, 1.5]/rho_mean;      %Thermal difusivity for component A B

kc0=mean(kt0)*1e-7;
kc20=mean(kt0)*1e-7*1;

% kc0=1e-12;        %Chemical diffusivity for SiO2
% kc20=1e-10;        %Chemical diffusivity for H2O

BC_type=[3,3];
BC_value=[0,0];
%% State equation data -% paper based huber
aa=[-112.528, 127.811, 112.04];
bb=[-0.381, -1.135, -0.411];
cc=[0.033, 0, 0];

%% Phase diagram parameters  
A1=-99.6;  C1=1160; B1=760-C1-A1;

a1=A1/(-A1-B1); b1=B1/(-A1-B1); c1=1; %-a1-b1;
% A1=50; B1=-360; C1=1433.15-273.15; %Matt paper
A2=0; B2=0; C2=0; ae=1;  
Ts=A1+B1+C1;
Tl=C1;

is_solid_solution=1; 
simplified_TS=1; % simplified solidus shape, then alpha =1;
if N_component==5
    simplified_5p=1;
else
    simplified_5p=0;
end
n_order=40;  % the order for T=(a2*(1-x).^n+(1-a2)*(1-x.^(1/n)))*(Tl-Ts)+Ts;
alpha_order=1;

% [Jacobians,Rhs,Constant_index ]=chemical_solver_initiator(is_solid_solution, 3, simplified_TS, simplified_5p);   %For initialize the Jacobion, must remove for compiling
chemical_solver_loader;
%% H2O dependancy for solidus and liquidus - % paper based 
%(0,0)    : (57,326)
%(1200,0) : (455,326)
%(0,8000) : (57,70)
% from 8000 bar to 0 bar, water 13% to 0 %
%           13  12  11  10    9     8     7    6     5    4    3    2   1     0
Data_point=[80, 83, 86, 91, 100,  115,  135, 164,  196, 230, 267, 307, 355, 401;... %70
             0, 82, 85,90, 98, 111.5,130.5,158,  189, 220, 255, 293, 343, 388;... %102
             0,  0, 84,87.5,95.5,108.5,126, 151.5, 180, 211, 243, 280, 329, 375;... %134
             0,  0, 0,   86, 93, 105,  121, 144,   170, 199, 230, 266, 316, 362;... %166
             0,  0, 0,    0,90.5,100,  115, 136,   159, 186, 217, 252, 302, 349;... %198
             0,  0, 0,    0,   0, 95,  109, 128.5, 150, 176, 205, 241, 291, 339;... %223
             0,  0, 0,    0,   0,  0,  103, 120,   141, 165, 194, 230, 280, 328;... %248
             0,  0, 0,    0,   0,  0,    0, 113,  133.5,155.5,185,220, 270, 320;... %268
             0,  0, 0,    0,   0,  0,    0,   0,   127,146.5,176, 211.5,261,313;... %285
             0,  0, 0,    0,   0,  0,    0,   0,     0, 141, 170, 206, 256, 309;... %295
             0,  0, 0,    0,   0,  0,    0,   0,     0,   0, 164, 200, 250, 305;... %305
             0,  0, 0,    0,   0,  0,   0,    0,     0,   0,   0, 194, 244, 301;... %316
             0,  0, 0,    0,   0,  0,   0,    0,     0,    0,  0,   0, 240, 298;... %322
             0,  0,  0,   0,   0,  0,    0,   0,     0,   0,   0,   0,   0, 297];   %326   
Data_y=[70, 102,134,166,198,223,248,268,285,295,305,316,322,326];
Data_y=8000-(Data_y-70)*8000/(326-70);
Data_point=(Data_point-57)*(1200-600)/(455-57)+600;
%
% (0,0) (61,738)
% (0,1700) (838,738)
% (0,40000) (61 33)
%             20   10    5   2    0
Data_point2=[358, 428,  512, 586, 699;...%120
             321, 397,  477, 547, 653;...%208
             288, 369,  443, 508, 604;...%297
             274, 351,  412, 473, 558;...%390
               0, 340,  391, 447, 517;...%473
               0, 332,  379, 427, 477;...%561
               0, 329,  376, 420, 462;...%596
               0,   0,  374, 414, 438;...%691
               0,   0,    0, 413, 436;...%728
               0,   0,    0,   0, 436];...%738  
Data_y2=[120,208,297,390,473,561,596,691,728,738];               
Data_y2=40000-(Data_y2-120)*40000/(738-120);
Data_point2=(Data_point2-61)*(1700-800)/(838-61)+800;
%% Domain
% CAB
Length01=4700;
Length0=9700;
Length1=14700;
Length2=8300;
% CAB end

Sill_length=50;
margin=0;
fine_ratio=3;  %fine meshed region, in comparison with the size of the initial sill region
Show_z=[Length2 Length2+Sill_length*fine_ratio];
Sill=[Length2+Sill_length*(fine_ratio-1)/2,Length2+Sill_length*(fine_ratio+1)/2];  %intruded sill position
min_show_range=Show_z(2)-Show_z(1);

% define the meshing and mesh adaptivity
Adaptive_mesh=1;     % set to one to turn on the adaptive meshing.
Adaptive_step_gap0=10; %frequency of mesh adaptivity done
Adaptive_step_time=1*Year; %minimal time gap before the new adaptation
max_adaptive_number=1;

min_dx=5;   % minimal cell length
max_dx=500; % maximum cell length
min_N=100;  % minimal number of allowed cells
max_N=8000; % maximum number of allowed cells
aspect_ratio=2.3; % maximum change ratio between adjacent cells,


dphi_max=3e-2; %targeting change of melt fraction in one cell.
dphi_min=4e-3; %when mesh needs to be coarse.=1e-2

phi_threshold=0.75;
phi_threshold2=0.15;

dS_max=1e-2;
dS_min=5e-3;

dMNV_min=20;    %change of MNV 
dMNV_max=30;


Last_adapted=0;
Adaptive_step_gap=Adaptive_step_gap0;
%%
% CAB
part1=linspace(0,Length2,max(round(min_N/8),5));
part2=linspace(0,Sill_length*fine_ratio,min_N);
part3=linspace(0,Length1,max(round(min_N/8),5));
part4=linspace(0,Sill_length*fine_ratio,min_N);
part5=linspace(0,Length0,max(round(min_N/8),5));
part6=linspace(0,Sill_length*fine_ratio,min_N);
part7=linspace(0,Length01,max(round(min_N/8),5));
nodez=[part1 part1(end)+part2(2:end) part1(end)+part2(end)+part3(2:end)...
    part1(end)+part2(end)+part3(end)+part4(2:end) part1(end)+part2(end)+part3(end)+part4(end)+part5(2:end)...
    part1(end)+part2(end)+part3(end)+part4(end)+part5(end)+part6(2:end) part1(end)+part2(end)+part3(end)+part4(end)+part5(end)+part6(end)+part7(2:end)];
% CAB end 
nodez=nodez';
% nodez=nodez-nodez(end); %Change the distance

N=length(nodez)-1;    

cellz=(nodez(1:end-1)+nodez(2:end))/2;  % cell center points

dz=nodez(2:end)-nodez(1:end-1);

Sill_index=zeros(2,1);
Sill_index(1)=find(cellz>Sill(1),1,'first');
Sill_index(2)=find(cellz<Sill(2),1,'last');

nodez_km=(nodez-nodez(end))/1000;
cellz_km=(cellz-nodez(end))/1000;
%% Calculate densities of three components in each phases from SiO2% and H2O%
N_number=N;

SiO2_crust=52;%SiO2% at bottom %  % percent
SiO2_crust2=52;%SiO2% at -15km
SiO2_crust3=55;%SiO2% at surface




H2O_crust=1.21;%1.21;%H2O% at bottom  %1.21;% for D0.2
H2O_crust2=1.21;%1.21;%H2O% at -15km
H2O_crust3=1.21;%1.21;%H2O% at surface


SiO2_sill=50;%0.09*(74-47)+47;  % percent  
H2O_sill=5;%1.21; % percent


CL_crust=0.5;
CU_crust=0.1;

CL_sill=2;
CU_sill=1;


Sillphi=0.935;% melt fraction of the sill %0.934: 50S1 0.876 50S2 ; 0.861 50S5
              %Geo24.5 phi 0.835 H2O 1.5:0.8

x=Calculate_Initial_state(PD_range, rho_mean, SiO2_crust, H2O_crust, CL_crust, N_component, Conservation_type, Precision);              
MM_crust=x(1);
NN_crust=x(2);
V_crust=x(3);
if N_component==5
    CL1=x(4);
end


x=Calculate_Initial_state(PD_range, rho_mean, SiO2_crust2, H2O_crust2, CL_crust, N_component, Conservation_type, Precision);              
MM_crust2=x(1);
NN_crust2=x(2);
V_crust2=x(3);
if N_component==5
    CL2=x(4);
end

x=Calculate_Initial_state(PD_range, rho_mean, SiO2_crust3, H2O_crust3, CL_crust, N_component, Conservation_type, Precision);              
MM_crust3=x(1);
NN_crust3=x(2);
V_crust3=x(3);
if N_component==5
    CL3=x(4);
end


x=Calculate_Initial_state(PD_range, rho_mean, SiO2_sill, H2O_sill, CL_sill, N_component, Conservation_type, Precision);              
MM_sill=x(1);
NN_sill=x(2);
V_sill=x(3);
if N_component==5
    CL_sill=x(4);
end

% For sill part, negative values are possible....so used the same procedure

N=N_number;


% MM=linspace(MM_crust,MM_crust2,N)'; 
index=find(cellz-nodez(end)>=-15e3-100,1,'first');
MM=interp1([cellz(1) cellz(index) cellz(end)], [MM_crust MM_crust2 MM_crust3], cellz,'linear' ); 
MM(Sill_index(1):Sill_index(2))=MM_sill;


% NN=linspace(NN_crust,NN_crust2,N)';
NN=interp1([cellz(1) cellz(index) cellz(end)], [NN_crust NN_crust2 NN_crust3], cellz ,'linear' ); 
NN(Sill_index(1):Sill_index(2))=NN_sill;


% V=zeros(N,1);
% V=linspace(V_crust,V_crust2,N)';
V=interp1([cellz(1) cellz(index) cellz(end)], [V_crust V_crust2 V_crust3], cellz ,'linear' ); 
V(Sill_index(1):Sill_index(2))=V_sill;


if N_component==5
    CL=interp1([cellz(1) cellz(index) cellz(end)], [CL1 CL2 CL3], cellz ,'linear' ); 
    CU=(MM+NN+V+CL)*CU_crust/100; %trace element
    
    CL(Sill_index(1):Sill_index(2))=CL_sill;
    CU(Sill_index(1):Sill_index(2))=(MM_sill+NN_sill+V_sill+CL_sill)*CU_sill/100;
end

if Add_CLCU==1
    CL=(MM+NN+V)*CL_crust/100; %trace element
    CU=(MM+NN+V)*CU_crust/100; %trace element

    CL(Sill_index(1):Sill_index(2))=(MM_sill+NN_sill+V_sill)*CL_sill/100;
    CU(Sill_index(1):Sill_index(2))=(MM_sill+NN_sill+V_sill)*CU_sill/100;
else
    CL=zeros(N,1);
    CU=zeros(N,1);
end

% make the lower part dryer
Dry_depth=Length2-300;  % the lowest part has much higher solidus
index_temp=find(cellz>Dry_depth,1,'first');
V(1:index_temp)=0.01*ones(index_temp,1); %linspace(0.001, V(index_temp), index_temp);
if Conservation_type==2
    if N_component==3
        temp=rho_mean-(MM(1:index_temp)+NN(1:index_temp)+V(1:index_temp));
        correction=MM(1:index_temp)./(MM(1:index_temp)+NN(1:index_temp));
        MM(1:index_temp)=MM(1:index_temp)+temp.*correction;
        NN(1:index_temp)=NN(1:index_temp)+temp.*(1-correction);
    else
        temp=rho_mean-(MM(1:index_temp)+NN(1:index_temp)+V(1:index_temp)+CL(1:index_temp));
        correction=MM(1:index_temp)./(MM(1:index_temp)+NN(1:index_temp)+CL(1:index_temp));
        correction2=NN(1:index_temp)./(MM(1:index_temp)+NN(1:index_temp)+CL(1:index_temp));
        MM(1:index_temp)=MM(1:index_temp)+temp.*correction;
        NN(1:index_temp)=NN(1:index_temp)+temp.*correction2;
        CL(1:index_temp)=CL(1:index_temp)+temp.*(1-correction2-correction);
    end
end
dry_index=index_temp;




if N_component==3
    T0=abs(cellz-nodez(end))/1000*25+10;   %thermal gradient 25C/km, surface at 10C

    H=(MM*cp0(1)+NN*cp0(2)+V*cp0(3)).*T0;  %background enthalpy, assuming no latent heat and 20C/km themal gradient
    H(Sill_index(1):Sill_index(2))=H(Sill_index(1):Sill_index(2))+mean(Lf)*(MM_sill+NN_sill+V_sill)*Sillphi*2.1;
    Sill_den=MM_sill+NN_sill+V_sill;
else
    T0=abs(cellz-nodez(end))/1000*25+10;   %thermal gradient 25C/km, surface at 10C

    % Cu is treated as trace element, not taken into enthaly calculation
    H=(MM*cp0(1)+NN*cp0(2)+(V)*cp0(3)).*T0;
    Lf_temp=cb*Lf(1)+(1-cb)*Lf(2);
    H(Sill_index(1):Sill_index(2))=H(Sill_index(1):Sill_index(2))+Lf_temp*(MM_sill+NN_sill+V_sill)*Sillphi*1.85;
    Sill_den=MM_sill+NN_sill+V_sill+CL_sill;
end

% % H(Sill_index(2)+1:end)=H_bot;
H0=H;

Show_z=[Length2 Length2+Sill_length*fine_ratio];
% Show_z=[cellz(1) cellz(end)];


Sum_M0=sum(MM.*dz);
Sum_N0=sum(NN.*dz);
Sum_V0=sum(V.*dz);
%% Time steps
fixed_dt=0;

Courant0=9.5e-1;%0.025;  %During the start of the sill intrusion
Courant1=9.5e-1;%0.15;  
min_phi_dt=5e-2;

Max_dt=200*Year;
Min_dt=5*Year;

if fixed_dt==1
    dt=1e-1*Year;    
else
    dt=1e-3*Year;
end
End_time=40e6*Year;   %250e3*Year;
% End_time=50e3*Year;   %250e3*Year;

Time=0;
Courant=Courant0;
%% Temperature
% T=931; %temperature in C
phi=0.0*ones(N,1);
S=0.0*ones(N,1);
% 0=kt0*ones(N,3);
kc=kc0*ones(N+1,3);
cp=cp0.*ones(N,1);
%% Water component and partition coefficient
K0=0.2;

% if N_component==5 
Kcl=0.05;  % Solid/melt partition coefficient of Cl 0-0.2
Kcu=0.5;   % Solid/melt partition coefficient of Cu 0.2-2.0
% end
%%


K=K0*ones(N,1);
Cl2=0.08*ones(N,1);
Cs2=Cl2*K0;

% Pressure at the top 
P_top=1; % unit bar

% Volitile saturation coefficient
% fitting the points at 0.5, 2, 4 kbar
T_cutmin=800;
T_cutmax=1200;
% coef1=[16,4, 1; 4 2,1; 0.25, 0.5, 1]\[9.2 ; 6;   2.5];
% coef2=[16,4, 1; 4 2,1; 0.25, 0.5, 1]\[7.1 ; 4.8; 2];

if is_solid_solution==0
    Sys_constant=[A1,B1,C1; aa; bb;cc; S_cap, 0; 0      ,          0,rho_mean];
else
    Sys_constant=[A1,B1,C1; aa; bb;cc; S_cap, 0; n_order,alpha_order,rho_mean];
end
rho=ones(N,3);
%%
% mu_f_dynamics=10^muf1*ones(N+1,1)/Dynamic_scaling;
mu_m_dynamics=ones(N+1,1)/Dynamic_scaling;
mu_g_dynamics=10^mug*ones(N+1,1)/Dynamic_scaling;


ZerosU1=[1,N+1]; 
ZerosU2=[1,N+1];
ZerosU3=[1,N+1];

P=zeros(N,1); %zeros(N,1)+P_top*1e5; %Pressure in Pa
if N_component==3
    P(end)=P_top*1e5+(MM(1)+NN(1)+V(1))/2*g*dz(end)/2;
    for i=N-1:-1:1
        P(i)=P(i+1)+((MM(i)+NN(i)+V(i))*dz(i+1)+(MM(i+1)+NN(i+1)+V(i+1))*dz(i))/2*g*Dynamic_scaling;
    end
else
    P(end)=P_top*1e5+(MM(1)+NN(1)+V(1)+CL(1))/2*g*dz(end)/2;
    for i=N-1:-1:1
        P(i)=P(i+1)+((MM(i)+NN(i)+V(i)+CL(i))*dz(i+1)+(MM(i+1)+NN(i+1)+V(i+1)+CL(i+1))*dz(i))/2*g*Dynamic_scaling;
    end
end

dP=zeros(N+1,1);

C_values=zeros(N+1,3);

U_new=zeros((N+1)*2,1);

T=ones(N,1);

%% Input data correction by the chemical model

Pg_real=P/1e8; % P in kbar
Pg_real0=Pg_real;
iter=zeros(N+1,1);
k2=zeros(N,1);
% Ts_new=zeros(N,1);
T_S_region=zeros(N,2);
Ts_local=zeros(N,1);
Tl_local=zeros(N,1);
S_cap=zeros(N,1);
Saturation=zeros(N,1);
if N_component==3 && Add_CLCU~=1
    Mass_data=zeros(N,7);
else
    Mass_data=zeros(N,13);
    Partition_CL_CU=zeros(N,2);
end
% pressure input here is kbar: 1e8 Pasca
% [phi(1),S(1),T(1),rho(1,:),iter(1),k2(1),Ts_new(1),Mass_data(1,:)]=Phase_component_solve_new5(H(1),MM(1),NN(1),V(1),Pg_real(1), 1000, Sys_constant,cp(1,:),  rho_constant, Lf, K, Ts,  Tl, min_por, dz(1),Data_point,Data_y,Data_point2,Data_y2,Precision,3, Conservation_type);
% for i=2:N
%     [phi(i),S(i),T(i),rho(i,:),iter(i),k2(i),Ts_new(i),Mass_data(i,:)]=Phase_component_solve_new5(H(i),MM(i),NN(i),V(i),Pg_real(i), T(i-1), Sys_constant,cp(i,:),   rho_constant, Lf, K, Ts, Tl, min_por, dz(i),Data_point,Data_y,Data_point2,Data_y2,Precision,3, Conservation_type);
% end


if is_solid_solution==0
    [phi(1),S(1),T(1),rho(1,:),iter(1),T_S_region(1,:),Ts_new(1),Mass_data(1,:),S_cap(1)]=Phase_component_updated2(H(1),MM(1),NN(1),V(1),Pg_real(1), 500, [1,1], Jacobians, Rhs, Constant_index, Sys_constant,cp(1,:),  rho_constant, Lf, K(1), Ts,  Tl, min_por, dz(1),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type);
    for i=2:N
        %     if i==40
        %         aaa=1;
        %     end
        [phi(i),S(i),T(i),rho(i,:),iter(i),T_S_region(i,:),Ts_new(i),Mass_data(i,:),S_cap(i)]=Phase_component_updated2(H(i),MM(i),NN(i),V(i),Pg_real(i-1),T(i-1), T_S_region(i-1,:), Jacobians, Rhs, Constant_index, Sys_constant,cp(i,:),   rho_constant, Lf, K(i), Ts, Tl, min_por, dz(i),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type);
    end
else
    if N_component==3
        if Add_CLCU==1
            input=[H(1),MM(1),NN(1),V(1),CL(1),CU(1)];
        else            
            input=[H(1),MM(1),NN(1),V(1)];
        end
        [phi(1),S(1),T(1),rho(1,:),iter(1),T_S_region(1,:),Ts_new(1),Ts_local(1),Mass_data(1,:),S_cap(1),Tl_local(1),~,~, Partition_CL_CU(1,:)]=Phase_component_updated_solid3(input,Pg_real(1), 500, [1,1], Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(1,:),  rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts,  Tl, min_por, max_melt, dz(1),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS,0);
        for i=2:N
            if Add_CLCU==1
                input=[H(i),MM(i),NN(i),V(i),CL(i),CU(i)];
            else
                input=[H(i),MM(i),NN(i),V(i)];
            end
            [phi(i),S(i),T(i),rho(i,:),iter(i),T_S_region(i,:),Ts_new(i),Ts_local(i),Mass_data(i,:),S_cap(i),Tl_local(i),~,~, Partition_CL_CU(i,:)]=Phase_component_updated_solid3(input,Pg_real(i-1),T(i-1), 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(i,:),   rho_constant, Lf, [K(i) Kcl Kcu], Ts, Tl, min_por, max_melt, dz(i),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS,0);
        end
    else
        [phi(1),S(1),T(1),rho(1,:),iter(1),T_S_region(1),Ts_new(1),Ts_local(1),Mass_data(1,:),S_cap(1),Partition_CL_CU(1,:)]=Phase_component_updated_solid_5p_simple(H(1),MM(1),NN(1),V(1),CL(1),CU(1), Pg_real(1), 1000, 1, Jacobians, Rhs, Constant_index,Extra_func, Sys_constant, cp(1,:),   rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts, Tl, min_por, dz, Data_point,Data_y,Data_point2,Data_y2, PD_range, Precision, Conservation_type, 1);
        for i=2:N
            [phi(i),S(i),T(i),rho(i,:),iter(i),T_S_region(i),Ts_local(i),Tl_local(i),Mass_data(i,:),S_cap(i),Partition_CL_CU(i,:)]=Phase_component_updated_solid_5p_simple(H(i),MM(i),NN(i),V(i),CL(i),CU(i), Pg_real(i), T(i-1), T_S_region(i-1), Jacobians, Rhs, Constant_index, Extra_func, Sys_constant, cp(i,:),   rho_constant, Lf, [K(i) Kcl Kcu Partition_CL_CU(i-1,:)], Ts, Tl, min_por, dz, Data_point,Data_y,Data_point2,Data_y2, PD_range, Precision, Conservation_type, 1);
        end
    end
end


% Mass_data=[m1,n1,v1,m2,n2,v2,v3];

Solid_saturation0=0.05;
Solid_saturation=zeros(N,1);
for i=1:N
    [~,~, temp]=Cal_liquidus(19,Pg_real0(i)*1000,Data_point2,Data_y2);
    Solid_saturation(i)=Solid_saturation0;
end
Cb=(PD_range(2)*sum(Mass_data(:,[1,4]),2)+PD_range(1)*sum(Mass_data(:,[2,5]),2))./sum(Mass_data(:,[1 2 4 5]),2);


Convective_cut=0.990;
Convective_cut2=0.990;

Calculate_capped_values2;
%% Injection
To_intrude=1;  %set to 1 to turn on the multi intrusion

% Injection depth method
% method 99: always at the same depth
% When failed to find lighter than sill crust, 
% method 0: inject at the crust with minimal density; 
% method 1: at the same location %not to use !!
% method 2: always intrude at the minimal density 
% method 3: find the deepest phi_crit position, otherwise intrude at the
% previous depth

% All methods are affected by the random shift 
Injection_depth_method=3; 
phi_crit=0.2; %the porosity where the initusion happens

%Number of nodes per sill
dzf=part2(2)-part2(1);
SillNodez=round(Sill_length/(dzf));

% injection_time=(5e3: 5e3: 3000e3)*Year; %
% injection_time=(5e3: 5e3: 2000e3)*Year; %
injection_time=(10e3: 10e3: 4000e3)*Year;
% injection_time=(15e3: 15e3: 3000e3)*Year;
% injection_time=(20e3: 20e3: 8000e3)*Year; %
% injection_time=(30e3: 30e3: 12000e3)*Year; %
% injection_time=(40e3: 40e3: 16000e3)*Year; %
% injection_time=(30e3: 30e3: 18000e3)*Year;
% injection_time=(50e3: 50e3: 20000e3)*Year;
% injection_time=(60e3: 60e3: 24000e3)*Year;
% injection_time=(75e3: 75e3: 30000e3)*Year;
% injection_time=[50e3: 50e3: 10000e3, (10000e3+40e3):40e3:18000e3, (18000e3+30e3):30e3:24000e3]*Year;
% injection_time=(10e3: 10e3: 3000e3)*Year;
injection_time=[injection_time End_time+100]; 

% injection_time=(30:50:500)*Year; %
Random_shift=10; %random shift of distance from the intended position

Injected=1;
% Tot_Sill=100; % Total magma intruded (in metres)

Injection_depth=30000;

%Intrusion depth - Nodez
depthN = find(nodez>Injection_depth, 1, 'first');
%Intrusion depth -cellz
depthCN = find(cellz==((nodez(depthN)+nodez(depthN-1))/2));
SillMass= Mass_data(Sill_index(1),:);
SillH=H(Sill_index(1));
SillT=T(Sill_index(1));
Sillphi=phi(Sill_index(1));
SillS=S(Sill_index(1));
Sillrho=rho(Sill_index(1),:);
if Conservation_type==1
    Sillden=sum(SillMass);
else
    i=Sill_index(1);
    m1=Mass_data(i,1); n1=Mass_data(i,2); v1=Mass_data(i,3); m2=Mass_data(i,4); n2=Mass_data(i,5); v2=Mass_data(i,6); v3=Mass_data(i,7);
    Sillden=((rhof_1*m1+rhof_1*n1)/max(m1+n1,1e-3)*(m1+n1+v1)+(rhom_1*m2+rhom_2*n2)/max(m2+n2,1e-3)*(m2+n2+v2))/(m1+n1+v1+m2+n2+v2+v3);
end

Last_intrude_time=0;


%% Evacuations - CAB
To_evacuate=1; %set to 1 to enable evacuation
% initialising arrays for layers
% buoyant layers which have porosity>0
buoy_phi_top=[];
buoy_phi_base=[];
% buoyant layers which have porosity>CMF
buoy_Hphi_top=[];
buoy_Hphi_base=[];
% layers which have porosity>0 and contain buoyant magma
all_phi_buoy_top=[];
all_phi_buoy_base=[];
% confined time for hrti
t2=[];
% hrti thickness
hrti=[];

%critical melt fraction to be classed as high melt fraction
CMF=0.6;

% define
%Thickness of magma to calculate the density of the rock above (m)
xabove = 2000;
%bulk density
rho_bulk=phi.*(1-S).*rho(:,1)+(1-phi).*rho(:,2)+phi.*S.*rho(:,3);

%Diameter of reservoir (m)
Diameter = 60000;
%Crust viscosity (Pa s)
crust_visc = 1e19;
%Rock strength (Pa)
crit_crust =10e6; %3e6;
%Initial RTI scale (-)
scale_RTI = 0.001;
%Elastic modulus (Pa)
E_layer = 1e10;
%Freezing parameter (-)
freeze=0.14;
% Thickness of crust above magma layer to take thermal gradient over (m)
dTabove = 200;

%Counter of evacuations
evacuation_counter = 0;

% Intrusion method, similar to sill intrusion
%1. density contrast based, or last location if no den_melt<den_crust point found.  
%2. fixed depth
Intrusion_Evac_method=1;  


% Intrusion depth(s) for evacuations
% Intrusion_Evac = [-15000 -5000];
Intrusion_Evac = [-5000];
Last_Intrusion_Evac=min(Intrusion_Evac);


Intrusion_Evac=sort(Intrusion_Evac);

%The distance of top of the system and the intended intrusion depth must be
%larger than the following value in order for the evacuation to happen. 
%It's a list for all the intended intrusion depths
Safe_intrusion_evac_dist=5e3*ones(1,length(Intrusion_Evac));  


% %Intrusion depth 1 for evacuations - Nodez
% initial_T_depth_N_1 = find(nodez>=Intrusion_Evac_1, 1, 'first');
% %Intrusion depth 1 for evacuations -cellz
% initial_T_Depth_1 = find(cellz==((nodez(initial_T_depth_N_1)+nodez(initial_T_depth_N_1-1))/2));
% 
% %Intrusion depth 2 for evacuations - Nodez
% initial_T_depth_N_2 = find(nodez>=Intrusion_Evac_2, 1, 'first');
% %Intrusion depth 2 for evacuations -cellz
% initial_T_Depth_2 = find(cellz==((nodez(initial_T_depth_N_2)+nodez(initial_T_depth_N_2-1))/2));


Finding_porosity;
% CAB end
%% outputs
%CAB
Out_time = 5000;
output_counter=1;
% Set as 1, to turn on txt outputs
Record_data=0;
Record_index=2;
Record_time=[0:Out_time*Year:End_time];
restart_No=5;

%CAB end
%% save the original component mass data
Mass_data_old=Mass_data;
component_A_old=sum(Mass_data_old(:,1).*dz+Mass_data_old(:,4).*dz);
component_B_old=sum(Mass_data_old(:,2).*dz+Mass_data_old(:,5).*dz);
component_V_old=sum(Mass_data_old(:,3).*dz+Mass_data_old(:,6).*dz+Mass_data_old(:,7).*dz);
%%
% Record intrusion marker
To_recored_intrusion_marker=0;
if To_recored_intrusion_marker==1
    marker_ratio=0.1;

    marker_nodes=find(cellz>(cellz(Sill_index(1))+marker_ratio*Sill_length),1,'first')-Sill_index(1);
    Intrusion_marker_static=zeros(N,1);
    Intrusion_marker_dynamics=zeros(N,1);
    Intrusion_marker_static(Sill_index(1):Sill_index(1)+marker_nodes)=1;
    Intrusion_marker_dynamics(Sill_index(1):Sill_index(1)+marker_nodes)=1;
end

%% Set up plottings 
% labels for parameters to be plot:
% 1. porosity
% 19. volatile fraction S
% 2. liquid velocity
% 3. solid velocity
% 20 volotile velocity
% 4. enthalpy
% 5. temperature
% 6. liquid composition
% 7. solid composition
% 8. bulk composition 
% 24. volatile compositions: Cl2, Cs2
% 25. Mass per unit volume of major component M2

% 10. force contribution from solid viscous force
% 11. Heat-chemical-compaction contribution, stepwise.
% 12. Heat-chemical-compaction contribution, total (abs).
% 13. Refined recored of bulk composition
% 14. Refined recored of HCC contribution, total.
% 15. Heat-chemical-compaction contribution, total
% 16. Chemical differenciation from compacton and reactive flow, total
% 17. Chemical differenciation from compacton and reactive flow, stepwise
% 18. Chemical differenciation from compacton and reactive flow, total (abs)
% 21. Pressures (liquid, solid, volatile)
% 22. Pressures gradient(liquid, solid, volatile)
% 23. Densities

% 50. melt viscosity
% 51. contributions
% 101. Intrusion marker static
% 102. Intrusion marker dynamics

% 99. total mass, for testing

% Plot_Label={'$$\phi$$','uf','um','H','T','Cf','Cm','Cb','St', 'Solid_vis(%)','$$\Delta\phi$$','$$\sum|\Delta\phi|$$'};
% when multiple variable in one figure, temperature will be scaled between solidus and liquidus
% option 11 and 12 must stand alone!

% Plot_configure=...
%     {[1],[2],[4],[5];...
%      [6],[7],[8],[11]};   % a 2*4 plot setting
% Plot_configure=...
%     {[1,5,8],[10],[11]};  
running_plot=1; %CAB added, if 1 it will plot whilst running
SiO2_range=[47 74]; %CAB added,
Cb2=(sum(Mass_data(:,[3,6,7]),2))./sum(Mass_data,2); %CAB added,

Display_type=2;  %set to 2 if the computer screen only support 1024*768..ffs!

%  Plot_configure=  {[1],[19],[8,6,7],[24], [5]};  
if N_component==3 && Add_CLCU~=1
    Plot_configure=  {[1,19],[8,6,7],[24], [5]};
    % Plot_configure=  {[1,19],[8,6,7],[24], [101],[102]};
else
    Plot_configure=  {[1 19],[8,6,7],[24], 30, 31,[5]};
end
 Contribution_type=2;  %1, instant 2 accumulative 
 % Plot_configure=  {[1],[8,6,7], [5};  
% Plot_configure=  {[1],[19],[5],[8,6,7],[24], [2,3,20],[23], [99]};  

Z_label_type=2;  % 1 in meter, 2 in km start from 0 on top

if Z_label_type==2
    Show_z=(Show_z-nodez(end))/1000;
end

Adaptive_show_range=1; % Set to one to let the show range track the active region

u_all=zeros(3*N+3,1);
if running_plot==1
    Plot_settings;
    Update_plot;
    drawnow;
    Show_refresh=10;  %the refresh rate in second
    Frame_num=1;
end

tic;
%%  Record videos (data)
Create_video=1;
Fixed_record_dt0=2000*Year;
Fixed_record_dt=Fixed_record_dt0;
Fixed_record=1;

video_data_index=0;
video_data_file_name=['Data' num2str(video_data_index)];
save(video_data_file_name,'Time','cellz','Mass_data','Pg_real','Ts_local','Tl_local','T', 'S_cap','phi')

%% Save simulation data
Create_save_data=1;
if Create_save_data==1
%     Save_data_times=[10, 20 50 100]*Year;  %define the time where data are saved.
    % Save_data_times=[0:10e4:2e6 (2e6+5e4):5e4:5e6]*Year;
    Save_data_times=[0:20e4:(End_time*1.1/Year)]*Year;
    % Save_data_times=[0:2e3:2e6 (2e6+5e4):5e4:5e6]*Year;
    % Save_data_times=0:2e4
    Save_data_times=[Save_data_times End_time+1000*Year];
    Save_data_index=1;
end
%%
%Extract volatile phases
To_extract_volatile=1; %When set to 1, extract the volatile phases created and recored the collected volatile　
if To_extract_volatile==1
    Record_gap=1000*Year;
    Record_vol_time=0;
    fileID = fopen('Volatile_leak.txt', 'w');
    if fileID == -1
        error('File could not be opened');
    end
    
    Leak_V=0;
    Leak_H=0;
    if N_component==5 || Add_CLCU==1
        Leak_CU=0;
        Leak_CL=0;
        fprintf(fileID, 'Time(y)\t V\t Cu\t Cl\t H(J)\n');
        fprintf(fileID, '%d\t %.4f\t %.4f\t %.4f\t %.4f\n', 0, 0, 0, 0,0);
    else
        fprintf(fileID, 'Time(y)\t V\t H(J)\n');
        fprintf(fileID, '%d\t %.4f\t %.4f\n', 0, 0,0);
    end
end
%% extra parameters
Invis_vol=zeros(N,1);
Top0=0;

Convergence_record_steps=10;
Convergence_record=zeros(Convergence_record_steps,2);
Num_non_convergence=0;
file_iter = fopen('Iteration_recored.txt', 'w');



%% Debug parameters

found_error=0;


function x=Calculate_Initial_state(PD_range, rho_mean, SiO2_crust, H2O_crust, CL_crust, N_component,Conservation_type,Precision)
    err=1;
    if N_component==3
        x=[1500, 1500, 20]';  %(M,N,V)
        M=x(1); N=x(2); V=x(3);
        J=zeros(3);
        J(1,:)=[PD_range(2)-SiO2_crust, PD_range(1)-SiO2_crust, 0];  %(74*M+N*47)/(M+N)=SiO2_crust
        J(2,:)=[-H2O_crust, -H2O_crust, 100-H2O_crust];              %100*V/(M+N+V)=H2O_crust
        while err>Precision
            if Conservation_type==1
                J(3,:)=[rhom_1 - 2*N - V - 2*M, rhom_2 - 2*N - V - 2*M, - M - N];    %rhom_1*N+rhom_2*M=(M+ N+V)*(M+N)
                rhs=[(PD_range(2)-SiO2_crust)*M+ (PD_range(1)-SiO2_crust)*N;...
                    (100-H2O_crust)*V-H2O_crust*(M+N);...
                    rhom_1*N+rhom_2*M-(M+N+V)*(M+N)];
            else
                J(3,:)=[1,1,1];    %M+N+V=rho_mean
                rhs=[(PD_range(2)-SiO2_crust)*M+ (PD_range(1)-SiO2_crust)*N;...
                    (100-H2O_crust)*V-H2O_crust*(M+N);...
                    M+N+V-rho_mean];
            end
            dx=J\rhs;
            x=x-dx;
            err=max(abs(dx(:)));
            M=x(1); N=x(2); V=x(3);
        end
        % MM_crust=x(1);
        % NN_crust=x(2);
        % V_crust=x(3);
    else
        x=[1500, 1500, 20, 3]';  %(M,N,V,CL)
        M=x(1); N=x(2); V=x(3); CL=x(4);
        J=zeros(4);
        J(1,:)=[PD_range(2)-SiO2_crust, PD_range(1)-SiO2_crust, 0, 0];  %(74*M+N*47)/(M+N)=SiO2_crust
        J(2,:)=[-H2O_crust, -H2O_crust, 100-H2O_crust, -H2O_crust];     %100*V/(M+N+V+CL)=H2O_crust
        J(3,:)=[-CL_crust, -CL_crust,-CL_crust, 100-CL_crust]; %100*CL/(M+N+V+CL)=CL_crust
        while err>Precision
            if Conservation_type==1
                J(3,:)=[rhom_1 - 2*N - V - 2*M, rhom_2 - 2*N - V - 2*M, - M - N];    %rhom_1*N+rhom_2*M=(M+ N+V)*(M+N)
                rhs=[(PD_range(2)-SiO2_crust)*M+ (PD_range(1)-SiO2_crust)*N;...
                    (100-H2O_crust)*V-H2O_crust*(M+N);...
                    rhom_1*N+rhom_2*M-(M+N+V)*(M+N)];
            else
                J(4,:)=[1,1,1,1];    %M+N+V+CL=rho_mean
                
                rhs=[(PD_range(2)-SiO2_crust)*M+ (PD_range(1)-SiO2_crust)*N;...
                    (100-H2O_crust)*V-H2O_crust*(M+N+CL);...
                    (100-CL_crust)*CL-CL_crust*(M+N+V)
                    M+N+V+CL-rho_mean];
            end
            dx=J\rhs;
            x=x-dx;
            err=max(abs(dx(:)));
            M=x(1); N=x(2); V=x(3); CL=x(4);
        end
    end
end


