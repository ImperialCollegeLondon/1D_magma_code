% function [av_rho_T, av_Mass_data_T, av_H_T, av_T_T, av_Ts_local_T, av_phi_T, av_S_T, av_cb_T, av_OG_Cb_T, av_S_cap,av_Tl_local_T,av_Partition_CL_CU ] = Averages_evacuated_updated(evacuation_counter, melt_base, melt_top, phi, S, H, Mass_data, ...
%                                                                                                     Pg_real, dz, MM, NN, V, CL, CU, Add_CLCU, Lf, cp, simplified_TS, Sys_constant, rho, Melt_density_type, ...
%                                                                                                         PD_range, Coef_melt_den, cellz, nodez, aa, bb,cc,rho_mean,g, Velocity_solver_type, Ts, Tl, ...
%                                                                                                         Conservation_type, rhom_1, rhom_2, Time, Year,OG_Cb)

%% HH 26/12/2025 An conservative averaging method updated from Averages_evacuated 

melt_base= buoy_Hphi_base(j);
melt_top=buoy_Hphi_top(j);
N_temp=abs(melt_base-melt_top)+1;  %number of cells after average
Length=sum(dz(melt_base:melt_top));
dze=Length/N_temp;
%% Calculate the total averaged mass and energy data of the cells to be averaged
% calculate the sum of each major element components
% Mass data = [m1 n1 v1 m2 n2 v2 v3] if no chlorine copper
% Mass data = [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3] if chlorine copper
Sum_mass=sum(Mass_data(melt_base:melt_top,:).*dz(melt_base:melt_top),1)/Length;
MM_temp=(Sum_mass(1)+Sum_mass(4));
NN_temp=(Sum_mass(2)+Sum_mass(5));
V_temp=(Sum_mass(3)+Sum_mass(6)+Sum_mass(7));
if Add_CLCU==1 
    CL_temp=(Sum_mass(8)+Sum_mass(9)+Sum_mass(10));
    CU_temp=(Sum_mass(11)+Sum_mass(12)+Sum_mass(13));
end


% total enthalpy
H_temp=(sum(H(melt_base:melt_top).*dz(melt_base:melt_top)))/Length;
av_H_T=H_temp; 
%% Call chemical model to recalculate all the system data
% Assume there is no pressure dependant change.
% To be most accurate, we need to fill in the pressure at the depth of the reinjection

if Add_CLCU==1
    input=[H_temp,MM_temp,NN_temp,V_temp,CL_temp,CU_temp];
else
    input=[H_temp,MM_temp,NN_temp,V_temp];
end
[av_phi_T,av_S_T, av_T_T, av_rho_T,~,~,~,av_Ts_local_T, av_Mass_data_T,av_S_cap, av_Tl_local_T,~,~, av_Partition_CL_CU]=Phase_component_updated_solid3(input,Pg_real(melt_top), 800, 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(melt_top,:),   rho_constant, Lf, [K(melt_top) Kcl Kcu], Ts, Tl, min_por, max_melt, dz(melt_top),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS,0);

av_mass_addition=sum(Mass_data_addition(melt_base:melt_top,:).*dz(melt_base:melt_top),1)/Length;
%% Average OG_cb  %still kept as the old way
av_OG_Cb_T = mean(OG_Cb(melt_base:melt_top));


%% Outputs
av_cb_T = (PD_range(2)*sum(av_Mass_data_T(:,[1,4]),2)+PD_range(1)*sum(av_Mass_data_T(:,[2,5]),2))./sum(av_Mass_data_T(:,[1 2 4 5]),2);


% Filename_intrude = ['magma_avintruded_',num2str(evacuation_counter),'.txt'];
% File_In = fopen(Filename_intrude, 'w');
% 
% fprintf(File_In,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );
% 
% % 1. av_H_T, 2. av_phi_T, 3. av_S_T, 4. av_T_T,  5. av_cb_T, 6. av_rho_T,
% % 7. av_Ts_local_T 
% 
% fprintf(File_In, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
%     'Enthalpy', 'Porosity', 'Vol frac', 'Temp', 'Bulk comp', 'Bulk rho', 'Ts_local');
% 
% % fprintf(File_In, '%3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
% %    [av_H_T; av_phi_T; av_S_T; av_T_T; av_cb_T; av_rho_T; av_Ts_local_T]); 
% 
% fclose(File_In);





% Filename_MD_intrude = ['mass_data_avintruded_',num2str(evacuation_counter),'.txt'];
% File_Out_Iem = fopen(Filename_MD_intrude, 'w');
% 
% fprintf(File_Out_Iem,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );
% 
% if Add_CLCU==1 
% 
%     fprintf(File_Out_Iem, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
%         'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3', 'l1', 'l2', 'l3', 'u1', 'u2', 'u3');
% 
% 
%     fprintf(File_Out_Iem, ' %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \n', ...
%             [ av_Mass_data_T(:,1);av_Mass_data_T(:,2);av_Mass_data_T(:,3);av_Mass_data_T(:,4);av_Mass_data_T(:,5);av_Mass_data_T(:,6);av_Mass_data_T(:,7); ...
%             av_Mass_data_T(:,8);av_Mass_data_T(:,9);av_Mass_data_T(:,10);av_Mass_data_T(:,11);av_Mass_data_T(:,12);av_Mass_data_T(:,13)]);
% 
% 
% 
% 
% else
% 
%     fprintf(File_Out_Iem, ' %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
%          'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3');
% 
% 
%     fprintf(File_Out_Iem, ' %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \n', ...
%             [av_Mass_data_T(:,1);av_Mass_data_T(:,2);av_Mass_data_T(:,3);av_Mass_data_T(:,4);av_Mass_data_T(:,5);av_Mass_data_T(:,6);av_Mass_data_T(:,7)]);
% 
% 
% end
% 
% 
% fclose(File_Out_Iem);


