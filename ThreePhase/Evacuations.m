% CAB 19/8/25
% Evacuation of magma and intrusion of the evacuated magma script.

melt_base = buoy_Hphi_base(j);
melt_top = buoy_Hphi_top(j);


% thickness that is evacuated
Thick_melt_evac = abs(melt_base-melt_top)+1;
Thick_melt_m = abs(cellz(melt_top-1) - cellz(melt_base)); %Thick_melt_evac;

%% Amount transferred and depth intruded
% finding the intrusion depth
deni = -5000;
while deni ==-5000
    for i=melt_top+1:N%-Thick_melt_evac
        
       % if rho_on_cellz(i)<=av_rho_T  % av_rho_magma will be calculated outside this script
        %    depth_int_N_c = i;
         %   deni=0;
        %end

        if i==N%-Thick_melt_evac
            if melt_top<initial_T_Depth_1
                depth_int_N_c = initial_T_depth_N_1;%-Thick_melt_evac;
            else
                depth_int_N_c = initial_T_depth_N_2;%-Thick_melt_evac;
            end
             % initial depth for transfer
            deni=0;
        end
        
        if deni==0
            break
        end
    
    end
end

% add focusing here if wanted.
% calculatig the thickness intruded

%Thick_melt_int_calc=0;
%Thick_melt_m_count=0;
%while Thick_melt_m_count<Thick_melt_m
%    for i=depth_int_N_c:N
%        Thick_melt_int_calc = Thick_melt_int_calc+1;
 %       Thick_melt_m_count = Thick_melt_m_count+(dz(i-1));
 %       if Thick_melt_m_count>=Thick_melt_m
  %          break
  %      end
  %  end

   % if i==N
   %     fprintf(File_echo, '%6s %5.5f %6s \n', 'ERROR IN EVACUATION');
   %     Thick_melt_m_count = Thick_melt_m;
   %     break
  %  end
        
% end

Thick_melt_int = Thick_melt_evac; %(depth_int_N_c + Thick_melt_int_calc) - depth_int_N_c;




%% Evacuation of magma (removing nodes)

if Z_label_type==1
    dist_cellz_removed = cellz(melt_base:melt_top);
    dist_nodez_removed = nodez(melt_base:melt_top);
else
    dist_cellz_removed = (cellz(melt_base:melt_top)-nodez(end))/1000;
   dist_nodez_removed = (nodez(melt_base:melt_top)-nodez(end))/1000;
end
    


% Nx1 arrays
dz_removed = dz(melt_base-1:melt_top-1);
dz(melt_base-1:melt_top-1) = [];



rho_removed = rho(melt_base:melt_top,:);
%av_rho_T = mean(rho_removed);
rho(melt_base:melt_top,:)=[];

OG_Cb_removed = OG_Cb(melt_base:melt_top);
OG_Cb(melt_base:melt_top)=[];


Mass_data_removed = Mass_data(melt_base:melt_top,:);
%av_Mass_data_T = mean(Mass_data_removed);
Mass_data(melt_base:melt_top,:)=[];

H_removed = H(melt_base:melt_top);
%av_H_T = mean(H_removed);
H(melt_base:melt_top)=[];

if Add_CLCU==1 
    Partition_CL_CU_removed = Partition_CL_CU(melt_base:melt_top,:);
    %av_Partition_CL_CI_T = mean(Partition_CL_CU_removed);
    Partition_CL_CU(melt_base:melt_top,:)=[];
else
    Partition_CL_CU_removed=zeros(melt_top-melt_base+1,2);    
    Partition_CL_CU=zeros(N,2);
end

T_removed = T(melt_base:melt_top);
%av_T_T = mean(T_removed);
T(melt_base:melt_top)=[];

Ts_local_removed = Ts_local(melt_base:melt_top);
%av_Ts_local_T = mean(Ts_local_removed);
Ts_local(melt_base:melt_top)=[];

Pg_real_removed = Pg_real(melt_base:melt_top);
%av_Pg_real_T = mean(Pg_real_removed);
Pg_real(melt_base:melt_top)=[];

S_cap_removed = S_cap(melt_base:melt_top);
%av_S_cap_T = mean(S_cap_removed);
S_cap(melt_base:melt_top)=[];

T_S_region_removed = T_S_region(melt_base:melt_top,:);
%av_T_S_region_T = mean(T_S_region_removed);
T_S_region(melt_base:melt_top,:)=[];

% 2Nx1 arrays

phi_removed = phi(melt_base:melt_top);
%av_phi_T = mean(phi_removed);
phi(melt_base:melt_top)=[];

S_removed = S(melt_base:melt_top);
%av_S_T = mean(S_removed);
S(melt_base:melt_top)=[];


% Caculate cb_removed for outputs only.
Cb_removed=(PD_range(2)*sum(Mass_data_removed(:,[1,4]),2)+PD_range(1)*sum(Mass_data_removed(:,[2,5]),2))./sum(Mass_data_removed(:,[1 2 4 5]),2);



%% Adapt intrusion depth

depth_int_N_c = depth_int_N_c - (melt_top-melt_base);
        

%% Intrusion of magma (intruding nodes)

dz = [dz(1:depth_int_N_c-2); dzf.*ones(Thick_melt_int,1); dz(depth_int_N_c-1:end)];

nodez = zeros(N+1+(Thick_melt_int - Thick_melt_evac),1);

for i=2:N+1+(Thick_melt_int - Thick_melt_evac)
    nodez(i) = nodez(i-1)+dz(i-1);
end

N=length(nodez)-1;

cellz = (nodez(1:end-1)+nodez(2:end))/2;

% Intrusion in Nx1 arrays
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Densities
rho=[rho(1:depth_int_N_c,:);  Sillrho.*ones(Thick_melt_int,1); rho(depth_int_N_c+1:end,:)];

%OG_Cb
OG_Cb=[OG_Cb(1:depth_int_N_c,:);  av_OG_Cb_T.*ones(Thick_melt_int,1); OG_Cb(depth_int_N_c+1:end,:)];

%Massdata
Mass_data=[Mass_data(1:depth_int_N_c,:);   av_Mass_data_T.*ones(Thick_melt_int,1); Mass_data(depth_int_N_c+1:end,:)];

%Enthalpy
H = [H(1:depth_int_N_c); av_H_T.*ones(Thick_melt_int,1);H(depth_int_N_c+1:end)];

if Add_CLCU==1
    Partition_CL_CU=[Partition_CL_CU(1:depth_int_N_c,:); nan(Thick_melt_int,2); Partition_CL_CU(depth_int_N_c+1:end,:)]; 
end

%Temperature
T = [T(1:depth_int_N_c); av_T_T.*ones(Thick_melt_int,1);T(depth_int_N_c+1:end)];

Ts_local = [Ts_local(1:depth_int_N_c); av_Ts_local_T.*ones(Thick_melt_int,1);Ts_local(depth_int_N_c+1:end)];  %just to fill the data, filled Ts_local is not correct value

Pg_real=[Pg_real(1:depth_int_N_c); Pg_real(depth_int_N_c).*ones(Thick_melt_int,1); Pg_real(depth_int_N_c+1:end)];
S_cap=[S_cap(1:depth_int_N_c); S_cap(depth_int_N_c).*ones(Thick_melt_int,1); S_cap(depth_int_N_c+1:end)];

T_S_region=[T_S_region(1:depth_int_N_c,:); T_S_region(depth_int_N_c,:).*ones(Thick_melt_int,1); T_S_region(depth_int_N_c+1:end,:)];


% do we need intrusion marker static and intrusion marker dynamics.

% Intrusion in 2Nx1 arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%Melt Fraction / Solid Fraction
phi = [phi(1:depth_int_N_c); av_phi_T.*ones(Thick_melt_int,1); phi(depth_int_N_c+1:end)]; 


%volatile saturation
S = [S(1:depth_int_N_c); av_S_T.*ones(Thick_melt_int,1); S(depth_int_N_c+1:end)]; 



%% Outputs - needs to be numbered


%list of data removed

%cellz, dz, rho
Filename_remove = ['magma_removed_',num2str(evacuation_counter),'.txt'];
File_OutR = fopen(Filename_remove, 'w');

fprintf(File_OutR,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );



% mass_data = [m1, n1, v1, m2, n2, v2, v3, l1, l2, l3, u1, u2, u3]

fprintf(File_OutR, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s\n', ...
    'Depth (km)', 'dz', 'phi', 'S', 'H', 'T', 'cb', 'rho_b', 'rho_l', 'rho_s', 'rho_v', 'Ts_local', 'Pg_real', 'S_cap', 'T_S_region', 'CL', 'CU');


for i=1:length(dz_removed)
    fprintf(File_OutR, '%10.4f \t %10.4f \t %10.4f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
        [(cellz(melt_base-1+i)-nodez(end))/1000; dz_removed(i); phi_removed(i); S_removed(i); H_removed(i); T_removed(i); Cb_removed(i); rho_on_cellz(melt_base-1+i); rho_removed(i,1); rho_removed(i,2); rho_removed(i,3); ...
        Ts_local_removed(i); Pg_real_removed(i); S_cap_removed(i); T_S_region_removed(i); Partition_CL_CU_removed(i,1); Partition_CL_CU_removed(i,2)]);
end



fclose(File_OutR)



Filename_MD_remove = ['mass_data_removed_',num2str(evacuation_counter),'.txt'];
File_Out_Rem = fopen(Filename_MD_remove, 'w');

fprintf(File_Out_Rem,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );

if Add_CLCU==1 

    fprintf(File_Out_Rem, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'Depth (km)', 'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3', 'l1', 'l2', 'l3', 'u1', 'u2', 'u3');

    for i=1:length(dz_removed)
        fprintf(File_Out_Rem, '%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f   \n', ...
            [(cellz(melt_base-1+i)-nodez(end))/1000; Mass_data_removed(i,1);Mass_data_removed(i,2);Mass_data_removed(i,3);Mass_data_removed(i,4);Mass_data_removed(i,5);Mass_data_removed(i,6);Mass_data_removed(i,7); ...
            Mass_data_removed(i,8);Mass_data_removed(i,9);Mass_data_removed(i,10);Mass_data_removed(i,11);Mass_data_removed(i,12);Mass_data_removed(i,13)]);
    end



else

    fprintf(File_Out_Rem, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'Depth (km)', 'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3');

    for i=1:length(dz_removed)
        fprintf(File_Out_Rem, '%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \n', ...
            [cellz(melt_base-1+i); Mass_data_removed(i,1);Mass_data_removed(i,2);Mass_data_removed(i,3);Mass_data_removed(i,4);Mass_data_removed(i,5);Mass_data_removed(i,6);Mass_data_removed(i,7)]);
    end

end


fclose(File_Out_Rem);

melt_evacuated_top = cellz(melt_base-1+length(dz_removed));







