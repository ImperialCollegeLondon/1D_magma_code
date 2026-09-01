%CAB 30/1/23
%Sill intrusion function file.
% (N_component,dz, Old_end, N, depthN, dzf, H, T, Ts, Ts_local,phi, S, Mass_data, rho, SillNodez, SillMass, SillH, SillT, SillTs, Sillphi, SillS, Sillrho, Pg_real, S_cap, T_S_region, Invis_vol)
%%
%Inserts nodez/cellz with the properties of sill intrusions into the
%variables.
%% Update mesh parameters

% dz
dz=[dz(1:depthN-1); dzf.*ones(SillNodez,1); dz(depthN:end)];

% nodez
nodez=zeros(N+1+SillNodez,1);
for i=2:N+1+SillNodez
    nodez(i)=nodez(i-1)+dz(i-1);
end

%N
N=length(nodez)-1;

%cellz
cellz = (nodez(1:end-1)+nodez(2:end))/2; % cell center points. 

Start_depthN_c=depthN-1;

%% Intrusion in Nx1 arrays

%Densities
rho=[rho(1:Start_depthN_c,:);   Sillrho.*ones(SillNodez,1); rho(Start_depthN_c+1:end,:)];

%Massdata
Mass_data=[Mass_data(1:Start_depthN_c,:);   SillMass.*ones(SillNodez,1); Mass_data(Start_depthN_c+1:end,:)];
Mass_data_addition=[Mass_data_addition(1:Start_depthN_c,:);   SillMass_addition.*ones(SillNodez,1); Mass_data_addition(Start_depthN_c+1:end,:)];

%Enthalpy
H = [H(1:Start_depthN_c); SillH.*ones(SillNodez,1);H(Start_depthN_c+1:end)];

if N_component==5
    Partition_CL_CU=[Partition_CL_CU(1:Start_depthN_c,:); nan(SillNodez,2); Partition_CL_CU(Start_depthN_c+1:end,:)]; 
end

%Temperature
T = [T(1:Start_depthN_c); SillT.*ones(SillNodez,1);T(Start_depthN_c+1:end)];
% Ts = [Ts(1:Start_depthN_c); SillTs.*ones(SillNodez,1);Ts(Start_depthN_c+1:end)];
Ts_local = [Ts_local(1:Start_depthN_c); SillTs.*ones(SillNodez,1);Ts_local(Start_depthN_c+1:end)];  %just to fill the data, filled Ts_local is not correct value

Pg_real=[Pg_real(1:Start_depthN_c); Pg_real(Start_depthN_c).*ones(SillNodez,1); Pg_real(Start_depthN_c+1:end)];
S_cap=[S_cap(1:Start_depthN_c); S_cap(Start_depthN_c).*ones(SillNodez,1); S_cap(Start_depthN_c+1:end)];

T_S_region=[T_S_region(1:Start_depthN_c,:); T_S_region(Start_depthN_c,:).*ones(SillNodez,1); T_S_region(Start_depthN_c+1:end,:)];


% Invis_vol=[Invis_vol(1:Start_depthN_c); zeros(SillNodez,1); Invis_vol(Start_depthN_c+1:end)];
if To_recored_intrusion_marker==1
    Intrusion_marker_static=[Intrusion_marker_static(1:Start_depthN_c); Intrusion_marker_static(Start_depthN_c).*zeros(SillNodez,1); Intrusion_marker_static(Start_depthN_c+1:end)];
    Intrusion_marker_static(Start_depthN_c:Start_depthN_c+marker_nodes)=1;

    Intrusion_marker_dynamics=[Intrusion_marker_dynamics(1:Start_depthN_c); Intrusion_marker_dynamics(Start_depthN_c).*zeros(SillNodez,1); Intrusion_marker_dynamics(Start_depthN_c+1:end)];
    Intrusion_marker_dynamics(Start_depthN_c:Start_depthN_c+marker_nodes)=1;
end
%% Intrusion in 2Nx1 arrays

%Melt Fraction / Solid Fractio
phi = [phi(1:Start_depthN_c); Sillphi.*ones(SillNodez,1); phi(Start_depthN_c+1:end)]; 


%volatile saturation
S = [S(1:Start_depthN_c); SillS.*ones(SillNodez,1); S(Start_depthN_c+1:end)]; 

%% Adapting evacuation 
% CAB added
if ~isempty(all_phi_buoy_top)

    for i = 1:length(all_phi_buoy_top)
        if all_phi_buoy_top(i)>Start_depthN_c
            all_phi_buoy_top(i) = all_phi_buoy_top(i) + SillNodez;
        end

        if all_phi_buoy_base(i)>Start_depthN_c
            all_phi_buoy_base(i) = all_phi_buoy_base(i) + SillNodez;
        end
    end

end
% CAB end

%% Calculating mass conservation
Sum_M0=Sum_M0+MM_sill*dzf*SillNodez;
Sum_N0=Sum_N0+NN_sill*dzf*SillNodez;
Sum_V0=Sum_V0 +V_sill*dzf*SillNodez;

Sillcomp =(PD_range(2)*sum(SillMass(:,[1,4]),2)+PD_range(1)*sum(SillMass(:,[2,5]),2))./sum(SillMass(:,[1 2 4 5]),2);
OG_Cb = [OG_Cb(1:Start_depthN_c); Sillcomp.*ones(SillNodez,1);OG_Cb(Start_depthN_c+1:end)];