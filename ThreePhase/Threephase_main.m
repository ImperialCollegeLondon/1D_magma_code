% HH 5/12/2021
% A 1D three phase solver for magma simulation.
% fraction change
%% 
clear;

To_Restart=0;  %set to 1 to restart the simulation from a saved state
if To_Restart==1
    Load_data_index=20;  % The restart step to load 
end



if To_Restart==0

    % Model_setting_threephase_sill_multi_intrustion2;    
    Model_setting_threephase_sill_multi_intrustion3;
    %%    
%     Plot_configure=  {[1],[19],[5],[8,6,7],[24], [2,3,20],[23], [99]}; 
    
    if running_plot==1
        Plot_settings
        Update_plot
        tic
        Frame_num=1;
    end


    
    if exist('Fixed_record','var')
        if  Fixed_record==1
            Last_record_time=0;
        end
    end

    counter=0;
    break_counter=0;
    break_percent=0;
    output_counter=1;
    iter_chemical=zeros(N+1,1);
    C_values=zeros(N+1,3);
    dP=zeros(N+1,1);

    %CAB - mass cons - original Cb. Only updated in sill intrusion and 
    OG_Cb = Cb;
    OG_composition=sum(Cb);

 
 %% CAB - outputs
 % % Overal model data % % %

    if Record_data ==1
        % CELLs
    
        FilenameO =  'output_0_CELLS.txt';
    
        File_Out = fopen(FilenameO, 'w');
    
        fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');
    
        fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \n', ...
            'Depth', 'phi', 'Vol. frac', 'Enth.', 'Temp (C)', 'Ts (C)', 'Tl (C)', 'cb', 'cl', 'cs', 'vb','vl', 'vs','vg') ;
    
        for i=1:1:N
            fprintf(File_Out, '%10.4f \t %10.4f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f\n ', ...
                [((cellz(i)-nodez(end))/1000); phi(i); S(i); H(i); T(i); ...
                Ts_local(i); Tl_local(i); Cb(i); (SiO2_range(2)*Mass_data(i,1) + SiO2_range(1)*Mass_data(i,2))/(sum(Mass_data(i,1:2),2)); (SiO2_range(2)*Mass_data(i,4) + SiO2_range(1)*Mass_data(i,5))/(sum(Mass_data(i,4:5),2)); ...
                Cb2(i); Mass_data(i,1)/sum(Mass_data(i,[1,2]),2); Mass_data(i,4)/sum(Mass_data(i,[4,5]),2); (Mass_data(i,7)+max(Mass_data(i,6)-Mass_solid(i)*S_cap(i),0))/rho_mean]);
    
    
        end
    
        fclose(File_Out);
    
        %NODES
    
        FilenameON = 'output_0_NODES.txt';
        
        File_OutN = fopen(FilenameON, 'w');
    
        fprintf(File_OutN, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');
    
        fprintf(File_OutN, '%10s \t %10s \t %10s \t %10s  \n', ...
            'Depth', 'Solid vel', 'Liq. vel', 'Vol. vel') ;
    
        for i=1:1:N
            fprintf(File_OutN, '%10.4f \t %4e \t %4e\t %4e  \n', ...
                [((nodez(i)-nodez(end))/1000); u_all(N+1+i); u_all(i); 0]);
        end

        fclose(File_OutN);
    end


    % ongoing text file for the averages evacuated 
    Filename_AVevac = 'Averages_evacuated.txt';
    File_AVevac = fopen(Filename_AVevac, 'w');
    fprintf(File_AVevac, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'Time(ka)', 'EvacFrom(km)', 'EvacTo(km)', 'Thick(km)', 'VolumeKm3', 'Comp', 'NodesEvac', 'NodesInt');


    % ongoing file to write information that comes up
    File_echo=fopen('echo_screen.txt','w');

    fprintf(File_echo, '%s %5.5f %s %5.5f %6s\n', 'Sill intruded at', ...
                            (cellz(Sill_index(2))-nodez(end))/1000, 'km at', Time/Year/1000, 'ka');


    % ongoing file to output when tolerance is not met
    File_precision  = fopen('In_Depth_break_data.txt', 'w');
    fprintf(File_precision, '%10s \t %10s \t %10s \t \n', 'Time (ka)', 'Precision', 'Residual');

    %outputted at every output of data
    File_breaks = fopen('Break_of_tolerance.txt','w');
    fprintf(File_breaks, '%10s \t %10s \t %10s \t %10s \n ', 'Time (ka)', 'Percentage', 'no. breaks', 'no.timesteps');

    % ongoing text file for tracking mass conservation
    File_comp=fopen('mass_cons.txt','w');
    fprintf(File_comp, '%10s \t %10s \t %10s \n', 'Time (ka)', 'Mass cons', 'Mass Res' );

% CAB END
%%

else
    eval(['load(''Record_' num2str(Load_data_index) '.mat'');']);
    % load('Record_error.mat');
    % Transport_method=2;
    
    File_echo=fopen('echo_screen.txt','w');
    File_precision  = fopen('In_Depth_break_data.txt', 'w');

    warning('off');
    if To_extract_volatile==1
        fileID = fopen('Volatile_leak.txt', 'a');
    end
    running_plot=1;
    Courant0=6.5e-1;
%     Plot_configure=  {[1],[19],[5],[8,6,7],[24], [2,3,20],[23], [99]}; 

    if running_plot==1
        Plot_settings
        Update_plot
        tic
        Frame_num=1;
    end
end
%%
Adjusted=0;
while Time<End_time
    iter=0;
    eps=1;
    
    Just_intruded=0;
    if Time>injection_time(Injected)  && To_intrude==1
        if Conservation_type==1 || Time==0
            sum_mass=sum(Mass_data,2);
        else
            m1=Mass_data(:,1); n1=Mass_data(:,2); v1=Mass_data(:,3); m2=Mass_data(:,4); n2=Mass_data(:,5); v2=Mass_data(:,6); v3=Mass_data(:,7);
            sum_mass=((rhof_1*m1+rhof_1*n1)./max(m1+n1,1e-3).*(m1+n1+v1)+(rhom_1*m2+rhom_2*n2)./max(m2+n2,1e-3).*(m2+n2+v2))./(m1+n1+v1+m2+n2+v2+v3);
        end
        index=find(cellz>Dry_depth,1,'first');
        depthN0=find(Sillden>sum_mass(index:end),1,'first');
        depthN0=depthN0+index;
        if Injection_depth_method==99            
            depthN=find(nodez>=(nodez(end)-Injection_depth)+(rand-0.5)*Random_shift*2, 1, 'first'); 
        elseif Injection_depth_method==2
            [ ~,depthN]=min(sum_mass);
        elseif Injection_depth_method==3
            temp=find(phi>phi_crit,1,'first');
            if isempty(temp)
                % intrude at the same depth
                depthN=find(nodez>=(nodez(end)-Injection_depth)+(rand-0.5)*Random_shift*2, 1, 'first');
            else
                depthN=find(nodez>=cellz(temp)+(rand-0.5)*Random_shift*2, 1, 'first'); 
                Injection_depth=nodez(end)-cellz(temp);
            end
        else
            if ~isempty(depthN0)
                Injection_depth=cellz(depthN0)+(rand-0.5)*Random_shift*2;
                depthN=find(nodez>=Injection_depth, 1, 'first');
            else
                if Injection_depth_method==0
                    [ ~,depthN]=min(sum_mass);
                else
                    depthN=find(nodez>=Injection_depth, 1, 'first'); %inject at the same depth, but need to find the nodes again
                end
            end
        end
        
        fprintf(File_echo, '%s %5.5f %s %5.5f %6s\n', 'Sill intruded at', ...
                            (nodez(depthN)-nodez(end))/1000, 'km at', Time/Year/1000, 'ka');

        SillTs=Cal_solidus(sum(SillMass([3,6,7]))/sum(SillMass)*100/K0, Pg_real(depthN)*1e3 ,Data_point,Data_y);
        if is_solid_solution==0
            [dz, nodez, N, cellz, H, T, Ts_new, phi, S, Mass_data, rho, Pg_real, S_cap] = sill_intrusion_break(dz, nodez(end), N, depthN, dzf, H, T, Ts_new, phi, S, Mass_data, rho, SillNodez, SillMass, SillH, SillT, SillTs, Sillphi, SillS, Sillrho, Pg_real, S_cap);
        else
                sill_intrusion_break_solid2;
        end
        % Update the bottom cell temperature and enthalpy
        T_botom=min(abs(cellz(1)-nodez(end))/1000*24.5+10,1100);
        H(1)=T_botom*cp(1,2)*rho_mean;
%         Injection_depth=Injection_depth+Sill_length*0.7;
        cp=ones(N,1)*cp(1,:);
        Just_intruded=1;

        Last_intrude_time=Time;
        Injected=Injected+1;

        Fixed_record_dt=min(Injected*Fixed_record_dt0, Fixed_record_dt0*2);

        if Time/Year>2e6
            Fixed_record_dt=2000*Year;
        end
        
        if Z_label_type==1
            Show_z(1)=Show_z(1)-Sill_length;
        else
            Show_z(1)=Show_z(1)-Sill_length/1000;
        end

        Adaptive_step_time=0.1*Year;
        Adaptive_step_gap=10;
        % Courant=Courant0;
    end

    if (Adaptive_mesh==1 && mod(counter,Adaptive_step_gap)==1 && Time-Last_adapted>Adaptive_step_time) || Just_intruded==1
        %         save('temp.mat')
        
        % old_check_buoy_top_km = cellz(all_phi_buoy_top)-nodez(end);
        % old_check_buoy_base_km = cellz(all_phi_buoy_base)-nodez(end);
        
        Adapt_mesh_3p;
       
        Calculate_capped_values2;

        Last_adapted=Time;
        if Time-Last_intrude_time>1000*Year
            Adaptive_step_time=1*Year;
            Adaptive_step_gap=Adaptive_step_gap0;
            % Courant=Courant1;
        end



        % CAB EDIT - needed for finding_porosity 
         rho_on_nodez=zeros(N+1,3);
         if To_evacuate==1
             if Conservation_type==1
                 rho_on_nodez(:,1)=interp1(cellz,rho(:,1),nodez,'linear','extrap');
                 rho_on_nodez(:,2)=interp1(cellz,rho(:,2),nodez,'linear','extrap');
                 rho_on_nodez(:,3)=interp1(cellz,rho(:,3),nodez,'linear','extrap');
                 rho_bulk=phi.*(1-S).*rho(:,1)+(1-phi).*rho(:,2)+phi.*S.*rho(:,3);
             else % for volume conservation, recalulate the dyanmic denties used for gravitational term
                 if Melt_density_type==2
                     SiO2=(Mass_data(:,1)*PD_range(2)+Mass_data(:,2)*PD_range(1))./max(sum(Mass_data(:,1:2),2),1e-5)/100;
                     H2O=Mass_data(:,3)./max(sum(Mass_data(:,1:3),2),1e-5);
                     VMX=H2O.*SiO2*Coef_melt_den(1)+SiO2*Coef_melt_den(2)+H2O*Coef_melt_den(3)+Coef_melt_den(4);
                     x_other=1-SiO2-H2O;
                     V_all=SiO2.*26.86e-6/0.06009+H2O.*26.27e-6/0.01802+x_other.*VMX;
                     rho_temp=1./V_all;
                     rho_on_nodez(2:end-1,1)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                     rho_on_nodez(1,1)=rho_temp(1);  rho_on_nodez(end,1)=rho_temp(end);
                 else
                     rho_temp=(Mass_data(:,1)*rhof_1+Mass_data(:,2)*rhof_2)./max(sum(Mass_data(:,1:2),2),1e-5);
                     rho_on_nodez(2:end-1,1)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                     rho_on_nodez(1,1)=rho_temp(1);  rho_on_nodez(end,1)=rho_temp(end);
                 end

                 rho_temp=(Mass_data(:,4)*rhom_1+Mass_data(:,5)*rhom_2)./max(sum(Mass_data(:,4:5),2),1e-5);
                 rho_on_nodez(2:end-1,2)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                 rho_on_nodez(1,2)=rho_temp(1);  rho_on_nodez(end,2)=rho_temp(end);

                 rho_on_nodez(rho_on_nodez(:,2)==0,2)=rhom_2; % patch the phi=1 case

                 P2=-rho_mean*g*(cellz-nodez(end))/1e5; %in bar
                 rho_temp=(aa(1)* max(T,300).^bb(1)+aa(2)*max(P2,400).^bb(2)+aa(3)*max(T,300).^bb(3).*max(P2,400).^cc(1))*1000;
                 rho_on_nodez(2:end-1,3)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                 rho_on_nodez(1,3)=rho_temp(1);  rho_on_nodez(end,3)=rho_temp(end);
                 % if Velocity_solver_type==2
                 %     rho_on_nodez(S_on_nodez>0,1)=rho_on_nodez(S_on_nodez>0,3);
                 % end
                 rho_on_cellz = interp1(nodez, rho_on_nodez, cellz, 'linear', 'extrap');
                 rho_bulk = phi.*(1-S).*rho_on_cellz(:,1)+(1-phi).*rho_on_cellz(:,2)+phi.*S.*rho_on_cellz(:,3);
             end


             % Finding_porosity_adaptmesh;
         end
        %CAB end
    end

    u_all_old=u_all;
    phi_old=phi;
    S_old=S;

    rho_old=rho;
    P_old=P;
    H_old=H;
    Mass_data_old=Mass_data;
    Ts_local_old=Ts_local;
    Tl_local_old=Tl_local;

    T_old=T;
    
    dt_old=dt;
    
    % CAB
    Mass_data_old_capped=Mass_data_capped;

    melt_fraction_old =melt_fraction;
    solid_fraction_old=solid_fraction;

    indexS_old=indexS;
   
    if To_recored_intrusion_marker==1
        Intrusion_marker_dynamics_old=Intrusion_marker_dynamics;
    end

    mu_g_dynamics=ones(N+1,1)*mu_g_dynamics(1); %%TO_FIX! 
    T_S_region_old=T_S_region;
    if Courant<Courant0
        Courant=min(Courant*1.03,Courant0);
    end


    
    while iter<min_iter || (eps>Precision && iter<Max_iter)
        
%%      
        Calculate_capped_values2;

        phi_on_nodez=interp1(cellz, melt_fraction_capped, nodez,'linear','extrap');   %here use 1-solid_fraction_capped data instead of phi 
        
        % index_margin=find(phi_on_nodez(1:end-1) ~= 0 & phi_on_nodez(2:end) == 0);
        % index_margin=[index_margin; index_margin(end)+1; index_margin(end)+2];
        % phi_on_nodez(index_margin)=0.5;
        
        melt_fraction_on_nodez=interp1(cellz,melt_fraction_capped,nodez,'linear','extrap');  %Here use capped data
        melt_fraction_on_nodez(melt_fraction_on_nodez<0)=0;

        index=find(melt_fraction_capped==0); 
        phi_on_nodez(index)=0;
        phi_on_nodez(index+1)=0;
        %
        S_on_nodez=interp1(cellz,S,nodez,'linear','extrap');
        % 
        S_on_nodez(S==0)=0;  %% how should I explain these two lines? fuck it, too many words.............now I forget it myself
        S_on_nodez(find(S==0)+1)=0;          
        

        %%
            sum_mass_FM=interp1(cellz,Mass_all,nodez,'linear','extrap');
            dP=zeros(N+1,1);
            % index=find(phi_on_nodez<=max(min_por,1e-2)); %this part only considers hydrostatic pressure
            % index2=setxor(1:N+1,index);  %this part considers dynamic pressure
            
            index=1:N;
            index2=[];
            

            dP(index)=-sum_mass_FM(index)*g;
            % dP(index2)=+C_values(index2,1).*(u_all(N+1+index2)-u_all(index2))./phi_on_nodez(index2)...
            %            -C_values(index2,3).*(u_all(2*N+2+index2)-u_all(N+1+index2))./phi_on_nodez(index2)...
            %     -((1-S_on_nodez(index2)).*rho_on_nodez(index2,1)+S_on_nodez(index2).*rho_on_nodez(index2,3))*g;
            
            % dP=max(dP,-rho_on_nodez(:,2)*g);
            % dP=min(dP,-rho_on_nodez(:,3)*g);
            dP=dP*Dynamic_scaling;
            P=zeros(N,1);
            P(N)=-dP(N)*dz(N)/2;
            for i=N-1:-1:1
                P(i)=P(i+1)-dP(i)*(dz(i)+dz(i+1))/2;
            end
            Pg_real=P/1e8+P_top/1000; %P in Pa and P_real in kbar
            %%
            
        %new melt viscosity model
        if melt_viscosity_type==0
            %old melt viscosity model
            muf_order=(muf2*Mass_data(:,1)+muf1*Mass_data(:,2))./max(sum(Mass_data(:,1:2),2),1e-3);
        else
            %new melt viscosity model
            SiO2=(Mass_data(:,1)*74+Mass_data(:,2)*47)./max(sum(Mass_data(:,1:2),2),1e-3); % SiO2%
            if Fix_H2O==1 || Fix_H2O==3
                H2O=H2O_crust*ones(N,1);
            else
                H2O=100*(Mass_data(:,3))./max(sum(Mass_data(:,1:3),2),1e-3); % H2O%  !!!
                for i=1:N
                    P3=Pg_real(i)*100;
                    Saturation=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5);
                    H2O(i)=H2O(i)/Saturation; %%should I do the scaling here????
                end
            end
            vis_B=vis_b1*SiO2+vis_b2*H2O+vis_b3*log(1+H2O);
            vis_C=vis_c1*SiO2+vis_c2*H2O+vis_c3*log(1+H2O);
            muf_order=-4.55+vis_B./(max(T,500)+273.15-vis_C);            
        end
        mu_f_dynamics=10.^muf_order;  
        mu_f_dynamics=interp1(cellz,mu_f_dynamics,nodez,'linear','extrap')/Dynamic_scaling;
        mu_m_dynamics=mu_all(min(ceil(melt_fraction_on_nodez(1:N+1)*N_vis)+1,N_vis))/Dynamic_scaling;  %viscosity dependant on melt fraction! Not porosity!


        C_values=zeros(N,3);
        % C_values=zeros(N,1);
        for i=1:N+1
            C_values(i,1)=c_gen2(mu_f_dynamics(i),A_type2,phi_on_nodez(i), 0, 1);
        end
        % C_values(phi_on_nodez<0.2)=C_values(phi_on_nodez<0.2)/100;
        
        iter_inner=1;
        rho_on_nodez=zeros(N+1,3);
        if Conservation_type==1
            rho_on_nodez(:,1)=interp1(cellz,rho(:,1),nodez,'linear','extrap');
            rho_on_nodez(:,2)=interp1(cellz,rho(:,2),nodez,'linear','extrap');
            rho_on_nodez(:,3)=interp1(cellz,rho(:,3),nodez,'linear','extrap');
        else % for volume conservation, recalulate the dyanmic denties used for gravitational term
            if Melt_density_type==2
                SiO2=(Mass_data(:,1)*PD_range(2)+Mass_data(:,2)*PD_range(1))./max(sum(Mass_data(:,1:2),2),1e-5)/100;
                H2O=Mass_data(:,3)./max(sum(Mass_data(:,1:3),2),1e-5);
                VMX=H2O.*SiO2*Coef_melt_den(1)+SiO2*Coef_melt_den(2)+H2O*Coef_melt_den(3)+Coef_melt_den(4);
                x_other=1-SiO2-H2O;
                V_all=SiO2.*26.86e-6/0.06009+H2O.*26.27e-6/0.01802+x_other.*VMX;
                rho_temp=1./V_all;
                rho_on_nodez(2:end-1,1)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                rho_on_nodez(1,1)=rho_temp(1);  rho_on_nodez(end,1)=rho_temp(end);
            else
                rho_temp=(Mass_data(:,1)*rhof_1+Mass_data(:,2)*rhof_2)./max(sum(Mass_data(:,1:2),2),1e-5);
                rho_on_nodez(2:end-1,1)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
                rho_on_nodez(1,1)=rho_temp(1);  rho_on_nodez(end,1)=rho_temp(end);
            end

            rho_temp=(Mass_data(:,4)*rhom_1+Mass_data(:,5)*rhom_2)./max(sum(Mass_data(:,4:5),2),1e-5);
            rho_on_nodez(2:end-1,2)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
            rho_on_nodez(1,2)=rho_temp(1);  rho_on_nodez(end,2)=rho_temp(end);

            rho_on_nodez(rho_on_nodez(:,2)==0,2)=rhom_2; % patch the phi=1 case 
            
            P2=-rho_mean*g*(cellz-nodez(end))/1e5; %in bar
            rho_temp=(aa(1)* max(T,300).^bb(1)+aa(2)*max(P2,400).^bb(2)+aa(3)*max(T,300).^bb(3).*max(P2,400).^cc(1))*1000;
            rho_on_nodez(2:end-1,3)=(rho_temp(1:end-1)+rho_temp(2:end))/2;
            rho_on_nodez(1,3)=rho_temp(1);  rho_on_nodez(end,3)=rho_temp(end);
            if Velocity_solver_type==2
                rho_on_nodez(S_on_nodez>0,1)=rho_on_nodez(S_on_nodez>0,3);
            end
        end
        %%
        
        
        if Velocity_solver_type==2
            % Two_phase solver
            u_all=Velocity_solve_portion_with_source(phi_on_nodez,u_all_old, Mass_data, zeros(1,N+1), C_values, mu_m_dynamics, mu_m_dynamics, rho_on_nodez(:,1), rho_on_nodez(:,2), cellz, nodez, dz, g, zeros(N+1,1));
        else
            % Three phase solver
            % pass the Mass_data_capped to the solver
            u_all=Threephase_velocity_solver_simple(C_values,mu_f_dynamics,mu_m_dynamics,mu_g_dynamics,rho_on_nodez, Precision,dz,g, Mass_data_capped, cellz, nodez, dt,0,  Conservation_type, indexS);
        end

        u_all(u_all==inf)=0;
        u_all(u_all==-inf)=0;
        u_all(isnan(u_all))=0;
        
%         u_all(:)=0;
            
        %rho_ref=(aa(1)* max(T(1),300).^bb(1)+aa(2)*max(Pg_real(1)*1000,700).^bb(2)+aa(3)*max(T(1),300).^bb(3).*max(Pg_real(1)*1000,700).^cc(1))*1000; %% CAB DENSITY OF VOLATILES - Haiyang's original.
        rho_ref=(aa(1)* max(T,300).^bb(1)+aa(2)*max(Pg_real*1000,700).^bb(2)+aa(3)*max(T,300).^bb(3).*max(Pg_real*1000,700).^cc(1))*1000; %% CAB DENSITY OF VOLATILES - use this one
        %% update step times
        if ~fixed_dt %&& iter==0
            if max(abs(u_all))~=0
                ind=u_all~=0;
                index=find(phi>min_phi_dt);
                index=intersect(index,find(phi<1-min_phi_dt));
                dt=min(abs([Courant*dz(index)./u_all(index); Courant*dz(index)./u_all(index+N+1)]));
            else
                dt=3*dt;
            end
            if isempty(dt)
                dt=50*Year;
            end
            if dt>1.7*dt_old
                dt=1.7*dt_old;
            end
            
            %                 if Time<35*Year
            %                     dt=min(dt,0.2*Year);
            %                 end
            
            dt=min(dt,Max_dt);
            dt=max(dt,Min_dt);
        end
%         dt=100;
        %% Component and Energy transport
        % Mass_data: [m1, n2, v1, m2, n2,v2, v3]
        if is_solid_solution==0
            T_above=interp1(cellz,T-Ts_new,nodez,'linear','extrap');
        else
            T_above=interp1(cellz,T-Ts_local,nodez,'linear','extrap');
        end
        kc=ones(N+1,1)*kc0.*(T_above>=-20);
        kc2=ones(N+1,1)*kc20.*(T_above>=-20);  % for water transport
     

        p1=ones(N+1,1); %phi_on_nodez.*(1-S_on_nodez);
        p2=ones(N+1,1); %(1-phi_on_nodez); 
        p3=ones(N+1,1); %phi_on_nodez.*S_on_nodez;      

        M1l=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,1), ones(N,1), ones(N,1), u_all(1:N+1), kc.*p1,dt, 0,[3,3],0);
        M1s=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,4), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc.*p2,dt, 0,[3,3],0);

        Nl=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,2), ones(N,1), ones(N,1), u_all(1:N+1), kc.*p1,dt, 0,[3,3],0);
        Ns=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,5), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc.*p2,dt, 0,[3,3],0);

        Vl=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,3), ones(N,1), ones(N,1), u_all(1:N+1), kc2.*p1,dt, 0,[3,3],0);
        Vs=CV_transport_scaled_v2(dz,min(Mass_data_old_capped(:,6), Mass_solid.*S_cap), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc2.*p2,dt, 0,[3,3],0);  %only those less than S_cap is transported as solid

        MM=M1s+M1l;
        % NN=rho_mean-MM-Vl-Vs;
        NN=Ns+Nl;


        if N_component==5 || Add_CLCU==1
            CLs=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,9), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc,dt, 0,[3,3],0);
            CLl=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,8), ones(N,1), ones(N,1), u_all(1:N+1)    , kc*10,dt, 0,[3,3],0);
            % CLg=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,10),ones(N,1), ones(N,1), u_all(2*N+3:3*N+3), kc,dt, 0,[1,2],[rho_ref*phi(1)*S(1)*CL_vol 0]);
            % CLg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,10), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0 0]);
            % CL=CLs+CLl+CLg;

            % CUs=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,12), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc,dt, 0,[3,3],0);
            % CUl=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,11), ones(N,1), ones(N,1), u_all(1:N+1)    , kc*10,dt, 0,[3,3],0);
            % no need to transport anymore for the new Cu and Sul calculation...just for a place holder
            CUs=Mass_data_old_capped(:,12);
            CUl=Mass_data_old_capped(:,11);
            % CUg=CV_transport_scaled_v2(dz,Mass_data_old_capped(:,13),ones(N,1), ones(N,1), u_all(2*N+3:3*N+3), kc,dt, 0,[1,2],[rho_ref*phi(1)*S(1)*CU_vol 0]);
            % CUg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,13), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0 0]);
            % CU=CUs+CUl+CUg;      
        end

        if N_component==3
            Temp=(Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0))/rho_mean;
        else 
            Vol_Mass=sum(Mass_data_old_capped(:,[7 10] ),2);
            Temp=(Vol_Mass+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0))/rho_mean;
        end

        low_kv=4;
        high_kv=6;
        low_kv_phi=0.01;
        high_kv_phi=0.02;
        Temp=max(Temp,low_kv_phi);
        Temp=min(Temp,high_kv_phi);
        Temp=(high_kv-low_kv)/(high_kv_phi-low_kv_phi)*(Temp-low_kv_phi)+low_kv;   %order -2 to 6 for 0.015-0.04
        kc3=kc20.*10.^Temp;
        % kc3(T>Ts_local-40)=0;
        if N_component==3 && Add_CLCU~=1
            if Vol_flux>0
                Vg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[1,3],[Vol_flux*rho_mean, 0]);
            else
                Vg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0, 0]);
            end
        else
            Vg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0, 0]);
            CLg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,10), ones(N,1), ones(N,1), zeros(N+1,1), kc3 ,dt, 0,[3,3],[0 0]);

            % CUg=CV_transport_scaled_one_sided(dz, Mass_data_old_capped(:,13), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0 0]);
            CUg=Mass_data_old_capped(:,13);
            % CL=CLs+CLl+Mass_data_old_capped(:,10);
        end
        % V=Vs+Vl+Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0);  
        V=Vs+Vl+Vg;
        if Add_CLCU==1
            CL=CLs+CLl+CLg;
            CU=CUs+CUl+CUg;
        else
            CL=zeros(N,1);
            CU=zeros(N,1);
        end

        if Conservation_type==2
            if Transport_method==1
                sum_mass=MM+NN+V;
                MM=MM*rho_mean./sum_mass;
                NN=NN*rho_mean./sum_mass;
                V=V*rho_mean./sum_mass;
            else
                % for j=1:2
                %     index=find(MM./(rho_mean-V)>0.98);
                %     extra=(MM(index)-(rho_mean-V(index))*0.98);
                %     for i=1:length(extra)
                %         MM(index(i)-1)=MM(index(i)-1)+extra(i)*dz(index(i))/dz(index(i)-1);
                %         % MM(index(i)+1)=MM(index(i)+1)+extra(i)/2*dz(index(i))/dz(index(i)+1);
                %         MM(index(i))=MM(index(i))-extra(i);
                %     end
                % end
                if N_component==3
                    % NN=rho_mean-MM-V;
                else
                    % NN=rho_mean-MM-V-CL;
                    CL=CLs+CLl+CLg;
                    CU=CUs+CUl+CUg;
                end

                % NN(NN<0)=0;
          
 
                % V=Vs+Vl+Vg;
                
                % extra=Vg-(Mass_data_old_capped(:,7)+max(Mass_data_old_capped(:,6)-Mass_solid.*S_cap,0));
                % 
                % index=find(extra>0);
                % if ~isempty(index)
                %     ratio=MM(index)./(MM(index)+NN(index));
                %     MM(index)=MM(index)-extra(index).*ratio;
                %     NN(index)=NN(index)-extra(index).*(1-ratio);
                %     SumM=sum(extra(index).*ratio.*dz(index));
                %     SumN=sum(extra(index).*(1-ratio).*dz(index));
                % 
                %     ratio=SumM/(SumM+SumN);
                %     index=find(extra<0);
                %     MM(index)=MM(index)-extra(index).*ratio;
                %     NN(index)=NN(index)-extra(index).*(1-ratio);
                % end
            end
        end
        H_source= zeros(N,1);
        
        kt=(kt0(1)*sum(Mass_data(:,[1 4]),2)+kt0(2)*sum(Mass_data(:,[2 5]),2))./max(sum(Mass_data(:,[1 2 5 5]),2),1e-3);
        if simplified_5p==1
            H=Threephase_enthalpy_solver3(Mass_data, Mass_data_old,T_old, H_old, u_all,cp, kt ,Lf, dz, cellz, nodez, dt, H_source, rho_ref*phi(1)*S(1), 0, 3);
        else
            H=Threephase_enthalpy_solver3(Mass_data, Mass_data_old,T_old, H_old, u_all,cp, kt ,Lf, dz, cellz, nodez, dt, H_source, rho_ref*phi(1)*S(1), 0, N_component);
        end
        

        if To_recored_intrusion_marker==1
            Intrusion_marker_dynamics=CV_transport_scaled_v2(dz,Intrusion_marker_dynamics_old, ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc*0,dt, 0,[3,3],0);
        end
        %% Solve chemical model
        phi_pre=phi;
        S_pre=S;

        i=1;
        iter_chemical=zeros(N,1);

        Tl_local=zeros(N,1);
        Saturation=zeros(N,1);
        Mass_data_pre=Mass_data;        
        Output_guess=[0,0,0, MM(1), NN(1),V(1),0, 500, S_cap(1)]; 
        K=K0*ones(N,1);
        Cases=zeros(N,1);
        if is_solid_solution==0
            [phi(1),S(1),T(1),rho(1,:),~,T_S_region(1,:),Ts_new(1),Mass_data(1,:),S_cap(1)]=Phase_component_updated2(H(1),MM(1),NN(1),V(1),Pg_real(1), Output_guess, [1,1], Jacobians, Rhs, Constant_index, Sys_constant,cp(1,:),  rho_constant, Lf, K0, Ts,  Tl, min_por, dz(1),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type);
            for i=2:N

                [phi(i),S(i),T(i),rho(i,:),~,T_S_region(i,:),Ts_new(i),Mass_data(i,:),S_cap(i)]=Phase_component_updated2(H(i),MM(i),NN(i),V(i),Pg_real(i-1),[Mass_data(i-1,:) T(i-1) S_cap(i-1)], T_S_region(i-1,:), Jacobians, Rhs, Constant_index, Sys_constant,cp(i,:),   rho_constant, Lf, K0, Ts, Tl, min_por, dz(i),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type);
            end
        else
            if N_component==3
                if Add_CLCU==1
                    input=[H(1),MM(1),NN(1),V(1),CL(1),CU(1)];
                else
                    input=[H(1),MM(1),NN(1),V(1)];
                end
                [phi(1),S(1),T(1),rho(1,:),~,T_S_region(1,:),Ts_new(1),Ts_local(1),Mass_data(1,:),S_cap(1),Tl_local(1),~,Saturation(1), Partition_CL_CU(1,:)]=Phase_component_updated_solid3(input,Pg_real(1), 500, [1,1], Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(1,:),  rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts,  Tl, min_por, max_melt, dz(1),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS,0);

                % [phi(1),S(1),T(1),rho(1,:),~,T_S_region(1,:),Ts_new(1),Ts_local(1),Mass_data(1,:),S_cap(1),Tl_local(1),Cases(1)]=Phase_component_updated_solid3(H(1),MM(1),NN(1),V(1),Pg_real(1), 500, 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(1,:),   rho_constant, Lf, K0, Ts, Tl, min_por, max_melt, dz(1),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type,simplified_TS, 1);
                phi(1:5)=0;

                if isempty(gcp('nocreate'))
                    parpool(Num_core_use); % Reopen the pool desired number of cores
                end
                try
                parfor i=2:N
                    warning('off', 'all');

                    if Add_CLCU==1
                        input=[H(i),MM(i),NN(i),V(i),CL(i),CU(i)];
                    else
                        input=[H(i),MM(i),NN(i),V(i)];
                    end
                    [phi(i),S(i),T(i),rho(i,:),~,T_S_region(i,:),Ts_new(i),Ts_local(i),Mass_data(i,:),S_cap(i),Tl_local(i),~,Saturation(i), Partition_CL_CU(i,:)]=Phase_component_updated_solid3(input,Pg_real(i-1),T_old(i), 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(i,:),   rho_constant, Lf, [K(i) Kcl Kcu], Ts, Tl, min_por, max_melt, dz(i),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS,0);
                    % [phi(i),S(i),T(i),rho(i,:),~,T_S_region(i,:),Ts_new(i),Ts_local(i),Mass_data(i,:),S_cap(i),Tl_local(i),Cases(i)]=Phase_component_updated_solid3(H(i),MM(i),NN(i),V(i),Pg_real(i-1),T_old(i), 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(i,:),   rho_constant, Lf, K0, Ts, Tl, min_por, max_melt, dz(i),Data_point,Data_y,Data_point2,Data_y2,Precision, Conservation_type,simplified_TS, force_solid);
                end
                catch
                    eval(['save(''Record_error'  '.mat''' ',' '''-regexp''' ',' '''^(?!Video_handel$).'');'])
                end
            else
                [phi(1),S(1),T(1),rho(1,:),~,T_S_region(1),Ts_local(1),Tl_local(1),Mass_data(1,:),S_cap(1),Partition_CL_CU(1,:)]=Phase_component_updated_solid_5p_simple(H(1),MM(1),NN(1),V(1),CL(1),CU(1), Pg_real(1), 1000, 1, Jacobians, Rhs, Constant_index,Extra_func, Sys_constant, cp(1,:),   rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts, Tl, min_por, dz, Data_point,Data_y,Data_point2,Data_y2, PD_range, Precision, Conservation_type, 1);
                phi(1:5)=0;
                if N<50 || Use_parrallel==0
                    for i=2:N
                        warning('off', 'all');
                        [phi(i),S(i),T(i),rho(i,:),~,T_S_region(i,:),Ts_local(i),Tl_local(i),Mass_data(i,:),S_cap(i),Partition_CL_CU(i,:)]=Phase_component_updated_solid_5p_simple(H(i),MM(i),NN(i),V(i),CL(i),CU(i),Pg_real(i-1), 1000, 1, Jacobians, Rhs, Constant_index, Extra_func, Sys_constant,cp(i,:),   rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts, Tl, min_por,  dz(i),Data_point,Data_y,Data_point2,Data_y2,PD_range, Precision,Conservation_type,1);
                    end
                else
                    if isempty(gcp('nocreate'))
                        parpool(Num_core_use); % Reopen the pool desired number of cores
                    end
                    parfor i=2:N
                        warning('off', 'all');
                        [phi(i),S(i),T(i),rho(i,:),~,T_S_region(i,:),Ts_local(i),Tl_local(i),Mass_data(i,:),S_cap(i),Partition_CL_CU(i,:)]=Phase_component_updated_solid_5p_simple(H(i),MM(i),NN(i),V(i),CL(i),CU(i),Pg_real(i-1), 1000, 1, Jacobians, Rhs, Constant_index, Extra_func, Sys_constant,cp(i,:),   rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts, Tl, min_por,  dz(i),Data_point,Data_y,Data_point2,Data_y2,PD_range, Precision,Conservation_type,1);
                    end
                end
            end
        end
        eps=max(abs(phi-phi_pre));
        % eps=max(max(abs(Mass_data-Mass_data_pre))/1e7);

                % Update_plot;
                % drawnow;
        iter=iter+1;
        % if any(Mass_data<0,'all')
        %     eval(['save(''Record_erro .mat''' ',' '''-regexp''' ',' '''^(?!Video_handel$).'');']);
        %     error('sorry, it fucked up..')
        % end

        if Enhanced_convergence==1 
            if iter==Max_iter && Courant>1e-3 && dt>1e-2*Year    
                Num_non_convergence=Num_non_convergence+1;
                Courant=Courant/2;
                iter=1;
                disp(['New Courant: ' num2str(Courant)] )

                phi=phi_old;
                S=S_old;

                rho=rho_old;
                P=P_old;
                H=H_old;
                Mass_data=Mass_data_old;
                Ts_local=Ts_local_old;
                Tl_local=Tl_local_old;

                T=T_old;

                dt=dt_old;
                
                %CAB
                Mass_data_capped=Mass_data_old_capped;

                melt_fraction =melt_fraction_old;
                solid_fraction=solid_fraction_old;                
            end
        end
    end
    %post-processings for the new Sulfur and Copper calculations
    % transport S and Cu
    CUs=CV_transport_scaled_v2(dz,Mass_data_addition(:,5), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc,dt, 0,[3,3],0);
    CUl=CV_transport_scaled_v2(dz,Mass_data_addition(:,4), ones(N,1), ones(N,1), u_all(1:N+1)    , kc*10,dt, 0,[3,3],0);
    CUg=CV_transport_scaled_one_sided(dz, Mass_data_addition(:,6), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0 0]); %may need to update this....
    CU=CUs+CUl+CUg;
    
    Suls=CV_transport_scaled_v2(dz,Mass_data_addition(:,2), ones(N,1), ones(N,1), u_all(N+2:2*N+2), kc,dt, 0,[3,3],0);
    Sull=CV_transport_scaled_v2(dz,Mass_data_addition(:,1), ones(N,1), ones(N,1), u_all(1:N+1)    , kc*10,dt, 0,[3,3],0);
    Sulg=CV_transport_scaled_one_sided(dz, Mass_data_addition(:,3), ones(N,1), ones(N,1), zeros(N+1,1), kc3,dt, 0,[3,3],[0 0]); %may need to update this....
    SUL=Suls+Sull+Sulg;
    
    S_CU_to_process=1:N;
    S_Cu_postprocessing;
    
    if Enhanced_convergence==1
        if iter<Max_iter/5
            Courant=min(Courant*1.03,Courant0);
        end
    end

    if To_extract_volatile==1

        Last_phi=find(phi>0,1,'last');
        Last_phi=find(cellz>cellz(Last_phi)+100,1,'first');
        
        temp=Mass_data(:,7)+max(Mass_data(:,6)-Mass_solid.*S_cap,0);
        vol_cells=find(temp>0)';
        vol_cells=vol_cells(vol_cells>=Last_phi);
        % if ~isempty(vol_cells)
        %     if vol_cells(end)==Last_phi
        %         d = diff(vol_cells);
        %         splitAt = find(d > 1);
        %         seriesStart = [1, splitAt + 1];
        %         seriesEnd = [splitAt, numel(vol_cells)];
        %         vol_cells = vol_cells(seriesStart(end):seriesEnd(end));
        %     else
        %         vol_cells=[];
        %     end
        % end

        % index=find(cellz-nodez(end)+Depth_to_extract>0);
        % vol_cells=intersect(vol_cells,index);  %only for above Depth_to_extract
        if ~isempty(vol_cells)
            release=temp(vol_cells)/2;
            % Invis_vol(vol_cells)=Invis_vol(vol_cells)+release;
            Leak_V=Leak_V+sum(release.*dz(vol_cells));
            Sum_V0=Sum_V0-sum(release.*dz(vol_cells)); % update the conservative V

            if N_component==5 || Add_CLCU==1
                temp=temp./max(V,1e-3)/2;
                Leak_CU=Leak_CU+sum(CL(vol_cells).*temp(vol_cells).*dz(vol_cells));
                Leak_CL=Leak_CL+sum(CU(vol_cells).*temp(vol_cells).*dz(vol_cells));
                Leak_SUL=Leak_SUL+sum(SUL(vol_cells).*temp(vol_cells).*dz(vol_cells));
                CU(vol_cells)=CU(vol_cells)-CU(vol_cells).*temp(vol_cells);
                CL(vol_cells)=CL(vol_cells)-CL(vol_cells).*temp(vol_cells);
                SUL(vol_cells)=SUL(vol_cells)-SUL(vol_cells).*temp(vol_cells);
            end

            Leak_H=Leak_H+sum(release.*dz(vol_cells).*T(vol_cells).*cp(vol_cells,3));            

            V(vol_cells)=V(vol_cells)-release;
            % scale=MM(vol_cells)./(MM(vol_cells)+NN(vol_cells));
            % MM(vol_cells)=MM(vol_cells)+scale.*release;
            % NN(vol_cells)=NN(vol_cells)+(1-scale).*release;
            for i=1:length(vol_cells)
                % try
                    if Add_CLCU==1
                        input=[H(vol_cells(i)),MM(vol_cells(i)),NN(vol_cells(i)),V(vol_cells(i)),CL(vol_cells(i)),CU(vol_cells(i))];
                    else
                        input=[H(vol_cells(i)),MM(vol_cells(i)),NN(vol_cells(i)),V(vol_cells(i))];
                    end
                    [phi(vol_cells(i)),S(vol_cells(i)),T(vol_cells(i)),rho(vol_cells(i),:),~,T_S_region(vol_cells(i),:),Ts_new(vol_cells(i)),Ts_local(vol_cells(i)),Mass_data(vol_cells(i),:),S_cap(vol_cells(i)),Tl_local(vol_cells(i)),~,~,Partition_CL_CU(vol_cells(i),:)]=Phase_component_updated_solid3(input,Pg_real(vol_cells(i)-1),T_old(vol_cells(i)), 1, Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(vol_cells(i),:),   rho_constant, Lf, [K(1) Kcl Kcu 3 50], Ts, Tl, min_por, max_melt, dz(vol_cells(i)),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS, 0);
                % catch
                %     eval(['save(''Record_error2'  '.mat''' ',' '''-regexp''' ',' '''^(?!Video_handel$).'');'])
                %     % error('error 2')
                % end
            end
            S_CU_to_process=vol_cells;
            S_Cu_postprocessing;
        end
        
        Record_vol_time=Record_vol_time+dt;
        if Record_vol_time>Record_gap
            Record_vol_time=0;
            if N_component==5 || Add_CLCU==1
                fprintf(fileID, '%d\t %.4f\t %.4f\t %.4f\t %.4f\n', round(Time/Year), Leak_V, Leak_CU, Leak_CL, Leak_H);
                Leak_CU=0;
                Leak_CL=0;
            else
                fprintf(fileID, '%d\t %.4f\t %.4f\n', round(Time/Year), Leak_V,  Leak_H);
            end
            Leak_V=0;
            Leak_H=0;
        end
    end
    


    Cb=(PD_range(2)*sum(Mass_data(:,[1,4]),2)+PD_range(1)*sum(Mass_data(:,[2,5]),2))./sum(Mass_data(:,[1 2 4 5]),2);

    % Sat=zeros(N,1);
    % PT=[-4e-4,-5e-4,-6e-4,-13e-4, -15.5e-4,-17e-4,-16e-4,-5e-4, 0, 26e-4,5e-3, 5e-3];
    % PTx=[0,   0.12    0.2  0.3    0.5       1       2       3   4  5,    11, 20];
    % P2=Pg_real*100;
    % for k=1:N
    %     ind=find(Pg_real(k)<=PTx,1,'first');
    %     dSdT=(PT(ind)*(Pg_real(k)-PTx(ind-1))+PT(ind-1)*(PTx(ind)-Pg_real(k)))/(PTx(ind)-PTx(ind-1));
    %     Sat(k)=((2.859e-2*P2(k)-1.495e-3*P2(k).^1.5+2.702e-5*P2(k).^2+0.257*P2(k).^0.5)+(T(k)-800)*dSdT)/100;
    % end
    % if N_component==3
    %     Cb2=V./(MM+NN+V);
    % else
    %     Cb2=V./(MM+NN+V+CL);
    % end

    %% CAB Evacuations
    if Conservation_type ==1
        rho_bulk=phi.*(1-S).*rho(:,1)+(1-phi).*rho(:,2)+phi.*S.*rho(:,3);
    else
        rho_on_cellz = interp1(nodez, rho_on_nodez, cellz, 'linear', 'extrap');
        rho_bulk = phi.*(1-S).*rho_on_cellz(:,1)+(1-phi).*rho_on_cellz(:,2)+phi.*S.*rho_on_cellz(:,3);        
    end

    


    if To_evacuate==1
        Finding_porosity;
    end
        


    if To_evacuate && ~isempty(buoy_phi_top)
        for j=1:length(buoy_phi_top)
            if crit_overpressure(j)<=overpressure_total(j)
                if buoy_Hphi_top(j)==-5000
                    % put in how this will mean that there is no evacuation
                    % because there is no high melt fraction present to
                    % evacuate
                     fprintf(File_echo, '%6s \n', 'No evacuation, as there is no high melt present in the layer which has reached critical overpressure.');
                else
                    evacuation_counter = evacuation_counter+1;
                    % [av_rho_T, av_Mass_data_T, av_H_T, av_T_T, av_Ts_local_T, av_phi_T, av_S_T, av_cb_T, av_OG_Cb_T,  av_S_cap,av_Tl_local_T,av_Partition_CL_CU]=Averages_evacuated_updated(evacuation_counter, buoy_Hphi_base(j),buoy_Hphi_top(j), phi, S, H, Mass_data, Pg_real, dz, MM, NN, V, CL, CU, Add_CLCU, Lf, cp, simplified_TS, Sys_constant, rho, Melt_density_type, PD_range, Coef_melt_den, cellz, nodez, aa, bb,cc,rho_mean,g, Velocity_solver_type, Ts, Tl, Conservation_type, rhom_1, rhom_2, Time, Year,OG_Cb); % calculates averages to evacuate and intrude. Also outputs averages that are intruded.
                    Averages_evacuated_updated;
                    Evacuations; % Evacuates nodes and then intrudes nodes. Also outputs the nodes which are removed.

                    % fprintf(File_echo, '%s %5.5f %s %5.5f %6s\n', 'Evacuated magma intruded at', ...
                    %         (nodez(depth_int_N_c)-nodez(end))/1000, 'km at', Time/Year/1000, 'ka');
                    % outputs details on the critical buoyancy and total
                    % buoyancy
                    % output: 1. hb, 2. hrti, 3. t2, 4. delta rho, 5. critical overpressue,
                    % 6. therm death overpessure, 7. crit rock ovepressure, 8. total
                    % ovepressure
                    % Filename_evac = ['Evac_details_', num2str(evacuation_counter), '.txt'];
                    % File_evac = fopen(Filename_evac, 'w');
                    % fprintf(File_evac, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');
                    % 
                    % fprintf(File_evac, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
                    %     'hb (km)', 'hrti (km)', 't2 (ka)', 'delta rho', 'Crit_over','Crit_therm', 'Crit_rock', ...
                    %     'OverTotal');
                    % 
                    % fprintf(File_evac, '%10.4f \t %10.4f \t %10.4f \t %10.4f \t %3e \t %3e \t %3e \t %3e \n ', ...
                    %     [buoyant_phi_m(j)./1000; hrti(j)./1000; t2(j)./Year./1000; ov_rho_bulk_buoy(j)-Av_bulk_density(j); ...
                    %     crit_overpressure(j); crit_thermdeath(j); crit_crust; overpressure_total(j)]);
                    % 
                    % fclose(File_evac);




                    %output the volume and composition evacuated
                    % 1. time 2. thickness of HMF 3. Volume 4. Composition
                    % fprintf(File_AVevac, '%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                    %     [Time/Year/1000; (melt_evacuated_top-nodez(end))/1000; (nodez(depth_int_N_c)-nodez(end))/1000;  buoyant_Hphi_m(j)/1000; (buoyant_Hphi_m(j)/1000)*((Diameter/1000)/2)^2*pi; av_cb_T;...
                    %     Thick_melt_evac; Thick_melt_int]);

                    
                    % reset
                    buoyant_phi_m(j) = -5000;
                    hrti(j) = scale_RTI*Diameter;
                    t2(j) = -5000;
                    
                 
                    %Adapt_mesh_3p;
                    Calculate_capped_values2;



                end
            end
        end
    end


    
% CAB end
    %% Outputs

    if Record_data == 1

        % CAB
        if (length(Record_time)>=Record_index)
            if (Time)>=Record_time(Record_index)+dt
                Record_index=Record_index+1;
                vv=['Output ', num2str(Time/Year/1000), ' ka'];
%                 disp(vv)

                % % % Overal model data % % %


                % CELLs

                FilenameO = [ 'output_', num2str(output_counter),'_CELLS.txt'];

                File_Out = fopen(FilenameO, 'w');

                fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

                fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
                    'Depth', 'phi', 'Vol. frac', 'Enth.', 'Temp (C)', 'Ts (C)', 'Tl (C)', 'cb', 'cl', 'cs', 'vb','vl', 'vs','vg', 'ov.density') ;

                for i=1:1:N
                    fprintf(File_Out, '%10.4f \t %10.4f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n ', ...
                        [((cellz(i)-nodez(end))/1000); phi(i); S(i); H(i); T(i); ...
                        Ts_local(i); Tl_local(i); Cb(i); (SiO2_range(2)*Mass_data(i,1) + SiO2_range(1)*Mass_data(i,2))/(sum(Mass_data(i,1:2),2)); (SiO2_range(2)*Mass_data(i,4) + SiO2_range(1)*Mass_data(i,5))/(sum(Mass_data(i,4:5),2)); ...
                        Cb2(i); Mass_data(i,1)/sum(Mass_data(i,[1,2]),2); Mass_data(i,4)/sum(Mass_data(i,[4,5]),2); (Mass_data(i,7)+max(Mass_data(i,6)-Mass_solid(i)*S_cap(i),0))/rho_mean; rho_bulk(i)]);


                end

                fclose(File_Out);

                %NODES

                FilenameON = ['output_',num2str(output_counter),'_NODES.txt'];
                
                File_OutN = fopen(FilenameON, 'w');

                fprintf(File_OutN, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

                fprintf(File_OutN, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
                    'Depth', 'Solid vel', 'Liq. vel', 'Vol. vel', 'Sol. dens', 'Liq. dens', 'Vol. dens') ;

                for i=1:1:N
                    fprintf(File_OutN, '%10.4f \t %4e \t %4e\t %4e \t %10.4f \t %10.4f \t %10.4f \n', ...
                        [((nodez(i)-nodez(end))/1000); u_all(N+1+i); u_all(i); 0; rho_on_nodez(i,2); rho_on_nodez(i,1); rho_on_nodez(i,3)]);
                end

                fclose(File_OutN);


                if mod(output_counter,restart_No)==0
                        save(['output_',num2str(output_counter),'.mat']');  %(FileNameO)   
                end

                % % % phi data from finding porosity % % %
                Filename_phi = ['phi_output_',num2str(output_counter),'.txt'];
                File_phi = fopen(Filename_phi, 'w');
                
                fprintf(File_phi, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n',...
                     'Time (ka)', 'No.', 'phi_top(-)', 'phi_base(-)', 'Hphi_t(-)', 'Hphi_b(-)', 'phi_t(km)', 'phi_b(km)', ...
                    'Hphi_t(km)', 'Hphi_b(km)', 'ov_rho_bulk(kgm-3/)');

                if ~isempty(phi_top)

                   for i=1:1:length(phi_top)

                        if Hphi_top(i) == - 5000
                            fprintf(File_phi, '%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \n ',...
                                [Time/Year/1000; i; phi_top(i); phi_base(i); Hphi_top(i); Hphi_base(i); (cellz(phi_top(i))-nodez(end))/1000;...
                                (cellz(phi_base(i))-nodez(end))/1000; 0; 0;...
                                ov_rho_bulk(i)]);
                        else
                            fprintf(File_phi, '%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \n ',...
                            [Time/Year/1000; i; phi_top(i); phi_base(i); Hphi_top(i); Hphi_base(i); (cellz(phi_top(i))-nodez(end))/1000;...
                            (cellz(phi_base(i))-nodez(end))/1000; (cellz(Hphi_top(i))-nodez(end))/1000; (cellz(Hphi_base(i))-nodez(end))/1000;...
                            ov_rho_bulk(i)]);
                        end

                    end
         

                else

                    fprintf(File_phi, '%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f\t %10.3f \t %10.3f \n ',...
                            [Time/Year/1000; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]);


                end

                fclose(File_phi);



                % % % buoy phi data from finding porosity % % %
                Filename_BUOYphi = ['BUOYphi_output_', num2str(output_counter), '.txt'];
                File_BUOYphi = fopen(Filename_BUOYphi,'w');

                

                fprintf(File_BUOYphi, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n ',...
                                'Time (ka)', 'No.', 'Prev.No', 'Phi top', 'Phi base', ...
                                'buoy top', 'buoy base','buoyCMFtop', 'buoyCMFbas', 'buoyPhi', ...
                                'buoy Hphi', 'PhiTop(km)', 'PhiBase(km)', 'buoyTop(km)', 'buoyBase(km)', ...
                                'buoyCMFT(km)', 'buoyCMFb(km)', 'buoyPhi(m)', 'buoyHphi(m)', 'CrustDens', ...
                                'layer dens', 'hrti_flag','t2(yr)', 'hrti_0(m)', 'hrti(m)', ...
                                'OverRTI', 'OverLayer', 'OverTotal', 'Cri_therm', 'Crit_over');

                if ~isempty(buoy_phi_top)

                    for i=1:1:length(buoy_phi_top)

                        if buoy_Hphi_top(i)==-5000
                        
                            fprintf(File_BUOYphi, ['%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f' ...
                                ' \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %3e \t %3e \t %3e \t %3e \t %3e  \n '], ...
                                [Time/Year/1000; i; prev_buoy_layer(i); all_phi_buoy_top(i); all_phi_buoy_base(i); ...
                                buoy_phi_top(i); buoy_phi_base(i); buoy_Hphi_top(i); buoy_Hphi_base(i); buoyant_phi(i);
                                buoyant_Hphi(i);(cellz(all_phi_buoy_top(i))-nodez(end))/1000; (cellz(all_phi_buoy_base(i))-nodez(end))/1000; (cellz(buoy_phi_top(i))-nodez(end))/1000; (cellz(buoy_phi_base(i))-nodez(end))/1000; ...
                                0;0; buoyant_phi_m(i); buoyant_Hphi_m(i); ov_rho_bulk_buoy(i); ...
                                Av_bulk_density(i); hrti_flag(i); t2(i)/Year;hrti_start(i); hrti(i); ...
                                overpressure_RTI(i); overpressure_layer(i); overpressure_total(i); crit_thermdeath(i); crit_overpressure(i)]);

                        else


                            fprintf(File_BUOYphi, ['%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f ' ...
                                '\t %10.3f \t %10.3f \t %10.3f \t %3e \t %3e \t %3e \t %3e  \n '], ...
                                [Time/Year/1000; i; prev_buoy_layer(i); all_phi_buoy_top(i); all_phi_buoy_base(i); buoy_phi_top(i); buoy_phi_base(i); ...
                                buoy_Hphi_top(i); buoy_Hphi_base(i); buoyant_phi(i);
                                buoyant_Hphi(i); (cellz(all_phi_buoy_top(i))-nodez(end))/1000;(cellz(all_phi_buoy_base(i))-nodez(end))/1000; ...
                                (cellz(buoy_phi_top(i))-nodez(end))/1000; (cellz(buoy_phi_base(i))-nodez(end))/1000; (cellz(buoy_Hphi_top(i))-nodez(end))/1000; ...
                                (cellz(buoy_Hphi_base(i))-nodez(end))/1000; buoyant_phi_m(i); buoyant_Hphi_m(i); ov_rho_bulk_buoy(i); Av_bulk_density(i); hrti_flag(i); t2(i)/Year;...
                                hrti_start(i); hrti(i); overpressure_RTI(i); overpressure_layer(i); overpressure_total(i); crit_thermdeath(i); crit_overpressure(i)]);





                        end


                    end

                else

                    frprintf(File_BUOYphi, ['%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t ...' ...
                            '%10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t...' ...
                            ' %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t %10.3f \t ...' ...
                            '%10.3f \t %10.3f \t %10.3f \n '], ...
                            [Time/Year/1000; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0])


                end

                fclose(File_BUOYphi);



                fprintf(File_breaks, '%10.5f \t %10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, break_percent, break_counter, counter);

                fprintf(File_comp, '%10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, Mass_cons, Mass_res);

                output_counter=output_counter+1;

            end
        end
    end



        
    counter=counter+1;

    % CAB - output when the residual is greater than precision
    if eps>Precision
        fprintf(File_echo, '%6s %5.5f %6s \n', 'Residual is greater than precision at ', Time/Year/1000, 'ka');
        fprintf(File_precision, '%10.5f \t %10.9f \t %10.9f \n', Time/Year/1000, Precision, eps);
        break_counter=break_counter+1;

    end
    break_percent = break_counter/counter*100;


    %calculate mass conservation
    composition = sum(Cb);
    OG_composition = sum(OG_Cb);
    Mass_cons=composition/OG_composition;
    Mass_res=abs(composition-OG_composition);
% CAB end


    % Adaptive mesh
    

    Refreshed=0;
    
%     Cal_water_contribution2;
    if running_plot==1
        if floor(toc/Show_refresh)>Frame_num
            Update_plot;
            drawnow
            Frame_num=Frame_num+1;
            Refreshed=1;
        end
    end
    % if To_extract_volatile==1
    %     V_loss=(sum(sum(Mass_data_old(:,[3,6,7])))-sum(sum(Mass_data(:,[3,6,7]))))/sum(sum(Mass_data_old(:,[3,6,7])));
    % end
    Sum_H=sum(H.*dz);
    % disp([num2str(counter), ' ', num2str(iter), ' ', num2str(max(sum(Mass_data(1:7,:),2))), ' ' num2str(Num_non_convergence)])  %, ' ', num2str((Sum_H-Sum_H0)/Sum_H0
    disp([num2str(counter), ' ', num2str(iter), ' ', num2str(max(sum(Mass_data(1:7,:),2))), ' ', num2str(min(sum(Mass_data(1:7,:),2))), ' ', num2str(Time/Year)  ' ' ])
    
    Conservation_check=[sum(sum(Mass_data(:,[1 4]),2).*dz)/Sum_M0, sum(sum(Mass_data(:,[2 5]),2).*dz)/Sum_N0, sum(sum(Mass_data(:,[3 6 7]),2).*dz)/Sum_V0]; 
    %Display the conservation
    %     disp([num2str(Conservation_check(1)), ' ', num2str(Conservation_check(2)), ' ', num2str(Conservation_check(2))])
   

    Time=Time+dt;
    
    temp=mod(counter,Convergence_record_steps);
    if temp==0
        Convergence_record(Convergence_record_steps,:)=[counter iter];
        writematrix(Convergence_record,'Iteration_recored.txt','WriteMode','append')
    else
        Convergence_record(temp,:)=[counter iter];
    end
    
    if exist('Create_save_data','var')
        if Create_save_data==1
            if Time>Save_data_times(Save_data_index)
                eval(['save(''Record_' num2str(Save_data_index) '.mat''' ',' '''-regexp''' ',' '''^(?!Video_handel$).'');'])
                Save_data_index=Save_data_index+1;
            end
        end
    end
    
    if Create_video==1 && Time>(video_data_index+1)*Fixed_record_dt0
        video_data_index=video_data_index+1;
        video_data_file_name=['Data' num2str(video_data_index)];
        save(video_data_file_name,'Time','cellz','Mass_data','Mass_data_addition', 'Pg_real','Ts_local','Tl_local','T', 'S_cap','phi')
    end
    if isempty(find(T>Ts_local-1,1)) && (Time>injection_time(end-1)|| To_intrude==0 )
        eval(['save(''Record_' num2str(Save_data_index) '.mat''' ',' '''-regexp''' ',' '''^(?!Video_handel$).'');'])
        break
    end
end

%% final txt overall output
% CAB 
if Record_data==1
    
    % CELLs

    FilenameO = [ 'output_', num2str(Time/Year/1000),'_CELLS.txt'];

    File_Out = fopen(FilenameO, 'w');

    fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

    fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'Depth', 'phi', 'Vol. frac', 'Enth.', 'Temp (C)', 'Ts (C)', 'Tl (C)', 'cb', 'cl', 'cs', 'vb','vl', 'vs','vg', 'ov.density') ;

    for i=1:1:N
        fprintf(File_Out, '%10.4f \t %10.4f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n ', ...
            [((cellz(i)-nodez(end))/1000); phi(i); S(i); H(i); T(i); ...
            Ts_local(i); Tl_local(i); Cb(i); (SiO2_range(2)*Mass_data(i,1) + SiO2_range(1)*Mass_data(i,2))/(sum(Mass_data(i,1:2),2)); (SiO2_range(2)*Mass_data(i,4) + SiO2_range(1)*Mass_data(i,5))/(sum(Mass_data(i,4:5),2)); ...
            Cb2(i); Mass_data(i,1)/sum(Mass_data(i,[1,2]),2); Mass_data(i,4)/sum(Mass_data(i,[4,5]),2); (Mass_data(i,7)+max(Mass_data(i,6)-Mass_solid(i)*S_cap(i),0))/rho_mean; rho_bulk(i)]);


    end

    fclose(File_Out);

    %NODES

    FilenameON = ['output_',num2str(Time/Year/1000),'_NODES.txt'];
    
    File_OutN = fopen(FilenameON, 'w');

    fprintf(File_OutN, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

    fprintf(File_OutN, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'Depth', 'Solid vel', 'Liq. vel', 'Vol. vel', 'Sol. dens', 'Liq. dens', 'Vol. dens') ;

    for i=1:1:N
        fprintf(File_OutN, '%10.4f \t %4e \t %4e\t %4e \t %10.4f \t %10.4f \t %10.4f \n', ...
            [((nodez(i)-nodez(end))/1000); u_all(N+1+i); u_all(i); 0; rho_on_nodez(i,2); rho_on_nodez(i,1); rho_on_nodez(i,3)]);
    end

    fclose(File_OutN);
end

fclose(File_AVevac);
fclose(File_echo);
fclose(File_precision);
fclose(File_breaks);
fclose(File_comp);

% CAB end


toc;
