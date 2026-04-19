%% 
% Haiyang Hu 2/2/2020
% A 1D control volumn two phase solver for magma simulation.
% update v2: Add sill intrusion
%            Add mesh adaptivity
% update new: use the new type of phase diagram which is defined as piecewise second order polynomial
% update v3: Change the definition of reactive flow and how it is
% update v4: Using new way to calculate the three contributions for melt fraction change
% update v5: Remove the old multi intrusion scheme and include the new way of applying three source terms. 
% update v6: Add mesh adaptivity
% update v7: use new viscosity scheme: see (Costa et al, 2009) 

% If code or comments are written by Catherine Booth, it is annotated CAB

%% 
%%
clear;
clc
warning('off')

% CAB - calling the m-file which contains the settings for the run.
Model_setting_v8_sill_exe_master; 

Inputs= readtable('Input_Files/1AA_2Phase_code_RESTART.txt');
[r,~] = size(Inputs);
names=string(Inputs.Var2);
Number=Inputs.Var3;
for i=1:r-1
    assignin('base',names(i),Number(i))
end
filename4restart = names(r);
filename4restart = string(filename4restart{1});

% turn to 1 to restart a simulation, otherwise load the initiator
if Restart==1
    path='';
    disp('Loading '+ filename4restart)
    load(filename4restart);  
    disp('Loaded')
   % Plot_settings;
    To_plot_CD_step=0;

    File_echo=fopen('echo_screen.txt','a');
    File_comp=fopen('mass_cons.txt','a');
    File_breaks = fopen('Break_of_tolerance.txt','a');
    File_precision = fopen('In_Depth_break_data.txt', 'a');
else


    kt_background=kt0;
    %%
    counter=0;
    break_counter=0;
    break_percent=0;
    output_counter=0;
    N=length(cellz); %CAB - cellz is defined in Model_setting line 170
    pos=1;
    Conservative_pos=[find(cellz>=Show_z(1),1),find(cellz<=Show_z(2),1,'last')];
    
    OG_composition=sum(Cb);
    
    dt_scale=1;
    
    dt=dt0;
    dt_old=1e10;
    dt_old_old=1e10;
    
    C_values=zeros(N+1,1);
    Cphi_all=[C_all(1:N).*phi(1:N); C_all(N+1:2*N).*phi(N+1:2*N)];
    
    gap_counter=0;
    
    T_old=T;
    phi_old=phi;
    C_all_old=C_all;
    iter=0;
    %rhof=0*C_all(1:N);
    %rhom=0*C_all(1:N);
    Counter_step_intrusion=0;
    Passing_time_since_intrude=0;
    
    iter=0;
    
    
    %% CAB - Sill intrusion
    SillRate = injection_rate*1e3*Year;
    SillNo = Tot_Sill/Sill_length;
    Initial_dN = depthN;
    Initial_dC = depthCN;
    SillDens = mean(injection_densB);
    
    % CAB - end
    
    %% CAB - initial output
    if Record_data==1
        
            %Save the variables as .txt
        
                FilenameO = 'output_0_CELLS.txt';
    
                File_Out = fopen(FilenameO, 'w');
            
                fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');
    
                if (HHJPet ==1 || SSPD==1)
                    if To_cal_mc==1
                        fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                            'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)','cb',  ...
                            'Cb SiO2','Cs SiO2', 'Cl SiO2', 'H_MF', 'C_MF', 'R_MF','H_MF_t','C_MF_t', 'R_MF_t', 'C_CB', 'R_CB', 'C_CB_t', 'R_CB_t');
                    
                        for i=1:1:N
                            fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f    \n', ...
                                [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                T(i); Ts(i); Tl(i); Cb(i);  Cb(i)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); ...
                                C_all(N+i)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1);0;0;0;0;0;0;0;0;0;0]);
                        end
                    else
                        fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', 'Depth (km)', ...
                            'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)','cb',  ...
                            'Cb SiO2','Cs SiO2', 'Cl SiO2');
                    
                        for i=1:1:N
                            fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f  \n', ...
                                [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                T(i); Ts(i); Tl(i); Cb(i);  Cb(i)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); ...
                                C_all(N+i)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)]);
                        end
                    end
    
    
    
                elseif FourMPD==1
                    if To_cal_mc==1
                        fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                                    'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)','cb',  ...
                                    'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld', 'H_MF', 'C_MF', 'R_MF','H_MF_t','C_MF_t', 'R_MF_t', 'C_CB', 'R_CB', 'C_CB_t', 'R_CB_t');
            
                        for i=1:1:N
                            fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                                [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                T(i); Ts(i); Tl(i); Cb(i);  Cb(i)*(MgO_range(2)-MgO_range(1))+MgO_range(1); ...
                                C_all(N+i)*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); ol(i); opx(i); cpx(i); feld(i); 0; 0; 0; 0; 0;0;0;0;0;0;]);
                               
                        end
            
                        
                    else
                        fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', 'Depth (km)', ...
                            'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)','cb',  ...
                            'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld');
                    
                        for i=1:1:N
                            fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                                [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                T(i); Ts(i); Tl(i); Cb(i);  Cb(i)*(MgO_range(2)-MgO_range(1))+MgO_range(1); ...
                                C_all(N+i)*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); ol(i); opx(i); cpx(i); feld(i)]);
                        end
    
                           
                    end
                end
            
                fclose(File_Out);
            
            
                FilenameON = 'output_time_0_NODES.txt';
            
                File_OutN = fopen(FilenameON, 'w');
            
                fprintf(File_OutN, '%4s  %10.10f %3s \n', 'Time', Time/Year/1000, 'kyr' );
            
                fprintf(File_OutN, '%10s \t %10s \t %10s \n', 'Depth (km)', ...
                    'Sol. Vel.', 'Liq. Vel');
            
                for i=1:1:N+1
                    fprintf(File_OutN, '%10.3f \t %4e \t %4e  \n', ...
                        [nodez(i)/1000-Base_crust; u_all(N+1+i); u_all(i)]);
                end
            
                fclose(File_OutN);
    
                output_counter=output_counter+1;
    
    
        
        
        if Output_Flag==0
        else
            % Save the variables in matlab.
            FileNameO = 'output_time_0.mat';
            save(FileNameO)                
        end
    
        Record_index=2;

    end

    if ~exist('Active_region','var')
        Active_region=[0,0];
    end
    
    if exist('Fixed_record','var')
       if  Fixed_record==1
          Last_record_time=0; 
       end
    end
    
    
    %% CAB - Intrude initial sill at initial intrusion depth
    File_echo=fopen('echo_screen.txt','w');
    File_comp=fopen('mass_cons.txt','w');
    File_breaks = fopen('Break_of_tolerance.txt','w');
    File_precision = fopen('In_Depth_break_data.txt', 'w');
    
    fprintf(File_comp, '%10s \t %10s \t %10s \n', ...
                               'Time (ka)', 'Mass cons', 'Mass Res' );
    
    fprintf(File_breaks, '%10s \t %10s \t %10s \t %10s \n ', 'Time (ka)', 'Percentage', 'no. breaks', 'no.timesteps');
    
    
    if res_type==0
        fprintf(File_precision, '%10s \t %10s \t %10s \t %10s \t %10s \t \n', 'Time (ka)', 'Precision', 'Residual', 'Residual1 cb', 'Residual2 h');
    else
        fprintf(File_precision, '%10s \t %10s \t %10s \t \n', 'Time (ka)', 'Precision', 'Residual');
    end
    
    
    fprintf(File_echo, '%s %5.5f %s %5.5f %6s\n', 'Sill intruded at', ...
                                cellz(Initial_dC)/1000-Base_crust, 'km by initial depth at', Time/Year/1000, 'ka');
    
    
    
    Sillcomp=injection_Cb;
    SillH=injection_H;
    SillTemp=injection_T;
    SillMF=injection_phi;
    SillCf=injection_Cl;
    SillCm=injection_Cs;
    
    sill_intrusion_master
    
    dt=0;
    Newton_solver2;
    
    Cb_all0=sum(Cb.*dz');
    % Update sillcount
    SillCount=1;



%%

%if Create_video==1
 %   open(Video_handel);
%end

end

%%


tic;
dt_history=zeros(1,10);
dt_history_index=1;
Easy_converge=0;

while Time<End_time
    %% CAB - INTRUDE SILLS         
        if SillCount<SillNo 
            if Time > SillCount*SillRate-dt
                if Time < SillCount*SillRate+dt

                    % Finding intrusion depth - either initial depth or
                    % density.

                    depthCN=0;
                    depthN=0;
                    

                    if DensF>0 && MFSill_flag>0

                        for i=1:1:N

                                if depthCN==0 && phi(i)>MFSill_int
                                    if RandF>0 % random intrusion of sills around depth
                                        RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                        depthCN=i+RF;
                                        depthN=i-1+RF;
                                        disp(Time)
                                        disp('Sill intruded randomly around melt fraction difference')
                                        fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                                            cellz(depthCN)/1000-Base_crust, 'km by melt fraction difference with a random shifting of', ...
                                            RF*dzF, 'm, at', Time/Year/1000, 'ka');
                
                                    else
                                        depthCN=i;
                                        depthN=i-1;
                                        disp(Time)
                                        disp('Sill intruded by melt fraction difference')
                                        fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                                            cellz(depthCN)/1000-Base_crust, 'km by melt fraction difference at', Time/Year/1000, 'ka');
                                    end
                                end

                                if SillDens>=rho_b(i) && depthCN==0 
                                    if RandF>0 %random intrusion of sills around depth
                                        RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                        depthCN=i+RF;
                                        depthN=i+RF-1;
                                        disp(Time)
                                        disp('Sill intruded randomly around density' )
                                        fprintf(File_echo, '%6s %5.5f %6s %5.5f %6s %5.5f %6s \n', 'Sill intruded at', ...
                                        cellz(depthCN)/1000-Base_crust, 'km by density with a random shifting of', ...
                                        RF*dzF, 'm, at', Time/Year/1000, 'ka');
            
                                    else %no random intrusion of sills around depth
                                        depthCN=i;
                                        depthN=i-1;
                                        disp(Time)
                                        disp('Sill intruded by density')
                                        fprintf(File_echo, '%6s %5.5f %6s %5.5f %6s \n', 'Sill intruded at', ...
                                        cellz(depthCN)/1000-Base_crust, 'km at', Time/Year/1000, 'ka');
                                    end 
                                
                                
                                
                                end
                        end

                         if depthCN==0 
                            if RandF>0 % random intrusion of sills around depth
                                RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                depthCN=Initial_dC+RF;
                                depthN=Initial_dN+RF;
                                disp(Time)
                                disp('Sill intruded randomly around initial depth (overaccretion) ')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion)  with a random shifting of', ...
                                    RF*dzF, 'm, at', Time/Year/1000, 'ka');
        
                            else
                                depthCN=Initial_dC;
                                depthN=Initial_dN;
                                disp(Time)
                                disp('Sill intruded by initial depth (overaccretion)')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion) at', Time/Year/1000, 'ka');
                            end
                         end
                 


                    
                    
                    elseif DensF>0 % density intrusion
                        for i = 1:1:N
                            if SillDens>=rho_b(i) && depthCN==0 
                                if RandF>0 %random intrusion of sills around depth
                                    RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                    depthCN=i+RF;
                                    depthN=i+RF-1;
                                    disp(Time)
                                    disp('Sill intruded randomly around density' )
                                    fprintf(File_echo, '%6s %5.5f %6s %5.5f %6s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by density with a random shifting of', ...
                                    RF*dzF, 'm, at', Time/Year/1000, 'ka');
        
                                else %no random intrusion of sills around depth
                                    depthCN=i;
                                    depthN=i-1;
                                    disp(Time)
                                    disp('Sill intruded by density')
                                    fprintf(File_echo, '%6s %5.5f %6s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km at', Time/Year/1000, 'ka');
                                end 
                                
                                
                                
                            end
                        end

                         if depthCN==0 
                            if RandF>0 % random intrusion of sills around depth
                                RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                depthCN=Initial_dC+RF;
                                depthN=Initial_dN+RF;
                                disp(Time)
                                disp('Sill intruded randomly around initial depth (overaccretion) ')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion)  with a random shifting of', ...
                                    RF*dzF, 'm, at', Time/Year/1000, 'ka');
        
                            else
                                depthCN=Initial_dC;
                                depthN=Initial_dN;
                                disp(Time)
                                disp('Sill intruded by initial depth (overaccretion)')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion) at', Time/Year/1000, 'ka');
                            end
                        end



        
                    


                    elseif MFSill_flag>0 % melt fraction                        
                        if RandF>0
                            RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                        else
                            RF=0;
                        end
                        index=find(phi(1:N)>MFSill_int,1,'first');
                        if ~isempty(index)
                            depthCN=index+RF;
                            depthN=depthN+1;
                        else
                            index=find(nodez>Last_intrusion_depth,1,'first');
                            depthCN=index+RF;
                            depthN=depthN+1;
                        end
                        % for i=1:1:N
                        % 
                        %     if depthCN==0 && phi(i)>MFSill_int
                        %         if RandF>0 % random intrusion of sills around depth
                        %             RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                        %             depthCN=i+RF;
                        %             depthN=i-1+RF;
                        %             disp(Time)
                        %             disp('Sill intruded randomly around melt fraction difference')
                        %             fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                        %                 cellz(depthCN)/1000-Base_crust, 'km by melt fraction difference with a random shifting of', ...
                        %                 RF*dzF, 'm, at', Time/Year/1000, 'ka');
                        % 
                        %         else
                        %             depthCN=i;
                        %             depthN=i-1;
                        %             disp(Time)
                        %             disp('Sill intruded by melt fraction difference')
                        %             fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                        %                 cellz(depthCN)/1000-Base_crust, 'km by melt fraction difference at', Time/Year/1000, 'ka');
                        %         end
                        %     end
                        % end
                        % 
                        % if depthCN==0 
                        %     if RandF>0 % random intrusion of sills around depth
                        %         RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                        %         depthCN=Initial_dC+RF;
                        %         depthN=Initial_dN+RF;
                        %         disp(Time)
                        %         disp('Sill intruded randomly around initial depth (overaccretion) ')
                        %         fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                        %             cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion)  with a random shifting of', ...
                        %             RF*dzF, 'm, at', Time/Year/1000, 'ka');
                        % 
                        %     else
                        %         depthCN=Initial_dC;
                        %         depthN=Initial_dN;
                        %         disp(Time)
                        %         disp('Sill intruded by initial depth (overaccretion)')
                        %         fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                        %             cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion) at', Time/Year/1000, 'ka');
                        %     end
                        % end

                    elseif OverF==0 % overaccretion

                        if depthCN==0 
                            if RandF>0 % random intrusion of sills around depth
                                RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                depthCN=Initial_dC+RF;
                                depthN=Initial_dN+RF;
                                disp(Time)
                                disp('Sill intruded randomly around initial depth (overaccretion) ')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (overaccretion)  with a random shifting of', ...
                                    RF*dzF, 'm, at', Time/Year/1000, 'ka');
        
                            else
                                depthCN=Initial_dC;
                                depthN=Initial_dN;
                                disp(Time)
                                disp('Sill intruded by initial depth (overaccretion)')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (overaccretion) at', Time/Year/1000, 'ka');
                            end
                        end


                    else  % underaccretion
                        if depthCN==0 
                            if RandF>0 % random intrusion of sills around depth
                                RF=round((2.*rand_injectN).*rand(1,1)-rand_injectN);
                                depthCN=Initial_dC+RF-SillCount*SillNodez;
                                depthN=Initial_dN+RF-SillCount*SillNodez;
                                disp(Time)
                                disp('Sill intruded randomly around initial depth (underaccretion) ')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s %5.5f %6s\n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion)  with a random shifting of', ...
                                    RF*dzF, 'm, at', Time/Year/1000, 'ka');
        
                            else
                                depthCN=Initial_dC-SillCount*SillNodez;
                                depthN=Initial_dN-SillCount*SillNodez;
                                disp(Time)
                                disp('Sill intruded by initial depth (underaccretion)')
                                fprintf(File_echo, '%s %5.5f %s %5.5f %6s \n', 'Sill intruded at', ...
                                    cellz(depthCN)/1000-Base_crust, 'km by initial depth (underaccretion) at', Time/Year/1000, 'ka');
                            end
                        end
                    end
                        
                        
                    
                    sill_intrusion_master;

                    Cb_all0=sum(Cb.*dz');
                    %Update sill count 
                    SillCount=SillCount+1;
                end
            end
        end
%%

        % increase timestep if needed, not for Newton's method
        if (fixed_dt==0 && iter<Min_iter) && ~Use_Newton
            [dt, ~] = dynamic_dt_master(1,iter,dt,Min_dt,Max_dt, N_dt, Courant, dz, u_all,N,inc_N);

            if dt==Max_dt
                    fprintf(File_echo, '%6s %5.5f %6s \n', 'Maximum time step reached', Time/Year/1000, 'ka');
            end
        end

        if (Adaptive_mesh==1 && mod(counter,Adaptive_step_gap)==1 && Time-Last_adapted>Adaptive_step_time) || Just_intruded
            Adapt_mesh_master
        end

        counter=counter+1;
        improve=1;
        iter=0;
        phi_old=phi;
        u_all_old=u_all;
        H_old=H;
        T_old=T;
        C_all_old=C_all;
        Cphi_all_old=Cphi_all;
        Cb_old=Cb;
        Gamma=zeros(2*length(cellz),1); % Melft
        Type_record=zeros(N,2);
        Type_fix=zeros(N,1);
        type_index=1;
        C_values=zeros(N+1,1);

        if Has_volatile==1
            S_old=S;
            Cs2_old=Cs2;
            Cl2_old=Cl2;            
        end
        
        %Injecting=0;
        %Sill_injection=0;
        %CAB - change this
        %if Sill_injection==1
         %   for i=1:size(injection_time_range,1)
          %      if Time>=injection_time_range(i,1) && Time<=injection_time_range(i,2)
           %         Injecting=1;
            %    end
            %end
        %end
       
%%
    if Use_Newton==1
        % Two_phase_Newton;
        if dt_intended<1e-5*Year
            dt_intended=1*Year;
        end
        dt=dt_intended;
        % dt_history(dt_history_index)=dt;
        % if dt_history_index<10
        %     dt_history_index=dt_history_index+1;
        % else
        %     dt_history_index=1;
        % end
        % if (max(dt_history)-min(dt_history))/max(min(dt_history),1e-6)<1.3 && min(dt_history)<1e-1*Year
        %     dt=Max_dtY*Year/10;
        % end
        
        Newton_solver2;
        
        if (iter==Max_Newton_iter)|| dt<Max_dtY*Year/1e4
            min_dx=min_dx/2;
            Adapt_mesh_master;
            dt=dt0Y*Year;
            Newton_solver2;
            min_dx=min_dx*2;

            if iter==Max_Newton_iter
            error('Newton solver not converged before max iteration');
            end
            Easy_converge=0;
        else
            Easy_converge=Easy_converge+1;
        end
        
        if Easy_converge>10
            min_dx=min_dx0;
            Easy_converge=0;
        end
        
        Cb_all=sum(Cb.*dz');
        disp(['conservation:' num2str(Cb_all/Cb_all0)])
    else
        kt=kt_background*ones(N,1);
        ind1=find(dz<max(dz)*0.9,1,'first');
        ind2=find(dz<max(dz)*0.9,1,'last');
        kt(ind1+1:ind2-1)=kt0;
        kc=kc0*ones(N,1);
        % Max_iter=30;

        Cphi_all_nonlinear=zeros(2*N,1);
        H=T.*cp+phi(1:N).*Lf;


        % CAB - get rid of this?
       % if Injecting==1
        %    [U_source,U_source_int,H_source,C_source,Index]=Cal_sources(injection_range,injection_rate,injection_T,cp,Lf,injection_phi,injection_Cl,injection_Cs,sink_location,nodez);
       % else
            U_source=zeros(N+1,1);
            U_source_int=zeros(N+1,1);
            H_source=zeros(N,1);
            C_source=zeros(N,4);
        %end
        %%
        %CAB - iteration for the timestep starts
        while improve>Precision && iter<Max_iter
            %%  Solving Momentum equation
            temp1=project_cell2node_master(nodez,dz,phi(1:N));     %project phi from cell center to nodes, fluid part;
            temp2=project_cell2node_master(nodez,dz,phi(N+1:2*N)); %project phi from cell center to nodes;
            Sum=temp1+temp2;
            phi_on_nodes=[temp1./Sum;temp2./Sum];
            phi_on_nodes=max(0,phi_on_nodes);
            phi_on_nodes=min(1,phi_on_nodes);

            % forcing 0 nodal phi values (10/01/2025)
            index=find(phi==0);
            phi_on_nodes(index)=0;
            phi_on_nodes(index+1)=0;

            if perm_min_vel_FLAG==1 && perm_min_FLAG==1
                ZerosU=union(find(abs(phi_on_nodes(1:N+1)-1)<Precision),find(abs(phi_on_nodes(1:N+1))<(Precision+perm_min_MF)));
            else
                ZerosU=union(find(abs(phi_on_nodes(1:N+1)-1)<Precision),find(abs(phi_on_nodes(1:N+1))<Precision));
            end

            ZerosU=union(ZerosU,[1,N+1]);
            
            %if Injecting==1
%                 for i=1:size(injection_range,1)+1
%                     ZerosU=setdiff(ZerosU,U_source_index(i,1):U_source_index(i,2));
%                 end

             %   ZerosU_values=zeros(length(ZerosU),1);
              %  ZerosU_values(1)=0;
               % for i=2:length(ZerosU)
                %    ZerosU_values(i)=U_source_int(ZerosU(i));
                %end
            %else
              ZerosU_values=zeros(length(ZerosU),1);
            %end
            Nonzero=setxor(1:(N+1),ZerosU);            
            Nonzero=reshape(Nonzero,[1,length(Nonzero)]);   %MUST HAVE THIS ffs!
            
            %% update the viscosities based on the composition

            if (HHJPet==1||SSPD==1)
                Ssio2_scaled = C_all(1:N);

            elseif FourMPD==1
                Ssio2_dim = 0*C_all(1:N);
                % calculating melt SiO2 from melt MgO
                %dim MgO
                C_all_dim = C_all(1:N)*(MgO_range(2)-MgO_range(1))+MgO_range(1);
                for i=1:N
                    if C_all_dim(i)>=crit_mg_si_melt
                        Ssio2_dim(i) = min_sio2_visc;
                    else
                        Ssio2_dim(i) = max_sio2_m*C_all_dim(i)+max_sio2_c;
                    end
                end
                    
                Ssio2_scaled=(Ssio2_dim- SiO2_range(1))/(SiO2_range(2)-SiO2_range(1));


            end
            

            %
            
            mu_f_dynamics=10.^((mu_f2-mu_f1)*Ssio2_scaled+mu_f1);  %traditional            
            mu_f_dynamics=project_cell2node_master(nodez,dz,mu_f_dynamics)/scaling_factor;
            
            mu_m_dynamics=mu_all(min(ceil((phi_on_nodes(1:N+1))*N_vis)+1,N_vis))/scaling_factor;

    %% densities       

            if (HHJPet==1 || SSPD==1)

                rhof=rhof_1*(1-C_all(1:N))+(C_all(1:N))*rhof_2;
            elseif FourMPD==1
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
            rhom=rhom_1*(1-C_all(N+1:2*N))+(C_all(N+1:2*N))*rhom_2;
            
            
            rho_b=rhof.*phi(1:N)+rhom.*phi(N+1:2*N);  % Bulk density
            
            
            rhof=project_cell2node_master(nodez,dz,rhof);
            rhom=project_cell2node_master(nodez,dz,rhom);
            
            for i=Nonzero
                C_values(i)=c_gen_master(mu_f_dynamics(i),grain_size,phi_on_nodes(i),C_value_A,C_value_B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF);
            end

            if MeltSeg_Flag>0
                [u_all, Portion]=Velocity_solve_portion_with_source_master(phi_on_nodes,u_all_old,Nonzero,ZerosU,ZerosU_values, C_values,mu_f_dynamics,mu_m_dynamics,0.5,rhof,rhom, Precision,length(nodez),dz,g,U_source_int);
            end
            %% Estimate dt
            if fixed_dt==0 && iter==Max_iter-1
                [dt,iter] = dynamic_dt_master(0,iter,dt,Min_dt,Max_dt,N_dt,Courant, dz, u_all, N,inc_N);

                if dt==Min_dt
                    fprintf(File_echo, '%6s %5.5f %6s \n', 'Minimum time step reached', Time/Year/1000, 'ka');
                end

            end
       %% Enthalpy and components transport
            H_nonlinear=H;
            u_bar=u_all(1:N+1).*phi_on_nodes(1:N+1)+u_all(N+2:2*N+2).*phi_on_nodes(N+2:2*N+2);
            
            H=CV_enthalpy_solve_source_master(H_old,phi,u_all,u_bar,dz,dt,kt,cp,Lf,BC_T_type,BC_T,H_source);
            
    %           H=H_old;
            % Solve for composition
            Cb_nonlinear=Cb;

            Cphi_all=CV_composition_solve_source_master(C_all,phi,phi_old,Cphi_all_old, u_all,dz,dt,kc,BC_C_type, BC_C_value,C_transport_method,C_source);
            % Update bulk composition
            Cb=Cphi_all(N+1:2*N)+Cphi_all(1:N);
      %% Update melt fraction and composition due to phase_diagram

            phi_old_nonlinear=phi;

     
            if HHJPet == 1

                if N<50 || Use_parallel==0
                    for i=1:N
                        [phi(i), T(i), C_all(i), C_all(i+N),TYPE(i), Tl(i), Ts(i)]=poro_component_solve_JPET_master(Cb(i), H(i),A1, B1, C1, A2, B2, C2, ae, Lf, cp,Precision_PD);
                    end
                else
                   if isempty(gcp('nocreate'))
                       parpool(Num_core_use);
                   end

                   parfor i=1:N
                       warning('off', 'all');
                       %set(0,'DefaultFigureVisible','off')
                       [phi(i), T(i), C_all(i), C_all_dummy(i),TYPE(i), Tl(i), Ts(i)]=poro_component_solve_JPET_master(Cb(i), H(i),A1, B1, C1, A2, B2, C2, ae, Lf, cp,Precision_PD);
                   end

                   C_all(N+1:2*N) = C_all_dummy;

                end


            end


            if SSPD ==1

                if N<50 || Use_parallel==0
                    for i=1:N
                        [phi(i), T(i), C_all(i), C_all(i+N),TYPE(i), Tl(i), Ts(i)]=poro_component_solve_solid_master(Cb(i), H(i), A1, B1, C1, alpha, n_PD ,Lf, cp, Precision_PD,step_size);
                    end
                else
                   if isempty(gcp('nocreate'))
                       parpool(Num_core_use);
                   end

                   parfor i=1:N
                       warning('off', 'all');
                       %set(0,'DefaultFigureVisible','off')
                       [phi(i), T(i), C_all(i), C_all_dummy(i),TYPE(i), Tl(i), Ts(i)]=poro_component_solve_solid_master(Cb(i), H(i), A1, B1, C1, alpha, n_PD ,Lf, cp, Precision_PD,step_size);
                   end


                   C_all(N+1:2*N) = C_all_dummy;
                end
                    

            end


            if FourMPD==1

                if N<50 || Use_parallel==0
                    for i=1:N
                        [phi(i),T(i),C_all(i), C_all(i+N), TYPE(i), Tl(i), Ts(i)]=poro_component_solve_4_components_v2_master(Cb(i), H(i), Lf, cp, min_mgo, max_mgo, crit_T1, crit_T2,...
                                                                 crit_cb2, crit_cb3, crit_cb4,crit_cb4_OG, crit_cb5,liq_k1, liq_a1, liq_b1,liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,...
                                                                lm1, lc1, lm3, lc3, lm8, lc8, mk, ck, BT, BC, DT, DC,JC, olpx_C,olpxm, olpxc, s5T, sol_k2,liq_P2_C,s1m,s1c, ...
                                                                lowest_end_T, liq_P3_C,l13_k,l13_a,l13_b,lowest_k,s8m,s8c,crit_T3,liq_P4_C, Precision_PD,step_size);
                    end
                else
                   if isempty(gcp('nocreate'))
                       parpool(Num_core_use);
                   end

                   parfor i=1:N
                       warning('off', 'all');
                       %set(0,'DefaultFigureVisible','off')
                       [phi(i),T(i),C_all(i), C_all_dummy(i), TYPE(i), Tl(i), Ts(i)]=poro_component_solve_4_components_v2_master(Cb(i), H(i), Lf, cp, min_mgo, max_mgo, crit_T1, crit_T2,...
                                                                 crit_cb2, crit_cb3, crit_cb4,crit_cb4_OG, crit_cb5,liq_k1, liq_a1, liq_b1,liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,...
                                                                lm1, lc1, lm3, lc3, lm8, lc8, mk, ck, BT, BC, DT, DC,JC, olpx_C,olpxm, olpxc, s5T, sol_k2,liq_P2_C,s1m,s1c, ...
                                                                lowest_end_T, liq_P3_C,l13_k,l13_a,l13_b,lowest_k,s8m,s8c,crit_T3,liq_P4_C, Precision_PD,step_size);
                    end


                   C_all(N+1:2*N) = C_all_dummy;
                end
            end

            
        
            phi(N+1:N*2)=1-phi(1:N);
            
            Cphi_all_nonlinear(:)=Cphi_all(:);
            Cphi_all=[C_all(1:N).*phi(1:N); C_all(N+1:2*N).*phi(N+1:2*N)];
            Cb=Cphi_all(N+1:2*N)+Cphi_all(1:N);

          
            if res_type == 0
                improve1=max(abs(Cb(1:N)-Cb_nonlinear(1:N)));
                improve2=max(abs(H(1:N)-H_nonlinear(1:N)))/Lf;
                improve=max(improve1,improve2);
            else
                improve=max(abs(phi(1:N)-phi_old_nonlinear(1:N)));
            end

            num_bad=100;   
            %if iter>100
            %    figure(8)
             %   plot(iter,improve/Precision,'o')
              %  hold on
            %end
%             Check_convergence; % Enable to show the convergence on the Cm plot
            iter=iter+1;            
        end
    end

        %%%%%% HH added for calculating the contribution of mel fraction/Cb
        if To_cal_mc==1&& iter>0
            use_new=0;  %0: old, 1: new

            phi_temp=CV_saturation_solve_explicit_master(u_all,phi_old,phi,dz,dt,use_new);
            %         phi_temp=CV_saturation_solve(u_all,phi_old,dz,dt,zeros(2*length(cellz),1),Precision); %Gamma= 0
            Con_M=phi_temp(1:N)-phi_old(1:N);

            Weight_new=0.5;
            Dif=(C_all(1:N)-C_all(N+1:2*N))*Weight_new+(C_all_old(1:N)-C_all_old(N+1:2*N))*(1-Weight_new);
            Index=Dif~=0;

            Temp=(C_all-C_all_old).*(phi*Weight_new+phi_old*(1-Weight_new));
            Temp=Temp(1:N)+Temp(N+1:2*N);
            Con_T=zeros(N,1);
            Con_T(Index)=-Temp(Index)./Dif(Index);

            if react_try==1

                Temp_C = CV_reactiveflow(dt,dz,C_all,phi,u_all,phi_old,u_all_old,Weight_new);
                %Temp_C_old = CV_reactiveflow(dt,dz,C_all_old,phi_old,u_all_old);
    
                %Temp_c_diff = Temp_C*Weight_new + Temp_C_old*(1-Weight_new)
                Con_C=zeros(N,1);
                Con_C(Index)=-Temp_C(Index)./Dif(Index);
                Con_C_dummy=phi(1:N)-phi_old(1:N)-Con_M-Con_T;

            else
            
                Con_C=phi(1:N)-phi_old(1:N)-Con_M-Con_T;
                Con_C_dummy=phi(1:N)-phi_old(1:N)-Con_M-Con_T;

            end


            % Calculate the places the contribution needs to be recoreded
            Record_cut=ones(N,1);
            if Contribution_cut<1
                Record_cut=(phi(1:N)<Contribution_cut);
            end

            Only_record_decrease=0;
            if Only_record_decrease==1
                Record_cut=Record_cut & (phi(1:N)<phi_old(1:N));
            end

            Con_T=Con_T.*Record_cut;
            Con_C=Con_C.*Record_cut;
            Con_C_dummy=Con_C_dummy.*Record_cut;
            Con_M=Con_M.*Record_cut;


            ACon_M=ACon_M+Con_M;
            ACon_T=ACon_T+Con_T;
            ACon_C=ACon_C+Con_C;
            ACon_C_dummy=ACon_C_dummy+Con_C_dummy;

           

            %     Index=phi(1:N)>1e-10;
            %     disp(max(abs((Con_T(Index)+Con_C(Index)+Con_M(Index)+phi_old(Index)-phi(Index))./phi(Index))));

            CD_M=CD_M+Con_M.*(C_all(1:N)-C_all(N+1:2*N)); %CD caused by compaction
            CD_R=CD_R+Con_C.*(C_all(1:N)-C_all(N+1:2*N));
            CD_R_dummy=CD_R_dummy+Con_C_dummy.*(C_all(1:N)-C_all(N+1:2*N));

            CD_Ms=Con_M.*(C_all(1:N)-C_all(N+1:2*N)); %CD caused by compaction
            CD_Rs=Con_C.*(C_all(1:N)-C_all(N+1:2*N));
            CD_Rs_dummy=Con_C_dummy.*(C_all(1:N)-C_all(N+1:2*N));
        end

        if FourMPD==1
            for i=1:N
                [ol(i),opx(i),cpx(i),feld(i),ol_mg(i),px_mg(i),opx_mg(i),cpx_mg(i),feld_mg(i),ol_mn(i),sol_mn(i)] = A_4_mineral_proportions_v2_master(...
                                    ol(i),opx(i),cpx(i), feld(i), ...
                                    T(i),phi(i),phi_old(i),Cb(i),C_all(N+i),C_all(i),crit_T1,crit_T2,crit_cb1,crit_cb2,crit_cb3,crit_cb4,crit_cb5,...
                                    s1m,s1c,s3m,s3c,s4m,s4c,s6m,s6c,s7m,s7c,s8m,s8c,s11m,s11c, ...
                                    DC, DT, EC, ET, NC, NT, copx_T_m, copx_T_c, copx_C_m, copx_C_c, ...
                                    px_k, olpx_C, ol_i_m0, ol_i_c0, ol_i_m1, ol_i_c1, ol_i_m2, ol_i_c2, ol_i_m3, ol_i_c3,...
                                    ol_e_m1, ol_e_c1, ol_e_m2, ol_e_c2, ol_e_m3, ol_e_c3, ol_e_m4, ol_e_c4, ...
                                    ol_T_high_m, ol_T_high_c, ol_T_low_m, ol_T_low_c, ol_T_lowest_m, ol_T_lowest_c, ...
                                    crit_sol_mn_cb, ol_mn_bridge,sol_mn_k_hCb, sol_mn_k_lCb, mn_k_m, mn_k_c, ...
                                    a_mn, b_mn, c_mn,sol_mni_m, sol_mni_c, ol_T_crit_cb2, ol_i_m, ol_i_c, ol_mg_m,...
                                    cb_low_high, copx_cb_m_1h, copx_cb_c_1h, copx_cb_m_1l, copx_cb_c_1l, ...
                                    copx_T_m_1h, copx_T_c_1h, copx_T_m_1l, copx_T_c_1l, ...
                                    copx_cb_m_2h, copx_cb_c_2h, copx_cb_m_2l, copx_cb_c_2l, ...
                                    copx_T_m_2h, copx_T_c_2h, copx_T_m_2l, copx_T_c_2l,...
                                    px_cb_m, px_cb_c, liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,liq_P2_C,Ts(i),...
                                    max_mgo, min_mgo);
                
            end
        end
        %if iter>100
         %   figure(8)
          %  clf;
        %end
    %CAB - iteration for the timestep ends
    % CAB - added improve which highlights the precision level
    %disp([counter, iter, improve/Precision])
    
    if improve/Precision>1 
        fprintf(File_echo, '%6s %5.5f %6s \n', ['Difference is a magnitude ' ...
            'greater than precision at'], Time/Year/1000, 'ka');

        if res_type==0
            fprintf(File_precision, '%10.5f \t %10.5f \t %10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, Precision, improve, improve1, improve2);
        else
            fprintf(File_precision, '%10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, Precision, improve);
        end

        break_counter=break_counter+1;
    end
     break_percent = break_counter/counter*100;

   % if HighTsF==1
        %for i=1:1:N
            %if Phase_Type==0
              %  if i<=depthTsCN
                  %  Ts(i)=TsD;
                  %  Tl(i)=A1*cb(i)^2+B1*cb(i)+C1;
                %else
                   % Ts(i)=TsU;
                   % Tl(i)=A1*cb(i)^2+B1*cb(i)+C1;
               % end
            %else
               % if i<=depthTsCN
               %     Ts(i)=TsD;
                %    [~,Tl(i),~,~]=solid_state_Ts_Tl_local(Cb(i),alpha_order,n_order,Lf,cp,A1,B1,C1);
               % else
                %    [Ts(i),Tl(i),~,~]=solid_state_Ts_Tl_local(Cb(i),alpha_order,n_order,Lf,cp,A1,B1,C1);
           %     end
            %end
        %end
  %  else
      %  for i=1:1:N
            %if Phase_Type==0
              %  Ts(i)=TsU;
              %  Tl(i)=A1*cb(i)^2+B1*cb(i)+C1;
            %else
             %   [Ts(i),Tl(i),~,~]=solid_state_Ts_Tl_local(Cb(i),alpha_order,n_order,Lf,cp,A1,B1,C1);
           % end
        %end
   % end


     
    composition=sum(Cb);
    Mass_cons=composition/OG_composition;
    Mass_res=abs(composition-OG_composition);

   

%%
    Time=Time+dt;
    Passing_time_since_intrude=Passing_time_since_intrude+dt;
    dt_old_old=dt;
    dt_old=dt;
    Step_counts=Step_counts+1;
    Counter_step_intrusion=Counter_step_intrusion+1;

    

   
  

%% CAB - OUTPUT

 if Record_data==1

         if (length(Record_time)>=Record_index)

         if (Time)>=Record_time(Record_index)-dt 
             if (Time)<Record_time(Record_index)+dt

                Record_index=Record_index+1;
                vv=['Output ',num2str(Time/Year/1000), ' ka'];
                    disp(vv)

                
                    %Save the variables as .txt
                    
    
                    FilenameO = ['output_',num2str(output_counter),'_CELLS.txt'];
    
                    File_Out = fopen(FilenameO, 'w');
                
                    fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

                    if (HHJPet==1 || SSPD==1)

                        if To_cal_mc==1

                            fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb SiO2','Cs SiO2', 'Cl SiO2 ', 'Dens', 'H_MF', 'C_MF', 'R_MF','H_MF_t','C_MF_t', 'R_MF_t', 'C_CB', 'R_CB', 'C_CB_t', 'R_CB_t');
                        
                            for i=1:1:N
                                fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                                    [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                    T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(N+i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1) ...
                                    ; (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); rho_b(i); ...
                                    ACon_T(i); ACon_M(i); ACon_C(i); Con_T(i); Con_M(i); Con_C(i); ...
                                    CD_M(i)*(SiO2_range(2)-SiO2_range(1));CD_R(i)*(SiO2_range(2)-SiO2_range(1));CD_Ms(i)*(SiO2_range(2)-SiO2_range(1));CD_Rs(i)*(SiO2_range(2)-SiO2_range(1))]);
                            end

                        else

                            fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb SiO2','Cs SiO2', 'Cl SiO2 ', 'Dens');
                        
                            for i=1:1:N
                                fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \n', ...
                                    [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                    T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(N+i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1) ...
                                    ; (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); rho_b(i)]);
                            end

                        end

                    elseif FourMPD==1

                        if To_cal_mc==1

                            fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld', 'Dens', 'H_MF', 'C_MF', 'R_MF','R_MF_dum','H_MF_t','C_MF_t', 'R_MF_t','R_MF_t_dum', 'C_CB', 'R_CB', 'R_CB_dum', 'C_CB_t', 'R_CB_t', 'R_CB_t_dum');
                        
                            for i=1:1:N
                                fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                                    [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                    T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(N+i))* ...
                                    (MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1);ol(i);opx(i);cpx(i);feld(i); rho_b(i); ...
                                    ACon_T(i); ACon_M(i); ACon_C(i); ACon_C_dummy(i); Con_T(i); Con_M(i); Con_C(i);Con_C_dummy(i); ...
                                    CD_M(i)*(MgO_range(2)-MgO_range(1));CD_R(i)*(MgO_range(2)-MgO_range(1));CD_R_dummy(i)*(MgO_range(2)-MgO_range(1));CD_Ms(i)*(MgO_range(2)-MgO_range(1));CD_Rs(i)*(MgO_range(2)-MgO_range(1));CD_Rs_dummy(i)*(MgO_range(2)-MgO_range(1))]);
                            end
                        

                        else
                
                            fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s \t %10s \t %10s \n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld', 'Dens');
                        
                            for i=1:1:N
                                fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
                                    [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                                    T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(N+i))* ...
                                    (MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1);ol(i);opx(i);cpx(i);feld(i); rho_b(i)]);
                            end
                        end

                    end
                
                    fclose(File_Out);
                
                
                    FilenameON = ['output_',num2str(output_counter),'_NODES.txt'];
                
                    File_OutN = fopen(FilenameON, 'w');
                
                    fprintf(File_OutN, '%4s  %10.10f %3s \n', 'Time', Time/Year/1000, 'kyr' );
                
                    fprintf(File_OutN, '%10s \t %10s \t %10s \t %10s \t %10s \n', 'Depth (km)', ...
                        'Sol. Vel.', 'Liq. Vel', 'Sol. Dens', 'Liq. Dens');
                
                    for i=1:1:N+1
                        fprintf(File_OutN, '%10.3f \t %4e \t %4e \t %10.3f \t %10.3f  \n', ...
                            [nodez(i)/1000-Base_crust; u_all(N+1+i); u_all(i); rhom(i); rhof(i)]);
                    end
                
                    fclose(File_OutN);


                
                if Output_Flag==0
                else
                    % Save the variables in matlab.
                   % FileNameO = ['output_time_',num2str(Time/Year/1000),'_kyr.mat'];
                   if mod(output_counter,restart_No)==0
                        save(['output_',num2str(output_counter),'.mat']');  %(FileNameO)   
                   end
                end
               
                 %% CAB ADDED

                    %figure(9)
                    %subplot(2,1,1)
                    %plot(Time/Year/1000, Mass_cons, 'b.')
                    %hold on
                    %xlabel('Time (ka)')
                    %ylabel('Mass conservation error (out of 1)')
                    %if (min_cons>Mass_cons)
                    %    min_cons=Mass_cons;
                    %end
                    %if (max_cons<Mass_cons)
                    %    max_cons=Mass_cons;
                    %end
                    %ylim([min_cons,max_cons])

                    %subplot(2,1,2)
                    %plot(Time/Year/1000, Mass_res, 'r.')
                    %hold on
                    %xlabel('Time (ka)')
                    %ylabel('Residual=|Comp now - Original comp|')
                    %if (min_res>Mass_res)
                     %   min_res=Mass_res;
                    %end
                    %if (max_res<Mass_res)
                     %   max_res=Mass_res;
                    %end
                    %ylim([min_res,max_res])
                    
                    fprintf(File_comp, '%10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, Mass_cons, Mass_res);
                    
                    fprintf(File_breaks, '%10.5f \t %10.5f \t %10.5f \t %10.5f \n', Time/Year/1000, break_percent, break_counter, counter);

                    output_counter=output_counter+1;
            end
        end
        end

 end
%CAB end        
%%
    % if isempty(find(H>=ESolidus(0,cp,rhof),1))
    Index1=find(phi(1:N)>0,1,'first');
    Index2=find(phi(1:N)>0,1,'last');
    if ~isempty(Index1)
        %         Active_region=max(Active_region,abs((cellz(Index1)-cellz(Index2))));
        Active_region=[min(Active_region(1),cellz(Index1)+833) max(Active_region(2),cellz(Index2)-10)];
    end
    
    if max(phi(1:N))<1e-6 && SillCount>=SillNo 
        break;
    end

    if With_monitor==1
        % if toc>Monitor_frame*Update_frequency
            Update_plot_master
            Monitor_frame=Monitor_frame+1;
        % end
    end



     
    Just_intruded=0;
end
toc;
fclose(File_echo);
fclose(File_comp);
fclose(File_breaks);
fclose(File_precision);

if Record_data==1    
 %   eval(['save(''mu' num2str(mu_m) '_cut' num2str(Contribution_cut)  'T=' num2str(round(Time/Year)) '.mat'')']);
  %  Record_index=Record_index+1;    
%end
    
        %Save the variables as .txt


        FilenameO = ['output_time_',num2str(Time/Year/1000),'_kyr_CELLS.txt'];

        File_Out = fopen(FilenameO, 'w');
    
        fprintf(File_Out, '%4s %6.4f %3s \n', 'Time', Time/Year/1000, 'kyr');

        if (HHJPet==1 || SSPD==1)

            if To_cal_mc==1
                fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                            'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                            'Cb MgO','Cs MgO', 'Cl MgO',  'Dens', 'H_MF', 'C_MF', 'R_MF','H_MF_t','C_MF_t', 'R_MF_t', 'C_CB', 'R_CB', 'C_CB_t', 'R_CB_t');
                    
                for i=1:1:N
                    fprintf(File_Out, ['%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f ' ...
                        '\t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \n'], ...
                        [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                        T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(N+i))* ...
                        (SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); rho_b(i); ...
                        ACon_T(i); ACon_M(i); ACon_C(i); Con_T(i); Con_M(i); Con_C(i); ...
                        CD_M(i)*(SiO2_range(2)-SiO2_range(1));CD_R(i)*(SiO2_range(2)-SiO2_range(1));CD_Ms(i)*(SiO2_range(2)-SiO2_range(1));CD_Rs(i)*(SiO2_range(2)-SiO2_range(1))]);
                end
            else
    
                fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb MgO','Cs MgO', 'Cl MgO',  'Dens');
                        
                for i=1:1:N
                    fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \n', ...
                        [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                        T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(N+i))* ...
                        (SiO2_range(2)-SiO2_range(1))+SiO2_range(1); (C_all(i))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1); rho_b(i)]);
                end
            end

        elseif FourMPD==1

            if To_cal_mc==1
                    fprintf(File_Out, ['%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s' ...
                        ' \t %10s \t %10s\t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s  \t %10s \t %10s\t %10s\n'], 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld', 'Dens', 'H_MF', 'C_MF', 'R_MF','R_MF_dum','H_MF_t','C_MF_t', 'R_MF_t','R_MF_t_dum', 'C_CB', 'R_CB', 'R_CB_dum', 'C_CB_t', 'R_CB_t', 'R_CB_t_dum');
                        
                for i=1:1:N
                    fprintf(File_Out, ['%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f ' ...
                        ' \t %10.4f  \t %10.4f \t %10.4f \t %10.4f ' ...
                        '\t %10.4f  \t %10.4f  \t %10.4f \t %10.4f \t %10.4f  \t %10.4f  \t %10.4f \t %10.4f  \n'], ...
                        [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                        T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(N+i))* ...
                        (MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1);...
                        ol(i);opx(i);cpx(i);feld(i); rho_b(i); ...
                        ACon_T(i); ACon_M(i); ACon_C(i); ACon_C_dummy(i); Con_T(i); Con_M(i); Con_C(i); Con_C_dummy(i);...
                        CD_M(i)*(MgO_range(2)-MgO_range(1));CD_R(i)*(MgO_range(2)-MgO_range(1));CD_R_dummy(i)*(MgO_range(2)-MgO_range(1));CD_Ms(i)*(MgO_range(2)-MgO_range(1));CD_Rs(i)*(MgO_range(2)-MgO_range(1));CD_Rs_dummy(i)*(MgO_range(2)-MgO_range(1))]);
                end

            else

                fprintf(File_Out, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s  \t %10s \t %10s\t %10s \t %10s \t %10s \t %10s \t %10s\n', 'Depth (km)', ...
                                'Melt frac.', 'Enth.', 'Temp. (C)', 'Ts (C)', 'Tl (C)', 'Cb', ...
                                'Cb MgO','Cs MgO', 'Cl MgO', 'ol', 'opx', 'cpx', 'feld', 'Dens');
                        
                for i=1:1:N
                    fprintf(File_Out, '%10.3f \t %10.4f \t %3e \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \t %10.4f \t %10.4f  \t %10.4f  \t %10.4f \t %10.4f  \n', ...
                        [cellz(i)/1000-Base_crust; phi(i); H(i); ...
                        T(i); Ts(i); Tl(i); Cb(i);  (Cb(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(N+i))* ...
                        (MgO_range(2)-MgO_range(1))+MgO_range(1); (C_all(i))*(MgO_range(2)-MgO_range(1))+MgO_range(1);ol(i);opx(i);cpx(i);feld(i); rho_b(i)]);
                end

            end

        end


        fclose(File_Out);
    
    
        FilenameON = ['output_time_',num2str(Time/Year/1000),'_kyr_NODES.txt'];
    
        File_OutN = fopen(FilenameON, 'w');
    
        fprintf(File_OutN, '%4s  %10.10f %3s \n', 'Time', Time/Year/1000, 'kyr' );
    
        fprintf(File_OutN, '%10s \t %10s \t %10s \t %10s \t %10s \n', 'Depth (km)', ...
            'Sol. Vel.', 'Liq. Vel', 'Sol. Dens', 'Liq. Dens');
    
        for i=1:1:N+1
            fprintf(File_OutN, '%10.3f \t %4e \t %4e \t %10.3f \t %10.3f  \n', ...
                [nodez(i)/1000-Base_crust; u_all(N+1+i); u_all(i); rhom(i); rhof(i)]);
        end
    
        fclose(File_OutN);



    if Output_Flag==0
    else
        % Save the variables in matlab.
        FileNameO = ['output_time_',num2str(Time/Year/1000),'_kyr.mat'];
        save(FileNameO)                
    end
end