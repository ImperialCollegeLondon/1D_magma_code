%% A mesh adaptivity procudure, based on melt fracton and Mass_data changes
% Haiyang Hu
% 2022-Sep-26
% both conservative and  non-conservative schemes are available
% 2023-Mar-12 modified for three phase model

% Adaption by CAB to add OG_Cb so that we can track mass conservation


%Mass_data=[m1,n1,v1,m2,n2,v2,v3];
num_adaptive=0;
to_adapt=1;
Mass_data_FM=zeros(N+1, 7);
OG_Cb_FM = zeros(N+1);

Mass_data_to_adapt=[1,2,3,4,5,6,7]; 

index=find(phi>0,1,'last');
if nodez(index+1)>Top0
    Top0=nodez(index+1);
end

while num_adaptive<=max_adaptive_number && to_adapt==1

    node_new=nodez;

    % cellz_save=cellz;
    % nodez_save=nodez;
    % u_save=u_all;
    % phi_save=phi;
    % T_save=T;
    % H_save=H;
    % Cb_save=Cb;
    % C_all_save=C_all;


    deleted=0;
    pass=0;

    phi_FM=project_cell2node(nodez,dz,phi(1:N));
    S_FM=project_cell2node(nodez,dz,S(1:N));
    Mass_data_FM=zeros(N+1,7);
    for i=Mass_data_to_adapt
        Mass_data_FM(:,i)=project_cell2node(nodez,dz,Mass_data(:,i));
    end
    OG_Cb_FM = project_cell2node(nodez,dz,OG_Cb);

    if N>min_N
        Top=find(phi>1e-3,1,'last');
        Bottom=find(phi>1e-3,1,'first');
        Top=cellz(Top);
        Bottom=cellz(Bottom);
        Reserve_depth=Bottom+(Top-Bottom)*0.85;
        for i=3:N-1
            if pass==1
                pass=0;
                continue;
            end

            newL=nodez(i+1)-nodez(i-1);
            Changes=Mass_data_FM(i+1, Mass_data_to_adapt)-Mass_data_FM(i-1, Mass_data_to_adapt);
            
            Changes=Changes/dMNV_min*0.05;
            if To_recored_intrusion_marker==1
                Changes2=[Changes Intrusion_marker_static(i-deleted)-Intrusion_marker_static(i-deleted-1) Intrusion_marker_dynamics(i-deleted)-Intrusion_marker_dynamics(i-deleted-1)];
            else
                Changes2=Changes;
            end

            if abs(cellz(i)< Reserve_depth &&...
                    phi_FM(i+1)-phi_FM(i-1))<dphi_min &&...
                    (S_FM(i+1)-S_FM(i-1))<dS_min &&...
                    prod(Changes2<0.05) &&...
                    newL<max_dx &&...
                    newL<(node_new(i+2-deleted)-node_new(i+1-deleted))*aspect_ratio &&...
                    newL<(node_new(i-1-deleted)-node_new(i-2-deleted))*aspect_ratio &&...
                    node_new(i-deleted)<Top0
               

                node_new(i-deleted)=[];
                %                 u_all([i-deleted,i-deleted*2+N, i-deleted*3+N*2])=[];


                H(i-deleted-1)=(H(i-deleted-1)*dz(i-deleted-1)+H(i-deleted)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
                H(i-deleted)=[];

                Mass_data(i-deleted-1,:)=(Mass_data(i-deleted-1,:)*dz(i-deleted-1)+Mass_data(i-deleted,:)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                Mass_data_addition(i-deleted-1,:)=(Mass_data_addition(i-deleted-1,:)*dz(i-deleted-1)+Mass_data_addition(i-deleted,:)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                OG_Cb(i-deleted-1,:)=(OG_Cb(i-deleted-1,:)*dz(i-deleted-1)+OG_Cb(i-deleted,:)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                %                 [phi(i-deleted-1),T(i-deleted-1),C_all(i-deleted-1), C_all(i-deleted*2-1+N), ~]=poro_component_solve(Cb(i-deleted-1), H(i-deleted-1),A1, B1, C1, A2, B2, C2, ae, Lf, rhof, cp,Precision);
                %                 Fix_TS=0;
                if is_solid_solution==0
                    [phi(i-deleted-1),S(i-deleted-1),T(i-deleted-1),rho(i-deleted-1,:),~,T_S_region(i-deleted-1,:),Ts_new(i-deleted-1),Mass_data(i-deleted-1,:),S_cap(i-deleted-1)]                                =Phase_component_updated2      (H(i-deleted-1),sum(Mass_data(i-deleted-1,[1,4]),2),sum(Mass_data(i-deleted-1,[2,5]),2), sum(Mass_data(i-deleted-1,[3,6,7]),2),Pg_real(i-deleted-1), T(i-deleted-1), [1,3], Jacobians, Rhs, Constant_index, Sys_constant,cp(i-deleted-1,:),   rho_constant, Lf, K0, Ts, Tl, min_por, dz(i-deleted-1),Data_point,Data_y,Data_point2,Data_y2,Precision,Conservation_type);
                else
                    if N_component==3
                        if Add_CLCU==1
                            input=[H(i-deleted-1),sum(Mass_data(i-deleted-1,[1,4]),2),sum(Mass_data(i-deleted-1,[2,5]),2), sum(Mass_data(i-deleted-1,[3,6,7]),2), sum(Mass_data(i-deleted-1,[8,9,10]),2),sum(Mass_data(i-deleted-1,[11,12,13]),2)];
                        else
                            input=[H(i-deleted-1),sum(Mass_data(i-deleted-1,[1,4]),2),sum(Mass_data(i-deleted-1,[2,5]),2), sum(Mass_data(i-deleted-1,[3,6,7]),2)];
                        end
                        [phi(i-deleted-1),S(i-deleted-1),T(i-deleted-1),rho(i-deleted-1,:),~,T_S_region(i-deleted-1,:),Ts_new(i-deleted-1),Ts_local(i-deleted-1),Mass_data(i-deleted-1,:),S_cap(i-deleted-1),Tl_local(i-deleted-1),~,~,Partition_CL_CU(i-deleted-1,:)]                               =Phase_component_updated_solid3         (input,Pg_real(i-deleted-1), T(i-deleted-1), [1,3], Jacobians, Rhs, Extra_func, Constant_index, Sys_constant,cp(i-deleted-1,:),   rho_constant, Lf, [K0 Kcl Kcu], Ts, Tl, min_por, max_melt, dz(1),Data_point,Data_y,Data_point2,Data_y2,PD_range,Precision, Conservation_type,simplified_TS, 0);
                    else
                        [phi(i-deleted-1),S(i-deleted-1),T(i-deleted-1),rho(i-deleted-1,:),~,T_S_region(i-deleted-1),Ts_local(i-deleted-1),Tl_local(i-deleted-1),Mass_data(i-deleted-1,:),S_cap(i-deleted-1),Partition_CL_CU(i-deleted-1,:)]=Phase_component_updated_solid_5p_simple(input, Pg_real(i-deleted-1), T(i-deleted-1), 1, Jacobians, Rhs, Constant_index, Extra_func, Sys_constant, cp(i-deleted-1,:),   rho_constant, Lf, [K0 Kcl Kcu Partition_CL_CU(i-deleted-2,:)], Ts, Tl, min_por, dz, Data_point,Data_y,Data_point2,Data_y2, PD_range, Precision, Conservation_type, 1);
                    end
                end

                phi(i-deleted)=[];   %Careful! delete item needs to be done in one-go
                S(i-deleted)=[];
                T(i-deleted)=[];
                    
                % Invis_vol(i-deleted)=[];

                if To_recored_intrusion_marker==1
                    Intrusion_marker_static(i-deleted)=[];
                    Intrusion_marker_dynamics(i-deleted)=[];
                end

                rho(i-deleted,:)=[];

                % Ts_new(i-deleted)=[];
                if is_solid_solution==1
                    Ts_local(i-deleted)=[];
                    S_cap(i-deleted)=[];
                else
                    % k2(i-deleted)=[];
                    S_cap(i-deleted)=[];
                end
                Mass_data(i-deleted,:)=[];
                Mass_data_addition(i-deleted,:)=[];
                OG_Cb(i-deleted,:)=[];
                % if N_component==5 || Add_CLCU==1
                %     Partition_CL_CU(i-deleted,:)=[];
                % end

                Pg_real(i-deleted)=[];
                cp(i-deleted,:)=[];

                deleted=deleted+1;
                pass=1;
            end
        end
    end


    if N<max_N
        added=0;
        temp=length(node_new);
        dz=node_new(2:end)-node_new(1:end-1);
        phi_FM=project_cell2node(node_new,dz,phi(1:temp-1));
%         S_FM=project_cell2node(nodez,dz,S(1:N));
%         for i=Mass_data_to_adapt
%             Mass_data_FM(:,i)=project_cell2node(nodez,dz,Mass_data(:,i));
%         end
        nodez_temp=node_new;
        %     H_FM=project_cell2node(node_new,dz,H);
        %     Cb_FM=project_cell2node(node_new,dz,Cb);
        %     C_all_FM=zeros(temp*2,1);
        %     C_all_FM(1:temp)=project_cell2node(node_new,dz,C_all(1:temp-1));
        %     C_all_FM(temp+1:end)=project_cell2node(node_new,dz,C_all(temp:end));
        %     phi_new=phi_FM;

        phi_cut=phi(1:temp-1);
        S_cut=S(1:temp-1);
        for i=3:temp-1
            olddphi=abs(phi_FM(i)-phi_FM(i-1));
            oldL=nodez_temp(i)-nodez_temp(i-1);
            % if (olddphi>dphi_max) && oldL>min_dx*2     %currently refining only considers phi
            if (olddphi >  dphi_max || (node_new(i+added)>Top0+100 && node_new(i+added-1)<Top0+100)|| (node_new(i+added)>Top0 && node_new(i+added)<Top0+100))      && (oldL > min_dx * 2)  
                
                temp=ceil(olddphi/dphi_max);
                if temp==0
                    insert_node=5;
                else
                    insert_node=min(temp,10);
                end
                % insert_node=min(ceil(olddphi/dphi_max), floor(oldL/min_dx));

                node_new=[node_new(1:i-2+added); linspace(node_new(i-1+added),node_new(i+added), insert_node+1)'; node_new(i+1+added:end)];
%                 u_all=[u_all(1:i-2+added);                          linspace(u_all(i-1+added),u_all(i+added), insert_node+1)';                    u_all(i+1+added: temp+added);...
%                        u_all(temp+1+added:temp+i-2+added*2);        linspace(u_all(temp+i-1+added*2),   u_all(temp+i+added*2), insert_node+1)';   u_all(temp+  i+1+added*2:temp*2+added*2);...
%                        u_all(temp*2+added*2+1: temp*2+added*3+i-2); linspace(u_all(temp*2+added*3+i-1), u_all(temp*2+added*3+i), insert_node+1)'; u_all(temp*2+i+1+added*3:end)];

                % conservative scheme,
                %             grad=(H(i+added)-H(i-2+added))/(dz(i-1+added)+dz(i+added)/2+dz(i-2+added)/2);
                %             X=(H(i-1+added)*insert_node-grad*insert_node*(insert_node-1)/2)/insert_node;
                %             H=[H(1:i-2+added);  linspace(X,X+(insert_node-1)*grad, insert_node)';   H(i+added:end)];
                H=[H(1:i-2+added);  ones(insert_node,1)*H(i-1+added);   H(i+added:end)];
                T=[T(1:i-2+added);  ones(insert_node,1)*T(i-1+added);   T(i+added:end)]; % this is just the gueseed T, need to recalculated with phase diagram
                % Invis_vol=[Invis_vol(1:i-2+added); ones(insert_node,1)*Invis_vol(i-1+added); Invis_vol(i+added:end)];

                if To_recored_intrusion_marker==1
                    Intrusion_marker_static=[Intrusion_marker_static(1:i-2+added);  ones(insert_node,1)*Intrusion_marker_static(i-1+added);   Intrusion_marker_static(i+added:end)];
                    Intrusion_marker_dynamics=[Intrusion_marker_dynamics(1:i-2+added);  ones(insert_node,1)*Intrusion_marker_dynamics(i-1+added);   Intrusion_marker_dynamics(i+added:end)];
                end
                
                %             grad=(Cb(i+added)-Cb(i-2+added))/(dz(i-1+added)+dz(i+added)/2+dz(i-2+added)/2);
                %             X=(Cb(i-1+added)*insert_node-grad*insert_node*(insert_node-1)/2)/insert_node;
                %             Cb=[Cb(1:i-2+added);  linspace(X,X+(insert_node-1)*grad, insert_node)';   Cb(i+added:end)];
                Mass_data=[Mass_data(1:i-2+added,:);  ones(insert_node,1)*Mass_data(i-1+added,:);   Mass_data(i+added:end,:)];
                Mass_data_addition=[Mass_data_addition(1:i-2+added,:);  ones(insert_node,1)*Mass_data_addition(i-1+added,:);   Mass_data_addition(i+added:end,:)];
                OG_Cb=[OG_Cb(1:i-2+added,:);  ones(insert_node,1)*OG_Cb(i-1+added,:);   OG_Cb(i+added:end,:)];
                if N_component==5 || Add_CLCU==1
                    % Partition_CL_CU=[Partition_CL_CU(1:i-2+added,:); ones(insert_node,1)*Partition_CL_CU(i-1+added,:); Partition_CL_CU(i+added:end,:)];
                    % Partition_CL_CU_temp=zeros(insert_node,2);
                    Mass_data_temp=zeros(insert_node,13);
                else
                    Mass_data_temp=zeros(insert_node,7);
                end
                Pg_real=[Pg_real(1:i-2+added);  ones(insert_node,1)*Pg_real(i-1+added);   Pg_real(i+added:end)];
                cp=[cp(1:i-2+added,:);  ones(insert_node,1)*cp(i-1+added,:);   cp(i+added:end,:)];
               
                
                phi_cut=[phi_cut(1:i-2+added);  phi_cut(i-1+added)*ones(insert_node,1);   phi_cut(i+added:end)];
                S_cut=[S_cut(1:i-2+added);  S_cut(i-1+added)*ones(insert_node,1);   S_cut(i+added:end)];
                rho=[rho(1:i-2+added,:); rho(i-1+added,:).*ones(insert_node,1); rho(i+added:end,:)];

                % Ts_new=[Ts_new(1:i-2+added); Ts_new(i-1+added)*ones(insert_node,1); Ts_new(i+added:end)];
                

                if is_solid_solution==1
                    Ts_local=[Ts_local(1:i-2+added); Ts_local(i-1+added)*ones(insert_node,1); Ts_local(i+added:end)];
                    % Tl_local=[Tl_local(1:i-2+added); Tl_local(i-1+added)*ones(insert_node,1); Tl_local(i+added:end)];

                    % T_S_region=[T_S_region(1:i-2+added,:); T_S_region(i-1+added,:).*ones(insert_node,1); T_S_region(i+added:end,:)];
                    S_cap=[S_cap(1:i-2+added); S_cap(i-1+added)*ones(insert_node,1); S_cap(i+added:end)];
                else
                    S_cap=[S_cap(1:i-2+added); S_cap(i-1+added)*ones(insert_node,1); S_cap(i+added:end)];
                end

                temp_N=length(Ts_new);
                K=K0*ones(temp_N,1);
                added=added+insert_node-1;
            end
        end
        
        

        phi=phi_cut;
        S=S_cut;
    end

    





    cellz_new=(node_new(1:end-1)+node_new(2:end))/2;  % cell center points

    N_old=N;
    N=length(node_new)-1;
    nodez=node_new;
    cellz=cellz_new;  % cell center points
    dz=nodez(2:end)-nodez(1:end-1);

    % if N>max_N
    %     min_dx=min_dx*1.1;
    %     max_dx=max_dx*1.1;
    % end

    % if N<min_N
    %     min_dx=min_dx*0.9;
    % %     max_dx=max_dx*0.9;
    % end
    disp(['Number of nodes: ', num2str(N)])
    if abs(N_old-N)/N_old<0.05
        to_adapt=0;
    end
    num_adaptive=num_adaptive+1;
end
% Plot_settings;