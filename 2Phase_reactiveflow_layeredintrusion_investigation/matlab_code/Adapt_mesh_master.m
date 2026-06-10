%% A mesh adaptivity procudure, based on melt fracton and bulk composition changes
% Haiyang Hu
% 2022-Sep-26
% both conservative and  non-conservative schemes are available

num_adaptive=0;
to_adapt=1;
while num_adaptive<=max_adaptive_number && to_adapt==1
phi_FM=project_cell2node(nodez,dz,phi(1:N));
Cb_FM=project_cell2node(nodez,dz,Cb(1:N));
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
if N>min_N
    for i=3:N-1
        if pass==1
            pass=0;
            continue;
        end
          
        newL=nodez(i+1)-nodez(i-1);
        if abs(phi_FM(i+1)-phi_FM(i-1))<dphi_min && abs(Cb_FM(i+1)-Cb_FM(i-1))<dCb_min &&newL<max_dx && newL<(node_new(i+2-deleted)-node_new(i+1-deleted))*aspect_ratio && newL<(node_new(i-1-deleted)-node_new(i-2-deleted))*aspect_ratio
            node_new(i-deleted)=[];
            u_all([i-deleted,i-deleted*2+N])=[];


            H(i-deleted-1)=(H(i-deleted-1)*dz(i-deleted-1)+H(i-deleted)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
            H(i-deleted)=[];
            Cb(i-deleted-1)=(Cb(i-deleted-1)*dz(i-deleted-1)+Cb(i-deleted)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
            Cb(i-deleted)=[];
            
            [phi(i-deleted-1),T(i-deleted-1),C_all(i-deleted-1), C_all(i-deleted*2-1+N), ~]=poro_component_solve(Cb(i-deleted-1), H(i-deleted-1),A1, B1, C1, A2, B2, C2, ae, Lf, rhof, cp,Precision);
            phi(i-deleted*2-1+N)=1-phi(i-deleted-1);

            phi([i-deleted, i-deleted*2+N])=[];   %Careful! delete item needs to be done in one-go
            T(i-deleted)=[];
            C_all([i-deleted, i-deleted*2+N] )=[];

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
    nodez_temp=node_new;
%     H_FM=project_cell2node(node_new,dz,H);
%     Cb_FM=project_cell2node(node_new,dz,Cb);
%     C_all_FM=zeros(temp*2,1);
%     C_all_FM(1:temp)=project_cell2node(node_new,dz,C_all(1:temp-1));
%     C_all_FM(temp+1:end)=project_cell2node(node_new,dz,C_all(temp:end));
%     phi_new=phi_FM;
    Cl=C_all(1:temp-1);
    Cs=C_all(temp:end);
    phi_cut=phi(1:temp-1);
    for i=3:temp-1
        olddphi=abs(phi_FM(i)-phi_FM(i-1));
        oldL=nodez_temp(i)-nodez_temp(i-1);
        if olddphi>dphi_max && oldL>min_dx*2
            insert_node=min(ceil(olddphi/dphi_max), floor(oldL/min_dx));

            node_new=[node_new(1:i-2+added) linspace(node_new(i-1+added),node_new(i+added), insert_node+1) node_new(i+1+added:end)];
            u_all=[u_all(1:i-2+added); linspace(u_all(i-1+added),u_all(i+added), insert_node+1)'; u_all(i+1+added: temp+added); u_all(temp+1+added:temp+i-2+added*2); linspace(u_all(temp+i-1+added*2),u_all(temp+i+added*2)', insert_node+1)'; u_all(temp+i+1+added*2:end)];

% conservative scheme,   
%             grad=(H(i+added)-H(i-2+added))/(dz(i-1+added)+dz(i+added)/2+dz(i-2+added)/2);
%             X=(H(i-1+added)*insert_node-grad*insert_node*(insert_node-1)/2)/insert_node;
%             H=[H(1:i-2+added);  linspace(X,X+(insert_node-1)*grad, insert_node)';   H(i+added:end)];
              H=[H(1:i-2+added);  ones(insert_node,1)*H(i-1+added);   H(i+added:end)];
%             grad=(Cb(i+added)-Cb(i-2+added))/(dz(i-1+added)+dz(i+added)/2+dz(i-2+added)/2);
%             X=(Cb(i-1+added)*insert_node-grad*insert_node*(insert_node-1)/2)/insert_node;
%             Cb=[Cb(1:i-2+added);  linspace(X,X+(insert_node-1)*grad, insert_node)';   Cb(i+added:end)];
              Cb=[Cb(1:i-2+added);  ones(insert_node,1)*Cb(i-1+added);   Cb(i+added:end)];
% non-conservative scheme 
%             h1=(H(i-2+added)*dz(i-1+added)+H(i-1+added)*dz(i-2+added))/(dz(i-1+added)+dz(i-2+added));   
%             h2=(H(i-1+added)*dz(i+added)+H(i+added)*dz(i-1+added))/(dz(i+added)+dz(i-1+added)); 
%             grad=(h2-h1)/dz(i-1+added);
%             H=[H(1:i-2+added);  linspace(h1+grad*dz(i-1+added)/insert_node*0.5,h2-grad*dz(i-1+added)/insert_node*0.5, insert_node)'; H(i+added:end)]; 
% 
%             h1=(Cb(i-2+added)*dz(i-1+added)+Cb(i-1+added)*dz(i-2+added))/(dz(i-1+added)+dz(i-2+added));   
%             h2=(Cb(i-1+added)*dz(i+added)+Cb(i+added)*dz(i-1+added))/(dz(i+added)+dz(i-1+added)); 
%             grad=(h2-h1)/dz(i-1+added);
%             Cb=[Cb(1:i-2+added);  linspace(h1+grad*dz(i-1+added)/insert_node*0.5,h2-grad*dz(i-1+added)/insert_node*0.5, insert_node)'; Cb(i+added:end)]; 
            phi_temp=zeros(insert_node,1);
            C_all_temp=zeros(insert_node,2);
            T_temp=zeros(insert_node,1);
            for j=1: insert_node
                [phi_temp(j),T_temp(j), C_all_temp(j,1), C_all_temp(j,2),  ~]=poro_component_solve(Cb(i-2+added+j), H(i-2+added+j),A1, B1, C1, A2, B2, C2, ae, Lf, rhof, cp,Precision);
            end
            phi_cut=[phi_cut(1:i-2+added);  phi_temp;   phi_cut(i+added:end)];
            Cl=[Cl(1:i-2+added);  C_all_temp(:,1);   Cl(i+added:end)];
            Cs=[Cs(1:i-2+added);  C_all_temp(:,2);   Cs(i+added:end)];
%             C_all=[C_all(1:i-2+added);  C_all_temp(:,1);   C_all(i+added:temp+added); C_all(temp+1+added:temp+i-2+added*2);  C_all_temp(:,2);   C_all(temp+i+added*2:end)];
            T=[T(1:i-2+added); T_temp;   T(i+added:end)];
            added=added+insert_node-1;
        end
    end
    C_all=[Cl;Cs];
    phi=[phi_cut;1-phi_cut];
end
Cphi_all=C_all.*phi;

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