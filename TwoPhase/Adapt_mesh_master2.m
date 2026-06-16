%% A mesh adaptivity procudure, based on melt fracton and Mass_data changes
% Haiyang Hu
% 2022-Sep-26
% both conservative and  non-conservative schemes are available
% 2023-Mar-12 modified for three phase model

% Adaption by CAB to add OG_Cb so that we can track mass conservation


%Mass_data=[m1,n1,v1,m2,n2,v2,v3];
num_adaptive=1;
to_adapt=1;



index=find(phi(1:N)>0,1,'last');
if nodez(index+1)>Top0
    Top0=nodez(index+1);
end
fine_margin=500; %meter

% Cb_sum=sum((C_all(1:N).*phi(1:N)+C_all(N+1:end).*phi(N+1:end)).*dz');
% H_sum=sum((cp*T+phi(1:N)*Lf).*dz');
while num_adaptive<=max_adaptive_number && to_adapt==1
    phi=phi(1:N);
    Cl=C_all(1:N);
    Cs=C_all(N+1:2*N);
    
    node_new=nodez;

    deleted=0;
    pass=0;
    
    ul=u_all(1:N+1);
    us=u_all(N+2:end);
    % OG_Cb_FM = project_cell2node(nodez,dz,OG_Cb);

    Change = abs(diff(phi(1:N)));

    Change = max(Change, abs(diff(Cl)));
    Change = max(Change, abs(diff(Cs)));

    if Has_volatile == 1
        Change = max(Change, abs(diff(S)));
        Change = max(Change, abs(diff(Cl2)));
        Change = max(Change, abs(diff(Cs2)));
    end


    newL=nodez(3:N+1)-nodez(1:N-1);
    if N>min_N

        for i=3:N-1 %cell number
            if pass==1
                pass=0;
                continue;
            end

           

            if Change(i-1)<min_change &&...
                    newL(i-1)<max_dx &&...
                    node_new(i-deleted)<Top0 && phi(i-deleted-1)<0.5 && phi(i-deleted)<0.5

                node_new(i-deleted)=[];

                
                ul(i-deleted)=[];
                us(i-deleted)=[];
                %conservative coarsen
                Cb_all=((phi(i-deleted-1).*Cl(i-deleted-1)+(1-phi(i-deleted-1))*Cs(i-deleted-1))*dz(i-deleted-1)+...
                       (phi(i-deleted)  .*Cl(i-deleted)  +(1-phi(i-deleted))  *Cs(i-deleted))  *dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
                
                Cb2_all=((phi(i-deleted-1).*Cl2(i-deleted-1)+(1-phi(i-deleted-1))*Cs2(i-deleted-1)+ S(i-deleted-1))*dz(i-deleted-1)+...
                         (phi(i-deleted)  .*Cl2(i-deleted)  +(1-phi(i-deleted))  *Cs2(i-deleted)  + S(i-deleted))  *dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                H_all=((cp*T(i-deleted-1)+phi(i-deleted-1)*Lf)*dz(i-deleted-1)+(cp*T(i-deleted)+phi(i-deleted)*Lf)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                phi(i-deleted-1)=(phi(i-deleted-1)+phi(i-deleted))/2;
                T(i-deleted-1)=(H_all-phi(i-deleted-1)*Lf)/cp;
                
                % if i-deleted-1==161
                %     aaa=1
                % end
%     Pressure=(nodez(end)-(node_new(i-deleted-1)+node_new(i-deleted))/2)*g*rho_mean/1e5+1; %in bar
%     index_pressure_nts=max(min(floor(Pressure/8000*N_Ts)+1,N_Ts),1);    
%     index_pressure_ntl=max(min(floor(Pressure/40e3* N_Tl) + 1, N_Tl),1);
%                 index_cb2_nts=max(min(floor(Cb2_all/0.13/Par_v * N_Ts) + 1, N_Ts),1);
%                 index_cb2_ntl=max(min(floor(Cb2_all/0.2/Par_v * N_Tl) + 1, N_Tl),1);
% 
%             lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts, index_cb2_nts);
%             lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl, index_cb2_ntl);
% 
%             coef_ts = Ts0_coefficient(lin_ts, :);   % (nI × 4)
%             coef_tl = Tl0_coefficient(lin_tl, :);   % (nI × 4)
% 
%             % unpack
%             a_ts = coef_ts(:,1); b_ts = coef_ts(:,2);
%             c_ts = coef_ts(:,3); d_ts = coef_ts(:,4);
% 
%             a_tl = coef_tl(:,1); b_tl = coef_tl(:,2);
%             c_tl = coef_tl(:,3); d_tl = coef_tl(:,4);
% 
% 
% 
%                     ts= a_ts*Pressure/8000+b_ts*Cb2_all/0.13/Par_v+c_ts*Pressure/8000*Cb2_all/0.13/Par_v+d_ts;
%                     tl= a_tl*Pressure/40e3+b_tl*Cb2_all/0.2/Par_v +c_tl*Pressure/40e3*Cb2_all/0.2/Par_v +d_tl;
%                     C1= tl;
%                     B1= ts-A1-tl;
% eps=0.1;
% eps2=1e-6;
% 
%                     MINT=((T(i-deleted-1)+C1+50)-sqrt((T(i-deleted-1)-C1-50)^2+eps))/2;
%                     Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
%                     SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps2))/2;
%                     Cl(i-deleted-1)=(Cb_all+SMIN+sqrt((SMIN-Cb_all)^2+eps2))/2;
                    
                    Cl(i-deleted-1)=(Cl(i-deleted-1)+Cl(i-deleted))/2;
                    Cs(i-deleted-1)=(Cb_all-phi(i-deleted-1)*Cl(i-deleted-1))/max((1-phi(i-deleted-1)),1e-6);

                % end
                
                dz(i-deleted)=[]; 
                phi(i-deleted)=[];   %Careful! deleting item needs to be done in one-go
                T(i-deleted)=[];
                Cs(i-deleted)=[];
                Cl(i-deleted)=[];
                if Has_volatile==1
                    S(i-deleted-1)=(S(i-deleted-1)+S(i-deleted))/2;

                    Cl2(i-deleted-1)=(Cb2_all-S(i-deleted-1))/(phi(i-deleted-1)+Par_v*(1-phi(i-deleted-1)));
                    Cs2(i-deleted-1)=Cl2(i-deleted-1)*Par_v;

                    
                    Cs2(i-deleted)=[];
                    Cl2(i-deleted)=[];
                    S(i-deleted)=[];
                    % Lsaturation(i-deleted)=[];
                    % Ssaturation(i-deleted)=[];
                end
                Ts(i-deleted)=[];
                Tl(i-deleted)=[];

                % OG_Cb(i-deleted-1,:)=(OG_Cb(i-deleted-1,:)*dz(i-deleted-1)+OG_Cb(i-deleted,:)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
                % OG_Cb(i-deleted,:)=0;                    

                % if To_recored_intrusion_marker==1
                %     Intrusion_marker_static(i-deleted)=[];
                %     Intrusion_marker_dynamics(i-deleted)=[];
                % end

                % rho(i-deleted,:)=[];





                % Pg_real(i-deleted)=[];
                % cp(i-deleted,:)=[];

                deleted=deleted+1;
                pass=1;
            end
        end
    end

    nodez=node_new;
    if N<max_N
        dz=nodez(2:end)-nodez(1:end-1);
        dphi=abs(phi(2:end)-phi(1:end-1));

        to_refine=zeros(length(phi)-1,1);
        index=find(phi>1e-2,1,'last');
        if ~isempty(index)
            index2=find(nodez<nodez(index)+fine_margin,1,'last');

            if index2-index<fine_margin/min_dx/2
                to_refine(index:index2)=1;
            end
        end
        to_refine=((dphi>max_change)| to_refine) & dz(1:end-1)'>min_dx*2  & dz(2:end)'>min_dx*2;
        index=find(to_refine);

        if ~isempty(index)
            out = index(1); % always keep the first element
            for i = 2:length(index)
                if index(i) ~= out(end) + 1
                    out = [out; index(i)];
                end
            end
            to_refine=out;
        else
            to_refine=[];
        end

        offset = 0; 

        phi_new=phi;
        T_new=T;
        Cl_new=Cl;
        Cs_new=Cs;
        if Has_volatile==1
            S_new=S;
            Cl2_new=Cl2;
            Cs2_new=Cs2;
        end

        for k=1:numel(to_refine)
            i = to_refine(k) + offset;

            Cb_all=((phi_new(i).*Cl(i)    +(1-phi(i))    *Cs(i))  *dz(i)+...
                  (phi_new(i+1).*Cl(i+1)  +(1-phi(i+1))  *Cs(i+1))*dz(i+1))/(dz(i)+dz(i+1));

            Cb2_all=((phi_new(i).*Cl2(i)    +(1-phi(i))    *Cs2(i))  *dz(i)+...
                  (phi_new(i+1).*Cl2(i+1)  +(1-phi(i+1))  *Cs2(i+1))*dz(i+1))/(dz(i)+dz(i+1));

            H_all=((cp*T(i)+phi(i)*Lf)*dz(i)+(cp*T(i+1)+phi(i+1)*Lf)*dz(i+1))/(dz(i)+dz(i+1));

            

            phi_gap=sum(phi_new(i:i+1))/2;
            phi_new = [phi_new(1:i); phi_gap; phi_new(i+1:end)];

            T_gap=(H_all-phi_gap*Lf)/cp;
            T_new = [T_new(1:i); T_gap; T_new(i+1:end)];

%     Pressure=(sum(node_new(i-1:i+1))+nodez_new(i))/4*g*rho_mean/1e5+1; %in bar
%     index_pressure_nts=max(min(floor(Pressure/8000*N_Ts)+1,N_Ts),1);    
%     index_pressure_ntl=max(min(floor(Pressure/40e3* N_Tl) + 1, N_Tl),1);
%                 index_cb2_nts=max(min(floor(Cb2_all/0.13/Par_v * N_Ts) + 1, N_Ts),1);
%                 index_cb2_ntl=max(min(floor(Cb2_all/0.2/Par_v * N_Tl) + 1, N_Tl),1);
% 
%             lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts, index_cb2_nts);
%             lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl, index_cb2_ntl);
% 
%             coef_ts = Ts0_coefficient(lin_ts, :);   % (nI × 4)
%             coef_tl = Tl0_coefficient(lin_tl, :);   % (nI × 4)
% 
%             % unpack
%             a_ts = coef_ts(:,1); b_ts = coef_ts(:,2);
%             c_ts = coef_ts(:,3); d_ts = coef_ts(:,4);
% 
%             a_tl = coef_tl(:,1); b_tl = coef_tl(:,2);
%             c_tl = coef_tl(:,3); d_tl = coef_tl(:,4);
% 
% 
% 
%                     ts= a_ts*Pressure/8000+b_ts*Cb2_all/0.13/Par_v+c_ts*Pressure/8000*Cb2_all/0.13/Par_v+d_ts;
%                     tl= a_tl*Pressure/40e3+b_tl*Cb2_all/0.2/Par_v +c_tl*Pressure/40e3*Cb2_all/0.2/Par_v +d_tl;
%                     C1= tl;
%                     B1= ts-A1-tl;
% eps=0.1;
% eps2=1e-6;
% 
%                     MINT=((T(i-deleted-1)+C1+50)-sqrt((T(i-deleted-1)-C1-50)^2+eps))/2;
%                     Cond5=(-B1-sqrt(B1^2-4*A1*(C1-MINT)))/2/A1;
%                     SMIN=(Cond5+1-sqrt((Cond5-1)^2+eps2))/2;
%                     Cl(i-deleted-1)=(Cb_all+SMIN+sqrt((SMIN-Cb_all)^2+eps2))/2;
% 
%                     Cl(i-deleted-1)=(Cl(i-deleted-1)+Cl(i-deleted))/2;
%                     Cs(i-deleted-1)=(Cb_all-phi(i-deleted-1)*Cl(i-deleted-1))/(1-phi(i-deleted-1));
%+++++
            Cl_gap=(Cl_new(i)+Cl_new(i+1))/2;
            Cl_new = [Cl_new(1:i); Cl_gap; Cl_new(i+1:end)];
            Cs_gap=(Cb_all-phi_gap*Cl_gap)/max((1-phi_gap),1e-6);
            Cs_new = [Cs_new(1:i); Cs_gap; Cs_new(i+1:end)];


            ul=[ul(1:i); ul(i+1);  ul(i+1:end) ];
            us=[us(1:i); us(i+1);  us(i+1:end) ];
            if Has_volatile==1
                S_gap=mean(S_new(i:i+1));
                S_new = [S_new(1:i); S_gap; S_new(i+1:end)];

                Cl2_gap=(Cb2_all-S_gap)/(phi_gap+(1-phi_gap)*Par_v);
                Cl2_new = [Cl2_new(1:i); Cl2_gap; Cl2_new(i+1:end)];
                Cs2_gap=Cl2_gap*Par_v;
                Cs2_new = [Cs2_new(1:i); Cs2_gap; Cs2_new(i+1:end)];
            end
            offset=offset+1;
        end


        node_new = nodez;
        shift = 0;

        for k = 1:numel(to_refine)
            i = to_refine(k) + shift+1;

            left_avg  = (node_new(i-1) + node_new(i)) / 2;
            right_avg = (node_new(i) + node_new(i+1)) / 2;

            node_new = [node_new(1:i-1), left_avg, right_avg, node_new(i+1:end)];
            shift = shift + 1;   
        end

        phi=phi_new;
        Cl=Cl_new;
        Cs=Cs_new;
        T=T_new;
        if Has_volatile==1
            S=S_new;
            Cl2=Cl2_new;
            Cs2=Cs2_new;
        end
    end
    cellz_new=(node_new(1:end-1)+node_new(2:end))/2;  % cell center points

    N_old=N;
    N=length(node_new)-1;
    nodez=node_new;
    cellz=(nodez(1:end-1)+nodez(2:end))/2;  % cell center points
    dz=nodez(2:end)-nodez(1:end-1);

    disp(['Number of nodes: ', num2str(N)])
    if abs(N_old-N)/N_old<0.05
        to_adapt=0;
    end
    num_adaptive=num_adaptive+1;

    phi=[phi; 1-phi];
    C_all=[Cl;Cs];
    u_all=[ul;us];

    H=phi(1:N)*Lf+cp*T;
    Cb=phi(1:N).*Cl+(1-phi(1:N)).*Cs;
end



% Cb_sum2=sum((Cl.*phi(1:N)+Cs.*(1-phi(1:N))).*dz');
% H_sum2=sum((cp*T+phi(1:N)*Lf).*dz');

% disp(['component conservation:' num2str((Cb_sum2-Cb_sum)/Cb_sum)]);
% disp(['enthalpy conservation:' num2str((H_sum2-H_sum)/H_sum)]);


phi_old=phi;
u_all_old=u_all;
H_old=H;
T_old=T;
C_all_old=C_all;
Cb_old=Cb;
if Has_volatile==1
    S_old=S;
    Cs2_old=Cs2;
    Cl2_old=Cl2;
end
Last_adapted=Time;


Advance_time=0;
Newton_solver3;