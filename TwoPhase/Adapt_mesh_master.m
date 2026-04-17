%% A mesh adaptivity procudure, based on melt fracton and Mass_data changes
% Haiyang Hu
% 2022-Sep-26
% both conservative and  non-conservative schemes are available
% 2023-Mar-12 modified for three phase model

% Adaption by CAB to add OG_Cb so that we can track mass conservation


%Mass_data=[m1,n1,v1,m2,n2,v2,v3];
num_adaptive=0;
to_adapt=1;

index=find(phi(1:N)>0,1,'last');
if nodez(index+1)>Top0
    Top0=nodez(index+1);
end


while num_adaptive<=max_adaptive_number && to_adapt==1
    phi=phi(1:N);
    Cl=C_all(1:N);
    Cs=C_all(N+1:2*N);

    node_new=nodez;

    deleted=0;
    pass=0;

    % OG_Cb_FM = project_cell2node(nodez,dz,OG_Cb);

    Change = abs(diff(phi(1:N)));

    Change = max(Change, abs(diff(Cl)));
    Change = max(Change, abs(diff(Cs)));

    if Has_volatile == 1
        Change = max(Change, abs(diff(S)));
        Change = max(Change, abs(diff(Cl2)));
        Change = max(Change, abs(diff(Cs2)));
    end


    newL=nodez(4:N)-nodez(2:N-2);
    if N>min_N

        for i=3:N-1
            if pass==1
                pass=0;
                continue;
            end

           

            if Change(i-deleted)<min_change &&...
                    newL(i-deleted)<max_dx &&...
                    node_new(i-deleted)<Top0

                node_new(i-deleted)=[];

                u_all([i-deleted,i-deleted*2+N])=[];

                %conservative coarsen
                Cb_all=((phi(i-deleted-1).*Cl(i-deleted-1)+(1-phi(i-deleted-1))*Cs(i-deleted-1))*dz(i-deleted-1)+...
                       (phi(i-deleted)  .*Cl(i-deleted)  +(1-phi(i-deleted))  *Cs(i-deleted))  *dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                H_all=((cp*T(i-deleted-1)+phi(i-deleted-1)*Lf)*dz(i-deleted-1)+(cp*T(i-deleted)+phi(i-deleted)*Lf)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));

                phi(i-deleted-1)=(phi(i-deleted-1)+phi(i-deleted))/2;
                T(i-deleted-1)=(H_all-phi(i-deleted-1)*Lf)/cp;

                if phi(i-deleted-1)<1e-8
                    Cs(i-deleted-1)=Cb_all;
                    Cl(i-deleted-1)=1;
                else
                    Cs(i-deleted-1)=(Cs(i-deleted-1)+Cs(i-deleted))/2;
                    Cl(i-deleted-1)=(Cb_all-(1-phi(i-deleted-1))*Cs(i-deleted-1))/phi(i-deleted-1);
                end
                
                phi(i-deleted)=[];   %Careful! delete item needs to be done in one-go
                T(i-deleted)=[];
                Cs(i-deleted)=[];
                Cl(i-deleted)=[];
                if Has_volatile==1
                    Cs2(i-deleted-1)=(Cs2(i-deleted-1)+Cs2(i-deleted))/2;
                    Cl2(i-deleted-1)=(Cl2(i-deleted-1)+Cl2(i-deleted))/2;
                    S(i-deleted-1)=(S(i-deleted-1)+S(i-deleted))/2;
                    Cs2(i-deleted)=[];
                    Cl2(i-deleted)=[];
                    S(i-deleted)=[];
                end

                % OG_Cb(i-deleted-1,:)=(OG_Cb(i-deleted-1,:)*dz(i-deleted-1)+OG_Cb(i-deleted,:)*dz(i-deleted))/(dz(i-deleted-1)+dz(i-deleted));
                % OG_Cb(i-deleted,:)=0;                    

                % if To_recored_intrusion_marker==1
                %     Intrusion_marker_static(i-deleted)=[];
                %     Intrusion_marker_dynamics(i-deleted)=[];
                % end

                % rho(i-deleted,:)=[];


                Ts(i-deleted)=[];
                Tl(i-deleted)=[];
                Lsaturation(i-deleted)=[];
                Ssaturation(i-deleted)=[];

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
        
        fine_margin=500; %meter
        index=find(phi>0,1,'last');
        index2=find(nodez>nodez(index)+fine_margin,1,'first');
        to_refine=zeros(length(phi)-1,1);
        if index2-index<fine_margin/min_dx/2
            to_refine(index:index2)=1;
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
            phi_new = [phi_new(1:i); phi_new(i:i+1); phi_new(i+1:end)];
            T_new = [T_new(1:i); T_new(i:i+1); T_new(i+1:end)];
            Cl_new = [Cl_new(1:i); Cl_new(i:i+1); Cl_new(i+1:end)];
            Cs_new = [Cs_new(1:i); Cs_new(i:i+1); Cs_new(i+1:end)];
            if Has_volatile==1
                S_new = [S_new(1:i); S_new(i:i+1); S_new(i+1:end)];
                Cl2_new = [Cl2_new(1:i); Cl2_new(i:i+1); Cl2_new(i+1:end)];
                Cs2_new = [Cs2_new(1:i); Cs2_new(i:i+1); Cs2_new(i+1:end)];
            end
            offset=offset+2;
        end
        

        node_new = nodez;
        shift = 0;

        for k = 1:numel(to_refine)
            i = to_refine(k) + shift+1;

            left_avg  = (node_new(i-1) + node_new(i)) / 2;
            right_avg = (node_new(i) + node_new(i+1)) / 2;

            node_new = [node_new(1:i-1), left_avg, node_new(i), right_avg, node_new(i+1:end)];
            shift = shift + 2;   
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

end

