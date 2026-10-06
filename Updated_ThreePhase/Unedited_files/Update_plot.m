if Z_label_type==1
    Ytick1=cellz;
    Ytick2=nodez;
else
    Ytick1=(cellz-nodez(end))/1000; % originally (cellz-nodez(end))/1000,but did not understand why 
    Ytick2=(nodez-nodez(end))/1000;
end

if Adaptive_show_range==1
    % if is_solid_solution==0 
    %     Active=find(T>=Ts_new);
    % else
    %     Active=find(T>=Ts_local);
    % end
    Active=find(phi>0);
    if ~isempty(Active)
        if Z_label_type==1
            Margin=[-100 100]';  %in meters
            Active_range=cellz([Active(1) Active(end)]);
        else
            Margin=[-100 100]';  %in meters
            Active_range=(cellz([Active(1) Active(end)])+Margin-nodez(end))/1000;
        end
        if Active_range(1)<Show_z(1) && Vol_flux==0
            % Show_z(1)=Active_range(1)+Margin(1);
            Show_z(1)=Active_range(1)-(Active_range(2)-Active_range(1))*0.13;
        end
        if Active_range(2)>Show_z(2)
            % Show_z(2)=Active_range(2)+Margin(2);
            Show_z(2)=Active_range(2)+(Active_range(2)-Active_range(1))*0.2;
        end
        for i=1:(Plot_size(1)*Plot_size(2))
            set(ha(i),'ylim',Show_z)
        end
        if exist('ax_volatile','var')
            set(ax_volatile,'ylim',Show_z)
        end
    end
end

if N_component==5 || simplified_5p==1 || Add_CLCU==1
    is_vol=find(sum(Mass_data(:,[7 10 13]),2)>1e-5);
    not_vol=setxor(1:N,is_vol);
end
for i=1:(Plot_size(1)*Plot_size(2))
    if isscalar(Plot_configure{i})
        j=Plot_configure{i};
        switch j
            case 1
                set(Handel_all{i}(j),'xdata',phi(1:N),'ydata',Ytick1);
            case 19
                set(Handel_all{i}(j),'xdata',S(1:N),'ydata',Ytick1);
            case 2
                set(Handel_all{i}(j),'xdata',u_all(1:N+1),'ydata',Ytick2);
            case 20
                set(Handel_all{i}(j),'xdata',u_all(2*N+3: 3*N+3),'ydata',Ytick1);
            case 4
                set(Handel_all{i}(j),'xdata',H,'ydata',Ytick1);
            case 5
                set(Handel_all{i},'xdata',T,'ydata',Ytick1);
                if is_solid_solution==0
                    set(Handel_ts,'xdata',Ts_new,'ydata',Ytick1);
                else
                    set(Handel_ts,'xdata',Ts_local,'ydata',Ytick1);
                    set(Handel_tl,'xdata',Tl_local,'ydata',Ytick1);
                end
            case 6
                set(Handel_all{i}(j),'xdata',C_all(1:N),'ydata',Ytick1);
            case 7
                set(Handel_all{i}(j),'xdata',C_all(N+1:2*N),'ydata',Ytick1);
            case 8
                set(Handel_all{i}(j),'xdata',Cb,'ydata',Ytick1);                
            case 25
                set(Handel_all{i}(j),'xdata',M2,'ydata',Ytick1);    
            case 9
                set(Handel_all{i}(j),'xdata',TYPE(1:N),'ydata',Ytick1);
            case 3
                set(Handel_all{i}(j),'xdata',u_all(N+2:2*N+2),'ydata',Ytick2);
            case 10
                set(Handel_all{i}(1),'xdata',Portion2(:,3),'ydata',Ytick1);
                set(Handel_all{i}(2),'xdata',Portion2(:,2),'ydata',Ytick1);
                set(Handel_all{i}(3),'xdata',Portion2(:,1),'ydata',Ytick1);
            case 11
                set(Handel_all{i}(3), 'xdata',Con_C,'ydata',Ytick1)
                set(Handel_all{i}(2), 'xdata',Con_M,'ydata',Ytick1)
                set(Handel_all{i}(1), 'xdata',Con_T,'ydata',Ytick1)
            case 12
                %                 dphi=phi(1:N)-phi_old(1:N);
                set(Handel_all{i}(3), 'xdata',ACon_C,'ydata',Ytick1)
                set(Handel_all{i}(2), 'xdata',ACon_M,'ydata',Ytick1)
                set(Handel_all{i}(1), 'xdata',ACon_T,'ydata',Ytick1)
            case 13
                set(Handel_all{i}(j),'xdata',refined_Cb);
            case 14
                set(Handel_all{i}(3), 'xdata',refined_ACon_C)
                set(Handel_all{i}(2), 'xdata',refined_ACon_M)
                set(Handel_all{i}(1), 'xdata',refined_ACon_T)
            case 15
                %                 dphi=phi(1:N)-phi_old(1:N);
                set(Handel_all{i}(3), 'xdata',ACon_C2,'ydata',Ytick1)
                set(Handel_all{i}(2), 'xdata',ACon_M2,'ydata',Ytick1)
                set(Handel_all{i}(1), 'xdata',ACon_T2,'ydata',Ytick1)
            case 16
                if S_or_Cb==1
                    set(Handel_all{i}(1), 'xdata',CD_M*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',CD_R*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                else
                    set(Handel_all{i}(1), 'xdata',CD_M,'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',CD_R,'ydata',Ytick1)
                end
            case 17
                if S_or_Cb==1
                    set(Handel_all{i}(1), 'xdata',CD_Ms*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',-CD_Rs*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                else
                    set(Handel_all{i}(1), 'xdata',CD_Ms,'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',CD_Rs,'ydata',Ytick1)
                end
            case 18
                if S_or_Cb==1
                    set(Handel_all{i}(1), 'xdata',CD_Mabs*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',CD_Rabs*(SiO2_range(2)-SiO2_range(1)),'ydata',Ytick1)
                else
                    set(Handel_all{i}(1), 'xdata',CD_Mabs,'ydata',Ytick1)
                    set(Handel_all{i}(2), 'xdata',CD_Rabs,'ydata',Ytick1)
                end
            case 21
                set(Handel_all{i}(3), 'xdata',(P(:,1)*scaling),'ydata',Ytick2)
                set(Handel_all{i}(2), 'xdata',(P(:,2)*scaling),'ydata',Ytick2)
                set(Handel_all{i}(1), 'xdata',(P(:,3)*scaling),'ydata',Ytick2)
            case 22
                set(Handel_all{i}(4), 'xdata',-rho(:,3)*g*Dynamic_scaling,'ydata',Ytick2)
                set(Handel_all{i}(1), 'xdata',dP(:,1)*scaling,'ydata',Ytick2,'linewidth',1)
                set(Handel_all{i}(2), 'xdata',dP(:,2)*scaling,'ydata',Ytick2,'linewidth',1.2)
                set(Handel_all{i}(3), 'xdata',dP(:,3)*scaling,'ydata',Ytick2,'linewidth',1.5)
            case 23
                set(Handel_all{i}(1), 'xdata',rho(phi.*(1-S)>0,1),'ydata',Ytick1(phi.*(1-S)>0)) %phi.*(1-S)>0
                set(Handel_all{i}(2), 'xdata',rho(:,2),'ydata',Ytick1)
                set(Handel_all{i}(3), 'xdata',rho(phi.*S>0,3),'ydata',Ytick1(phi.*S>0))
            case 24
                Sat=zeros(N,1);
                PT=[-4e-4,-5e-4,-6e-4,-13e-4, -15.5e-4,-17e-4,-16e-4,-5e-4, 0, 26e-4,5e-3, 5e-3];
                PTx=[0,   0.12    0.2  0.3    0.5       1       2       3   4  5,    11, 20];
                P2=Pg_real*100;
                for k=1:N
                    ind=find(Pg_real(k)<=PTx,1,'first');
                    dSdT=(PT(ind)*(Pg_real(k)-PTx(ind-1))+PT(ind-1)*(PTx(ind)-Pg_real(k)))/(PTx(ind)-PTx(ind-1));
                    Sat(k)=((2.859e-2*P2(k)-1.495e-3*P2(k).^1.5+2.702e-5*P2(k).^2+0.257*P2(k).^0.5)+(T(k)-800)*dSdT)/100;
                end
                Cb2=(sum(Mass_data(:,[3,6,7]),2))./sum(Mass_data,2);

                Solid_saturation=ones(N,1)*Solid_saturation0;
                %                 for k=1:N
                % %                     [~,~, temp]=Cal_liquidus(19,Pg_real(k)*1000,Data_point2,Data_y2);
                %                     Solid_saturation(k)=Solid_saturation0;%(temp+2)/100;
                %                 end
                if N_component==3
                    Cl2=(100*Mass_data(:,3))./(sum(Mass_data(:,1:3),2));
                    Cs2=(100*Mass_data(:,6))./(sum(Mass_data(:,4:6),2));
                else
                    Cl2=(100*Mass_data(:,3))./(sum(Mass_data(:,[1:3 8]),2));
                    Cs2=(100*Mass_data(:,6))./(sum(Mass_data(:,[4:6 9]),2));
                end
                set(handle_saturation,'xdata',Sat*100,'ydata',Ytick1);
                %                 if is_solid_solution==0
%                     set(handle_saturation2,'xdata',S_cap*100,'ydata',Ytick1);
%                 else
                    set(handle_saturation2,'xdata',(S_cap)*100,'ydata',Ytick1);
%                 end
                set(Handel_all{i}(1), 'xdata',Cl2.*(phi>1e-2),'ydata',Ytick1)
                set(Handel_all{i}(2), 'xdata',Cs2,'ydata',Ytick1)    
                set(Handel_all{i}(3), 'xdata',Cb2*100,'ydata',Ytick1)
                set(Handel_all{i}(4), 'xdata',(Mass_data(:,7)+max(Mass_data(:,6)-Mass_solid.*S_cap,0))/rho_mean*100,'ydata',Ytick1)
            case 99
                set(Handel_all{i}(1), 'xdata',sum(Mass_data,2),'ydata',Ytick1)
            case 98
                set(Handel_all{i}(1), 'xdata',phi.*(1-S),'ydata',Ytick1)
            case 97
                set(Handel_all{i}(1), 'xdata',phi.*S.*(rho(:,3)>1),'ydata',Ytick1)
            case 30
                Cb_cl=CL./(MM+NN+V+CL+CU)*100;
                Cll=Mass_data(:,8)./sum(Mass_data(:,[1 2 3 8 11]),2)*100;
                Cls=Mass_data(:,9)./sum(Mass_data(:,[4 5 6 9 12]),2)*100;
                Clg=Mass_data(:,10)./sum(Mass_data(:,[7 10 13]),2)*100;
                set(Handel_all{i}(1), 'xdata',Cb_cl,'ydata',Ytick1);
                % set(Handel_all{i}(2), 'xdata',Cll,'ydata',Ytick1);
                % set(Handel_all{i}(3), 'xdata',Cls,'ydata',Ytick1);
                % set(Handel_all{i}(4), 'xdata',Clg,'ydata',Ytick1);
                % x=Partition_CL_CU(:,1);
                % x(not_vol)=nan;
                % set(Handel_all{i}(5), 'xdata',x,'ydata',Ytick1);      
            case 31
                Cb_cu=CU./(MM+NN+V+CL+CU)*100;
                % Cul=Mass_data(:,11)./sum(Mass_data(:,[1 2 3 8 11]),2)*100;
                % Cus=Mass_data(:,12)./sum(Mass_data(:,[4 5 6 9 12]),2)*100;
                % Cug=Mass_data(:,13)./sum(Mass_data(:,[7 10 13]),2)*100;
                set(Handel_all{i}(1), 'xdata',Cb_cu,'ydata',Ytick1);
                % set(Handel_all{i}(2), 'xdata',Cul,'ydata',Ytick1);
                % set(Handel_all{i}(3), 'xdata',Cus,'ydata',Ytick1);
                % set(Handel_all{i}(4), 'xdata',Cug,'ydata',Ytick1);
                % x=Partition_CL_CU(:,2);
                % x(not_vol)=nan;
                % set(Handel_all{i}(5), 'xdata',x,'ydata',Ytick1);   
            case 32
                Cb_sul=SUL./(MM+NN+V+CL+CU)*100;
                % Cul=Mass_data(:,11)./sum(Mass_data(:,[1 2 3 8 11]),2)*100;
                % Cus=Mass_data(:,12)./sum(Mass_data(:,[4 5 6 9 12]),2)*100;
                % Cug=Mass_data(:,13)./sum(Mass_data(:,[7 10 13]),2)*100;
                set(Handel_all{i}(1), 'xdata',Cb_sul,'ydata',Ytick1);
                % set(Handel_all{i}(2), 'xdata',Cul,'ydata',Ytick1);
                % set(Handel_all{i}(3), 'xdata',Cus,'ydata',Ytick1);
                % set(Handel_all{i}(4), 'xdata',Cug,'ydata',Ytick1);
                % x=Partition_CL_CU(:,2);
                % x(not_vol)=nan;
                % set(Handel_all{i}(5), 'xdata',x,'ydata',Ytick1);   
            case 50
                if melt_viscosity_type==0
                    %old melt viscosity model
                    muf_order=(muf2*Mass_data(:,1)+muf1*Mass_data(:,2))./sum(Mass_data(:,1:2),2);
                else
                    %new melt viscosity model
                    SiO2=(Mass_data(:,1)*74+Mass_data(:,2)*47)./max(sum(Mass_data(:,1:2),2),1e-3); % SiO2%
                    if Fix_H2O==1 || Fix_H2O==3
                        H2O=H2O_crust*ones(N,1);
                    else
                        H2O=100*Mass_data(:,3)./sum(Mass_data(:,1:3),2); % H2O%
                    end
                    vis_B=vis_b1*SiO2+vis_b2*H2O+vis_b3*log(1+H2O);
                    vis_C=vis_c1*SiO2+vis_c2*H2O+vis_c3*log(1+H2O);
                    muf_order=-4.55+vis_B./(T-vis_C);
                end                
            case 51
                if Conservation_type==1
                    set(Handel_all{i}(1), 'xdata',dphidx(:,1),'ydata',Ytick1);
                    set(Handel_all{i}(2), 'xdata',dphidx(:,2),'ydata',Ytick1);
                    set(Handel_all{i}(3), 'xdata',dphidx(:,3),'ydata',Ytick1);
                else
                    set(Handel_all{i}(1), 'xdata',acc_dphidx(:,1),'ydata',Ytick1);
                    set(Handel_all{i}(2), 'xdata',acc_dphidx(:,2),'ydata',Ytick1);
                    set(Handel_all{i}(3), 'xdata',acc_dphidx(:,3),'ydata',Ytick1);
                end
            case 101
                set(Handel_all{i}(1), 'xdata',Intrusion_marker_static,'ydata',Ytick1);
            case 102
                set(Handel_all{i}(1), 'xdata',Intrusion_marker_dynamics,'ydata',Ytick1);
        end
    else
        for j=Plot_configure{i}
            switch j
                case 1
                    set(Handel_all{i}(1),'xdata',phi(1:N),'ydata',Ytick1);
                    set(Handel_all{i}(2),'xdata',zeros(phi_top_i,1),'ydata',(cellz(phi_top)-nodez(end))/1000);
                    set(Handel_all{i}(3),'xdata',zeros(phi_top_i,1),'ydata',(cellz(phi_base)-nodez(end))/1000);
                   % set(Handel_all{i}(4),'xdata',ones(Hphi_top_i,1),'ydata',(cellz(Hphi_top)-nodez(end))/1000);
                   % set(Handel_all{i}(5),'xdata',ones(Hphi_top_i,1),'ydata',(cellz(Hphi_base)-nodez(end))/1000);
                case 19
                    set(Handel_all{i}(j),'xdata',S(1:N),'ydata',Ytick1);
                case 2
                    set(Handel_all{i}(j),'xdata',u_all(1:N+1)/UM_0,'ydata',Ytick2);
                case 4
                    set(Handel_all{i}(j),'xdata',H/(HLiquidus-HSolidus),'ydata',Ytick1);
                case 5
                    set(Handel_all{i}(j),'xdata',(T-Solidus)/(Max_Liquidus-Solidus),'ydata',Ytick2);
                case 6
                    Cl1l=(74*Mass_data(:,1)+47*Mass_data(:,2))./(sum(Mass_data(:,1:2),2));
                    set(Handel_all{i}(j),'xdata',Cl1l,'ydata',Ytick1);
                case 7
                    Cs1S=(74*Mass_data(:,4)+47*Mass_data(:,5))./(sum(Mass_data(:,4:5),2));
                    set(Handel_all{i}(j),'xdata',Cs1S,'ydata',Ytick1);
                case 8
                    set(Handel_all{i}(j),'xdata',Cb,'ydata',Ytick1);
                case 9
                    set(Handel_all{i}(j),'xdata',TYPE(1:N)/3,'ydata',Ytick1);
                case 3
                    set(Handel_all{i}(j),'xdata',u_all(N+2:2*N+2)/UM_0,'ydata',Ytick2);
                case 20
                    set(Handel_all{i}(j),'xdata',u_all(2*N+3: 3*N+3)/UM_0,'ydata',Ytick2);
            end
        end
    end    
end



if exist('Handel_T','var')
    eval("set(Handel_T, 'string', string('Time= ')+num2str(round(Time/Year))+string('y'))")
end
if exist('Handel_dt','var')
    eval("set(Handel_dt, 'string', string('dt= ')+num2str(round(dt/Year,2))+string('y'))")
end

if exist('Handel_N','var')
    eval("set(Handel_N, 'string', string('Nodes= ')+num2str(N))")
end

if exist('Gca1','var')
    set(Gca1,'xlim',[round(min(M1),1)-1,   round(max(M1),1)+1]);
end

if exist('Gca2','var')
    set(Gca2,'xlim',[round(min(M2),1)-1,   round(max(M2),1)+1]);
end
