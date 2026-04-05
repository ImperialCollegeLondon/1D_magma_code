if Step_counts>=Update_frame || Time>=End_time
    Step_counts=0;

    for i=1:(Plot_size(1)*Plot_size(2))
        if isscalar(Plot_configure{i})
            j=Plot_configure{i};
            switch j
                case 1
                    set(Handel_all{i}(j),'xdata',phi(1:N),'ydata',cellz/1000-Base_crust);
                    if Has_volatile==1
                        set(Handel_all{i}(2),'xdata',S,'ydata',cellz/1000-Base_crust);
                    end
                case 2
                    set(Handel_all{i}(j),'xdata',u_all(1:N+1),'ydata',nodez);
                case 4
                    set(Handel_all{i}(j),'xdata',H,'ydata',cellz/1000-Base_crust);
                case 5
                    set(Handel_all{i}(1),'xdata',T,'ydata',cellz/1000-Base_crust);
                    set(Handel_all{i}(2),'xdata',Ts,'ydata',cellz/1000-Base_crust);
                    set(Handel_all{i}(3),'xdata',Tl,'ydata',cellz/1000-Base_crust);
                case 6
                    temp=C_all(1:N);
                    temp(phi<1e-6)=nan;
                    set(Handel_all{i}(j),'xdata',temp(1:N),'ydata',cellz/1000-Base_crust);
                case 7
                    set(Handel_all{i}(j),'xdata',C_all(N+1:2*N),'ydata',cellz/1000-Base_crust);
                case 8
                    set(Handel_all{i}(j),'xdata',Cb,'ydata',cellz/1000-Base_crust);
                case 9
                    set(Handel_all{i}(j),'xdata',TYPE(1:N),'ydata',cellz/1000-Base_crust);
                case 3
                    set(Handel_all{i}(j),'xdata',u_all(N+2:2*N+2),'ydata',nodez);
                case 10
                    set(Handel_all{i}(1),'xdata',Portion2(:,3),'ydata',cellz/1000-Base_crust);
                    set(Handel_all{i}(2),'xdata',Portion2(:,2),'ydata',cellz/1000-Base_crust);
                    set(Handel_all{i}(3),'xdata',Portion2(:,1),'ydata',cellz/1000-Base_crust);
                case 11
                    set(Handel_all{i}(3), 'xdata',Con_C,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(2), 'xdata',Con_M,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(1), 'xdata',Con_T,'ydata',cellz/1000-Base_crust)
                case 12
                    %                 dphi=phi(1:N)-phi_old(1:N);
                    set(Handel_all{i}(3), 'xdata',ACon_C,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(2), 'xdata',ACon_M,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(1), 'xdata',ACon_T,'ydata',cellz/1000-Base_crust)
                case 13
                    set(Handel_all{i}(j),'xdata',refined_Cb);
                case 14
                    set(Handel_all{i}(3), 'xdata',refined_ACon_C)
                    set(Handel_all{i}(2), 'xdata',refined_ACon_M)
                    set(Handel_all{i}(1), 'xdata',refined_ACon_T)
                case 15
                    %                 dphi=phi(1:N)-phi_old(1:N);
                    set(Handel_all{i}(3), 'xdata',ACon_C,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(2), 'xdata',ACon_M,'ydata',cellz/1000-Base_crust)
                    set(Handel_all{i}(1), 'xdata',ACon_T,'ydata',cellz/1000-Base_crust)
                case 16
                    if S_or_Cb==1
                        if FourMPD==1
                            set(Handel_all{i}(1), 'xdata',CD_M*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',CD_R*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                        else
                            set(Handel_all{i}(1), 'xdata',CD_M*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',CD_R*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                        end
                    else
                        set(Handel_all{i}(1), 'xdata',CD_M,'ydata',cellz/1000-Base_crust)
                        set(Handel_all{i}(2), 'xdata',CD_R,'ydata',cellz/1000-Base_crust)
                    end
                case 17
                    if S_or_Cb==1
                        if FourMPD==1
                            set(Handel_all{i}(1), 'xdata',CD_Ms*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',-CD_Rs*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                        else
                            set(Handel_all{i}(1), 'xdata',CD_Ms*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',-CD_Rs*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                        end
                    else
                        set(Handel_all{i}(1), 'xdata',CD_Ms,'ydata',cellz/1000-Base_crust)
                        set(Handel_all{i}(2), 'xdata',CD_Rs,'ydata',cellz/1000-Base_crust)
                    end
                case 18
                    if S_or_Cb==1
                        if FourMPD==1
                            set(Handel_all{i}(1), 'xdata',CD_Mabs*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',CD_Rabs*(MgO_range(2)-MgO_range(1)),'ydata',cellz/1000-Base_crust)
                        else
                            set(Handel_all{i}(1), 'xdata',CD_Mabs*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                            set(Handel_all{i}(2), 'xdata',CD_Rabs*(SiO2_range(2)-SiO2_range(1)),'ydata',cellz/1000-Base_crust)
                        end
                    else
                        set(Handel_all{i}(1), 'xdata',CD_Mabs,'ydata',cellz/1000-Base_crust)
                        set(Handel_all{i}(2), 'xdata',CD_Rabs,'ydata',cellz/1000-Base_crust)
                    end
                case 19
                    phi_n=project_cell2node(nodez,dz,phi(1:N));
                    dP=u_all(N+2:2*N+2)*C_value_A.*mu_f_dynamics/grain_size^2./max(phi_n,1e-6).^C_value_B-rhof*g;
                    set(Handel_all{i}(1), 'xdata',-dP*scaling_factor,'ydata',nodez);
                    set(Handel_all{i}(2), 'xdata',rhof*g*scaling_factor,'ydata',nodez);
                    set(Handel_all{i}(3), 'xdata',rhom*g*scaling_factor,'ydata',nodez);
                    %                         xlim([rhof_2*g*scaling_factor max(-dP*scaling_factor)*1.1])
                case 20
                    dP=dP+rhof*g;
                    rhof=rhof_2*C_all(1:N)+(1-C_all(1:N))*rhof_1;
                    rhom=rhom_2*C_all(N+1:2*N)+(1-C_all(N+1:2*N))*rhom_1;
                    drho=(rhom-rhof);
                    temp=drho.*dz'*g;
                    dhydro=zeros(N,1);
                    dP_cell=(dP(1:end-1)+dP(2:end))/2.*dz';
                    for k=1:N
                        dhydro(k)=sum(temp(k:end));
                        Pressure(k)=sum(dP_cell(k:end)); %#ok<SAGROW>
                    end
                    set(Handel_all{i}(1), 'xdata',-Pressure*scaling_factor,'ydata',nodez);
                    set(Handel_all{i}(2), 'xdata', dhydro*scaling_factor,'ydata',cellz);
                    %
                case 1001
                    set(Handel_all{i}(1), 'ydata',cellz/1000-Base_crust,'xdata',Cb2);
                    set(Handel_all{i}(2), 'ydata',cellz/1000-Base_crust,'xdata',Cl2);
                    set(Handel_all{i}(3), 'ydata',cellz/1000-Base_crust,'xdata',Cs2);
                    set(Handel_all{i}(4), 'ydata',cellz/1000-Base_crust,'xdata',Lsaturation);
                    set(Handel_all{i}(5), 'ydata',cellz/1000-Base_crust,'xdata',Ssaturation);

            end
        else
            for j=Plot_configure{i}
                switch j
                    case 1
                        set(Handel_all{i}(j),'xdata',phi(1:N),'ydata',cellz/1000-Base_crust);
                    case 2
                        set(Handel_all{i}(j),'xdata',u_all(1:N+1),'ydata',nodez/1000-Base_crust); %/UM_0
                    case 4
                        set(Handel_all{i}(j),'xdata',H/(HLiquidus-HSolidus),'ydata',cellz/1000-Base_crust);
                    case 5
                        set(Handel_all{i}(j),'xdata',(T-Solidus)/(Max_Liquidus-Solidus),'ydata',cellz/1000-Base_crust);
                    case 6
                        if FourMPD==1
                            temp=((C_all(1:N))*(MgO_range(2)-MgO_range(1))+MgO_range(1));
                        else
                            temp=((C_all(1:N))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1));
                        end
                        temp(phi(1:N)<1e-6)=nan;
                        set(Handel_all{i}(j),'xdata',temp,'ydata',cellz/1000-Base_crust);
                    case 7
                        if FourMPD==1
                            set(Handel_all{i}(j),'xdata',((C_all(N+1:2*N))*(MgO_range(2)-MgO_range(1))+MgO_range(1)),'ydata',cellz/1000-Base_crust);
                        else
                            set(Handel_all{i}(j),'xdata',((C_all(N+1:2*N))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)),'ydata',cellz/1000-Base_crust);
                        end

                    case 8
                        if FourMPD==1
                            set(Handel_all{i}(j),'xdata',((Cb)*(MgO_range(2)-MgO_range(1))+MgO_range(1)),'ydata',cellz/1000-Base_crust);
                        else
                            set(Handel_all{i}(j),'xdata',((Cb)*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)),'ydata',cellz/1000-Base_crust);
                        end
                    case 9
                        set(Handel_all{i}(j),'xdata',TYPE(1:N)/3,'ydata',cellz/1000-Base_crust);
                    case 3
                        set(Handel_all{i}(j),'xdata',u_all(N+2:2*N+2),'ydata',nodez/1000-Base_crust);
                end
            end
        end
    end

    if exist('Handel_com_conserve','var')
        set(Handel_com_conserve,'string',num2str(sum(Cb(Conservative_pos(1):Conservative_pos(2)).*dz(Conservative_pos(1):Conservative_pos(2))')/Total_composition))
    end
    if exist('Handel_ent_conserve','var')
        set(Handel_ent_conserve,'string',num2str(sum(H.*dz')/Total_enthalpy))
    end
    if exist('Handel_T','var')
        set(Handel_T, 'string', ['T= ' num2str(round(Time/Year,2)) 'y'])
    end
    if exist('Handel_dt','var')
        set(Handel_dt, 'string', ['dT= ' num2str(dt/Year)])
    end
    if exist('Handel_CDN','var')
        CDN=(sum(Cb(CDN_index2).*(cellz(CDN_index2)-Lrange1(1))'.*dz(CDN_index2)')-CD0)/(CD1-CD0);
        set(Handel_CDN, 'string', string('CDN=')+num2str(round(CDN,3)))
    end
    set(TITLE,'string', "Time="+ num2str(Time/Year,4) +'y   '+ "N\_nodes="+ num2str(N))
end

drawnow;