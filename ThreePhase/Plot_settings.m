% HH 9/30/2021
% Initializing all the plots for three phase simulation

Plot_Label={'Porosity, $$\phi (-)$$','uf','um','H','Temperature($$^\circ$$C)','$$C_{l1}$$','$$C_{s1}$$','$$\bar{C}_{SiO_2}$$','St', 'Solid_vis(%)','$$\Delta\phi$$','$$\sum|\Delta\phi|$$', 'Cb', '$$\sum|\Delta\phi|$$', '$$\sum\Delta\phi$$','$$\sum\Delta\bar{C}$$','$$\Delta\bar{C}$$','$$\sum|\Delta\bar{C}|$$', 'Volatile fraction, $$S (-)$$', 'ug'};


marker_to_use=1; 
Markers={'none','o'};

Step_counts=0;
% Initialize all the plottings
Plot_size=size(Plot_configure);
Handel_all=Plot_configure;

UM_0=1;

% Cb=Cm.*phi(N+1:2*N)+Cf.*phi(1:N);

S_or_Cb=1;  %1 for showing SiO2
SiO2_range=[47 74];
Font_Size=11;

if Display_type==1
    gap=[0.02 0.02];
    marg_h=[0.06, 0.07];
    marg_w=[0.05,0.01];
else
    gap=[0.02 0.02];
    marg_h=[0.072, 0.075];
    marg_w=[0.08,0.01];    
end


if Z_label_type==1
    Ytick1=cellz;
    Ytick2=nodez;
else
    Ytick1=(cellz-nodez(end))/1000;
    Ytick2=(nodez-nodez(end))/1000;
%     Show_z=(Show_z-nodez(end))/1000;
end

figure(13)
clf; set(gcf,'color','w');
set(gcf,'Units','normalized')
set(gcf,'Position',[0         0    1.0000    0.9325]);

Plot_configure=reshape(Plot_configure',[1,Plot_size(1)*Plot_size(2)]);

[ha, pos] = tight_subplot(Plot_size(1), Plot_size(2), gap, marg_h, marg_w);

if N_component==5
    is_vol=find(sum(Mass_data(:,[7 10 13]),2)>1e-5);
    not_vol=setxor(1:N,is_vol);
end
for i=1:(Plot_size(1)*Plot_size(2))
    if isscalar(Plot_configure{i})
        j=1;
        set(gcf,'CurrentAxes',ha(i)); hold on;set(gca,'tickdir','out') ; set(gca,'fontsize', Font_Size);    
        switch Plot_configure{i}
            case 1
                Handel_all{i}(j)=plot(phi(1:N),Ytick1,'Marker',Markers{marker_to_use});

                %Handel_all{i}(2)=plot(0,(cellz(phi_top)-nodez(end))/1000, 'kx', 'MarkerSize', 12);
                %Handel_all{i}(3)=plot(0,(cellz(phi_base)-nodez(end))/1000, 'kx', 'MarkerSize', 12);
                %Handel_all{i}(4)=plot(1,(cellz(Hphi_top)-nodez(end))/1000, 'k+', 'MarkerSize', 12);
                %Handel_all{i}(5)=plot(1,(cellz(Hphi_base)-nodez(end))/1000, 'k+', 'MarkerSize', 12);
                
                ylim(Show_z)
                
                xlim([0,1])
                xlabel('Porosity, $$\phi$$(-)','interpreter','Latex','FontSize',Font_Size)
            case 2
                Handel_all{i}(j)=plot(u_all(1:N+1),Ytick2);
                ylim(Show_z)
            case 4
                Handel_all{i}(j)=plot(H,Ytick1);

                ylim(Show_z);
%                 if Phase_Type==1
%                     %                 xlim([0,max(H)*1.5]);
%                     gapH=max([max(H)-min(H),abs(max(H)-ESolidus(0,cp,rhof_1)),100]);
%                     xlim([min(min(H)-1e-5-gapH/20,ESolidus(0,cp,rhof_1)-gapH/20-1e-5),max(max(H)+1e-5+gapH/2,ESolidus(0,cp,rhof_1)+gapH/3+1e-5)])
%                 else
%                     gapH=max([max(H)-min(H),abs(max(H)-ESolidus(0,cp,rhof_1)),100]);
%                     xlim([min(H)-gapH,max(H)+gapH])
%                 end
            case 5
                Handel_all{i}(j)=plot(T,Ytick1);
%                 if Phase_Type==1
                if is_solid_solution==0
                    Handel_ts=plot(Ts_new,Ytick1,'--','color','red','linewidth',1.5);
                else
                    Handel_ts=plot(Ts_local,Ytick1,'--','color','red','linewidth',1.5);
                    Handel_tl=plot(Tl_local,Ytick1,'--','color','black','linewidth',1.5);
                end
                    
                    % gapT=max(max(T)-min(T),abs(max(T)-Ts));
                    % if (max(T)-Ts)*(min(T)-Ts)<0
                    %     xlim([min([T; Ts_new])-10 max([T; Ts_new])+40])
                    % else
                    %     %                     xlim([Solidus(0)-gapT/3-1e-5,max(T)+gapT/3+1e-5])
                    %     xlim([min(min(T)-1e-5-gapT/10,Ts-gapT/10-1e-5),max(max(T)+1e-5+gapT/2,Ts+gapT/3+1e-5)])
                    % end
                    xlim([500 max(Tl_local)+50])

                ylim(Show_z);
                xlabel('Temperature ($$^\circ$$C)','interpreter','Latex','FontSize',Font_Size)
            case 6
                Handel_all{i}(j)=plot(C_all(1:N),Ytick1);
                ylim(Show_z)
                xlim([0,1]);
            case 7
                Handel_all{i}(j)=plot(C_all(N+1:2*N),Ytick1);
                ylim(Show_z)
                xlim([0,1]);
            case 8
                Gca1=gca;
                Handel_all{i}(j)=plot(M1,Ytick1);
                xlabel('Major component mass per unit,$$M_1$$','interpreter','Latex','FontSize',Font_Size)
                xlim([round(min(M1))-1,   round(max(M1))+1]);
                ylim(Show_z)
            case 25
                Gca2=gca;
                Handel_all{i}(j)=plot(M2,Ytick1);
                xlabel('Volatile component mass per unit,$$M_2$$','interpreter','Latex','FontSize',Font_Size)
                xlim([round(min(M2))-1,   round(max(M2))+1]);
                ylim(Show_z)
            case 9
                Handel_all{i}(j)=plot(TYPE(1:N),Ytick1);
                ylim(Show_z)
                xlim([-0.3,3.3]);
            case 3
                Handel_all{i}(j)=plot(u_all(N+2:2*N+2),Ytick2);
                ylim(Show_z)  
            case 10
                Handel_all{i}(1)=plot(Force_coupling,Ytick2,'r');
                set(Handel_all{i}(1),'visible','off');
                Handel_all{i}(2)=plot(Force_solid,Ytick2,'g');
                Handel_all{i}(3)=plot(Force_liquid,Ytick2,'b');
                ylim(Show_z)
                xlim([0,100]);
            case 11
                Handel_all{i}(1)=plot(Con_T,Ytick1,'b');
                Handel_all{i}(2)=plot(Con_M,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(Con_C,Ytick1,'g');
                xlabel('$$\Delta\phi$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Heat','Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 12
                Handel_all{i}(1)=plot(ACon_T,Ytick1,'b');
                Handel_all{i}(2)=plot(ACon_M,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(ACon_C,Ytick1,'g');
                xlabel('$$\sum|\Delta\phi|$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Heat','Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 13      
                refined_Cb  =interp1(Ytick1,Cb,cell_refined,'linear','extrap');
                Handel_com_conserve=text(0.7,0.8*Length,"1",'fontsize',Font_Size);
                Handel_all{i}(j)=plot(refined_Cb,cell_refined);
                xlim([0,1]);
            case 14
                refined_ACon_T=interp1(Ytick1,ACon_T,cell_refined,'linear','extrap');
                refined_ACon_M=interp1(Ytick1,ACon_M,cell_refined,'linear','extrap');
                refined_ACon_C=interp1(Ytick1,ACon_C,cell_refined,'linear','extrap');                
                Handel_all{i}(1)=plot(refined_ACon_T,cell_refined,'b');
                Handel_all{i}(2)=plot(refined_ACon_M,cell_refined,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(refined_ACon_C,cell_refined,'g');
                xlabel('$$\sum|\Delta\phi|$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Heat','Comp','React'},'FontSize',Font_Size)    
                ylim(Show_z)
            case 15
                Handel_all{i}(1)=plot(ACon_T2,Ytick1,'b');
                Handel_all{i}(2)=plot(ACon_M2,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(ACon_C2,Ytick1,'g');
                xlabel('$$\sum\Delta\phi$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Heat','Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 16
                Handel_all{i}(1)=plot(CD_M,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_R,Ytick1,'g');
                if S_or_Cb==1
                    xlabel('$$\sum\Delta\bar{C}$$ (SiO$$_2 \%$$)','interpreter','Latex','FontSize',Font_Size)                    
                else
                    xlabel('$$\sum\Delta\bar{C}$$ (-)','interpreter','Latex','FontSize',Font_Size)
                end                
                legend({'Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 17
                Handel_all{i}(1)=plot(CD_Ms,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_Rs,Ytick1,'g');
                xlabel('$$\Delta\bar{C}$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 18
                Handel_all{i}(1)=plot(CD_Mabs,Ytick1,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_Rabs,Ytick1,'g');
                xlabel('$$\sum|\Delta\bar{C}|$$','interpreter','Latex','FontSize',Font_Size)
                legend({'Comp','React'},'FontSize',Font_Size)
                ylim(Show_z)
            case 19
                Handel_all{i}(j)=plot(S(1:N),Ytick1);
                xlabel('volatile fraction, $$S$$ (-)','interpreter','Latex','FontSize',Font_Size)
                ylim(Show_z)
                xlim([0,1])
            case 20
                Handel_all{i}(j)=plot(u_all(2*N+3:3*N+3),Ytick2,'o-');
                ylim(Show_z)
            case 21
                Handel_all{i}(1)=plot(P(:,1),Ytick2,'b','linewidth',1.4);
                Handel_all{i}(2)=plot(P(:,2),Ytick2,'r','linewidth',1.2);
                Handel_all{i}(3)=plot(P(:,3),Ytick2,'g','linewidth',1.0);
                xlabel('Pressures (Pa)','interpreter','Latex','FontSize',Font_Size)
                legend({'Melt','Solid','Volatile'},'FontSize',Font_Size)
                ylim(Show_z)
            case 22
                Handel_all{i}(4)=plot(-rho(:,3)*g*Dynamic_scaling,Ytick1,'--','color','black');
                plot(-rho(:,1)*g*Dynamic_scaling,Ytick1,'--','color','black')
                plot(-rho(:,2)*g*Dynamic_scaling,Ytick1,'--','color','black')
                Handel_all{i}(1)=plot(dP*Dynamic_scaling,Ytick2,'b','linewidth',1.4);
%                 Handel_all{i}(2)=plot(dP*Dynamic_scaling,Ytick2,'r','linewidth',1.2);
%                 Handel_all{i}(3)=plot(dP*Dynamic_scaling,Ytick2,'g','linewidth',1.0);
                xlabel('Pressures gradient(bar/m)','interpreter','Latex','FontSize',Font_Size)
                ylim(Show_z)
                xlim([-rho(1,2)*g*Dynamic_scaling*1.1 rho(1,2)*g*Dynamic_scaling/15])
            case 23
                Handel_all{i}(1)=plot(rho(:,1).*(phi.*(1-S)>1e-3),Ytick1,'b','linewidth',1.8);
                Handel_all{i}(2)=plot(rho(:,2).*((1-phi)>1e-3),Ytick1,'r','linewidth',1.4);
                Handel_all{i}(3)=plot(rho(:,3).*(phi.*S>1e-3),Ytick1,'g','linewidth',1.2);
                legend({'Melt','Solid','Volatile'},'FontSize',Font_Size)
                xlabel('Density(kg/m$$^3$$)','interpreter','Latex','FontSize',Font_Size)
                ylim(Show_z)
%                 xlim([-rho(1,2)*g*scaling*1.1 rho(1,2)*g*scaling/15])
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
                if N_component==3
                    Cb2=V./(MM+NN+V);
                else
                    Cb2=V./(MM+NN+V+CL);
                end
                handle_saturation=plot(Sat*100,Ytick1,'--','color','black','linewidth',1.8);

                handle_saturation2=plot((S_cap)*100,Ytick1,'--','color',[0.5,0.2,0.7],'linewidth',1.8);

                Handel_all{i}(3)=plot(Cb2*100,Ytick1,'b','linewidth',2,'Marker',Markers{marker_to_use});
                Handel_all{i}(1)=plot(Mass_data(:,1)./sum(Mass_data(:,[1 2]),2),Ytick1,'r','linewidth',1.2,'Marker',Markers{marker_to_use});
                Handel_all{i}(2)=plot(Mass_data(:,4)./sum(Mass_data(:,[4 5]),2),Ytick1,'g','linewidth',1.2,'Marker',Markers{marker_to_use});
                Handel_all{i}(4)=plot((Mass_data(:,7)+max(Mass_data(:,6)-Mass_solid.*S_cap,0))/rho_mean,Ytick1,'color',[255, 200, 100]/255 ,'linewidth',1.2,'Marker',Markers{marker_to_use});

                xlabel('Volatile concentration (H$$_2$$O$$\%$$)','interpreter','Latex','FontSize',Font_Size) %(SiO$$_2 \%$$)
                ylim(Show_z)
                legend([Handel_all{i}(3) Handel_all{i}(1),Handel_all{i}(2),Handel_all{i}(4)],{'$$\bar{C}_{H_2O}$$','$$C_{l2}$$','$$C_{s2}$$','$$C_{v}$$'},'interpreter','Latex','fontsize',Font_Size)
                %                 xlim([-rho(1,2)*g*scaling*1.1 rho(1,2)*g*scaling/15])
            case 99
                Handel_all{i}(1)=plot(sum(Mass_data,2), Ytick1, 'linewidth',2);
                ylim(Show_z)
            case 98
                Handel_all{i}(1)=plot(phi.*(1-S), Ytick1, 'linewidth',2);
                xlim([0 0.5])
                ylim(Show_z)
                xlabel('Melt volume fraction (-)','interpreter','Latex','FontSize',Font_Size)
            case 97
                Handel_all{i}(1)=plot(phi.*S.*(rho(:,3)>1), Ytick1, 'linewidth',2);
                xlabel('volatile volume fraction (-)','interpreter','Latex','FontSize',Font_Size)
                xlim([0 0.2])
                ylim(Show_z)
            case 30
                Cb_cl=CL./(MM+NN+V+CL+CU);
                Handel_all{i}(1)=plot(Cb_cl*100,cellz,'linewidth',2.5);
                Cll=Mass_data(:,8)./sum(Mass_data(:,[1 2 3 8 11]),2)*100;
                Cls=Mass_data(:,9)./sum(Mass_data(:,[4 5 6 9 12]),2)*100;
                Clg=Mass_data(:,10)./sum(Mass_data(:,[7 10 13]),2)*100;
                Handel_all{i}(2)=plot(Cll,Ytick1,'linewidth',2,'visible','off'); %melt
                Handel_all{i}(3)=plot(Cls,Ytick1,'linewidth',1.5,'visible','off'); %solid
                Handel_all{i}(4)=plot(Clg,Ytick1,'linewidth',1.5,'visible','off'); %volatile
                xlabel('Cl (%)','fontsize',Font_Size)
                yticklabels({})
                xlim([0 20])
                ylim(Show_z)    
                ax_extra=axes('Position',get(ha(i),'Position'),'XAxisLocation','top','Color','none','XColor','b','YColor','none');
                hold on;
                x=ones(N,1);%Partition_CL_CU(:,1);
%                 x(not_vol)=nan;
                ylim(Show_z)
                Handel_all{i}(5)=plot(ax_extra,x,Ytick1,'--','linewidth',1.5,'Parent',ax_extra,'color','b'); %D_cl
                legend(ax_extra,[Handel_all{i}],{'$$\bar{C}_{Cl}$$','$$C_{Cl_l}$$','$$C_{Cl_s}$$','$$C_{Cl_g}$$','$$D^{g/l}_{Cl}$$'},'interpreter','Latex','fontsize',Font_Size)
                xlabel('D_{Cl}^{g/l} (-)','fontsize',Font_Size)
                xlim([0 110])
            case 31
                Cb_cu=CU./(MM+NN+V+CL+CU);
                Handel_all{i}(1)=plot(Cb_cu*100,cellz,'linewidth',2.5);
                Cul=Mass_data(:,11)./sum(Mass_data(:,[1 2 3 8 11]),2)*100;
                Cus=Mass_data(:,12)./sum(Mass_data(:,[4 5 6 9 12]),2)*100;
                Cug=Mass_data(:,13)./sum(Mass_data(:,[7 10 13]),2)*100;
                Handel_all{i}(2)=plot(Cul,Ytick1,'linewidth',2,'visible','off'); %melt
                Handel_all{i}(3)=plot(Cus,Ytick1,'linewidth',1.5,'visible','off'); %solid
                Handel_all{i}(4)=plot(Cug,Ytick1,'linewidth',1.5,'visible','off'); %volatile
                xlabel('Cu (%)','fontsize',Font_Size)
                yticklabels({})
                xlim([0 13])
                ylim(Show_z)
                ax_extra=axes('Position',get(ha(i),'Position'),'XAxisLocation','top','Color','none','XColor','b','YColor','none');
                hold on;
                ylim(Show_z)
                x=ones(N,1);%Partition_CL_CU(:,2);
%                 x(not_vol)=nan;
                Handel_all{i}(5)=plot(ax_extra,x,Ytick1,'--','linewidth',1.5,'Parent',ax_extra,'color','b'); %D_cl
                legend(ax_extra,[Handel_all{i}],{'$$\bar{C}_{Cu}$$','$$C_{Cu_l}$$','$$C_{Cu_s}$$','$$C_{Cu_g}$$','$$D^{g/l}_{Cu}$$'},'interpreter','Latex','fontsize',Font_Size)
                xlabel('D_{Cu}^{g/l} (-)','fontsize',Font_Size)
                xlim([20 250])
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
                Handel_all{i}(1)=plot(muf_order,cellz,'linewidth',1.5);
                xlim([0,8])
                grid on
                xlabel('Melt viscosity order (-)','fontsize',Font_Size)
            case 51
                dphidx=zeros(N,3);
                if Contribution_type==2
                    acc_dphidx=zeros(N,3);
                end

                Handel_all{i}(1)=plot(zeros(N,1),cellz,'linewidth',1.6,'color','b');
                Handel_all{i}(2)=plot(zeros(N,1),cellz,'linewidth',1.4,'color','g');
                Handel_all{i}(3)=plot(zeros(N,1),cellz,'linewidth',1.2,'color','r');
                legend({'H','C','T'})
                xlabel('Contribution on \phi (-)','fontsize',Font_Size)
            case 101
                Handel_all{i}(j)=plot(Intrusion_marker_static,Ytick1,'linewidth',2);
                ylim(Show_z)

                xlim([0,1])
                xlabel('Static intrusion marker (-)','interpreter','Latex','FontSize',Font_Size)
            case 102
                Handel_all{i}(j)=plot(Intrusion_marker_dynamics,Ytick1,'linewidth',2);
                ylim(Show_z)

                xlim([0,1])
                xlabel('Dynamic intrusion marker (-)','interpreter','Latex','FontSize',Font_Size)
%                 xlim([-0.5, 0.5])
        end
        if i~=1
            set(gca,'YTickLabel',[]);
        else
            if Z_label_type==1
                ylabel('Distance (m)','fontsize',Font_Size+2)
            else
                ylabel('Depth (km)','fontsize',Font_Size+2)
            end            
        end
    else
        Linewidth=2;
        set(gcf,'CurrentAxes',ha(i)); hold on;set(gca,'tickdir','out') ; set(gca,'fontsize', Font_Size); 
        for j=Plot_configure{i}
            Linewidth=Linewidth-0.5;
            Linewidth=max(Linewidth,0.5);
            switch j
                case 1
                    Handel_all{i}(1)=plot(phi(1:N),Ytick1,'linewidth',Linewidth);
                    Handel_all{i}(2)=plot(0,(cellz(phi_top)-nodez(end))/1000, 'kx', 'MarkerSize', 12);
                    Handel_all{i}(3)=plot(0,(cellz(phi_base)-nodez(end))/1000, 'kx', 'MarkerSize', 12);
                   % Handel_all{i}(4)=plot(1,(cellz(Hphi_top)-nodez(end))/1000, 'k+', 'MarkerSize', 12);
                   % Handel_all{i}(5)=plot(1,(cellz(Hphi_base)-nodez(end))/1000, 'k+', 'MarkerSize', 12);

                    (cellz-nodez(end))/1000;
                    xlabel('Porosity (-)','fontsize',Font_Size,'color','b')
                case 19
                    ax_volatile=axes('Position',get(ha(i),'Position'),'XAxisLocation','top','Color','none','XColor','r','YColor','none');
                    hold on;
                    ylim(Show_z)
                    xlim([0 1])
                    Handel_all{i}(j)=plot(ax_volatile,S(1:N),Ytick1,'linewidth',Linewidth,'Parent',ax_volatile,'color','r');
                    xlabel('Volatile fraction(-)','fontsize',Font_Size,'color','r')
                    
                case 2
                    Handel_all{i}(j)=plot(u_all(1:N+1)/UM_0,Ytick2,'linewidth',Linewidth,'Marker',Markers{marker_to_use});
                case 4
                    HSolidus=ESolidus(0.5,cp);
                    HLiquidus=max(Eliquidus(1,cp,1, Lf,A1, B1, C1, A2, B2, C2, ae),Eliquidus(0,cp,1, Lf,A1, B1, C1, A2, B2, C2, ae));
                    %                 Handel_ent_conserve=text(max(H),0.8*Show_z(2)+0.2*Show_z(1),"1",'fontsize',12);
                    Handel_all{i}(j)=plot(H/(HLiquidus-HSolidus),Ytick1,'linewidth',Linewidth);
                case 5
                    Max_Liquidus=max(Liquidus(1,A1, B1, C1, A2, B2, C2, ae),Liquidus(0,A1, B1, C1, A2, B2, C2, ae));
                    Handel_all{i}(j)=plot((T-Solidus)/(Max_Liquidus-Solidus),Ytick1,'linewidth',Linewidth);
                case 6
                    Cl1l=(74*Mass_data(:,1)+47*Mass_data(:,2))./(sum(Mass_data(:,1:2),2));
                    Handel_all{i}(j)=plot(Cl1l,Ytick1,'linewidth',Linewidth,'color','r','Marker',Markers{marker_to_use});
                    if S_or_Cb==1
%                         TICK=[46, 53, 60, 67, 74];
%                         set(gca,'XTick',TICK);
%                         TICK=[-0.1, 0.2, 0.5, 0.8  1.1];
%                         xticks((TICK-SiO2_range(1))/(SiO2_range(2)-SiO2_range(1)))
%                         TICK=TICK*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1);
%                         set(gca,'xticklabel',{num2str(TICK(1)),num2str(TICK(2)),num2str(TICK(3)),num2str(TICK(4)),num2str(TICK(5))})
                        xlim([42,74])
%                         xticks([min(SiO2_range), SiO2_range(1)*2/3+SiO2_range(2)*1/3, SiO2_range(1)*1/3+SiO2_range(2)*2/3, max(SiO2_range)] )
                    end
                    xlabel(ha(i),'Compositions (SiO$$_2 \%$$)','Interpreter','Latex','FontSize',Font_Size)
                case 7
                    Cs1S=(74*Mass_data(:,4)+47*Mass_data(:,5))./(sum(Mass_data(:,4:5),2));
                    Handel_all{i}(j)=plot(Cs1S,Ytick1,'linewidth',Linewidth,'color','g','Marker',Markers{marker_to_use});
                case 8
                if S_or_Cb==1
%                     Handel_all{i}(j)=plot(M1,Ytick1,'linewidth',Linewidth);
                    Handel_all{i}(j)=plot(Cb,Ytick1,'linewidth',Linewidth*1.3,'color','b','Marker',Markers{marker_to_use});
                else
                    Handel_all{i}(j)=plot(Cb,Ytick1,'linewidth',Linewidth);
                    hold(A_handel, 'on')
                    ax2=axes('position',A_handel.Position,"Color",'none');
                    ax2.XAxisLocation = 'top';
                    ax2.YAxis.TickValues = [];
                    set(ax2,'xcolor','r')
                    xlim(ax2,[43,79])
                    xlabel(ax2,'Bulk $$SiO_2$$ ($$\%$$)','Interpreter','Latex','FontSize',Font_Size)
                    xlabel(A_handel,'Melt fraction (-)','Interpreter','Latex','FontSize',Font_Size)
                    ylim(ax2,Show_z)
                end
                case 9
                    Handel_all{i}(j)=plot(TYPE(1:N)/3,Ytick1,'linewidth',Linewidth);
                case 3
                    Handel_all{i}(j)=plot(u_all(N+2:2*N+2),Ytick2,'linewidth',Linewidth,'Marker',Markers{marker_to_use});
                case 20
                    Handel_all{i}(j)=plot(u_all(2*N+3:3*N+3),Ytick2,'color','g','Marker',Markers{marker_to_use});
                    xlabel(ha(i),'Velocities (m/s)','Interpreter','Latex','FontSize',Font_Size)
            end
            if j==1
                if Z_label_type==1
                    ylabel('Distance (m)','fontsize',Font_Size+2)
                else
                    ylabel('Depth (km)','fontsize',Font_Size+2)
                end
            end
        end
        ylim(ha(i),Show_z)
        if ~(ismembertol(2,Plot_configure{i},1e-5)||ismembertol(6,Plot_configure{i},1e-5))
            xlim(ha(i),[0,1])
        end
        if i~=1
            legend(ha(i),Plot_Label{Plot_configure{i}},'interpreter','Latex','fontsize',Font_Size)
            legend('Location', 'northwest')
        end
    end
    if i~=1
        set(gca,'YTickLabel',[]);
    else
        if Z_label_type==1
            ylabel('Distance (m)','fontsize',Font_Size+2)
        else
            ylabel('Depth (km)','fontsize',Font_Size+2)
        end
    end
end
%%
Background = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',Background)

text(Background,0.17,0.987, Model_name,'fontsize',Font_Size);
Handel_T=text(Background,0.47,0.98, "Time="+num2str(Time/Year),'fontsize',Font_Size);
Handel_dt=text(Background, 0.422,0.95, "dT="+num2str(dt/Year),'fontsize',Font_Size);
if N_component==3
    Handel_N=text(Background, 0.572,0.95, "Nodes="+num2str(N),'fontsize',Font_Size);
else
    Handel_N=text(Background, 0.642,0.96, "Nodes="+num2str(N),'fontsize',Font_Size);
end

for i=1:(Plot_size(1)*Plot_size(2))
    set(ha(i),'YMinorTick','on')
end


