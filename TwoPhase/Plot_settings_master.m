% HH 2/2/2020
% Initializing all the plots
%Plot_Label={'Melt fraction, $$\phi (-)$$','uf','um','Enthalpy','Temperature($$^\circ$$C)','Cl','Cs','$$\bar{C}$$','St', 'Solid_vis(%)','$$\Delta\phi$$','$$\sum|\Delta\phi|$$', 'Cb', '$$\sum|\Delta\phi|$$', '$$\sum\Delta\phi$$','$$\sum\Delta\bar{C}$$','$$\Delta\bar{C}$$','$$\sum|\Delta\bar{C}|$$'};
Font_Size=13;
% YLIMITS=[crust_bottom crust_top];
YLIMITS=[-Injection_depth-Sill_length/1000-0.1, -Injection_depth+0.1];
%YLIMITS=[-35 -20];

Step_counts=0;
% Initialize all the plottings
Plot_size=size(Plot_configure);
Handel_all=Plot_configure;
% % this part is to scaling the velocity.
PHI0=0.5;  
K0=PHI0^C_value_B*grain_size^2/C_value_A;
UM_0=K0/10^((mu_f1+mu_f2)/2)*(1-PHI0)*(rhom_1-rhof_2)*g*scaling_factor;   %dimensionless speed scaling

% Cb=Cm.*phi(N+1:2*N)+Cf.*phi(1:N);

S_or_Cb=1;  %1 for showing SiO2

figure(13)
clf; set(gcf,'color','w','position',[10,10,1650,1200]);
% title(['Compaction Length=',num2str(Compaction_length)])
Plot_configure=reshape(Plot_configure',[1,Plot_size(1)*Plot_size(2)]);

for i=1:(Plot_size(1)*Plot_size(2))
    A_handel=subplot(Plot_size(1),Plot_size(2),i);hold on
    if i==1 
%         Handel_T=text(0.65,0.72*Show_z(2)+0.28*Show_z(1), "T="+num2str(Time/Year),'fontsize',Font_Size-3);%
%         aaa=title(["Time=" num2str(Time/Year) 'y'],'fontsize',Font_Size-3);
%         Handel_dt=text(0.65,0.76*Show_z(2)+0.24*Show_z(1), "dT="+num2str(dt/Year),'fontsize',Font_Size-3);
%         Handel_CDN=text(0.65,0.70*Show_z(2)+0.30*Show_z(1), "CDN="+num2str(round(CDN,3)),'fontsize',Font_Size-3);
    end
    if i==2 
        title([ ],'fontsize',Font_Size-3);
%         Handel_dt=text(0.65,0.76*Show_z(2)+0.24*Show_z(1), "dT="+num2str(dt/Year),'fontsize',Font_Size-3);
%         Handel_CDN=text(0.65,0.70*Show_z(2)+0.30*Show_z(1), "CDN="+num2str(round(CDN,3)),'fontsize',Font_Size-3);
    end
    if length(Plot_configure{i})==1
        j=Plot_configure{i};
        switch j
            case 1
                box off
                Handel_all{i}(j)=plot(phi(1:N),cellz/1000-Base_crust,'b','linewidth',2.5);
%                 Handel_all{i}(j)=plot(phi(1:N),cellz,'-o','markersize',3);
                %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                 xlim([0,1])
                ylim(YLIMITS);
                xlabel('Melt fraction (-)')
                if Has_volatile==1
                    ax_volatile=axes('Position',get(gca,'Position'),'XAxisLocation','top','Color','none','XColor','r','YColor','none');
                    hold on;
                    ylim(YLIMITS)
                    xlim([0 0.15])
                    Handel_all{i}(2)=plot(ax_volatile,S,cellz/1000-Base_crust,'linewidth',1 ,'Parent',ax_volatile,'color','r');
                    xlabel('Volatile fraction(-)','fontsize',Font_Size,'color','r')
                end

                
                
               
            case 2
                Handel_all{i}(j)=plot(u_all(1:N+1),nodez/1000-Base_crust);
                ylim(YLIMITS);
            case 4
                HSolidus=Ts.*cp;
                Handel_all{i}(j)=plot(H,cellz/1000-Base_crust);
                if Phase_Type==0
                    plot(HSolidus*ones(N,1),cellz,'--','color','red')
                else
                    for i=1:N
                        [~,~,HSolidus(i),HLiquidus(i)]=solid_state_Ts_Tl_local(Cb(i),alpha_order,n_order,Lf,cp,A1,B1,C1);
                    end
                    plot(HSolidus*ones(N,1),cellz/1000-Base_crust,'--','color','red')
                    plot(HLiquidus*ones(N,1),cellz/1000-Base_crust,'--','color','red')

                end
                Handel_ent_conserve=text(max(H),0.8*Show_z(2)+0.2*Show_z(1),"1",'fontsize',Font_Size);  
%                 ylim([0,max(nodez)]);
                ylim(YLIMITS)
                if Phase_Type==0
                    %                 xlim([0,max(H)*1.5]);
                    gapH=max([max(H)-min(H),abs(max(H)-Ts*cp),100]);
                    xlim([min(min(H)-1e-5-gapH/20,Ts*cp-gapH/20-1e-5),max(max(H)+1e-5+gapH/2,Ts*cp+gapH/3+1e-5)])
                else
                    [Max_Cb_0_S,Max_Cb_0_L,~,~]=solid_state_Ts_Tl_local(0,alpha_order,n_order,Lf,cp,A1,B1,C1);
                    gapH=max([max(H)-min(H),abs(max(H)-Max_Cb_0_L),100]);
                    xlim([min(H)-gapH,max(H)+gapH])
                end

              case 5
                  box on
                Handel_all{i}(1)=plot(T,cellz/1000-Base_crust,'b','linewidth',2.5);
                Handel_all{i}(2)=plot(Ts0,cellz/1000-Base_crust,'r:','linewidth',2.5);
               Handel_all{i}(3)=plot(Tl0,cellz/1000-Base_crust,'r','linewidth',2.5);
     
                
                %if Phase_Type==0
                    
                   % gapT=max(max(T)-min(T),abs(max(T)-Ts));
%         
                  %  xlim([min(T)-10 max(max(T)+10,C1)])

                %else       

                    %[Max_Cb_1_S,Max_Cb_1_L,~,~]=solid_state_Ts_Tl_local(0,min_cb_ts,ssm1,ssc1,min_ts,liq_k,liq_a,liq_b,cp,Lf);
                    %[Max_Cb_0_S,Max_Cb_0_L,~,~]=solid_state_Ts_Tl_local(1,min_cb_ts,ssm1,ssc1,min_ts,liq_k,liq_a,liq_b,cp,Lf);

                    %MAX=max(Max_Cb_0_L,Max_Cb_1_L);
                    %MIN=min(Max_Cb_0_S,Max_Cb_1_S);

                   % Gap=MAX-MIN;
                   % xlim([0,MAX+Gap])
                %end

                
                ylim(YLIMITS)
                xlabel('Temperature (^oC)')



            case 6
                Handel_all{i}(j)=plot(C_all(1:N),cellz/1000-Base_crust,'color','red');
                %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS);
                xlim([0,1]);
            case 7
                Handel_all{i}(j)=plot(C_all(N+1:2*N),cellz/1000-Base_crust,'color','green');
                %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS);
                xlim([0,1]);
            %% catherine added
            case 71 % olivine
                Handel_all{i}(j)=plot(ol,cellz/1000-Base_crust,'g','linewidth',2.5);
                ylim(YLIMTS);
                xlim([0,1]);
            case 72 %pyroxene
                Handel_all{i}(j)=plot(px,cellz/1000-Base_crust,'r','linewidth',2.5);
                ylim(YLIMITS);
                xlim([0,1])
            
            case 8                
                %Handel_com_conserve=text(0.7,0.8*Show_z(2)+0.2*Show_z(1),"1",'fontsize',Font_Size);
                if S_or_Cb==0
                    Handel_all{i}(j)=plot(Cb*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1),cellz/1000-Base_crust,'color','blue');
                else
                    
                    Handel_all{i}(j)=plot(Cb*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1),cellz,'color','blue');
                end
                %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS);
               % xlim([0,1]);
            case 9
                Handel_all{i}(j)=plot(TYPE(1:N),cellz/1000-Base_crust);
                %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS);
                xlim([-0.3,3.3]);
            case 3
                Handel_all{i}(j)=plot(u_all(N+2:2*N+2),nodez/1000-Base_crust);
               % yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS); 
            case 10
                Handel_all{i}(1)=plot(Force_coupling,nodez/1000-Base_crust,'r');
                set(Handel_all{i}(1),'visible','off');
                Handel_all{i}(2)=plot(Force_solid,nodez/1000-Base_crust,'g');
                Handel_all{i}(3)=plot(Force_liquid,nodez/1000-Base_crust,'b');
                ylim(YLIMITS);
                xlim([0,100]);
            case 11
                Handel_all{i}(1)=plot(Con_T,cellz/1000-Base_crust,'b');
                Handel_all{i}(2)=plot(Con_M,cellz/1000-Base_crust,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(Con_C,cellz/1000-Base_crust,'g');
                xlabel('$$\Delta\phi$$','interpreter','Latex')
                ylabel('z (m)')
                legend({'Heat','Comp','React'})
                ylim(YLIMITS);
            case 12
                Handel_all{i}(1)=plot(ACon_T,cellz,'b');
                Handel_all{i}(2)=plot(ACon_M,ne,cellz,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(ACon_C,cellz,'g');
                xlabel('$$\sum|\Delta\phi|$$','interpreter','Latex')
                ylabel('z (m)')
                legend({'Heat','Comp','React'})
                ylim(YLIMITS);
            case 13                
                refined_Cb  =interp1(cellz,Cb,cell_refined,'linear','extrap');
                Handel_com_conserve=text(0.7,0.8*Length,"1",'fontsize',Font_Size);
                Handel_all{i}(j)=plot(refined_Cb,cell_refined);
                xlim([0,1]);
            case 14
                refined_ACon_T=interp1(cellz,ACon_T,cell_refined,'linear','extrap');
                refined_ACon_M=interp1(cellz,ACon_M,cell_refined,'linear','extrap');
                refined_ACon_C=interp1(cellz,ACon_C,cell_refined,'linear','extrap');                
                Handel_all{i}(1)=plot(refined_ACon_T,cell_refined,'b');
                Handel_all{i}(2)=plot(refined_ACon_M,cell_refined,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(refined_ACon_C,cell_refined,'g');
                xlabel('$$\sum|\Delta\phi|$$','interpreter','Latex')
                ylabel('z (km)')
                legend({'Heat','Comp','React'})    
                ylim(YLIMITS);
            case 15
                Handel_all{i}(1)=plot(ACon_T,cellz/1000-Base_crust,'b');
                Handel_all{i}(2)=plot(ACon_M,cellz/1000-Base_crust,'r','linewidth',1.5);
                Handel_all{i}(3)=plot(ACon_C,cellz/1000-Base_crust,'g');
                xlabel('$$\sum\Delta\phi$$','interpreter','Latex')
                ylabel('z (km)')
                legend({'Heat','Comp','React'})
                ylim(YLIMITS);
            case 16
                Handel_all{i}(1)=plot(CD_M,cellz,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_R,cellz,'g');
                if S_or_Cb==1
                    if FourMPD==1
                        xlabel('$$\sum\Delta\bar{C}$$ (MgO$$\%$$)','interpreter','Latex') 
                    else
                        xlabel('$$\sum\Delta\bar{C}$$ (SiO$$_2 \%$$)','interpreter','Latex') 
                    end
                else
                    xlabel('$$\sum\Delta\bar{C}$$ (-)','interpreter','Latex')
                end                
                ylabel('z (km)')
                legend({'Comp','React'})
                ylim(YLIMITS);
            case 17
                Handel_all{i}(1)=plot(CD_Ms,cellz,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_Rs,cellz,'g');
                xlabel('$$\Delta\bar{C}$$','interpreter','Latex')
                ylabel('z (km)')
                legend({'Comp','React'})
                ylim(YLIMITS);
            case 18
                Handel_all{i}(1)=plot(CD_Mabs,cellz,'r','linewidth',1.5);
                Handel_all{i}(2)=plot(CD_Rabs,cellz,'g');
                xlabel('$$\sum|\Delta\bar{C}|$$','interpreter','Latex')
                ylabel('z (km)')
                legend({'Comp','React'})
                ylim(YLIMITS);
            case 19
                rhof=rhof_2*C_all(1:N)+(1-C_all(1:N))*rhof_1;
                rhom=rhom_2*C_all(N+1:2*N)+(1-C_all(N+1:2*N))*rhom_1;
                Handel_all{i}(1)=plot(dP,nodez,'b','linewidth',1.5);
                Handel_all{i}(2)=plot(rhof*g*scaling_factor,cellz,'--','color','r','linewidth',1.5);
                Handel_all{i}(3)=plot(rhom*g*scaling_factor,cellz,'--','color','g','linewidth',1.5);
                xlabel('Pressure gradient (Pa/m)','interpreter','Latex')
                ylabel('z (km)')
                legend({'-dP','$$\rho_f g$$','$$\rho_m g$$'},'interpreter','Latex')
%                 xlim([rhof_2*g*scaling_factor rhom_1*g*scaling_factor])
                xlim([min(rhof)*0.9*g*scaling_factor max(rhom)*1.2*g*scaling_factor])
               % yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                ylim(YLIMITS);
            case 1001
                Handel_all{i}(1)=plot(Cb2,cellz/1000-Base_crust,'b','linewidth',1.5);
                Handel_all{i}(2)=plot(Cl2,cellz/1000-Base_crust,'r','linewidth',1.2);
                Handel_all{i}(3)=plot(Cs2,cellz/1000-Base_crust,'green','linewidth',1.2);
                
                Handel_all{i}(4)=plot(Lsaturation,cellz/1000-Base_crust,'--','color','r','linewidth',1.2);
                Handel_all{i}(5)=plot(Ssaturation,cellz/1000-Base_crust,'--','color','green','linewidth',1.2);
                ylim(YLIMITS);
                xlim([0 0.15])
        end
        if S_or_Cb~=10 && j~=19 && j~=20
            %xlabel(Plot_Label{j},'interpreter','Latex')
        end
        ylabel('z (km)')
    else
        Linewidth=2;
        for j=Plot_configure{i}
            Linewidth=Linewidth-0.5;
            Linewidth=max(Linewidth,0.5);
            switch j
                case 1
                    box on
                    Handel_all{i}(j)=plot(phi(1:N),cellz/1000-Base_crust,'linewidth',Linewidth);
                case 2
                    Handel_all{i}(j)=plot(u_all(1:N+1)/UM_0,nodez/1000-Base_crust,'linewidth',Linewidth,'color','r');
                   % yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                case 4
                    HSolidus=ESolidus(0.5,cp);
                    HLiquidus=max(Eliquidus(1,cp,1, Lf,A1, B1, C1, A2, B2, C2, ae),Eliquidus(0,cp,1, Lf,A1, B1, C1, A2, B2, C2, ae));
                    %                 Handel_ent_conserve=text(max(H),0.8*Show_z(2)+0.2*Show_z(1),"1",'fontsize',12);
                    Handel_all{i}(j)=plot(H/(HLiquidus-HSolidus),cellz/1000-Base_crust,'linewidth',Linewidth);


                case 6 %liquid sio2
                    box on
                    if FourMPD==1
                        
                        Handel_all{i}(j)=plot(((1-C_all(1:N))*(MgO_range(2)-MgO_range(1))+MgO_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','r');
                        if S_or_Cb==1
                            xlim([MgO_range(1),MgO_range(2)])
                        end
                        ylim(YLIMITS);
                        xlabel(A_handel,'Compositions (MgO$$\%$$)','Interpreter','Latex')

                    else

                        Handel_all{i}(j)=plot(((1-C_all(1:N))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','r');
                        if S_or_Cb==1
                            %TICK=[76, 69, 61, 53,45];
                            %[45, 53, 61, 69, 76];
    %                         TICK=[-0.1, 0.2, 0.5, 0.8  1.1];
                            %xticks(1-((TICK)-SiO2_range(1))/(SiO2_range(2)-SiO2_range(1)))
    %                         TICK=TICK*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1);
                            %set(gca,'xticklabel',{num2str(TICK(1)),num2str(TICK(2)),num2str(TICK(3)),num2str(TICK(4)),num2str(TICK(5))})
                            xlim([SiO2_range(1)-10,SiO2_range(2)+10])
    %                         xticks([min(SiO2_range), SiO2_range(1)*2/3+SiO2_range(2)*1/3, SiO2_range(1)*1/3+SiO2_range(2)*2/3, max(SiO2_range)] )
                        end
    %                     xlabel(A_handel,'Liquid and Solid composition (-)','Interpreter','Latex')
                        %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                        ylim(YLIMITS);
                        xlabel(A_handel,'Compositions (SiO$$_2 \%$$)','Interpreter','Latex')

                    end

                case 7 %solid
                    if FourMPD==1
                         Handel_all{i}(j)=plot(((1-C_all(N+1:2*N))*(MgO_range(2)-MgO_range(1))+MgO_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','g');
                    else
                         Handel_all{i}(j)=plot(((1-C_all(N+1:2*N))*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','g');
                    end
                    
                   
                    %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                case 8 %bulk sio2
                if S_or_Cb==1
                    if FourMPD==1
                        Handel_all{i}(j)=plot((Cb*(MgO_range(2)-MgO_range(1))+MgO_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','b');
                         ylim(YLIMITS);
                        xlabel('MgO (%)')
                        xlim([MgO_range(1),MgO_range(2)])
                    else
                        Handel_all{i}(j)=plot((Cb*(SiO2_range(2)-SiO2_range(1))+SiO2_range(1)),cellz/1000-Base_crust,'linewidth',2.5,'color','b');
                        ylim(YLIMITS);
                        xlabel('SiO2 (%)')
                        xlim([SiO2_range(1)-10,SiO2_range(2)+10])
                    end

                else
                    if FourMPD==1
                        Handel_all{i}(j)=plot(Cb,cellz/1000-Base_crust,'linewidth',Linewidth,'color','b');
                        hold(A_handel, 'on')
                        ax2=axes('position',A_handel.Position,"Color",'none');
                        ax2.XAxisLocation = 'top';
                        ax2.YAxis.TickValues = [];
                        set(ax2,'xcolor','r')
                        xlim(ax2,[43,79])
                        xlabel(ax2,'Bulk $$SiO_2$$ ($$\%$$)','Interpreter','Latex')
                        xlabel('Melt fraction (-)')
                       % yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                        ylim(YLIMITS);
                    else
                        Handel_all{i}(j)=plot(Cb,cellz/1000-Base_crust,'linewidth',Linewidth,'color','b');
                        hold(A_handel, 'on')
                        ax2=axes('position',A_handel.Position,"Color",'none');
                        ax2.XAxisLocation = 'top';
                        ax2.YAxis.TickValues = [];
                        set(ax2,'xcolor','r')
                        xlim(ax2,[43,79])
                        xlabel(ax2,'Bulk $$MgO$$ ($$\%$$)','Interpreter','Latex')
                        xlabel('Melt fraction (-)')
                       % yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
                        ylim(YLIMITS);
                    end

                end
              
               

                
                case 9
                    Handel_all{i}(j)=plot(TYPE(1:N)/3,cellz/1000-Base_crust,'linewidth',Linewidth);
                case 3
                    Handel_all{i}(j)=plot(u_all(N+2:2*N+2),nodez/1000-Base_crust,'linewidth',Linewidth,'color','g');
                    %yticks([-35 -32.5 -30 -27.5 -25 -22.5 -20])
            end
        end
        ylim(YLIMITS);
        %if ~(ismembertol(2,Plot_configure{i},1e-5)||ismembertol(6,Plot_configure{i},1e-5))
         %   xlim(A_handel,[0,1])
        %end
        %legend(A_handel,Plot_Label{Plot_configure{i}},'interpreter','Latex')
        ylabel('z (km)')
    end
end
%%
Background = axes('Position',[0 0.9 1 1],'Visible','off');
TITLE=text(Background,0.4,0.06, "Time="+ num2str(Time/Year/1000) +'ka   '+ "N\_nodes="+ num2str(N),'fontsize',Font_Size-2,'FontWeight','bold');

%if Create_video==1
%Current_frame=getframe(13);
%writeVideo(Video_handel,Current_frame);
%end


