function [av_rho_T, av_Mass_data_T, av_H_T, av_T_T, av_Ts_local_T, av_phi_T, av_S_T, av_cb_T, av_OG_Cb_T] = Averages_evacuated(evacuation_counter, melt_base, melt_top, phi, S, H, Mass_data, ...
                                                                                                    MM, NN, V, CL, CU, Add_CLCU, Lf, cp, simplified_TS, Sys_constant, rho, Melt_density_type, ...
                                                                                                        PD_range, Coef_melt_den, cellz, nodez, aa, bb,cc,rho_mean,g, Velocity_solver_type, Ts, Tl, ...
                                                                                                        Conservation_type, rhom_1, rhom_2, Time, Year,OG_Cb)

% CAB 17/9/25
% Calculation of the averages for evacuations. The averages are calculated
% over the the high melt fraction layer and only if transfer criteria is
% met

%% Average porosity
av_phi_T = mean(phi(melt_base:melt_top))
;


%% Average volatile fraction
av_S_T = mean(S(melt_base:melt_top));


%% Average enthalpy
av_H_T = mean(H(melt_base:melt_top));

%% Average OG_cb
av_OG_Cb_T = mean(OG_Cb(melt_base:melt_top));


%% Average mass data
% Take averages of each part of mass data
% Mass data = [m1 n1 v1 m2 n2 v2 v3] if no chlorine copper
% Mass data = [m1 n1 v1 m2 n2 v2 v3 l1 l2 l3 u1 u2 u3] if chlorine copper

m1_T = mean(Mass_data(melt_base:melt_top,1));
n1_T = mean(Mass_data(melt_base:melt_top,2));
v1_T = mean(Mass_data(melt_base:melt_top,3));
m2_T = mean(Mass_data(melt_base:melt_top,4));
n2_T = mean(Mass_data(melt_base:melt_top,5));
v2_T = mean(Mass_data(melt_base:melt_top,6));
v3_T = mean(Mass_data(melt_base:melt_top,7));

if Add_CLCU==1 
   l1_T = mean(Mass_data(melt_base:melt_top,8));
   l2_T = mean(Mass_data(melt_base:melt_top,9));
   l3_T = mean(Mass_data(melt_base:melt_top,10));
   u1_T = mean(Mass_data(melt_base:melt_top,11));
   u2_T = mean(Mass_data(melt_base:melt_top,12));
   u3_T = mean(Mass_data(melt_base:melt_top,13));
end



% Mass conservation correction
% (m1+m2-M)=0
% n1+n2-N = 0
% v1+v2+v3-V=0

av_M=mean(MM(melt_base:melt_top));
av_N = mean(NN(melt_base:melt_top));
av_V = mean(V(melt_base:melt_top));

if Add_CLCU==1  
    av_CL = mean(CL(melt_base:melt_top));
    av_CU = mean(CU(melt_base:melt_top));
end

% m
err = av_M - m1_T - m2_T;
fraction = m1_T/max(m1_T+m2_T,1e-5);
m1_T = m1_T + err*fraction;
m2_T = m2_T + err*(1-fraction);

% n
err = av_N - n1_T - n2_T;
fraction = n1_T/max(n1_T+n2_T,1e-5);
n1_T = n1_T + err*fraction;
n2_T = n2_T + err*(1-fraction);

%v 
err = av_V - v1_T - v2_T - v3_T;
fraction = v1_T/max(v1_T+v2_T+v3_T,1e-5);
fraction2 = v2_T/max(v1_T+v2_T+v3_T, 1e-5);
v1_T = v1_T + err*fraction;
v2_T = v2_T + err*fraction2;
v3_T = v3_T + err*(1-fraction-fraction2);

if Add_CLCU==1 
    % cl
    err = av_CL - l1_T -l2_T - l3_T;
    fraction = l1_T/max(l1_T+l2_T+l3_T,1e-5);
    fraction2 = l2_T/max(l1_T+l2_T+l3_T, 1e-5);
    l1_T = l1_T + err*fraction;
    l2_T = l2_T + err*fraction2;
    l3_T = l3_T + err*(1-fraction-fraction2);

    % cu
    err = av_CU - u1_T -u2_T - u3_T;
    fraction = u1_T/max(u1_T+u2_T+u3_T,1e-5);
    fraction2 = u2_T/max(u1_T+u2_T+u3_T, 1e-5);
    u1_T = u1_T + err*fraction;
    u2_T = u2_T + err*fraction2;
    u3_T = u3_T + err*(1-fraction-fraction2);

    if v3_T+l3_T+u3_T<0.1
        v2_T=v2_T+v3_T;
        v3_T=0;
        l2_T=l2_T+l3_T;
        l3_T=0;
        u2_T=u2_T+u3_T;
        u3_T=0;
    end


    av_Mass_data_T = max([m1_T,n1_T,v1_T,m2_T,n2_T,v2_T,v3_T,l1_T,l2_T,l3_T,u1_T,u2_T,u3_T],0);
else
    if v3_T<0.1
        v2_T=v2_T+v3_T;
        v3_T=0;
    end
    av_Mass_data_T = max([m1_T,n1_T,v1_T,m2_T,n2_T,v2_T,v3_T],0);
end



%% Average temperature 
% the following must be true: [(m1+n1+v1)*cp1+(m2+n2+v2)*cp2+v3*cp3]*T+(m1*Lf1+n1*Lf2+v1*Lf2)-H=0  
av_T_T = (av_H_T - (m1_T*Lf(1) + n1_T*Lf(2)+v1_T*Lf(2)))/((m1_T+n1_T+v1_T)*cp(1)+(m2_T+n2_T+v2_T)*cp(2)+v3_T*cp(3));



%% Average Ts local
cb_dummy = av_M/(av_M+av_N);

if simplified_TS==1
    av_Ts_local_T=Sys_constant(6,2)*(1-cb_dummy).^Sys_constant(6,1)*(Tl-Ts)+Ts;
else
    av_Ts_local_T=(Sys_constant(6,2)*(1-cb_dummy).^Sys_constant(6,1)+(1-Sys_constant(6,2))*(1-cb_dummy^(1/Sys_constant(6,1))))*(Tl-Ts)+Ts;
end



%% Average density
if Conservation_type==1
    rho_l = mean(rho(melt_base:melt_top),1);
    rho_s = mean(rho(melt_base:melt_top),2);
    rho_v = mean(rho(melt_base:melt_top),3);
 
else




    if Melt_density_type==2
        SiO2=(m1_T*PD_range(2)+n1_T*PD_range(1))./max((m1_T+n1_T),1e-5)/100;
        H2O=v1_T./max(m1_T+n1_T+v1_T,1e-5);
        VMX=H2O.*SiO2*Coef_melt_den(1)+SiO2*Coef_melt_den(2)+H2O*Coef_melt_den(3)+Coef_melt_den(4);
        x_other=1-SiO2-H2O;
        V_all=SiO2.*26.86e-6/0.06009+H2O.*26.27e-6/0.01802+x_other.*VMX;
        rho_l=1./V_all;
    else
        rho_l=(m1_T*rhof_1+n1_T*rhof_2)./max((m1_T+n1_T),1e-5);
        
    end

    rho_s=(m2_T*rhom_1+n2_T*rhom_2)./max(m2_T+n2_T,1e-5);
    

    
    P2=-rho_mean*g*(cellz(melt_top)-nodez(end))/1e5; %in bar
    rho_v=(aa(1)* max(av_T_T,300).^bb(1)+aa(2)*max(P2,400).^bb(2)+aa(3)*max(av_T_T,300).^bb(3).*max(P2,400).^cc(1))*1000;
    
    %if Velocity_solver_type==2
     %   rho_on_nodez(S_on_nodez>0,1)=rho_on_nodez(S_on_nodez>0,3);
    %end

    if Velocity_solver_type==2
        if av_S_T>0
            rho_l = rho_v;
        end
    end

end


av_rho_T = av_phi_T*(1-av_S_T)*rho_l + (1-av_phi_T)*rho_s + av_phi_T*av_S_T*rho_v;


%% Average S_Cap - put here in case we want to use

%av_S_Cap_T = (m2_T*Sys_Constant(5,1)+n2*Sys_Constant(5,2))/max(m2_T+n2_T,1e-2);



%% Outputs
av_cb_T = (PD_range(2)*sum(av_Mass_data_T(:,[1,4]),2)+PD_range(1)*sum(av_Mass_data_T(:,[2,5]),2))./sum(av_Mass_data_T(:,[1 2 4 5]),2);


Filename_intrude = ['magma_avintruded_',num2str(evacuation_counter),'.txt'];
File_In = fopen(Filename_intrude, 'w');

fprintf(File_In,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );

% 1. av_H_T, 2. av_phi_T, 3. av_S_T, 4. av_T_T,  5. av_cb_T, 6. av_rho_T,
% 7. av_Ts_local_T 

fprintf(File_In, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
    'Enthalpy', 'Porosity', 'Vol frac', 'Temp', 'Bulk comp', 'Bulk rho', 'Ts_local');

fprintf(File_In, '%3e \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n', ...
   [av_H_T; av_phi_T; av_S_T; av_T_T; av_cb_T; av_rho_T; av_Ts_local_T]); 

fclose(File_In);





Filename_MD_intrude = ['mass_data_avintruded_',num2str(evacuation_counter),'.txt'];
File_Out_Iem = fopen(Filename_MD_intrude, 'w');

fprintf(File_Out_Iem,'%4s %10.10f %3s \n', 'Time', Time/Year/1000, 'ka' );

if Add_CLCU==1 

    fprintf(File_Out_Iem, '%10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
        'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3', 'l1', 'l2', 'l3', 'u1', 'u2', 'u3');

    
    fprintf(File_Out_Iem, ' %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \n', ...
            [ av_Mass_data_T(:,1);av_Mass_data_T(:,2);av_Mass_data_T(:,3);av_Mass_data_T(:,4);av_Mass_data_T(:,5);av_Mass_data_T(:,6);av_Mass_data_T(:,7); ...
            av_Mass_data_T(:,8);av_Mass_data_T(:,9);av_Mass_data_T(:,10);av_Mass_data_T(:,11);av_Mass_data_T(:,12);av_Mass_data_T(:,13)]);
   



else

    fprintf(File_Out_Iem, ' %10s \t %10s \t %10s \t %10s \t %10s \t %10s \t %10s \n', ...
         'm1', 'n1', 'v1', 'm2', 'n2', 'v2', 'v3');

    
    fprintf(File_Out_Iem, ' %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f  \n', ...
            [av_Mass_data_T(:,1);av_Mass_data_T(:,2);av_Mass_data_T(:,3);av_Mass_data_T(:,4);av_Mass_data_T(:,5);av_Mass_data_T(:,6);av_Mass_data_T(:,7)]);
   

end


fclose(File_Out_Iem);


