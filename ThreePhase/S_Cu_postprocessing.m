


Mass_total=(MM+NN+V+CL);
% Molar masses (g/mol)
M_H2O = 18.015;
M_Cl = 35.45;    % elemental chlorine
M_S = 32.065;


CU_fraction=CU./Mass_total;
SUL_fraction=SUL./Mass_total;
for i=S_CU_to_process
    
    % if i==81
    %     aaa=1
    % end
    m1=Mass_data(i,1); n1=Mass_data(i,2); v1=Mass_data(i,3);
    m2=Mass_data(i,4); n2=Mass_data(i,5); v2=Mass_data(i,6); v3=Mass_data(i,7);
    l1=Mass_data(i,8); l2=Mass_data(i,9); l3=Mass_data(i,10);

    if m1+n1>0 %melt phase exist
        SiO2=(PD_range(2)-PD_range(1))*m1/(m1+n1)+PD_range(1);
        H2O=v1/(m1+n1+v1)*100;
    else
        % cl_ghost=(-B1-sqrt(B1^2-4*A1*(C1-Ts_local(i))))/2/A1;
        cl_ghost=MM(i)/(MM(i)+NN(i));
        SiO2=(PD_range(2)-PD_range(1))*cl_ghost+PD_range(1);
        H2O=min(V(i)/(MM(i)+NN(i)+V(i))/K0,Saturation(i))*100;
    end

    % 'SiO2' 'TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O' 'H2O'
    wt_percent=zeros(1,10);
    wt_percent(1)=SiO2;
    wt_percent(10)=H2O;
    for j=1:8
        wt_percent(1+j)    =polyval(Oxides_coefficients(j,:), SiO2);
    end
    % rescale back to 100% and keep SiO2 and H2O percent the same
    Target=100-SiO2-H2O;
    ratio=Target/sum(wt_percent(2:9));
    wt_percent(2:9)=wt_percent(2:9)*ratio;
    
    % calculate_S_redox(Na2O, MgO, Al2O3, SiO2, K2O, CaO, TiO2, MnO, FeOt, Fe3_over_Fe_total, T_input_C, S_ppm, mode, Pressure, Ni, Cu)
    % Ni assumed to be 100ppm Cu in this step assumed to be 80ppm, these only affects aFeS whose contribution is much smaller than other terms
    Fe3ratio=0.5;
    results=calculate_S_redox(wt_percent(8), wt_percent(6), wt_percent(3), wt_percent(1), wt_percent(9), wt_percent(7), wt_percent(2), wt_percent(5), wt_percent(4),    Fe3ratio, T(i), 0, 1, Pg_real(i)/10, 100, 80);

    dFMQ=results.DeltaQFM;
    D_s=40000/(238+1701/(1+10^(0.78-2*(dFMQ-1)))); %reference
    S6_fraction=results.S6_over_SumS; %ratio of S6+/S2-
    SCSS=exp(results.ln_SCSS)*1e-6;

   
    %
    p1=(m1+n1+v1+l1)/(MM(i)+NN(i)+V(i)+CL(i));
    p2=(m2+n2+v2+l2)/(MM(i)+NN(i)+V(i)+CL(i));
    p3=(v3+l3)/(MM(i)+NN(i)+V(i)+CL(i));

    if p1>0
        Sulfide_exist=0;
        S1=SUL_fraction(i)/(p3*D_s+p1);
        
        if S1*(1-S6_fraction)>SCSS
            S1=SCSS/(1-S6_fraction);           
            Sulfide_exist=1;
        end
    else
        S1=0;
    end
    S3=S1*D_s;
    S2=(SUL_fraction(i)-p1*S1-p3*S3)/max(p2, 1e-12);

    % if p1>0
    %     if S1*(1-S6_fraction)>SCSS
    %         S1=SCSS/(1-S6_fraction);
    %         S3=S1*D_s;
    % 
    % 
    %         Sulfide_exist=1;
    %     end
    % 
    % else
    %     S1=0;
    %     S2=SUL(i)/(m2+n2+v2+l2);
    %     S3=0;
    % end

    % obtain the mole fraction of H2O, Nacl_eq, and H2S for Copper partition coefficient
    if v3>0
        Total_m=(v3+l3)*(1+S3);
        w_H2O=v3/Total_m;
        w_Cl=l3/Total_m;
        W_S=S3;
    
        total_moles=w_H2O/M_H2O+w_Cl/M_Cl+S3/M_S;
        X_Cl=w_Cl/M_Cl/total_moles;
        X_H2O=w_H2O/M_H2O/total_moles;
        X_H2S=S3/M_S/total_moles*(1-S6_fraction);
    else
        X_Cl=0;X_H2S=0;X_H2O=0;
    end

    X_FeO=results.X_Fe2;
    T2=T(i)+273.15; %Temperature in K

    % Cu v/m partition from (Tattitch et al 2017)
    Dcu1=8e4*X_Cl^2*X_H2O^14*(1+180*X_H2S)+380*X_Cl+0.8;

    Dcu2=exp(1.13+0.39*1e4/T2+0.02*dFMQ-0.86*log(X_FeO));

    
    if p1>0
        if Sulfide_exist==0
            u2=0;
            u1=CU_fraction(i)/(p3*Dcu1+p1);
            u3=u1*Dcu1;
        else
            u1=CU_fraction(i)/(p3*Dcu1+p1+ p2*S2/0.35*Dcu2);
            u3=u1*Dcu1;
            u2=u1*Dcu2*S2/0.35;
        end
    else
        u1=0;
        u2=CU_fraction(i)/(p2+p3*Dcu1/Dcu2);
        u3=u2*Dcu1/Dcu2;
    end
    Mass_data_addition(i,:)=[S1*p1*Mass_total(i),S2*p2*Mass_total(i),S3*p3*Mass_total(i),u1*p1*Mass_total(i),u2*p2*Mass_total(i),u3*p3*Mass_total(i)];
end