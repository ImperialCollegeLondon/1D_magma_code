%% Finding_porosity
% The script contains the following sections:

    % Saving previous porosity

    % Initialising arrays

    % Find porosity

    % Find buoyant magma in the areas of porosity.
        % (in this section, we find the shallowest buoyant magma layer)

    % Calculate the thickness of buoyant magma layer
    
    % Compare current buoyant porosity to previous buoyant porosity

    % Calculate density of magma layer

    % Calculate density of high melt fraction portion of the magma layer

    % Calculate hRTI

    % Calculate overpressure

    % Calculate thermal death crtierion

    % Decide the critical overpressure 



%% Saving previous porosity

%old_buoy_phi_top = buoy_phi_top;
%old_buoy_phi_base = buoy_phi_base;
%old_buoy_Hphi_top = buoy_Hphi_top;
%old_buoy_Hphi_base = buoy_Hphi_base;
old_t2=t2;
old_hrti=hrti;

old_all_phi_buoy_top = all_phi_buoy_top;
old_all_phi_buoy_base = all_phi_buoy_base;

%% Initialising arrays
phi_top = []; % top of porosity
phi_base = []; % base of porosity

Hphi_top = []; % top of high (>CMF) porosity
Hphi_base =[]; % base of high (>CMF) porosity

buoy_phi_base = []; % base of buoyant porosity
buoy_phi_top = []; % top of buoyant porosity

buoy_Hphi_base=[]; % base of high (>CMF) buoyant layer of porosity
buoy_Hphi_top=[]; % top of high (>CMF) buoyant layer of porosity

ov_rho_bulk = []; % density of the crust overlying porosity
ov_rho_bulk_buoy=[]; % density of the crust overlying buoyant porosity 
Av_bulk_density  = []; % Average density of the buoyant porosity layer
Av_Hphi_bulk_density  = []; % Average density of the high melt fraction (>CMF) buoyant porosity layer.

all_phi_buoy_top =[];
all_phi_buoy_base=[];

hrti_flag = [];
hrti_start = []; % the previous hrti thickness. If the layer has just developed, this is equal to the initial starting hrti.
prev_buoy_layer = [];

t2  = []; %time since hrti has developed (confined time)
hrti  = []; % thickness of the rti
hb = []; %thickness of the buoyant magma layer

crit_thermdeath = []; %thermal death criterion value
kt_layer = []; % thermal conductivity of high porosity layer
latent_layer = []; % latent heat of high porosity layer
cp_layer = []; % specific heat capactity of high porosity layer
kappa_layer = []; % thermal diffusivity of high porosity layer
muf_order_layer = []; % melt viscosity of high porosity layer

crit_overpressure  = []; %critical overpressure 


%Indexing for array
phi_top_i=0; 
phi_base_i=0;
Hphi_top_i=0;
Hphi_base_i=0;
buoy_phi_top_i=0;
buoy_phi_base_i=0;
buoy_Hphi_top_i=0;
buoy_Hphi_base_i=0;



%% EDIT THIS OUT - this will become an input parameter 



%% Find porosity

% Top from base - finding porosity first
for i=N-1:-1:1

    if phi(i)>0 && phi(i+1)==0
        phi_top_i = phi_top_i+1;
        phi_top(phi_top_i) = i;
    end

    if phi(i)>0 && phi(i-1)==0
        phi_base_i = phi_base_i+1;
        phi_base(phi_base_i)=i;
    end

end


% For the areas of porosity, find the highest layer of high porosity 
for j=1:length(phi_top)

    Hphi_dummy=0;
    
    for i=phi_top(j):-1:phi_base(j)
        if phi(i)>=CMF && phi(i+1)<CMF
            Hphi_top_i = Hphi_top_i+1;
            %Hphi_top(Hphi_top_i) = i;
            Hphi_top(j) = i;
            Hphi_dummy=1;
            break
        end
    end

    if Hphi_dummy==1
        for i=Hphi_top(j):-1:phi_base(j)
            if phi(i)>=CMF && phi(i-1)<CMF
                Hphi_base_i = Hphi_base_i+1;
                %Hphi_base(Hphi_base_i)=i;
                Hphi_base(j)=i;
                break
            end
        end
    end

    if Hphi_dummy==0
        Hphi_top(j)=-5000;
        Hphi_base(j)=-5000;
    end

end



%% Find buoyant magma in the areas of porosity (the shallowest layer of magma which is buoyant)


if ~isempty(phi_top)
    ov_rho_bulk = zeros(length(phi_top),1);



    for j=1:length(phi_top)
        no_buoy=0; % reset everytime.
    
        %find density of the overlying crust.
        Tot_ov_rho_bulk=0;

        % calculate amount above
        
        for i=phi_top(j)+1:1:N
            find_end_dz = cellz(i)-cellz(phi_top(j));%sum(dz(phi_top(j)+1:i));
            if find_end_dz >= xabove
                xaboven = i;
                break
            end
        end




       % for i=phi_top(j)+1:1:xaboven   %phi_top(j)+xaboven
        %    Tot_ov_rho_bulk=Tot_ov_rho_bulk+rho_bulk(i);
        %end
    
        ov_rho_bulk(j) = mean(rho_bulk(phi_top(j)+1:xaboven)); %Tot_ov_rho_bulk/(abs(xaboven - (phi_top(j)+1))); %overlying crust density 
    
        for i=phi_top(j):-1:phi_base(j)
    
            if ov_rho_bulk(j)<=rho_bulk(i)
                %Tot_ov_rho_bulk=Tot_ov_rho_bulk+rho_bulk(i);
                ov_rho_bulk(j) = mean(rho_bulk(i:xaboven)); %Tot_ov_rho_bulk/(abs(xaboven-i));
            else
                buoy_phi_top_i = buoy_phi_top_i+1;
                buoy_phi_top(buoy_phi_top_i) = i;
                all_phi_buoy_top(buoy_phi_top_i) = phi_top(j);
                all_phi_buoy_base(buoy_phi_top_i) = phi_base(j);
                break
            end
    
            if i==phi_base(j) % this means there is no buoyant magma presented in this layer.
                no_buoy=1;
            end
    
        end
    
        if no_buoy == 0
    
            ov_rho_bulk_buoy(buoy_phi_top_i) = ov_rho_bulk(j);
            
            % find base of buoyant layer 
    
            for i = buoy_phi_top(buoy_phi_top_i): -1: phi_base(j)
                
                if rho_bulk(i)<ov_rho_bulk_buoy(buoy_phi_top_i) && rho_bulk(i-1)>= ov_rho_bulk_buoy(buoy_phi_top_i)
                    buoy_phi_base(buoy_phi_top_i) =i;
                    break
                end

                if i==phi_base(j) && rho_bulk(i)<ov_rho_bulk_buoy(buoy_phi_top_i)
                    buoy_phi_base(buoy_phi_top_i)=i;
                end
    
            end
    
    
    
            % find top of high melt fraction (set as -5000 if there is none)
    
            for i = buoy_phi_top(buoy_phi_top_i):-1:buoy_phi_base(buoy_phi_top_i)
    
                if phi(i)>=CMF && phi(i+1)<CMF
                    buoy_Hphi_top(buoy_phi_top_i)=i;
                    break
                end
    
                if i == buoy_phi_base(buoy_phi_top_i)
                    buoy_Hphi_top(buoy_phi_top_i)=-5000;
                end
    
            end 
    
    
            % find base of buoyant high porosity (set as -5000 if there is
            % none)

    
            if buoy_Hphi_top(buoy_phi_top_i)==-5000
                buoy_Hphi_base(buoy_phi_top_i)=-5000;
            else
                for i = buoy_Hphi_top(buoy_phi_top_i):-1:buoy_phi_base(buoy_phi_top_i)
                    if phi(i)>=CMF && phi(i-1)<CMF
                        buoy_Hphi_base(buoy_phi_top_i)=i;
                        break
                    end
                end
            end
            
        end
    end
end

% HH
index_error=find(buoy_Hphi_base==0);
buoy_Hphi_base(index_error)=-5000;
buoy_Hphi_top(index_error)=-5000;
    %% calculate thickness of buoyant magma 
if ~isempty(buoy_phi_top)
    buoyant_phi_m = abs(cellz(buoy_phi_top-1) - cellz(buoy_phi_base)); %thickness in m
    buoyant_phi = abs(buoy_phi_top - buoy_phi_base+1); % thickness in index
    %buoyant_Hphi_m  = abs(cellz(buoy_Hphi_top) - cellz(buoy_Hphi_base)); % thickness in m
    buoyant_Hphi = abs(buoy_Hphi_top - buoy_Hphi_base+1); %thickness in index

    for i=1:length(buoy_phi_top)
        if buoy_Hphi_top(i)==-5000
            buoyant_Hphi_m(i)=0;
        else
            if buoy_Hphi_top(i)-1>0  && buoy_Hphi_base(i)>0
                buoyant_Hphi_m(i)  = abs(cellz(buoy_Hphi_top(i)-1) - cellz(buoy_Hphi_base(i))); % thickness in m    
            else
                buoyant_Hphi_m(i)=0;
            end
        end
    end
    
    hb = buoyant_phi_m; 
    
    
    %% Compare current buoyant porosity to previous buoyant porosity.
    hrti_flag = zeros(length(buoy_phi_top),1);
    hrti_start = zeros(length(buoy_phi_top),1);
    prev_buoy_layer = zeros(length(buoy_phi_top),1);

    

    
    for i=1:length(buoy_phi_top)
        track_hrti=0;
        for j=1:length(old_all_phi_buoy_top)

            

            
                
                for k=all_phi_buoy_base(i):all_phi_buoy_top(i)%buoy_phi_base(i):buoy_phi_top(i)
                    for l = old_all_phi_buoy_base(j):old_all_phi_buoy_top(j) %old_buoy_phi_base(j):old_buoy_phi_top(j)

                        if k==l
                            hrti_flag(i) = 1;
                            hrti_start(i) = old_hrti(j);
                            prev_buoy_layer(i) = j;
                            track_hrti=1;
                            break
                        end

                        if k==l && old_t2(j)==-5000
                            hrti_flag(i) =0 ;
                            hrti_start(i) = scale_RTI*Diameter;
                            prev_buoy_layer(i)=-5000;
                            track_hrti=1;
                            break
                        end

                        
                        
                    end
                    
                    if track_hrti==1
                        break
                    end

                end

                if track_hrti==1
                    break 
                end

               
           
            
        end
        %if track_hrti==1
          %  break
        %end
        if track_hrti==0 
            hrti_flag(i) =0 ;
            hrti_start(i) = scale_RTI*Diameter;
            prev_buoy_layer(i)=-5000;
            track_hrti=1;
        end
    end
end


%% Calculate density of magma layer


if ~isempty(buoy_phi_top)
    Av_bulk_density  = zeros(length(buoy_phi_top),1);


    for i=1:length(buoy_phi_top)
        Av_bulk_density(i) = mean(rho_bulk(buoy_phi_base(i):buoy_phi_top(i)));
    end

end


%% Caculate density of the high melt fraction portion of the magma layer


if ~isempty(buoy_phi_top)
    Av_Hphi_bulk_density  = zeros(length(buoy_phi_top),1);



    for i=1:length(buoy_phi_top)
        if buoy_Hphi_top(i)==-5000
            Av_Hphi_bulk_density(i)=0;
        else
            try
            Av_Hphi_bulk_density(i) = mean(rho_bulk(buoy_Hphi_base(i):buoy_Hphi_top(i)));
            catch
                buoy_Hphi_base(i)
                buoy_Hphi_top(i)
                error('bouyancy index error')
            end
        end
    end





end



%% Calculate hRTI


if ~isempty(buoy_phi_top)
    t2  = zeros(length(buoy_phi_top),1);
    hrti  = zeros(length(buoy_phi_top),1);



    for i=1:length(buoy_phi_top)
        if hrti_flag(i)==0
            t2(i)=0;
        else
            t2(i) = old_t2(prev_buoy_layer(i))+dt;
        end

        tau = (6*pi*crust_visc)/(Diameter*9.81*(ov_rho_bulk_buoy(i)-Av_bulk_density(i)));
        hrti(i) = hrti_start(i)*exp(dt/tau);
    end






end






%% Calculate overpressure

ov_rho_bulk_buoy=ov_rho_bulk_buoy';

if ~isempty(buoy_phi_top)
    overpressure_RTI = 0.5.*9.81.*hrti.*(ov_rho_bulk_buoy-Av_bulk_density);
    overpressure_layer = 0.5.*9.81.*buoyant_phi_m.*(ov_rho_bulk_buoy - Av_bulk_density);
else
    overpressure_RTI=0;
    overpressure_layer=0;
end
   

overpressure_total = overpressure_RTI + overpressure_layer;



%% Calculate thermal death criterion
% Uses:
    % Evacuation specific:
        % freezing parameter
        % Elastic modulus 
        % Temperature gradient
    
    
    % Used previously in three - phase code:
        % specific heat capacity (cp) - cp is constant so no concern. 
        % thermal diffusivity kt/(rho*cp) - kt and rho is dependent on the
        % composition. rho is the density of the high melt fraction layer.
        % shear viscosity of magma entering the dike (we took it as the melt
        % viscosity of the melt entering the dike)
        % Latent heat (It is dependent on the composition)




if ~isempty(buoy_phi_top)
    crit_thermdeath = zeros(length(buoy_phi_top),1);
    kt_layer = zeros(length(buoy_phi_top),1);
    latent_layer = zeros(length(buoy_phi_top),1);
    cp_layer = zeros(length(buoy_phi_top),1);
    kappa_layer = zeros(length(buoy_phi_top),1);
    muf_order_layer = zeros(length(buoy_phi_top),1);
    

    for i=1:length(buoy_phi_top)

        if buoy_Hphi_top(i)==-5000
            crit_thermdeath(i) = 0;
        else


            %kt
            kt_layer(i) = mean((kt0(1)*sum(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),[1 4]),2)+kt0(2)*sum(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),[2 5]),2)) ./max(sum(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),[1 2 5 5]),2),1e-3));
    
            
            m1_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),1);
            n1_layer =  Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),2);
            v1_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),3);
            m2_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),4);
            n2_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),5);
            v2_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),6);
            v3_layer = Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),7);

            %Lf 
            latent_points = (m1_layer*Lf(1)+n1_layer*Lf(2)+v1_layer*Lf(2))./(m1_layer+n1_layer+v1_layer);
            latent_layer(i) = mean(latent_points);

            %cp
            cp_points = ((m1_layer+n1_layer+v1_layer).*cp(1)+(m2_layer+n2_layer+v2_layer).*cp(2)+v3_layer.*cp(3))./(m1_layer+n1_layer+v1_layer+m2_layer+n2_layer+v2_layer+v2_layer);
            cp_layer(i) = mean(cp_points);


            % thermal diffusivity, kappa
           
            kappa_layer(i) = kt_layer(i);%/(cp_layer(i)*Av_Hphi_bulk_density(i));


            %mu_f

            if melt_viscosity_type == 0 

                muf_order_layer(i) = (muf2*Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),1)+muf1*Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),2))./max(sum(Mass_data(:,1:2),2),1e-3);
           
            else

                SiO2_layer=(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),1)*74+Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),2)*47)./max(sum(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),1:2),2),1e-3); % SiO2%

                if Fix_H2O==1 || Fix_H2O==3

                    H2O_layer=H2O_crust*ones(buoyant_Hphi(i),1);

                else

                    H2O_layer=100*(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),3))./max(sum(Mass_data(buoy_Hphi_base(i):buoy_Hphi_top(i),1:3),2),1e-3); % H2O%  !!!
                    
                    for j=1:length(H2O_layer)%buoy_Hphi_base(i):buoy_Hphi_top(i)
                        P3_layer=Pg_real(j+buoy_Hphi_base(i)-1)*100;
                        Saturation_layer=(2.859e-2*P3_layer-1.495e-3*P3_layer.^1.5+2.702e-5*P3_layer.^2+0.257*P3_layer.^0.5);
                        H2O_layer(j)=H2O_layer(j)/Saturation_layer;
                    end

                   

                end

                vis_B_layer=vis_b1*SiO2_layer+vis_b2*H2O_layer+vis_b3*log(1+H2O_layer);
                vis_C_layer=vis_c1*SiO2_layer+vis_c2*H2O_layer+vis_c3*log(1+H2O_layer);
                for j=1:length(vis_B_layer)
                    muf_order_layer(j)=-4.55+vis_B_layer(j)./(max(T(j+buoy_Hphi_base(i)-1),500)+273.15-vis_C_layer(j)); 
                end
            end

            mu_f_layer(i) = 10.^mean(muf_order_layer);

            %critical overpressure
            

            crit_A = (2.*cp_layer(i).*E_layer^2)/(freeze.*latent_layer(i));
            crit_B = ((3.*kappa_layer(i) * mu_f_layer(i))./pi)^(1/2);

            for j=buoy_Hphi_top(i)+1:1:N
                find_dT_dz = cellz(j)- cellz(buoy_Hphi_top(i));%sum(dz(buoy_Hphi_top(i)+1:j));
                if find_dT_dz >= dTabove
                    dTaboven = j;
                    break
                end
            end

            crit_C = abs((T(dTaboven) - T(buoy_Hphi_top(i)))/(find_dT_dz));
            


            crit_thermdeath(i) = (crit_A*crit_B*crit_C)^(2/5);

           


        end

    end
end



%% Decide the critical overpressure

if ~isempty(buoy_phi_top)

    crit_overpressure  = zeros(length(buoy_phi_top),1);

    for i = 1:length(buoy_phi_top)

        if crit_thermdeath(i)<crit_crust
            crit_overpressure(i) = crit_crust;
        else
            crit_overpressure(i) = crit_thermdeath(i);
        end
    end

    crittt= crit_overpressure;
end

























