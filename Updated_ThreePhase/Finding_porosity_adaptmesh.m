%% This is a lite version of Finding Porosity which is only used when the mesh is adapted

% The script contains the following sections:

    %initialising arrays

    %Finding porosity

    % Find buoyant magma in the areas of porosity.
        % (in this section, we find the shallowest buoyant magma layer)

    % Check hRTI is still present 


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



%Indexing for array
phi_top_i=0; 
phi_base_i=0;
Hphi_top_i=0;
Hphi_base_i=0;
buoy_phi_top_i=0;
buoy_phi_base_i=0;
buoy_Hphi_top_i=0;
buoy_Hphi_base_i=0;      


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

% HH check
%index_error=find(buoy_Hphi_base==0);
%buoy_Hphi_base(index_error)=-5000;
%buoy_Hphi_top(index_error)=-5000;

%% Compare hrti

% check lengths of hrti and new all_phi_buoy_top

check_buoy_top_km = cellz(all_phi_buoy_top)-nodez(end);
check_buoy_base_km = cellz(all_phi_buoy_base)-nodez(end);

if length(all_phi_buoy_top)~=length(hrti)
    % 3 possible cases from here:
        % 1. all_phi_buoy_top is empty, so set hrti to be empty
        % 2. all_phi_buoy_top is larger than hrti, so give the value of
        % hrti to the shallowest layer and set any new layers to the
        % initial hrti
        % 3. all_phi_buoy_top is smaller than hrti, so figure out what hrti
        % layer remains. 

        if isempty(all_phi_buoy_top) % if all_phi_buoy is empty, then hrti is empty
            hrti = [];
            t2=[];

        elseif length(all_phi_buoy_top)>length(hrti) % if length of all_phi_buoy is bigger than length of hrti

            save_hrti = hrti; %save previous hrti
            save_t2 = t2;

            for j=1:1:length(all_phi_buoy_top)
                hrti(j) = scale_RTI*Diameter; % just set all to be initial hrti
                t2(j)=0;
            end


            % now find if previous buoyancy layer overlaps with new
            % buoyancy layer since the the mesh adapted. Gives the
            % shallowest overlap, the previous hrti.
            for i=1:1:length(save_hrti)
                check_hrti=0;
                for j=1:1:length(all_phi_buoy_top)
                    
                    overlap_top = min(old_check_buoy_top_km(i), check_buoy_top_km(j));
                    overlap_base = max(old_check_buoy_base_km(i), check_buoy_base_km(j));

                    if overlap_base <= overlap_top
                        hrti(j) = save_hrti(i);
                        t2(j) = save_t2(i);
                        check_hrti = 1;
                    end


                    if check_hrti==1
                        break
                    end
                end
            end


           

        else % if length of all_phi_buoy is smaller than length of hrti

            save_hrti = hrti; %save previous hrti
            save_t2 = t2;

            for j=1:1:length(all_phi_buoy_top)
                hrti(j) = scale_RTI*Diameter; % just set all to be initial hrti
                t2(j) = 0;
            end

            % now find if previous buoyancy layer overlaps with new
            % buoyancy layer since the the mesh adapted. Gives the
            % shallowest overlap, the previous hrti.
            for j=1:1:length(all_phi_buoy_top)
                check_hrti=0;
                for i=1:1:length(save_hrti)
                    
                    overlap_top = min(old_check_buoy_top_km(i), check_buoy_top_km(j));
                    overlap_base = max(old_check_buoy_base_km(i), check_buoy_base_km(j));

                    if overlap_base <= overlap_top
                        hrti(j) = save_hrti(i);
                        t2(j) = save_t2(i);
                        check_hrti = 1;
                    end

                    if check_hrti==1
                        break
                    end
                end
            end



        end


end