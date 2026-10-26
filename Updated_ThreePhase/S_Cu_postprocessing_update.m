%% This script is a post processing tool to calculate the mass of Sulfur and Copper within the model

%Total mass represented by the model components.
Mass_total=(MM+NN+V+CL); % remove CL from these calculations, as CL is a minor component 


% Molar masses (g/mol)
M_H2O = 18.015;
M_Cl = 35.45;    % elemental chlorine
M_S = 32.065; % elemental sulfur


% Convert total CU and S concentrations to fractions of the total mass
CU_fraction=CU./Mass_total;
SUL_fraction=SUL./Mass_total;

% Process only the cells where S and Cu are required 
parfor i=S_CU_to_process
    
    % Assign values from Mass_data
    m1=Mass_data(i,1); n1=Mass_data(i,2); v1=Mass_data(i,3); % Melt (evolved, refractory, fluid)
    m2=Mass_data(i,4); n2=Mass_data(i,5); v2=Mass_data(i,6); % Solid (evolved, refractory, fluid)
    v3=Mass_data(i,7); % Fluid component
    l1=Mass_data(i,8); l2=Mass_data(i,9); l3=Mass_data(i,10); % Chlorine (melt, solid, fluid) 
    
    % If m1 or n1 are greater than 0 then a melt phase exists
    if m1+n1>0
        % Melt phase exists

        %Calculate the SiO2 content of the melt
        SiO2=(PD_range(2)-PD_range(1))*m1/(m1+n1)+PD_range(1);

        % Calculate the H2O concentration in the melt (%)
        H2O=v1/(m1+n1+v1)*100;

    else % Water enter solid phase until saturation (disable current)
        % No melt phase present 
        % Code estimates a hypothetical melt composition from the bulk
        % solid+volatile composition. This is used for the subsequent S/Cu
        % chemistry calculations

        % Calculate the ghost melt composition (only using the evolved
        % parameters??
        cl_ghost=MM(i)/(MM(i)+NN(i)); % edit - it is melt 

        % Calculate the SiO2 of cl_ghost
        SiO2=(PD_range(2)-PD_range(1))*cl_ghost+PD_range(1);

        % Calculate the H2O that can partition into the ghost melt, limit it must be smaller than
        % saturation, other it has the value of saturation
        H2O=min(V(i)/(MM(i)+NN(i)+V(i))/K0,Saturation(i))*100;

    end
    

    %% Construct the melt major-element composition
    % Columns or wt_percent correspond to:
    % 1. 'SiO2'. 
    % 2. 'TiO2', 
    % 3. 'Al2O3',
    % 4. 'FeO',
    % 5. 'MnO',
    % 6. 'MgO',
    % 7. 'CaO',
    % 8. 'Na2O',
    % 9. 'K2O' 
    % 10. 'H2O'
    wt_percent=zeros(1,10);

    % SiO2 and H2O are calculated directly above
    wt_percent(1)=SiO2;
    wt_percent(10)=H2O;

    % The other eight oxide concentrations are a function of SiO2
    for j=1:8
        wt_percent(1+j)    =polyval(Oxides_coefficients(j,:), SiO2);
    end

    %% Rescale the other oxides so that the total composition sums to 100
    % wt%
    
    % SiO2 and H2O are deliberately kept at their calculated values.
    Target=100-SiO2-H2O;

    % Ratio required to make the remaining eight oxides sum to Target.
    ratio=Target/sum(wt_percent(2:9));

    % Apply the correction from TiO2 to K2O (2 to 9)
    wt_percent(2:9)=wt_percent(2:9)*ratio;
    
    

    %% Calculate sulphur redox/speciation and sulphide saturation.

    % Function input order:
    % Na2O. MgO, Al2O3, SiO2, K2O, CaO, TiO2, MnO, FeOt, Fe3+/Fe_total,
    % temperature, S concentration, mode, pressure, Ni, Cu, H2O

    % Ni (100ppm) and Cu (80ppm, in this step) are assumed values here because they only have a relatively
    % small effect on the calculated aFes(??) whose contribution is much
    % smaller than other terms
       
    %H2O content of the melt 
    H2O_melt=v1/max(m1+n1+v1,1e-12)*100;
    
    % update comment when calculate_s_redox has been gone through
    results=calculate_S_redox_update(wt_percent(8), wt_percent(6), wt_percent(3), wt_percent(1), wt_percent(9), wt_percent(7), wt_percent(2), wt_percent(5), wt_percent(4),    Fe3ratio, T(i), 0, 1, P_kbar(i)/10, 100, 80, H2O_melt);

    % Difference in oxygen fugacity relative to the FMQ buffer (?)
    dFMQ=results.DeltaQFM; % this gets DeltaFMQ from results
    
    % Partition coefficient of sulphur 
    D_s=40000/(238+1701/(1+10^(0.78-2*(dFMQ-1)))); %reference - Zajacz paper - figure 10. Does it matter that the numerator represents 4wt% water?

    % Fraction of total sulphur present as S6+
    S6_fraction=results.S6_over_SumS; %ratio of S6+/S2-

    % SCSS - Maximum solubility of reduced sulfur (S2-) in the melt.
    % Defined by the sulfide saturation limit (SCSS)
    % Returned in logarithmic form, convert to a fraction by multiplying
    % by taking the exponential and then multiplying by 1e-6 to convert
    % from ppm to mass fraction
    SCSS=exp(results.ln_SCSS)*1e-6; 

   
    % p1, p2 and p3 represent the fractions of the total system assigned to
    % phases 1, 2 and 3
    p1=(m1+n1+v1+l1)/(MM(i)+NN(i)+V(i)+CL(i)); % melt
    p2=(m2+n2+v2+l2)/(MM(i)+NN(i)+V(i)+CL(i)); % solid
    p3=(v3+l3)/(MM(i)+NN(i)+V(i)+CL(i)); % volatile


    %% Sulphur partitioning
    % Note that the sulphur partitions between the two phases, with the
    % concentrations related by:
    % S3 = D_s S1, where S3 is the sulphur concentration in the volatile phase (phase 3) and S1
    % is the sulphur concentration in the melt phase (phase 1)
    % 
    % The bulk sulphur mass balance is given by
    % SUL_fraction = p1*S1+p2*S2+p3*S3 where phase 2 is the solid phase
    %
    % In the following code, we initially assume that there is no sulphur
    % in the solid phase. 
    % 
    % Putting the two above equations and assumption together gives us the following equation:
    % SUL_fraction = p1*S1+p3*S1*D_s, which we rearrange to give us S1 if
    % the melt phase is present
   
    if p1>0
        % melt phase exists

        % Flag indicating whether a separate sulphide phase exists
        Sulfide_exist=0;

        % Calculate sulphur concentration in the melt phase (phase 1)
        % SUL_fraction = total sulphur concentration in the node
        % D_s controls the relationship between S in the melt and volatile phase
         S1=SUL_fraction(i)/(p3*D_s+p1); % is this wrong? % melt concentration* melt fraction + volatile coefficient*(1 - melt fraction)
        
        % Check whether the sulphur in the melt phase (phase 1) exceeds
        % sulphide saturation
        % (1- S6_fraction) converts total S to reduced sulphur (?)
        if S1*(1-S6_fraction)>SCSS

            % Limit melt phase sulphur to the saturation concentration.
            S1=SCSS/(1-S6_fraction); 

            % A separate sulphide phase must therefore exist
            Sulfide_exist=1;
        end
    else
        % No phase 1/melt exists, so assign no S to phase 1.
        S1=0;
    end

    % Calculate the sulphur content of the volatile phase based on the
    % partition coefficient
    S3=S1*D_s;
    % Calculate the sulphur content of the solid phase from mass balance
    S2=(SUL_fraction(i)-p1*S1-p3*S3)/max(p2, 1e-12);


    %% H2O - Cl - H2S COMPOSITION FOR Cu partitioning 

    if v3>0
        % Total mass of the volatile/Cl/ S-bearing phase    
        Total_m=(v3+l3)*(1+S3); %(1+S3) is used to convert the mass of the non-sulphur components into the mass including sulphur
        
        % Mass fraction of H2O and cl in the volatile phase
        w_H2O=v3/Total_m; 
        w_Cl=l3/Total_m;

        % Sulphur fraction in the volatile phase
        W_S=S3;
        
        % Convert H2O and Cl from mass fractions to mole fractions
        total_moles=w_H2O/M_H2O+w_Cl/M_Cl+S3/M_S;
        X_Cl=w_Cl/M_Cl/total_moles;
        X_H2O=w_H2O/M_H2O/total_moles;

        % Mole fractions of H2S.
        % Only the reduced sulphur fraction (1-S6_fraction) is included as
        % H2S
        X_H2S=S3/M_S/total_moles*(1-S6_fraction);
    else
        % No volatile phase -> no H2O, Cl or H2S for Cu partitioning.
        X_Cl=0;
        X_H2S=0;
        X_H2O=0;
    end
    
    % Fe2+ fraction/activity-related parameter returned by the redox model.
    X_FeO=results.X_Fe2;

    % Convert temperature from Celsius to Kelvin.
    T2=T(i)+273.15; %Temperature in K




    %% Copper partitioning

    % Dcu1 - Cu partition coefficient between melt and volatiles phases depends
    % strongly on:
    % X_Cl = chlorine mole fraction
    % X_H2O = water mole fraction
    % X_H2S = H2S mole fraction
    Dcu1=8e4*X_Cl^2*X_H2O^14*(1+180*X_H2S)+380*X_Cl+0.8; % Tattich and Blundy 2017 (there is an error range here for the numbers)
    % In tattich and blundy it is not X_Cl, but instead XNaCl?
    
    % Dcu2 - Cu partition coefficient between sulphide liquid and sulphide
    % melt (?). Sulphide liquid separates from the melt.

    % Dependence on temperature, oxygen fugacity (dFMQ) and FeO content. 
    Dcu2=exp(1.13+0.39*1e4/T2+0.02*dFMQ-0.86*log(X_FeO)); %Li and Audetat - equation 5. Constants are in Table 4, Cu.
    % Asusmptions for this partition coefficient: (1) calculated for a
    % magma differentiating from basaltic composition to granitic
    % composition. (2) Magma has a MORB major element composition (see Table 1 in
    % Li and Audetat). (3) It is assumed to have 1.5wt% water.

    % is the diffusion coefficient wrong here too?

    
    if p1>0
        % Melt phase exists
        if Sulfide_exist==0
            % Before sulphide saturation, set copper in the solid to 0
            u2=0;

            % Calculate copper concentration in the melt phase from the
            % bulk copper concentration, given that copper partitions
            % between the melt and volatile phases.
            u1=CU_fraction(i)/(p3*Dcu1+p1);

            % Calculate copper concentration in the volatile  phase
            u3=u1*Dcu1;
        else
            % MassCuCh is mass fraction of copper in the chalcopyrite. We
            % assume all the sulfur is chalcopyrite. Jon Blundy work. 

            % After sulphide saturation, allow copper into the solid phase.
            % Include the solid-phase contribution in the bulk copper mass
            % balance when calculating the copper in the melt phase.
            u1=CU_fraction(i)/(p3*Dcu1+p1+ p2*S2/MassCuCh*Dcu2); 
            % Calculate copper concentration in the volatile and solid
            % phase
            u3=u1*Dcu1;
            u2=u1*Dcu2*S2/MassCuCh;
        end
    else
        % Melt phase does not exist. Calculate the amount of copper in the
        % solid and volatile phase.
        u1=0;
        u2=CU_fraction(i)/(p2+p3*Dcu1/Dcu2);
        u3=u2*Dcu1/Dcu2;
    end
    % Store the mass contribution of sulphur and copper from each
    % component. Multiply each concentration by its phase proportion and
    % total cell mass. 
    Mass_data_addition(i,:)=[S1*p1*Mass_total(i),S2*p2*Mass_total(i),S3*p3*Mass_total(i),u1*p1*Mass_total(i),u2*p2*Mass_total(i),u3*p3*Mass_total(i)];
end