function results = calculate_S_redox(Na2O, MgO, Al2O3, SiO2, K2O, CaO, TiO2, MnO, FeOt, ...
                                      Fe3_over_Fe_total, T_input_C, S_ppm, mode, Pressure, Ni, Cu)
% calculate_S_redox - Replicates the Table S6 S redox calculator from:
% O'Neill & Mavrogenes (2022) GCA
%
% INPUTS:
%   Na2O, MgO, Al2O3, SiO2, K2O, CaO, TiO2, MnO, FeOt - weight percent oxides
%   Fe3_over_Fe_total - Fe3+/(Fe2+ + Fe3+) ratio (0 to 1)
%   T_input_C - temperature in °C (used if mode=1)
%   S_ppm - sulfur concentration in ppm
%   Pressure in Gpa
%   Ni and Cu in ppm
%   mode - 1 = use input T, 2 = calculate T from olivine-melt thermometer
%
% OUTPUTS (structure with fields):
%   T_C - temperature in °C (input or calculated)
%   T_K - temperature in Kelvin
%   logfO2 - log10(fO2) relative to 1 bar
%   DeltaQFM - log10(fO2) relative to QFM buffer
%   X_Na, X_Mg, X_Al, X_Si, X_K, X_Ca, X_Ti, X_Mn, X_Fe, X_Fe2 - cation mole fractions
%   Fe2_over_Fe_total - Fe2+/(Fe2+ + Fe3+) ratio
%   ln_CS2_minus - ln(sulfide capacity)
%   ln_CS6_plus - ln(sulfate capacity)
%   lnK_SO3_S2 - ln(equilibrium constant for SO3/S2)
%   LN_S6_over_S2 - ln(S6+/S2-)
%   S6_over_SumS - S6+/(S6+ + S2-) ratio
%   fS2 - fugacity of S2 in bars
%   fSO3 - fugacity of SO3 in bars
%   fSO2 - fugacity of SO2 in bars
%
% Example:
%   results = calculate_S_redox(2.63, 5.88, 13.36, 50.93, 0.17, 10.81, ...
%                               2.46, 0.247, 14.21, 0.112, 1400, 1763.2, 1);
    %% Constants
    R = 8.31441; % Gas constant (J/mol/K)
    LN10 = log(10); % Natural log of 10
    
    %% Molar masses of oxides (g/mol)
    M_Na2O = 30.99;   % Na2O
    M_MgO = 40.32;    % MgO
    M_Al2O3 = 50.98;  % Al2O3 (for Al, but in oxide fraction calculations we use 50.98 for Al2O3)
    M_SiO2 = 60.08;   % SiO2
    M_K2O = 47.10;    % K2O (actually 94.2/2 for cation basis)
    M_CaO = 56.08;    % CaO
    M_TiO2 = 79.90;   % TiO2
    M_MnO = 70.94;    % MnO
    M_FeO = 71.85;    % FeO(t)
    
    %% Step 1: Determine temperature
    if mode == 1
        % Use input temperature
        T_C = T_input_C;
    elseif mode == 2
        % Calculate T from modified Putirka (2008) olivine-melt thermometer
        % Modified to use Al2O3 instead of MgO for Al-bearing compositions
        T_C = 815.3 + 265.3 * (MgO/40.32) / (MgO/40.32 + FeOt/71.85) ...
              + 15.37 * MgO + 8.61 * FeOt + 6.646 * (Na2O + K2O);
    else
        error('mode must be 1 (input T) or 2 (calculate T)');
    end
    
    T_K = T_C + 273.15; % Convert to Kelvin
    
       
    %% Step 3: Calculate cation mole fractions (columns U to AD)
    % Single cation basis (NaO0.5, MgO, AlO1.5, SiO2, KO0.5, CaO, TiO2, MnO, FeO)
    % For Fe, we need to split into Fe2+ and Fe3+
    
    % First, calculate Fe2+ and Fe3+ from Fe3+/sumFe
    Fe2O3_wt = FeOt * (Fe3_over_Fe_total) * (159.69 / (2 * 71.85)); % Convert FeO to Fe2O3 for Fe3+
    FeO_wt = FeOt * (1 - Fe3_over_Fe_total); % Fe2+ as FeO
    
    % But for cation fractions, we use FeO(t) as Fe2+ initially, then correct
    % The spreadsheet uses FeO(t) as total Fe, then calculates Fe2+ fraction
    
    % Actually, the spreadsheet uses FeO(t) as is, then calculates X_Fe2 = X_Fe(t) * Fe2+/sumFe
    % So we use FeOt as the total FeO
    
    % Calculate oxide mole fractions (using oxide molar masses)
    n_Na2O = Na2O / M_Na2O;
    n_MgO = MgO / M_MgO;
    n_Al2O3 = Al2O3 / M_Al2O3;
    n_SiO2 = SiO2 / M_SiO2;
    n_K2O = K2O / M_K2O;
    n_CaO = CaO / M_CaO;
    n_TiO2 = TiO2 / M_TiO2;
    n_MnO = MnO / M_MnO;
    n_FeOt = FeOt / M_FeO; % Total Fe as FeO
    
    % For single cation basis, divide each oxide by its cation count
    % Na2O -> 2 Na per oxide, so Na cations = 2 * n_Na2O
    % But in spreadsheet, they use Na2O/30.99 (i.e., treating as NaO0.5)
    % So we use the cation equivalents directly
    cations_Na = Na2O / M_Na2O;  % 1 cation per NaO0.5
    cations_Mg = MgO / M_MgO;
    cations_Al = Al2O3 / M_Al2O3;
    cations_Si = SiO2 / M_SiO2;
    cations_K = K2O / M_K2O;
    cations_Ca = CaO / M_CaO;
    cations_Ti = TiO2 / M_TiO2;
    cations_Mn = MnO / M_MnO;
    cations_Fe = FeOt / M_FeO; % Total Fe cations
    
    sum_cations = cations_Na + cations_Mg + cations_Al + cations_Si + ...
                  cations_K + cations_Ca + cations_Ti + cations_Mn + cations_Fe;
    
    % Cation mole fractions (columns U to AD in spreadsheet)
    X_Na = cations_Na / sum_cations;
    X_Mg = cations_Mg / sum_cations;
    X_Al = cations_Al / sum_cations;
    X_Si = cations_Si / sum_cations;
    X_K = cations_K / sum_cations;
    X_Ca = cations_Ca / sum_cations;
    X_Ti = cations_Ti / sum_cations;
    X_Mn = cations_Mn / sum_cations;
    X_Fe = cations_Fe / sum_cations; % Total Fe
    
    % Now calculate Fe2+ fraction (column AE)
    Fe2_over_Fe_total = 1 - Fe3_over_Fe_total;
    X_Fe2 = X_Fe * Fe2_over_Fe_total; % Column AE = AC * AE (X_Fe * Fe2/Fe_total)
    %%
    DeltaQFM=4 * (log10(Fe3_over_Fe_total / Fe2_over_Fe_total) + 1.36 - 2*X_Na - 3.7*X_K - 2.4*X_Ca);
    logfO2=DeltaQFM-25050/T_K+8.58;



    %% Step 4: Calculate ln(CS2-) (column AG)
    % Eq. from O'Neill (2021) for sulfide capacity
    % ln(CS2-) = 8.77 - 23590/T + (1673/T)*(6.7*(X_Na+X_K)+4.9*X_Mg+8.1*X_Al+ ...
    %           8.9*(X_Fe2+X_Mn)+5*X_Ca+1.8*X_Ti-22.2*X_Ca*(X_Fe2+X_Mn)+ ...
    %           7.2*((X_Fe2+X_Mn)*X_Si)) - 2.06*ERF(-7.2*(X_Fe2+X_Mn))
    
    % Note: The ERF term in the spreadsheet is -2.06*ERF(-7.2*(X_Fe2+X_Mn))
    % But ERF(-x) = -ERF(x), so this is +2.06*ERF(7.2*(X_Fe2+X_Mn))
    
    term1 = 8.77 - 23590/T_K;
    term2 = (1673/T_K) * (6.7*(X_Na + X_K) + 4.9*X_Mg + 8.1*X_Ca + ...
            8.9*(X_Fe + X_Mn) + 5*X_Ti + 1.8*X_Al - ...
            22.2*X_Ti*(X_Fe + X_Mn) + 7.2*((X_Fe + X_Mn)*X_Si));
    term3 = -2.06 * erf(-7.2*(X_Fe + X_Mn));
    
    % Or equivalently: term3 = 2.06 * erf(7.2*(X_Fe2 + X_Mn));
    
    ln_CS2_minus = term1 + term2 + term3;
    
    %% Step 5: Calculate ln(CS6+) (column AH)
    % From Eq. 12a in the paper:
    % ln(CS6+) = -8.02 + (21100 + 44000*X_Na + 18700*X_Mg + 4300*X_Al + ...
    %            44200*X_K + 35600*X_Ca + 12600*X_Mn + 16500*X_Fe2) / T
    
    ln_CS6_plus = -8.02 + (21100 + 44000*X_Na + 18700*X_Mg + ...
                  4300*X_Al + 44200*X_K + 35600*X_Ca + ...
                  12600*X_Mn + 16500*X_Fe2) / T_K;
    
    %% Step 6: Calculate ln K(SO3/S2) (column AI)
    % From JANAF tables, Eq. 6b in paper:
    % ln K = 55921/T - 25.07 + 0.6465*ln(T)
    % But the spreadsheet stores -ln K, so we compute it as such
    lnK_SO3_S2 = -55921/T_K + 25.07 - 0.6465*log(T_K);
    % This equals -ln(K(6))
    
    %% Step 7: Calculate ln(S6+/S2-) (column AJ)
    % Eq. 7: ln(S6+/S2-) = ln(CS6+) - ln(CS2-) + ln(K(6)) + 2*ln(fO2)
    % But AI = -ln(K(6)), so:
    % AJ = ln(CS6+) - ln(CS2-) - AI + 2*ln(fO2)
    % Since ln(fO2) = LN10 * log10(fO2)
    
    LN_S6_over_S2 = ln_CS6_plus - ln_CS2_minus - lnK_SO3_S2 + 2*LN10*logfO2;
    
    %% Step 8: Calculate S6+/SumS (column AK)
    % S6+/(S6+ + S2-) = 1 - 1/(1 + exp(LN_S6_over_S2))
    S6_over_SumS = 1 - 1/(1 + exp(LN_S6_over_S2));
    
    %% Step 9: Calculate fS2, fSO3, fSO2 (columns AM, AN, AO)
    % fS2 = (S_ppm*(1-S6_over_SumS)/exp(ln_CS2_minus))^2 * 10^logfO2
    % fSO3 = S_ppm*S6_over_SumS / exp(ln_CS6_plus)
    % fSO2 = exp(43660/T - 9.8263 + 0.12917*ln(T)) * sqrt(fSO3) * 10^logfO2
    
    fS2 = (S_ppm * (1 - S6_over_SumS) / exp(ln_CS2_minus))^2 * 10^logfO2;
    fSO3 = S_ppm * S6_over_SumS / exp(ln_CS6_plus);
    fSO2 = exp(43660/T_K - 9.8263 + 0.12917*log(T_K)) * sqrt(fS2) * 10^logfO2;

    %% SCSS related calculations
    % free energy term  
    GP=(137778-91.666*T_K+8.474*T_K*log(T_K))/(R*T_K)+(-291*Pressure+351*erf(Pressure))/T_K;
    % activation FeS
    Fe_over_FENICU=1/(1+(Ni/(FeOt*Fe2_over_Fe_total))*0.013+(Cu/(FeOt*Fe2_over_Fe_total))*0.025);
    a_FeS=log(Fe_over_FENICU*(1-X_Fe2));
    % activation FeO
    a_FeO=log(X_Fe2)+(((1-X_Fe2)^2)*(28870-14710*X_Mg+1960*X_Ca+43300*X_Na+95380*X_K-76880*X_Ti)+(1-X_Fe2)*(-62190*X_Si+31520*X_Si^2))/(8.31441*T_K);
    ln_SCSS=ln_CS2_minus+GP+a_FeS-a_FeO;
    
    %% Package results
    results.T_C = T_C;
    results.T_K = T_K;
    results.logfO2 = logfO2;
    results.DeltaQFM = DeltaQFM;
    results.X_Na = X_Na;
    results.X_Mg = X_Mg;
    results.X_Al = X_Al;
    results.X_Si = X_Si;
    results.X_K = X_K;
    results.X_Ca = X_Ca;
    results.X_Ti = X_Ti;
    results.X_Mn = X_Mn;
    results.X_Fe = X_Fe;
    results.X_Fe2 = X_Fe2;
    results.Fe2_over_Fe_total = Fe2_over_Fe_total;
    results.ln_CS2_minus = ln_CS2_minus;
    results.ln_CS6_plus = ln_CS6_plus;
    results.lnK_SO3_S2 = lnK_SO3_S2;
    results.LN_S6_over_S2 = LN_S6_over_S2;
    results.S6_over_SumS = S6_over_SumS;
    results.fS2 = fS2;
    results.fSO3 = fSO3;
    results.fSO2 = fSO2;
    results.ln_SCSS=ln_SCSS;
end