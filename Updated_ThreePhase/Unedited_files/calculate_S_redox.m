function results = calculate_S_redox(Na2O, MgO, Al2O3, SiO2, K2O, CaO, TiO2, MnO, FeOt, ...
                                      Fe3_over_Fe_total, T_input_C, S_ppm, mode, Pressure, Ni, Cu, H2O_wt)
% calculate_S_redox - S redox + SCSS calculator
% Updated to Li & Zhang (2022) hydrosulfide capacity for H2O in the melt>4%:
%   - CHS term now uses 1.8*(X_Al + X_Fe3) instead of 1.8*X_Al
%   - Added Na+K-Al term inside the 1673/T bracket
%   - Added ln(X_OH + X_H2O) inside the erf() argument
% for melt H2O <4%, linear combination of Li & Zhang (2022), and ONeil (2021)
% All other calculations (fS2, fSO3, fSO2, CS6+, S6+/S2-, etc.) unchanged.

    if nargin < 17, H2O_wt = 0; end

    %% Constants
    R = 8.31441;
    LN10 = log(10);

    %% Molar masses (g/mol)
    M_Na2O = 30.99;   M_MgO = 40.32;   M_Al2O3 = 50.98;
    M_SiO2 = 60.08;   M_K2O = 47.10;   M_CaO = 56.08;
    M_TiO2 = 79.90;   M_MnO = 70.94;   M_FeO = 71.85;

    %% Temperature
    if mode == 1
        T_C = T_input_C;
    elseif mode == 2
        T_C = 815.3 + 265.3 * (MgO/40.32) / (MgO/40.32 + FeOt/71.85) ...
              + 15.37 * MgO + 8.61 * FeOt + 6.646 * (Na2O + K2O);
    else
        error('mode must be 1 or 2');
    end
    T_K = T_C + 273.15;

    %% Cation mole fractions (single-cation basis) -- unchanged
    cations_Na = Na2O / M_Na2O;
    cations_Mg = MgO / M_MgO;
    cations_Al = Al2O3 / M_Al2O3;
    cations_Si = SiO2 / M_SiO2;
    cations_K  = K2O / M_K2O;
    cations_Ca = CaO / M_CaO;
    cations_Ti = TiO2 / M_TiO2;
    cations_Mn = MnO / M_MnO;
    cations_Fe = FeOt / M_FeO;

    sum_cations = cations_Na + cations_Mg + cations_Al + cations_Si + ...
                  cations_K + cations_Ca + cations_Ti + cations_Mn + cations_Fe;

    X_Na = cations_Na / sum_cations;
    X_Mg = cations_Mg / sum_cations;
    X_Al = cations_Al / sum_cations;
    X_Si = cations_Si / sum_cations;
    X_K  = cations_K / sum_cations;
    X_Ca = cations_Ca / sum_cations;
    X_Ti = cations_Ti / sum_cations;
    X_Mn = cations_Mn / sum_cations;
    X_Fe = cations_Fe / sum_cations;

    Fe2_over_Fe_total = 1 - Fe3_over_Fe_total;
    X_Fe2 = X_Fe * Fe2_over_Fe_total;
    X_Fe3 = X_Fe * Fe3_over_Fe_total;   % NEW: needed for 1.8*(X_Al + X_Fe3)

    %% Redox state -- unchanged
    DeltaQFM = 4 * (log10(Fe3_over_Fe_total / Fe2_over_Fe_total) + 1.36 ...
                    - 2*X_Na - 3.7*X_K - 2.4*X_Ca);
    logfO2 = DeltaQFM - 25050/T_K + 8.58;

    %% Water speciation for ln(X_OH + X_H2O) -- NEW (Ni & Zhang, 2018)
    % Cation moles for the KOH expression use X_Si already computed above.
    % For XH2Ot, use moles per 100 g total (anhydrous + H2O).
    M_H2O = 18.01528;
    total_mass = 100 + H2O_wt;
    SiO2_t  = SiO2  * 100 / total_mass;
    TiO2_t  = TiO2  * 100 / total_mass;
    Al2O3_t = Al2O3 * 100 / total_mass;
    FeOt_t  = FeOt  * 100 / total_mass;
    MnO_t   = MnO   * 100 / total_mass;
    MgO_t   = MgO   * 100 / total_mass;
    CaO_t   = CaO   * 100 / total_mass;
    Na2O_t  = Na2O  * 100 / total_mass;
    K2O_t   = K2O   * 100 / total_mass;

    n_SiO2  = SiO2_t  / M_SiO2;
    n_TiO2  = TiO2_t  / M_TiO2;
    n_Al2O3 = Al2O3_t / M_Al2O3;
    n_FeO   = FeOt_t  / M_FeO;
    n_MnO   = MnO_t   / M_MnO;
    n_MgO   = MgO_t   / M_MgO;
    n_CaO   = CaO_t   / M_CaO;
    n_Na2O  = Na2O_t  / M_Na2O;
    n_K2O   = K2O_t   / M_K2O;

    O_moles = 2*n_SiO2 + 2*n_TiO2 + 3*n_Al2O3 + n_FeO + n_MnO + ...
              n_MgO + n_CaO + n_Na2O + n_K2O;
    
    % a linear combination of Li & Zhang (2022), and  ONeil (2021)
    if H2O_wt<=4
        correction=H2O_wt/4;
        H2O_wt=4;
    else
        correction=1;
    end
    XH2Ot = (H2O_wt/M_H2O) / (H2O_wt/M_H2O + O_moles);
    KOH = exp(2.6*X_Si - 4339*X_Si/T_K);
    tmp = (KOH - 4)/KOH * (XH2Ot - XH2Ot^2);
    tmp = max(0, min(0.25, tmp));
    XOH = (0.5 - sqrt(0.25 - tmp)) / (KOH - 4) * 2 * KOH;
    XH2Om = XH2Ot - 0.5*XOH;
    ln_XOH_XH2O = log(XOH + XH2Om)*correction;
    %% Na+K-Al term for Li & Zhang (2022) -- NEW
    Na_K_Al = X_Na + X_K - X_Al;
    if Na_K_Al > -0.0466 && Na_K_Al <= -0.015
        NKA_coeff = 26.365;   NKA_const = 0.9587;   % Clemente et al. range
    else
        NKA_coeff = 19.634;   NKA_const = 0.2542;   % Scaillet & MacDonald range
    end
    NKA_term = (1673/T_K) * (NKA_coeff * Na_K_Al + NKA_const);

    %% Step 4: ln(CS2-) -- modified to Li & Zhang (2022)
    % Changes vs original:
    %   - 1.8*X_Al  ->  1.8*(X_Al + X_Fe3)
    %   - added NKA_term
    %   - erf argument now includes ln_XOH_XH2O
    X_Fe_star = X_Fe + X_Mn;

    % term1 = 8.77 - 23590/T_K; %Oneil's original term
    term1 = (7.81-19748/T_K)*correction+(8.77 - 23590/T_K)*(1-correction);
    term2 = (1673/T_K) * (6.7*(X_Na + X_K) + 4.9*X_Mg + 8.1*X_Ca + ...
            8.9*X_Fe_star + 5*X_Ti + 1.8*(X_Al + X_Fe3) - ...
            22.2*X_Ti*X_Fe_star + 7.2*X_Fe_star*X_Si);
    term3 = -2.06 * erf(-7.2*X_Fe_star )+ ln_XOH_XH2O;%
    
    % NKA_term currently turned off! see the literature
    ln_CS2_minus = term1 + term2 + term3 + NKA_term*0;

    %% Step 5: ln(CS6+) -- unchanged
    ln_CS6_plus = -8.02 + (21100 + 44000*X_Na + 18700*X_Mg + ...
                  4300*X_Al + 44200*X_K + 35600*X_Ca + ...
                  12600*X_Mn + 16500*X_Fe2) / T_K;

    %% Step 6: ln K(SO3/S2) -- unchanged
    lnK_SO3_S2 = -55921/T_K + 25.07 - 0.6465*log(T_K);

    %% Step 7: ln(S6+/S2-) -- unchanged
    LN_S6_over_S2 = ln_CS6_plus - ln_CS2_minus - lnK_SO3_S2 + 2*LN10*logfO2;

    %% Step 8: S6+/SumS -- unchanged
    S6_over_SumS = 1 - 1/(1 + exp(LN_S6_over_S2));

    %% Step 9: fS2, fSO3, fSO2 -- unchanged
    fS2  = (S_ppm * (1 - S6_over_SumS) / exp(ln_CS2_minus))^2 * 10^logfO2;
    fSO3 = S_ppm * S6_over_SumS / exp(ln_CS6_plus);
    fSO2 = exp(43660/T_K - 9.8263 + 0.12917*log(T_K)) * sqrt(fS2) * 10^logfO2;

    %% SCSS related calculations -- unchanged
    GP = (137778 - 91.666*T_K + 8.474*T_K*log(T_K))/(R*T_K) + ...
         (-291*Pressure + 351*erf(Pressure))/T_K;

    Fe_over_FENICU = 1/(1 + (Ni/(FeOt*Fe2_over_Fe_total))*0.013 + ...
                           (Cu/(FeOt*Fe2_over_Fe_total))*0.025);
    a_FeS = log(Fe_over_FENICU * (1 - X_Fe2));

    a_FeO = log(X_Fe2) + (((1 - X_Fe2)^2) * (28870 - 14710*X_Mg + ...
            1960*X_Ca + 43300*X_Na + 95380*X_K - 76880*X_Ti) + ...
            (1 - X_Fe2)*(-62190*X_Si + 31520*X_Si^2)) / (8.31441*T_K);

    ln_SCSS = ln_CS2_minus + GP + a_FeS - a_FeO;   % now the Li&Zhang SCSS

    %% Package results -- unchanged + a few new fields
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
    results.X_Fe3 = X_Fe3;                 % NEW
    results.Fe2_over_Fe_total = Fe2_over_Fe_total;
    results.ln_CS2_minus = ln_CS2_minus;
    results.ln_CS6_plus = ln_CS6_plus;
    results.lnK_SO3_S2 = lnK_SO3_S2;
    results.LN_S6_over_S2 = LN_S6_over_S2;
    results.S6_over_SumS = S6_over_SumS;
    results.fS2 = fS2;
    results.fSO3 = fSO3;
    results.fSO2 = fSO2;
    results.ln_SCSS = ln_SCSS;
    results.ln_XOH_XH2O = ln_XOH_XH2O;     % NEW
    results.Na_K_Al = Na_K_Al;             % NEW
    results.NKA_term = NKA_term;           % NEW
end