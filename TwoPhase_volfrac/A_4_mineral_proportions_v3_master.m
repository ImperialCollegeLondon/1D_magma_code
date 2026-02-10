function [ol, opx, cpx, feld,ol_mg,px_mg,opx_mg,cpx_mg,feld_mg,ol_mn,sol_mn] = A_4_mineral_proportions_v3_master(ol,opx,cpx, feld, ...
                                T,phi,phi_old,cb,cs,cl,crit_T1,crit_T2,crit_cb1,crit_cb2,crit_cb3,crit_cb4,crit_cb5,...
                                s1m,s1c,s3m,s3c,s4m,s4c,s6m,s6c,s7m,s7c,s8m,s8c,s11m,s11c, ...
                                DC, DT, EC, ET, NC, NT, copx_T_m, copx_T_c, copx_C_m, copx_C_c, ...
                                px_k, olpx_C, ol_i_m0, ol_i_c0, ol_i_m1, ol_i_c1, ol_i_m2, ol_i_c2, ol_i_m3, ol_i_c3,...
                                ol_e_m1, ol_e_c1, ol_e_m2, ol_e_c2, ol_e_m3, ol_e_c3, ol_e_m4, ol_e_c4, ...
                                ol_T_high_m, ol_T_high_c, ol_T_low_m, ol_T_low_c, ol_T_lowest_m, ol_T_lowest_c, ...
                                crit_sol_mn_cb, ol_mn_bridge,sol_mn_k_hCb, sol_mn_k_lCb, mn_k_m, mn_k_c, ...
                                a_mn, b_mn, c_mn,sol_mni_m, sol_mni_c, ol_T_crit_cb2, ol_i_m, ol_i_c, ol_mg_m,...
                                cb_low_high, copx_cb_m_1h, copx_cb_c_1h, copx_cb_m_1l, copx_cb_c_1l, ...
                                copx_T_m_1h, copx_T_c_1h, copx_T_m_1l, copx_T_c_1l, ...
                                copx_cb_m_2h, copx_cb_c_2h, copx_cb_m_2l, copx_cb_c_2l, ...
                                copx_T_m_2h, copx_T_c_2h, copx_T_m_2l, copx_T_c_2l,...
                                px_cb_m, px_cb_c, liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,liq_P2_C,Ts_local,...
                                max_mg, min_mg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% ADD phi_old

%% Calculate magnesium content of minerals for the temperature

ol_mg_T_p = s1m*crit_T2 +s1c;
s12m = (ol_mg_T_p-NC)/(crit_T2-NT);
s12c = ol_mg_T_p - s12m*crit_T2;

ol_mg=0;
px_mg=0;
opx_mg=0;
cpx_mg=0;
feld_mg=0;

cb1=cb;%+0.005;


DC_dim = DC*(max_mg-min_mg)+min_mg;

crit_cb2_dim = crit_cb2*(max_mg-min_mg)+min_mg;
cb_dim = cb*(max_mg-min_mg)+min_mg;

D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);

if cb>crit_cb2

    D2C = (D2C_dim - min_mg)/(max_mg-min_mg);
else
    D2C = (DC_dim - min_mg)/(max_mg-min_mg);
end

% olivine mg 

%% Olivine

% Olivine MnO relationship

ol_mg_T_P = s1m*crit_T2+s1c;
s12_m = (ol_mg_T_P - NC)/(crit_T2 - NT);
s12_c = ol_mg_T_P - s12_m*crit_T2;
cb_dim = cb*(max_mg-min_mg)+min_mg;

% Ol MnO initial comps
if cb>olpx_C
    ol_mn_initial = ol_i_m0*cb_dim+ol_i_c0;
elseif cb>crit_cb2
    ol_mn_initial = ol_i_m1*cb_dim +ol_i_c1;
elseif cb>crit_cb5
    ol_mn_initial = ol_i_m2*cb_dim+ol_i_c2;
else
    ol_mn_initial = ol_i_m3*cb_dim + ol_i_c3;
end

if ol_mn_initial <0 
    ol_mn_initial =0;
end


% Ol MnO end comps
if cb>crit_cb2
    ol_mn_end = ol_e_m1*cb_dim+ol_e_c1;
elseif cb>crit_cb4
    ol_mn_end = ol_e_m2*cb_dim+ol_e_c2;
elseif cb>crit_cb5
    ol_mn_end = ol_e_m3*cb_dim+ol_e_c3;
else
    ol_mn_end = ol_e_m4*cb_dim+ol_e_c4;
end

if ol_mn_end<0
    ol_mn_end =0;
end

% Ol MnO start temp
if cb > crit_cb2
    ol_mn_T1 = ol_T_high_m*cb_dim+ol_T_high_c;
elseif cb > crit_sol_mn_cb
    ol_mn_T1 = ol_T_low_m*cb_dim+ol_T_low_c;
elseif cb>crit_cb5
    ol_mn_T1 = ol_mn_bridge;
else
    ol_mn_T1 = ol_T_lowest_m*cb_dim + ol_T_lowest_c;
end

if ol_mn_T1 > crit_T2
    ol_mn_T1 = crit_T2;
end

% Ol MnO end temp
if cb>=s1m*crit_T2+s1c
    ol_mn_T2 = (cb - s1c)/s1m;
elseif cb>=DC
    ol_mn_T2 = crit_T2;
elseif cb>=crit_cb4
    ol_mn_T2 = 872;
else
    ol_mn_T2 =800;
end


% Calculate ol MnO gradient and constant
ol_mn_m = (ol_mn_end - ol_mn_initial)/(ol_mn_T2 - ol_mn_T1);
ol_mn_c = ol_mn_initial - ol_mn_m*ol_mn_T1;

% Solid MnO calculations
if cb>crit_cb4
    sol_mn_k = sol_mn_k_hCb;
elseif cb<=crit_cb5
    sol_mn_k = sol_mn_k_lCb;
else
    sol_mn_k = mn_k_m*cb_dim + mn_k_c;
    if sol_mn_k<1100
        sol_mn_k = 1100;
    end
end

sol_mn_T1 = ol_mn_T1;
ol_mg_end = s1m*DT+s1c;
opx_mg_end = s3m*DT+s3c;
phi_end = (D2C - cb)/(D2C - crit_cb1);
ol_end = (D2C - opx_mg_end)/(ol_mg_end - opx_mg_end)*(1-phi_end);

if cb>crit_cb5
    sol_mn_c1 = (ol_mn_initial*ol_end)/(1-phi_end);
else
    sol_mn_c1 = sol_mni_m*cb_dim + sol_mni_c;
    if sol_mn_c1>1.7
        sol_mn_c1 = 1.7;
    end
end

if sol_mn_c1 < 0 
    sol_mn_c1 =0;
end

sol_mn_T2 = ol_mn_T2;
sol_mn_c2 = a_mn * cb_dim^2 + b_mn * cb_dim + c_mn;

% calculating solid MnO reciprocal relationship. 
sol_mn_b = ((sol_mn_k - sol_mn_T1 )*sol_mn_c1 - (sol_mn_k - sol_mn_T2)*sol_mn_c2)/(sol_mn_T2 - sol_mn_T1);
sol_mn_a  = (sol_mn_k - sol_mn_T1)*sol_mn_c1 - sol_mn_b*(sol_mn_k - sol_mn_T1);


if cb>DC
    ol_mg = s1m*T+s1c;
    if T<Ts_local
        ol_mg = s1m*Ts_local+s1c;
    end
    ol_mn = 0; 
    sol_mn = 0;
elseif cb>crit_cb2
    if T >= crit_T2
        ol_mg = s1m*T+s1c;
        ol_mn = 0;
        sol_mn = 0;

    elseif T<copx_T_m_2h*cb + copx_T_c_2h
        
        ol_mn = ol_mn_m*T + ol_mn_c;
        sol_mn = sol_mn_b + sol_mn_a/(sol_mn_k - T);
        
        ol_mg_m_lowT = 800;
        ol_mg_m_highT = copx_T_m_2h*cb + copx_T_c_2h;
        ol_mg_m_highc = s12m*ol_mg_m_highT+s12c;
        ol_mg_m_lowc = ((33-5*(18.9-cb_dim))-min_mg)/(max_mg-min_mg);
        ol_mg_m_dum = (ol_mg_m_highc-ol_mg_m_lowc)/(ol_mg_m_highT-ol_mg_m_lowT);
        ol_mg_c_dum = ol_mg_m_highc - ol_mg_m_dum*ol_mg_m_highT;



        ol_mg = ol_mg_m_dum*T + ol_mg_c_dum;

        if T<Ts_local
            ol_mg_m_lowT = 800;
            ol_mg_m_highT = copx_T_m_2h*cb + copx_T_c_2h;
            ol_mg_m_highc = s12m*ol_mg_m_highT+s12c;
            ol_mg_m_lowc = ((33-5*(18.9-cb_dim))-min_mg)/(max_mg-min_mg);
            ol_mg_m_dum = (ol_mg_m_highc-ol_mg_m_lowc)/(ol_mg_m_highT-ol_mg_m_lowT);
            ol_mg_c_dum = ol_mg_m_highc - ol_mg_m_dum*ol_mg_m_highT;
            ol_mg = ol_mg_m_dum*Ts_local + ol_mg_c_dum;
        end

    else
        ol_mg = s12m*T+s12c;
        ol_mn = ol_mn_m*T +ol_mn_c;
        sol_mn = sol_mn_b + (sol_mn_a)/(sol_mn_k - T);
    end
    
elseif cb>crit_cb4

    if T>ol_T_crit_cb2
        ol_mg=0;
        ol_mn=0;
        sol_mn = sol_mn_b + sol_mn_a/(sol_mn_k-T);

    elseif T<copx_T_m_2h*cb + copx_T_c_2h
        ol_mg_i = ol_i_m*cb + ol_i_c;
        ol_mg_c = ol_mg_i - ol_T_crit_cb2*ol_mg_m;
        ol_mg = ol_mg_m*T + ol_mg_c;

        if T<Ts_local
             ol_mg = ol_mg_m*Ts_local + ol_mg_c;
        end


        ol_mn = ol_mn_m*T + ol_mn_c;
        sol_mn = sol_mn_b + sol_mn_a/(sol_mn_k - T);

    else
       
        ol_mg_i = ol_i_m*cb + ol_i_c;
        ol_mg_c = ol_mg_i - ol_T_crit_cb2*ol_mg_m;
        meet_mg = ol_mg_m*(copx_T_m_2h*cb + copx_T_c_2h) + ol_mg_c;

        s12_up = meet_mg - s12m*((copx_T_m_2h*cb + copx_T_c_2h));

        ol_mg = s12m*T+s12_up;



        ol_mn = ol_mn_m*T +ol_mn_c;
        sol_mn = sol_mn_b + (sol_mn_a)/(sol_mn_k - T);


    end

else

    high_T_critcb4 = 1050;
    low_T_critcb4 = 800;
    crit_cb4_dim = crit_cb4*(max_mg-min_mg)+min_mg;

    max_cb_critcb4 = 22-(crit_cb4_dim - cb_dim)*2;
    min_cb_critcb4 = 12-(crit_cb4_dim - cb_dim)*2;

    k_critcb4 = 40.57971*cb_dim + 1025.362;

    b_crit_cb4 = ((k_critcb4 - high_T_critcb4)*max_cb_critcb4-(k_critcb4 -low_T_critcb4)*min_cb_critcb4)/(high_T_critcb4-low_T_critcb4);

    a_critcb4 = (k_critcb4 - high_T_critcb4)*max_cb_critcb4-(k_critcb4 - high_T_critcb4)*b_crit_cb4;

    ol_mg = ((a_critcb4/(k_critcb4-T)+b_crit_cb4)-min_mg)/(max_mg-min_mg);

    if T<Ts_local
        ol_mg = ((a_critcb4/(k_critcb4-Ts_local)+b_crit_cb4)-min_mg)/(max_mg-min_mg);
    end

    ol_mn = ol_mn_m*T +ol_mn_c;
    sol_mn = sol_mn_b + (sol_mn_a)/(sol_mn_k - T);

end

if sol_mn<0
    sol_mn=0;
end

if ol_mn<0
    ol_mn=0;
end
%% Orthopyroxene Mg content. 

if cb>cb_low_high

    if T< copx_T_m_2h*cb+copx_T_c_2h && cb< EC

        if copx_T_m_2h*cb + copx_T_c_2h >= ET
            opx_mg = s3m*(copx_T_m_2h*cb + copx_T_c_2h) + s3c; 
            
        else
            opx_mg = s4m*(copx_T_m_2h*cb+copx_T_c_2h)+s4c;
        end

    elseif T>=ET
        opx_mg = s3m*T + s3c;
        if T<Ts_local
            opx_mg = s3m*Ts_local + s3c;
        end
        %if copx_T_m_2h*cb+copx_T_c_2h >=ET && T<=crit_T2
        %        opx_mg=s3m*crit_T2+s3c;
       % end
    else
        opx_mg = s4m*T + s4c;
        if T<Ts_local
            opx_mg =  s4m*Ts_local + s4c;
        end
        %if copx_T_m_2h*cb+copx_T_c_2h >=ET
         %       opx_mg=s3m*crit_T2+s3c;
        
        %end
    end 
else

    if T< copx_T_m_2l*cb+copx_T_c_2l && cb< EC
        if copx_T_m_2l*cb + copx_T_c_2l >= ET
            opx_mg = s3m * (copx_T_m_2l*cb+copx_T_c_2l)+s3c;
        else
            opx_mg = s4m*(copx_T_m_2l*cb+copx_T_c_2l)+s4c;
        end
    elseif T>= ET
        opx_mg = s3m*T+s3c;
        if T<Ts_local
            opx_mg = s3m*Ts_local + s3c;
        end

    else
        opx_mg = s4m*T+s4c;
        if T<Ts_local
            opx_mg =  s4m*Ts_local + s4c;
        end
    end

end

%% Pyroxene Mg content 

px_k = 1230;

if cb> cb_low_high
    px_T = copx_T_m_1h*cb + copx_T_c_1h;
else
    px_T = copx_T_m_1l*cb + copx_T_c_1l;
end

px_cb = px_cb_m*px_T + px_cb_c;


if cb>liq_P2_C
    px_k = 1235;
    px_mg_b = ((px_k - px_T)*px_cb - (px_k - DT)*DC)/(DT - px_T);
    px_mg_a = (px_k - DT)*DC - px_mg_b*(px_k - DT);
else
    if cb<=crit_cb5
        px_k = 1100;
        lowest_start_T = liq_k3 - liq_a3/(cb - liq_b3);
        px_cs = s7m*T + s7c;
    else
        px_k_m = (1235 - 1150)/(liq_P2_C-crit_cb5);
        px_k_c = 1235 - px_k_m*liq_P2_C;
        px_k = px_k_m*cb + px_k_c;%1235;
        lowest_start_T = liq_k2 - liq_a2/(cb - liq_b2);
        px_cs = s7m*T + s7c;
    end

    px_mg_b = ((px_k - px_T)*px_cb - (px_k- lowest_start_T)*px_cs)/(lowest_start_T - px_T);
    px_mg_a = (px_k - lowest_start_T)*px_cs - px_mg_b*(px_k - lowest_start_T);

end


if T>=crit_T2 
    px_mg = s3m*T + s3c;
    if T<Ts_local
        px_mg = s3m*Ts_local + s3c;
    end
elseif cb>cb_low_high && T>=copx_T_m_1h*cb + copx_T_c_1h && T<crit_T2
    px_mg = px_mg_a /(px_k - T) + px_mg_b;

    if T<Ts_local
        px_mg = px_mg_a /(px_k - Ts_local) + px_mg_b;
    end

    if copx_T_m_1h*cb+copx_T_c_1h>=crit_T2
        px_mg = s3m*(crit_T2)+s3c;
    end

elseif cb<=cb_low_high && T>=copx_T_m_1l*cb+copx_T_c_1l && T<crit_T2
    px_mg = px_mg_a/(px_k - T) + px_mg_b;
    
    if T<Ts_local
        px_mg = px_mg_a /(px_k - Ts_local) + px_mg_b;
    end

    if copx_T_m_1l*cb+copx_T_c_1l<770 % && T<975
        px_mg = px_mg_a/(px_k - (copx_T_m_1l*cb + copx_T_c_1l)) + px_mg_b;
        %if px_mg<cs
            %px_mg = px_mg_a/(px_k - T) + px_mg_b;
        %end
    end
elseif cb>cb_low_high 
    px_mg = px_mg_a/(px_k - (copx_T_m_1h*cb+copx_T_c_1h)) + px_mg_b;
    if copx_T_m_1h*cb+copx_T_c_1h>=crit_T2
        px_mg = s3m*(crit_T2)+s3c;
        %if px_mg<cs
         %   px_mg = px_mg_a/(px_k - T) + px_mg_b;
        %end
    end
elseif cb<=cb_low_high
    px_mg = px_mg_a/(px_k - (copx_T_m_1l*cb + copx_T_c_1l)) + px_mg_b;
end


if px_mg <0
    px_mg =0;
end

%if px_mg < px_cb
   % px_mg=px_cb;
%end


if cb>cb_low_high
    temp_line10_cpx1 = copx_T_m_1h*cb + copx_T_c_1h;
    mgo_line10_cpx1 = copx_cb_m_1h*cb + copx_cb_c_1h;
else
    temp_line10_cpx1 = copx_T_m_1l*cb + copx_T_c_1l;
    mgo_line10_cpx1 = copx_cb_m_1l*cb+copx_cb_c_1l;
end

if mgo_line10_cpx1>s6m*temp_line10_cpx1+s6c
    mgo_line10_cpx1 = s6m*temp_line10_cpx1 +s6c; 
end


if cb>cb_low_high
    temp_line10_cpx2 = copx_T_m_2h*cb+copx_T_c_2h;
    mgo_line10_cpx2  = copx_cb_m_2h*cb + copx_cb_c_2h;
else
    temp_line10_cpx2 = copx_T_m_2l*cb+copx_T_c_2l;
    mgo_line10_cpx2 = copx_cb_m_2l*cb+copx_cb_c_2l;
end

s10m = (mgo_line10_cpx1 - mgo_line10_cpx2)/(temp_line10_cpx1 - temp_line10_cpx2);
s10c = mgo_line10_cpx1 - s10m*temp_line10_cpx1;

cpx_mg = s6m*T +s6c;
if T<Ts_local
    cpx_mg = s6m*Ts_local +s6c;
end

if px_mg <= s6m*T+s6c && T>=temp_line10_cpx1
    cpx_mg = px_mg;
else
    if T>=temp_line10_cpx1
        cpx_mg = s6m*T+s6c;
        if T<Ts_local
            cpx_mg = s6m*Ts_local+s6c;
        end
    elseif T<temp_line10_cpx1 && T>temp_line10_cpx2
        cpx_mg = s10m*T+s10c;
        if T<Ts_local
            cpx_mg = s10m*Ts_local+s10c;
        end
        
    else
        cpx_mg = s10m*temp_line10_cpx2+s10c;
    end

end


if opx_mg<px_mg
    opx_mg=px_mg;
end

if cpx_mg>px_mg
    cpx_mg=px_mg;
end

%% feldpsar
feld_mg = (0 - min_mg)/(max_mg-min_mg);


%% Minerals

%
if phi==1 
    ol = 0;
    opx = 0;
    cpx = 0;
    feld = 0;
    ol_mg=0;
    opx_mg=0;
    cpx_mg=0;
    feld_mg=0;
elseif phi_old>0 || (phi<1 && phi>0)

    if T>=crit_T1
        ol = 1 - phi;
        opx=0;
        cpx=0;
        feld=0;

    elseif T<crit_T1 && T>crit_T2

        if cb>crit_cb3
            ol = 1 - phi;
            opx=0;
            cpx=0;
            feld=0;
        elseif cb<crit_cb2 && cb>crit_cb4
            ol=0;
            opx = 1 - phi;
            cpx=0;
            feld=0;
        else
            ol = (cs - opx_mg)/(ol_mg-opx_mg)*(1-phi);
            opx = (ol_mg - cs)/(ol_mg-opx_mg)*(1-phi);
            if ol<0
                ol=0;
                opx=1-phi;
            end
        end

    else

        opx_mg_end = s3m*DT +s3c; 
        ol_mg_end = s1m*DT +s1c;
        phi_end = (D2C - cb)/(D2C - crit_cb1); 

        cb_olmn_start = (1/(max_mg-min_mg))*(-min_mg - ol_i_c0/ol_i_m0);

       % if cb>crit_cb3
        %    ol = 1 - phi ; 

        

        %else


            if ol_mn > 0 && T<=ol_mn_T1 && cb<cb_olmn_start
                ol = (0-sol_mn)/(0-ol_mn)*(1-phi);
            else
                if cb>crit_cb4
                    ol=(D2C - opx_mg_end)/(ol_mg_end-opx_mg_end)*(1-phi_end);%ol_end;
                else
                    ol=0;
                end

            end

            if cb >= s1m*DT +s1c
                ol=1-phi;
                opx=0;
                cpx=0;
                feld=0;

            elseif cb>=DC 
                opx = (ol_mg_end - cs)/(ol_mg_end - opx_mg_end)*(1-phi);
                ol = (cs - opx_mg_end)/(ol_mg_end - opx_mg_end)*(1-phi);
                if ol>1
                    ol=1-phi;
                    opx=0;
                end
                if ol<0
                    ol=0;
                    opx=1-phi;
                end
            else

                if ol<0
                    ol=0;
                end


                feld = ((px_mg - cs)/(px_mg - feld_mg))*(1-phi-ol);
                px = ((cs - feld_mg)/(px_mg - feld_mg))*(1-phi-ol);
    
                opx = (px_mg*px - cpx_mg*px)/(opx_mg - cpx_mg);
                cpx = (px_mg*px - opx_mg*px)/(cpx_mg - opx_mg);
                

            end

            if feld<0
                opx = opx + feld; 
                feld = 0;
            end

            if opx<0
                if cpx<=0 && feld>0
                    feld=feld+opx;
                    opx=0;
                elseif feld<=0 && cpx>0
                    cpx=cpx+opx;
                    opx=0;
                end
            end

            if cpx<0
                if opx<=0 && feld>0
                    feld=feld+cpx;
                    cpx=0;
                elseif feld<=0 && opx>0
                    opx=opx+cpx;
                    cpx=0;
                end
            end
                

            


       % end

       if opx_mg==0 && cpx_mg==0 && cb<DC
           opx=0;
           cpx=0;
           feld=1-phi-ol;
       end

       




      
    end

       if T <= Ts_local && cb<crit_cb5

           
           phi_end =  ((cb+0.01) - cb)/((cb+0.01) - cl); 
           
           ol = (0-sol_mn)/(0-ol_mn)*(1-phi_end);
           
           
        

           
           px = (((cb+0.01) - feld_mg)/(px_mg - feld_mg))*(1-phi_end-ol);
    
           opx = (px_mg*px - cpx_mg*px)/(opx_mg - cpx_mg);
           cpx = (px_mg*px - opx_mg*px)/(cpx_mg - opx_mg);

           feld = 1 - opx - px - ol - phi;
         
       end

       if cb==0
           ol=0;
           opx=0;
           cpx=0;
           feld=1-phi;
       end








end