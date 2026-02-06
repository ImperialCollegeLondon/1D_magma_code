function [Ts_local,Tl_local,Hs_local,Hl_local] = solid_state_Ts_Tl_local_master(HHJPet,SSPD,FourMPD, cb,cp,Lf,DT,DC,lm1,lc1, crit_cb4,crit_cb5,s5T,...
                                                    liq_P2_C, liq_P3_C, liq_k1,liq_a1,liq_b1,liq_k2,liq_a2,liq_b2,liq_k3,liq_a3,liq_b3,A1,B1,C1,alpha, n_PD)


%% JPet local solidus and liquidus
if HHJPet==1
    %local solidus
    Ts_local = A1+B1+C1;
       
    %local liquidus 
    Tl_local = A1*cb^2 + B1*cb +C1;
end


%% Solid solution local solidus and liquidus
if SSPD==1
    Ts_max=C1;
    Ts_min=A1+B1+C1;
    %local solidus 
    Ts_local = (alpha*(1-cb)^n_PD+(1-alpha)*(1-cb^(1/n_PD)))*(Ts_max-Ts_min)+Ts_min;

    %local liquidus 
    Tl_local = A1*cb^2+B1*cb+C1;
end


%% Layered intrusion phase diagram

if FourMPD==1
    % local solidus
    if cb>= (DT-lc1)/lm1 
        Ts_local=lm1*cb+lc1;
    elseif cb>=DC
        Ts_local=DT;
    elseif cb>=crit_cb4
        Ts_local=s5T;
    elseif cb>=crit_cb5
        Ts_local=800;
    else
        Ts_local=800;
    end
    
    
    %local liquidus
    if cb>=liq_P2_C
        Tl_local=liq_k1-(liq_a1/(cb-liq_b1)); 
    elseif cb>=liq_P3_C
        Tl_local=liq_k2-(liq_a2/(cb-liq_b2));
    else
        Tl_local=liq_k3-(liq_a3/(cb-liq_b3));
    end
end


Hs_local=Tl_local*cp+Lf;
Hl_local=Ts_local*cp;



end