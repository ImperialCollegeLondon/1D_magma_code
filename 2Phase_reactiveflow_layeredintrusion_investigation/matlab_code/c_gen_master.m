function [C]=c_gen_master(mu,d,phi,A,B,C_type,C_type2,KC_perm_MF_a, KC_perm_MF_b, KC_perm_MF_c, perm_min_FLAG, perm_min_MF)
%%
a=KC_perm_MF_a;%0.2;
b=KC_perm_MF_b;%0.6;
c=KC_perm_MF_c;%-5;

suspension_scale=1;

if perm_min_FLAG==1
    phi_f = max((phi-perm_min_MF)/(1-perm_min_MF),0);
else
    phi_f=phi;
end

if C_type2==1  %KC permeability
    if C_type==1

        if phi_f<1e-8
            C=A*mu./d^2*1e-8.^(2-B);
        elseif phi_f<=a
            C=A*mu./d^2*phi_f.^(2-B);
        elseif phi_f>=b
            C=(mu/d^2)*(1-phi_f).*phi_f.^c*suspension_scale;
        else
            Sl=-2;
            h2=exp(Sl./((phi_f-a)/(b-a)))./(exp(Sl./((phi_f-a)/(b-a)))+exp(Sl./(1-(phi_f-a)/(b-a))));
            C=(A*mu./d^2*phi_f.^(2-B)).*(1-h2)+(mu/d^2)*(1-phi_f).*phi_f.^c.*h2*suspension_scale;
        end
    else  %pure porous media type
        if phi_f<1e-8
            C=A*mu./d^2*1e-8.^(2-B);
        else
            C=A*mu./d^2*phi_f.^(2-B);
        end
    end
else   %RG permeability,   A=d^2/beta
    if phi_f<1e-8
        C=(1-1e-8)^2*mu/A/1e-8;
    else
        C=(1-phi_f)^2*mu/A/phi_f;
    end
end
end