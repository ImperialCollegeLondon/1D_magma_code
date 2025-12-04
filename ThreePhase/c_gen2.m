function [C]=c_gen2(mu,A,phi, S, type)
%%
tol=1e-2;

a=0.4;
b=0.6;
c=-5;
suspension_scale=1;
if phi<tol
    phi=tol;
end
if phi>0.99
    phi=0.99;
end
if type==1 % solid/liquid
    if S>1-tol
        S=1-tol;
    end
   
    if phi<=a
        C=(1-phi)^2*mu/A/phi/(1-S^2);
    elseif phi>=b
        C=(mu/3.5e-3^2)*(1-phi).*(phi.*(1-S)).^c*suspension_scale;
    else        
        Sl=-2;
        h2=exp(Sl./((phi-a)/(b-a)))./(exp(Sl./((phi-a)/(b-a)))+exp(Sl./(1-(phi-a)/(b-a))));
        C=(1-phi)^2*mu/A/phi/(1-S^2).*(1-h2)+(mu/3.5e-3^2)*(1-phi).*(phi.*(1-S)).^c.*h2*suspension_scale;
    end

else  %liquid/volotile
    if S<tol
        S=tol;
    end
    C=(1-phi)^2*mu/A/phi/S^2;
end
end


% h2=tan((x2-a)./(b-a).*pi/4);
% h2=linspace(0,1,length(x2));
% c3=c1(end)*(1-h2)+c2(1)*h2;




