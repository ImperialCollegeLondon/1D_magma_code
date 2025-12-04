function [phi,T, cl, cs, iter,Tl_local,Ts_local]=poro_component_solve_solid_master(cb, ent, A1, B1, C1, alpha, n_PD, Lf, cp,tol,step_size)
%% A solid solution 2-phase phase diagram
% Haiyang Hu    
% 12-07-2023
% liquidus is defined as a quadratic function T=A1*Cl^2+B1*Cl+C1
% solidus is defined as a function T=(alpha*(1-x).^n+(1-alpha)*(1-x.^(1/n)))*(Ts_max-Ts_min)+Ts_min where Ts Tl are the min and max of solidus/liquidus
% cb: bulk composition
% ent: enthalpy per unit kg
% solve for melt fraction, temperature, liquid and solid component 

    iter=0;

    %step_size=0.3;
    eps=0;
    cb=min(cb,1);
    cb=max(0,cb);
    
    Ts_max=C1;
    Ts_min=A1+B1+C1;
    Ts_local=(alpha*(1-cb)^n_PD+(1-alpha)*(1-cb^(1/n_PD)))*(Ts_max-Ts_min)+Ts_min; %local solidus    
    Tl_local=A1*cb^2+B1*cb+C1; %local liquidus
    
    ETl=Tl_local*cp+Lf;
    Ets=Ts_local*cp;
    if(ent>=ETl)
        region=3;
        phi=1;        
        cs=0;
        cl=cb;
        T=(ent-1*Lf)/cp;
    elseif(ent<=Ets)
        region=2;
        phi=0;
        cs=cb;
        cl=1;
        T=(ent-0*Lf)/cp;
    else
        region=2;
        err=1;        
        %initial guess of [cl cs phi T];
        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
        x=[cl cs phi T]'; 
        h=ent;

        while (err>tol)
            LHS=reshape([phi,0.0,B1+A1.*cl.*2.0,0.0,-phi+1.0,0.0,0.0,(Ts_min-Ts_max).*(alpha.*n_PD.*(-max(cs,eps)+1.0).^(n_PD-1.0)-(max(cs,eps).^(1.0./n_PD-1.0).*(alpha-1.0))./n_PD),cl-cs,Lf,0.0,0.0,0.0,cp,-1.0,-1.0],[4,4]);
            rhs=[-cb+cl.*phi-cs.*(phi-1.0);-h+T.*cp+Lf.*phi;C1-T+B1.*cl+A1.*cl.^2;-T+Ts_min-((alpha-1.0).*(max(cs,0).^(1.0./n_PD)-1.0)+alpha.*(-cs+1.0).^n_PD).*(Ts_min-Ts_max)];
            dx=LHS\rhs;
            err=max([abs(dx(1:3));  abs(dx(4))/100]);
            x=x-dx*step_size;
            cl=x(1); cs=x(2); phi=x(3); T=x(4);
            iter=iter+1;
        end
    end
end
