% HH 2/2/2020
% solving melt fraction, temperature and composition based on the phase diagram of (Hu et al, 2022)
function [por1,T,FC, SC,region,Tl,Ts]=poro_component_solve_JPET_master(bc, ent,A1, B1, C1, A2, B2, C2, ae, Lf, cp,tol)
    % solve for melt fraction, temperature, liquid and solid component based on a
    % eutectic phase diagram
  
    iter=0;
    region=0;
    bc=min(bc,1);
    bc=max(0,bc);

    Ts=A1+B1+C1; %solidus
    Tl=A1*bc^2+B1*bc+C1; %liquidus
    
    ETl=Tl*cp+Lf;
    Ets=Ts*cp;
    if(ent>ETl)
        test=1;
        region=3;
        if(bc>ae)
            SC=1;
        else
            SC=0;
        end
        FC=bc;
        T=(ent-test*Lf)/cp;
    elseif(ent<Ets)
        test=0;
        SC=bc;
        FC=ae;
        T=(ent-test*Lf)/cp;
    else      
        if(bc<=ae)
            He = Ets + Lf*bc/ae;
            if (ent<=He)
                test=(ent-Ets)/Lf;
                region=1;
                FC=ae;
                SC=(bc-test*ae)/(1-test); 
                SC=min(1,SC);
                SC=max(0,SC);
                T=(ent-test*Lf)/cp;
            else
                region=2;
                test=1;
                test_pre=0;
                while (abs(test_pre-test)>tol)
                   fx=(Lf/cp)*test^3+(C1-ent/cp)*test^2+B1*bc*test+A1*bc^2;
                   fdashx=3*Lf/cp*test^2+2*(C1-ent/cp)*test+B1*bc;
                   test_pre=test;
                   test=test-fx/fdashx;    
                   iter=iter+1;
                end
                SC=0;
                T=(ent-test*Lf)/cp;
                FC=(-B1-sqrt(B1^2-4*A1*(C1-T)))/2/A1; 
            end
        elseif(bc>ae)
               He = Ets + Lf*(1-bc)/(1-ae);
               if (ent<=He)
                 test=(1-bc)*(ent - Ets)/(He - Ets)/(1-ae);
                 region=1;
               else
                   region=2;
                   test=1;
                   test_pre=0;
                   while(abs(test_pre-test)>tol)
                       fx=(Lf/cp)*test^3+(C2-ent/cp+A2+B2)*test^2+(2*A2*bc-2*A2+B2*bc-B2)*test+A2*(bc-1)^2;
                       fdash=3*(Lf/cp)*test^2+2*(C2-ent/cp+A2+B2)*test+(2*A2*bc-2*A2+B2*bc-B2);
                       test_pre=test;
                       test=test-fx/fdash;
                       iter=iter+1;
                   end
                   SC=1;
                   T=(ent-test*Lf)/cp;
                   FC=(-B2+sqrt(B2^2-4*A2*(C2-T)))/2/A2;   
               end
        end   
    end   
    por1=test;   
por1=real(por1);
por1=min(por1,1);
por1=max(por1,0);


% if por1<Precision
%     por1=0;
% end
% 
% if por1>1-Precision
%     por1=1;
% end
end
