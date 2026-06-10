function [phi,T, cl, cs, type,Tl_local,Ts_local]=poro_component_solve_4_components_v2_master(cb, ent, Lf, cp, min_mgo, max_mgo, crit_T1, crit_T2,...
                                                             crit_cb2, crit_cb3, crit_cb4,crit_cb4_OG, crit_cb5,liq_k1, liq_a1, liq_b1,liq_k2, liq_a2, liq_b2, liq_k3, liq_a3, liq_b3,...
                                                            lm1, lc1, lm3, lc3, lm8, lc8, mk, ck, BT, BC, DT, DC,JC, olpx_C,olpxm, olpxc, s5T, sol_k2,liq_P2_C,sm1,sc1, ...
                                                            lowest_end_T, liq_P3_C,l13_k,l13_a,l13_b,lowest_k,s8m,s8c,crit_T3,liq_P4_C,tol,step_size)
%% A solid solution 2-phase phase diagram - four components
% Catherine Booth    
% 12-03-2024


% liquidus is defined
% solidus is defined 
% cb: bulk composition
% ent: enthalpy per unit kg
% solve for melt fraction, temperature, liquid and solid component 
    

    Option5=0;

    iter=0;
    type=0;

   % step_size=0.35;

    eps=0;
    cb=min(cb,1);
    cb=max(0,cb);

    cb_err=0;
    cb_err=0;
    cl_err=0;
    cs_err=0;
    t_err=0;
    phi_err=0;
   
    cl_error=0;
    cs_error=0;
    phi_error=0;
    T_error=0;
    Option=0;
    Option1=0;

    ts_m_loc = (s5T - (lm8*crit_cb4+lc8))/(0.05);
    ts_c_loc = (lm8*crit_cb4+lc8)-crit_cb4*ts_m_loc;
    
    % local solidus
    if cb>= (DT-lc1)/lm1 
        Ts_local=lm1*cb+lc1;
    elseif cb>=DC
        Ts_local=DT;
    elseif cb>=crit_cb4
              
        Ts_local=s5T;
    elseif cb>=crit_cb5
        Ts_local=lowest_end_T;

    else

        Ts_local=lowest_end_T;
    end
   


    %local liquidus
    if cb>=liq_P2_C
        Tl_local=liq_k1-(liq_a1/(cb-liq_b1)); 
    elseif cb>=liq_P3_C
        Tl_local=liq_k2-(liq_a2/(cb-liq_b2));
    else
        Tl_local=liq_k3-(liq_a3/(cb-liq_b3));
    end

    

    min_cl = liq_P4_C;
    
    ETl=Tl_local*cp+Lf;
   
    Ets=Ts_local*cp;

    % Eutetic

    if cb>=DC && cb<(DT-lc1)/lm1  % between these comps there is a eutetic

        if cb>crit_cb3
            cs_temp = sm1*Ts_local +sc1;
        else
            DC_dim = DC*(max_mgo-min_mgo)+min_mgo;
    
            crit_cb2_dim = crit_cb2*(max_mgo-min_mgo)+min_mgo;
            cb_dim = cb*(max_mgo-min_mgo)+min_mgo;
    
            D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);
            D2C = (D2C_dim - min_mgo)/(max_mgo-min_mgo);
    
            sol_k1 = mk*cb + ck;
            sol_b1=(BC*(sol_k1-BT)-D2C*(sol_k1-DT))/(DT-BT);
            sol_a1= (BC - sol_b1)*(sol_k1-BT);
        
            cs_temp=sol_a1/(sol_k1-Ts_local)+sol_b1;
        end

        
        cl_temp=liq_a1/(liq_k1-Ts_local)+liq_b1;
        phi_temp=(cb-cs_temp)/(cl_temp-cs_temp);
        if phi_temp<0 %new
            phi_temp=0;
        elseif phi_temp>1
            phi_temp=1;
        end
        ETe=Ts_local*cp+Lf*(phi_temp);
%% NEW
    elseif cb>=crit_cb4 && cb<DC

        DC_dim = DC*(max_mgo-min_mgo)+min_mgo;

        crit_cb2_dim = crit_cb2*(max_mgo-min_mgo)+min_mgo;
        cb_dim = cb*(max_mgo-min_mgo)+min_mgo;

        D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);

        if cb>crit_cb2

            D2C = (D2C_dim - min_mgo)/(max_mgo-min_mgo);
        else
            D2C = (DC_dim - min_mgo)/(max_mgo-min_mgo);
        end

        sol_b2 = ((sol_k2-DT)*D2C-(sol_k2-s5T)*(cb+0.01))/(s5T-DT);
        sol_a2 = (sol_k2-DT)*(D2C - sol_b2);

        cs_temp = sol_a2/(sol_k2-Ts_local)+sol_b2;
        cl_temp = liq_a3/(liq_k3-Ts_local)+liq_b3;
        %if cl_temp<0
        %    cl_temp=0;
       % end
        %disp('cb')
        
        %disp(cb*(56.7173-0.0591)+0.0591)
        %disp(liq_a3)
        
        %disp(liq_k3)
        %disp(ent)
        %disp(Ts_local)
        %disp(liq_b3)
        %disp('end')
        phi_temp = (cb-cs_temp)/(cl_temp-cs_temp);
        if phi_temp<0 %new
            phi_temp=0;
        elseif phi_temp>1
            phi_temp=1;
        end
        ETe=Ts_local*cp+Lf*(phi_temp);

    elseif cb<crit_cb4 && cb>=crit_cb5
        
        lowest_start_T =liq_k2 - (liq_a2/(cb - liq_b2));
        if cb>=crit_cb4_OG
            lowest_m = (JC -DC)/(crit_cb4_OG - crit_cb4);
            lowest_c = DC - lowest_m*crit_cb4;
            lowest_start_cs = lowest_m*cb+lowest_c;
        else
            lowest_start_cs = s8m*lowest_start_T +s8c;
        end

        low_k=sol_k2;
        lowest_end_cs = cb + 0.01;
                    
        lowest_end_T1 = 850;
        lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T1)*lowest_end_cs)/(lowest_end_T1 - lowest_start_T);
        lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);

        cs_temp = lowest_a/(low_k-Ts_local)+lowest_b;
        cl_temp = liq_a3/(liq_k3-Ts_local)+liq_b3;
        %if cl_temp<0
          %  cl_temp=0;
       % end

        phi_temp = (cb-cs_temp)/(cl_temp-cs_temp);
        if phi_temp<0 %new
            phi_temp=0;
        elseif phi_temp>1
            phi_temp=1;
        end
        ETe=Ts_local*cp+Lf*(phi_temp);
        
                  


    elseif cb<crit_cb5
       
        lowest_start_T =liq_k3 - (liq_a3/(cb - liq_b3));
        lowest_start_cs = l13_a/(l13_k - lowest_start_T) +l13_b;
        if lowest_start_cs<=0.01
            lowest_start_cs=0.01;
        end
        low_k=lowest_k;
        lowest_end_cs = cb + 0.01;
        lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T)*lowest_end_cs)/(lowest_end_T - lowest_start_T);
        lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);
        
        if lowest_start_T==lowest_end_T
            cs_temp = cb+0.01;
        else
         cs_temp = lowest_a/(low_k-Ts_local)+lowest_b;
        end
        cl_temp = liq_a3/(liq_k3-Ts_local)+liq_b3;
        %if cl_temp<0
          %  cl_temp=0;
        %end

        phi_temp = (cb-cs_temp)/(cl_temp-cs_temp);
        if phi_temp<0 %new
            phi_temp=0;
        elseif phi_temp>1
            phi_temp=1;
        end
        ETe=Ts_local*cp+Lf*(phi_temp);
        
        
        

%% end NEW
    else
        ETe=0;
    end
    
    if Ets>=ETl
        phi=0;
        cs=cb;
        cl=0;
        T=(ent-0*Lf)/cp;
        type=7;
        cb_err=0;
        cl_err=0;
        cs_err=0;
        t_err=0;
        phi_error=0;
        phi_err=0;
        cl_error=0;
       cs_error=0;
        phi_error=0;
        T_error=0;

    elseif(ent>=ETl)
        phi=1;        
        cs=1;
        cl=cb;
        T=(ent-1*Lf)/cp;
        type=0;
        
        
    elseif(ent<=Ets) 

        phi=0;
        cs=cb;
        cl=0;
        T=(ent-0*Lf)/cp;
        type=1;
        cb_err=0;
        cl_err=0;
        cs_err=0;
        t_err=0;
        phi_error=0;
        phi_err=0;
        cl_error=0;
        cs_error=0;
        phi_error=0;
        T_error=0;

    

    elseif (ent>Ets) && (ent<=ETe)
        
        phi=phi_temp*(ent-Ets)/(ETe-Ets);
        cl=cl_temp;%0;
        if cl<0
            cl=0;
        end

        cs=(cb-cl*phi)/(1-phi);
        

        %(cb-phi*phi_temp)/(1-phi);
        T=(ent-phi*Lf)/cp;
      
        
        type=4;
        
        cb_err=0;
        cl_err=0;
        cs_err=0;
        t_err=0;
        phi_error=0;
        phi_err=0;
        cl_error=0;
        cs_error=0;
        phi_error=0;
        T_error=0;

    else
        err=1;        
        %initial guess of [cl cs phi T];
        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
        x=[cl cs phi T]';
        h=ent;
        i=0;
        
        Option=1; %for cb>=crit_cb1 && cb<=crit_cb2, try first region
        if cb<crit_cb5
            Option=5;
            cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
            x=[cl cs phi T]';
        elseif cb<crit_cb4
            Option=4;
            %step_size=0.15;
            %if cb>=0.135
             %   step_size=0.075;
            %end

        end

        if cb<DC && cb>=olpx_C
            %step_size=0.2;
        end
       
        Option1=0;
        Y=0;
        sub_iter=0;
        while (err>tol)
            if  Option==1 %T>=crit_T1 
                LHS = reshape([phi, 0, liq_a1/(cl-liq_b1)^2, 0, 1-phi, 0, 0, lm1, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k1-(liq_a1/(cl-liq_b1))-T; lm1*cs+lc1-T];
            elseif Option==2 || Option==3  %^crit_T2<=T && T<crit_T1
                if cb<crit_cb2
                    LHS = reshape([phi, 0, liq_a1/(cl-liq_b1)^2, 0, 1-phi, 0, 0, lm3, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k1-(liq_a1/(cl-liq_b1))-T; lm3*cs+lc3-T]; 
                elseif cb>crit_cb3
                    LHS = reshape([phi, 0, liq_a1/(cl-liq_b1)^2, 0, 1-phi, 0, 0, lm1, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k1-(liq_a1/(cl-liq_b1))-T; lm1*cs+lc1-T];
                else%if cb>=crit_cb2 && cb<=crit_cb3
                    
                    if Option==2
                        DC_dim = DC*(max_mgo-min_mgo)+min_mgo;

                        crit_cb2_dim = crit_cb2*(max_mgo-min_mgo)+min_mgo;
                        cb_dim = cb*(max_mgo-min_mgo)+min_mgo;

                        D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);
                        D2C = (D2C_dim - min_mgo)/(max_mgo-min_mgo);

                        sol_k1 = mk*cb + ck;
                        sol_b1=(BC*(sol_k1-BT)-D2C*(sol_k1-DT))/(DT-BT);
                        sol_a1= (BC - sol_b1)*(sol_k1-BT);

                        LHS = reshape([phi, 0, liq_a1/(cl-liq_b1)^2, 0, 1-phi, 0, 0, sol_a1/(cs-sol_b1)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                        RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k1-(liq_a1/(cl-liq_b1))-T; sol_k1-(sol_a1/(cs-sol_b1))-T];
                    end
                end

                if Option==3 %Option = 3
                    %                         if cb<=olpx_C && T<=bm*cb+bc
                    LHS = reshape([phi, 0, liq_a1/(cl-liq_b1)^2, 0, 1-phi, 0, 0, lm3, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k1-(liq_a1/(cl-liq_b1))-T; lm3*cs+lc3-T];
                    %                         end
                end
                
            elseif Option==4  % T<crit_T2
               % if cb>=DC %this might never be reached but keep in case inputs change
                %    LHS = reshape([phi, 0, liq_a2/(cl-liq_b2)^2, 0, 1-phi, 0, 0, lm3, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                 %   RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k2-(liq_a2/(cl-liq_b2))-T; lm3*cs+lc3-T]; 
               % else
                
                 if cb>=DC

                    LHS = reshape([phi, 0, liq_a2/(cl-liq_b2)^2, 0, 1-phi, 0, 0, lm3, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k2-(liq_a2/(cl-liq_b2))-T; lm3*cs+lc3-T];


                 elseif cb>=crit_cb4  
                  

                    DC_dim = DC*(max_mgo-min_mgo)+min_mgo;

                    crit_cb2_dim = crit_cb2*(max_mgo-min_mgo)+min_mgo;
                    cb_dim = cb*(max_mgo-min_mgo)+min_mgo;

                    D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);

                    if cb>crit_cb2

                        D2C = (D2C_dim - min_mgo)/(max_mgo-min_mgo);
                    else
                        D2C = (DC_dim - min_mgo)/(max_mgo-min_mgo);
                    end

                    sol_b2 = ((sol_k2-DT)*D2C-(sol_k2-s5T)*(cb+0.01))/(s5T-DT);
                    sol_a2 = (sol_k2-DT)*(D2C - sol_b2);

                    LHS = reshape([phi, 0, liq_a2/(cl-liq_b2)^2, 0, 1-phi, 0, 0, sol_a2/(cs-sol_b2)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k2-(liq_a2/(cl-liq_b2))-T; sol_k2-(sol_a2/(cs-sol_b2))-T];

                 elseif cb>=crit_cb5

                    
                    
                    lowest_start_T =liq_k2 - (liq_a2/(cb - liq_b2));
                    if cb>=crit_cb4_OG
                        lowest_m = (JC -DC)/(crit_cb4_OG - crit_cb4);
                        lowest_c = DC - lowest_m*crit_cb4;
                        lowest_start_cs = lowest_m*cb+lowest_c;
                    else
                        lowest_start_cs = s8m*lowest_start_T +s8c;
                    end
                    low_k=sol_k2;
                    lowest_end_cs = cb + 0.01;
                    lowest_end_T1 = 850;
                    

                    lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T1)*lowest_end_cs)/(lowest_end_T1 - lowest_start_T);
                    lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);
                    LHS = reshape([phi, 0, liq_a2/(cl-liq_b2)^2, 0, 1-phi, 0, 0, lowest_a/(cs-lowest_b)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k2-(liq_a2/(cl-liq_b2))-T; low_k-(lowest_a/(cs-lowest_b))-T];

                 end

                    %lowest_start_T =liq_k3 - (liq_a3/(cb - liq_b3));
                    %lowest_start_cs = l13_a/(l13_k - lowest_start_T) +l13_b;
                    %low_k=lowest_k;
                    %lowest_end_cs = cb + 0.01;

                    %lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T)*lowest_end_cs)/(lowest_end_T - lowest_start_T);
                    %lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);
                   % LHS = reshape([phi, 0, liq_a2/(cl-liq_b2)^2, 0, 1-phi, 0, 0, lowest_a/(cs-lowest_b)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                  %  RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k2-(liq_a2/(cl-liq_b2))-T; low_k-(lowest_a/(cs-lowest_b))-T];
                 %end

            elseif Option==5
                if cb>=DC

                    LHS = reshape([phi, 0, liq_a3/(cl-liq_b3)^2, 0, 1-phi, 0, 0, lm3, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k3-(liq_a3/(cl-liq_b3))-T; lm3*cs+lc3-T];


                 elseif cb>=crit_cb4  
                  

                    DC_dim = DC*(max_mgo-min_mgo)+min_mgo;

                    crit_cb2_dim = crit_cb2*(max_mgo-min_mgo)+min_mgo;
                    cb_dim = cb*(max_mgo-min_mgo)+min_mgo;

                    D2C_dim = DC_dim + (cb_dim-crit_cb2_dim);

                    if cb>crit_cb2

                        D2C = (D2C_dim - min_mgo)/(max_mgo-min_mgo);
                    else
                        D2C = (DC_dim - min_mgo)/(max_mgo-min_mgo);
                    end

                    sol_b2 = ((sol_k2-DT)*D2C-(sol_k2-s5T)*(cb+0.01))/(s5T-DT);
                    sol_a2 = (sol_k2-DT)*(D2C - sol_b2);

                    LHS = reshape([phi, 0, liq_a3/(cl-liq_b3)^2, 0, 1-phi, 0, 0, sol_a2/(cs-sol_b2)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k3-(liq_a3/(cl-liq_b3))-T; sol_k2-(sol_a2/(cs-sol_b2))-T];
                    
                 elseif cb>=crit_cb5

                    
                    
                    lowest_start_T =liq_k2 - (liq_a2/(cb - liq_b2));
                    if cb>=crit_cb4_OG
                        lowest_m = (JC -DC)/(crit_cb4_OG - crit_cb4);
                        lowest_c = DC - lowest_m*crit_cb4;
                        lowest_start_cs = lowest_m*cb+lowest_c;
                    else
                        lowest_start_cs = s8m*lowest_start_T +s8c;
                    end
                    low_k=sol_k2;
                    lowest_end_cs = cb + 0.01;
                    
                    lowest_end_T1 = 850;
                    lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T1)*lowest_end_cs)/(lowest_end_T1 - lowest_start_T);
                    lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);
                    LHS = reshape([phi, 0, liq_a3/(cl-liq_b3)^2, 0, 1-phi, 0, 0, lowest_a/(cs-lowest_b)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k3-(liq_a3/(cl-liq_b3))-T; low_k-(lowest_a/(cs-lowest_b))-T];
                   
                 else

                    lowest_start_T =liq_k3 - (liq_a3/(cb - liq_b3));
                    lowest_start_cs = l13_a/(l13_k - lowest_start_T) +l13_b;
                    if lowest_start_cs<=0.01
                        lowest_start_cs=0.01;
                    end
                    low_k=lowest_k;
                    lowest_end_cs = cb + 0.01;

                    lowest_b = ((low_k-lowest_start_T)*lowest_start_cs - (low_k - lowest_end_T)*lowest_end_cs)/(lowest_end_T - lowest_start_T);
                    lowest_a = (low_k - lowest_start_T)*lowest_start_cs - lowest_b*(low_k - lowest_start_T);
                    LHS = reshape([phi, 0, liq_a3/(cl-liq_b3)^2, 0, 1-phi, 0, 0, lowest_a/(cs-lowest_b)^2, cl-cs, Lf, 0, 0, 0, cp, -1, -1],[4,4]);
                    RHS = [cl*phi+cs*(1-phi)-cb; cp*T+Lf*phi-h; liq_k3-(liq_a3/(cl-liq_b3))-T; low_k-(lowest_a/(cs-lowest_b))-T];
                 end

                

            end
            i=i+1;
           % disp(LHS)
            %disp(RHS)
            %disp(cb)
           % disp(phi)
            %disp(T)
            %disp(cs)
            %disp(cl)
            %disp(Option)
            %% 
            dx=LHS\RHS;
            
            %if Y==0
                x=x-dx*step_size;
            %else
              %  w=0.05;
             %   x=(1-w)*x-dx*step_size*w;
            %end
            
            cl=x(1); cs=x(2); phi=x(3); T=x(4);

            

             err=max([abs(dx(1:3));  abs(dx(4))/100]);

            
             
             cl_oh=isnan(cl);
             cs_oh=isnan(cs);
             phi_oh=isnan(phi);
             T_oh=isnan(T);

            %% Going through the options, dependent on the bulk comp.
            %check if they're used by Option1

             if cb<(DT-lc1)/lm1 && cb>= DC
                 if err<tol && Option==1 && T<crit_T1
                     Option=2;
                     Option1=1;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

             elseif cb>=olpx_C

                 



                if err<tol && Option==1 && T<crit_T1
                     Option=2;
                     Option1=3;
                     
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if err<tol && Option==2 && (T<crit_T2 || phi<0)
                     Option1=4;
                    if phi<0 
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                    end
                    if cl<min_cl ||cl>cs || cs>cb+0.01  
                        Option1=5;
                        cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                    end
                     Option =4; 
                     
                     err=1e5;
                     sub_iter=0;
                     continue
                 end


                 

                 if err<tol && Option==4 && (T<crit_T3 )

                     if phi<0 
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                    end
                    if cl<min_cl ||cl>cs || cs<cb+0.01  
                        Option1=5;
                        cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                    end
                     Option1=6;
                    
                     Option =5; 
                     
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if Option==4 &&(cl>cs || cs>1) && phi<0 %(cl>cs) %&& T<crit_T3%(cl<min_cl)%||cs<cb+0.01) %&& T<crit_T3  cs<cb+0.01||
                     Option1 = 60;
                     if phi<0 
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                     end
                     if cl<min_cl || cs<cb+0.01  
                       
                        cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     end
                     Option =5;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if  Option==2 && (cl>cs||cl<min_cl) 
                     
                      if  cs>cb+0.01  && cl<=min_cl %&&&& cs>1
                        cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                        Option1=8;
                      else
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                        Option1=9;
                      end
                    

                     Option=4;
                     
                     err=1e5;
                     sub_iter=0;
                     continue
                 end


                 if Option5==0 &&  Option==5 && (phi<0 ||cl<0)
                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=40;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     Option5=1;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end  

                 

                 

                
                
                 


                 

               
                         
                
             elseif cb>=crit_cb2

                if err<tol && Option==1 && T<crit_T1
                     Option=2;
                     Option1=12;
                     err=1e5;
                     sub_iter=0;
                     continue
                end
                

           

                 if err<tol && Option==2 && (T<olpxm*cb+olpxc )%|| cl>cs)
                     if cl>cs
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                     end

                     Option=3;
                     Option1=15;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end


                


                 if err<tol && Option==3 && T<crit_T2 

                     Option=4;
                     Option1=16;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end


                 if  Option==2 && (cl>cs||cl<min_cl||cs>1) 
                     
                      if cl<=min_cl && cs>cb+0.01
                        cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                        Option1=13;
                      else
                        cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                        Option1=14;
                      end

                     Option=3;
                     
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if err<tol && Option==4 && T<crit_T3 

                     Option=5;
                     Option1=19;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 
                 if Option5==0 &&  Option==5 && (phi<0 ||cl<0)
                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=40;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     Option5=1;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end    

                











             elseif cb>=crit_cb4

                 

                 if err<tol && Option==1 && T<crit_T1
                     Option=2;
                     Option1=19;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 

                 if err<tol && Option==2 && T<crit_T2 
                     Option =4; 
                     Option1=20;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if err<tol && Option==4 && T<crit_T3
                     Option =5; 
                     Option1=22;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end
                
                if Option5==0 && Option==5 && (phi<0 ||cl<0)
                     

                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=40;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     Option5=1;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end
            
                

                 

                 

             elseif cb>=crit_cb5
                
                % if err<tol && Option==1 && T<crit_T1
                   %  Option=2;
                  %   Option1=23;
                   %  err=1e5;
                  %   sub_iter=0;
                   %  continue
               %  end

                % if  Option==2 && (cl>cs||cl<min_cl||cs>1) 
                     
                    %  if cl<=min_cl && cs>cb+0.01
                    %    cl=0.5; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     %   Option1=13;
                    %  else
                    %    cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                    %    Option1=14;
                    %  end

                   %  Option=4;
                     
                   %  err=1e5;
                  %   sub_iter=0;
                  %   continue
                % end


                % if err<tol && Option==2 && T<crit_T2 

                    % cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                    % Option =4;
                   %  err=1e5;
                   %  sub_iter=0;
                    % continue
                 %end

                 if  err<tol && Option==4 && (T<crit_T3 ) % 

                     Option1=25;

                     %if cl<=min_cl || phi <0 || cs<cb+0.01
                        %cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                        %Option1=14;
                      %else
                       % cl=0.5; cs=0.5; phi=0.5; T=(Ts_local+Tl_local)/2;
                        %Option1=14;
                      %end
                                          
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     continue
                 end

                 if Option==4 && (isnan(phi)||isnan(cs)||isnan(cl)||isnan(T))
                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=34;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end

                 if Option5==0  && Option==5 && (cl<0||phi<0||cs<0||isnan(phi))
                     
                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=40;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     Option5=1;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end


                
                 if err<tol && Option==4 && (cl<0||phi<0||cs<0)
                     cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                     x=[cl cs phi T]';
                     h=ent;
                     Option1=38;
                     Option=5;
                     err=1e5;
                     sub_iter=0;
                     %step_size=0.05;
                     %showme=1;
                     continue
                 end



                %if  Option== &&  phi<0 && (cl>cs || cs>1) 
                    
                    %Option1 = 34;
                  %  phi=0.1;
                   % cs=cb+0.01;
                  %  cl=0;
                   % T=(Ts_local+Tl_local)/2;
                   
                     %Option=5;
                    % 
                     %err=1e5;
                   %  sub_iter=0;
                     %continue
                 %end
                
                 

               

                

             else
                 %Option=5;
                 %err=1e5;
                 %sub_iter=0;
             end


            
            
            
       

          


             

             
             
              %disp(err)
             
            iter=iter+1;
            sub_iter=sub_iter+1;
            type=3;
            if sub_iter>100 && Option~=5 
                Option1=27;
                if cb>=olpx_C
                    if Option==1
                        Option=2;
                    elseif Option==2
                        Option=4;
                    elseif Option==4
                        Option=5;
                        
                    end
                else
                    Option=Option+1;
                    
                end
                
                

                sub_iter=0;
                
            elseif sub_iter>250 && Option ==5 && Y==0
                Y=1;
                sub_iter=0;
                showme=1;
               % step_size=0.1;
                %disp('HERE')
                
                disp(cb*(max_mgo-min_mgo)+min_mgo)
                cl=min_cl; cs=cb+0.01; phi=0.5; T=(Ts_local+Tl_local)/2;
                x=[cl cs phi T]';
                err=1e5;
                %disp('moved')
             
            elseif sub_iter>250 && Option==5 && Y==1
                Option=Option+1;
                if Option==6
                    disp(err)
                    error('Chemical solution failure!') %Always raise this flag if it fails to converge
                end
            end

            




            if i>1000
                aaa=i;
                bbb=err;
                cb_err=cb;
                cl_err=cl;
                cs_err=cs;
                t_err=T;
                phi_err=phi;
                cl_error=dx(1);
                cs_error=dx(2);
                phi_error=dx(3);
                T_error=dx(4)/100;
                Option_err = Option;
                disp(cb*(max_mgo-min_mgo)+min_mgo)
                disp('Past iteration cut off ')
                disp(cb)
                h=h+1;
                break
            else
                cb_err=0;
                cl_err=0;
                cs_err=0;
                t_err=0;
                phi_err=0;
               
                cl_error=0;
                cs_error=0;
                phi_error=0;
                T_error=0;
            end
        end
        
    end


                





    

end
