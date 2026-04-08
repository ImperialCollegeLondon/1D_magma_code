
i=3054;
            a0=dSdT(i)/100;
            b0=b0_all(i);
D1=0.2;
cap_A=S_cap(1);
cap_B=S_cap(2);
eps=1e-15;

            column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3, i+(N+1)*2+N*4,i+(N+1)*2+N*5, i+(N+1)*2+N*6];
            in=X(column_index);

            phi_1=in(1); 
            T_1=in(2); 
            cs_1=in(3);
            cl_1=in(4);
            S_1=in(5); 
            cs2_1=in(6); 
            cl2_1=in(7);
        
Sats=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1/(1-phi_1+eps);
    T_capped=(T_1+500+sqrt((T_1-500)^2+eps))/2;
    Satl=a0*T_capped+b0-cl2_1/(phi_1+eps);

    % Par_l=cl2_1/(phi_1+eps)*D1;
    % Par=cs2_1/(1-phi_1+eps)-((cap_A*cs_1+cap_B*(1-cs_1))+Par_l-sqrt(((cap_A*cs_1+cap_B*(1-cs_1))-Par_l+eps)^2))/2;

    R=(phi_1/(phi_1+eps)*Satl+(1-phi_1/(phi_1+eps))*Sats);

    constraint=sqrt(S_1^2+R^2+eps)-(S_1+R)


            %%
            phi=linspace(0,1e-8,1e3);
            y=phi./(phi+1e-12);
            figure(9)
            clf; plot(phi,y)

            %%
i=2940;
column_index=[i+(N+1)*2, i+(N+1)*2+N*2, i+(N+1)*2+N*5, i+(N+1)*2+N*6 ];
            in=X(column_index);

            phi_1=in(1); 
            cs_1=in(2); 
            cs2_1=in(3); 
            cl2_1=in(4);

D1=0.2;
cap_A=S_cap(1);
cap_B=S_cap(2);
eps=1e-12;

beta=1e4;
            Sat=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1+0.02;
    Par=cl2_1*D1-cs2_1;
    Par=-1/beta*log(exp(-beta*(phi_1)^2)+exp(-beta*Par^2)); %log softmin
    % Par=(Par+sqrt(Par^2+eps))/2;
    condition=(sqrt(Sat^2+Par^2+eps)-(Sat+Par))/1e4