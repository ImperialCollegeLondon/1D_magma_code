
i=2153;
            a0=dSdT(i)/100;
            b0=b0_all(i);
D1=0.2;
cap_A=S_cap(1);
cap_B=S_cap(2);
eps=1e-12;

            column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2,  i+(N+1)*2+N*4,i+(N+1)*2+N*5, i+(N+1)*2+N*6];
            in=X(column_index);

            phi_1=in(1); 
            T_1=in(2); 
            cs_1=in(3);
            S_1=in(4); 
            cs2_1=in(5); 
            cl2_1=in(6);
        
    T_capped=(T_1+500+sqrt((T_1-500)^2+eps))/2;
    Satl=a0*T_capped+b0-cl2_1;
    
    Sats=cap_A*cs_1+cap_B*(1-cs_1)-cs2_1+0.02;
    Par=cl2_1*D1-cs2_1;

    % condition=(Satl+S_1-(sqrt(Satl-S_1)^2+eps))/1e4;
    condition=(Satl+S_1-sqrt((Satl-S_1)^2+eps))/1e4



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


    %%

 a=linspace(-1,1,1e3);   
b=(1-tanh(a/3e-2))/2;
figure(9)
clf;
plot(a,b)


%%
i=49;
eps=3e-2;
        lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts(I), index_cb2_nts(I));
        lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl(I), index_cb2_ntl(I));

        coef_ts = Ts0_coefficient(lin_ts, :);   % (nI × 4)
        coef_tl = Tl0_coefficient(lin_tl, :);   % (nI × 4)

        % unpack
        a_ts = coef_ts(:,1); b_ts = coef_ts(:,2);
        c_ts = coef_ts(:,3); d_ts = coef_ts(:,4);

        a_tl = coef_tl(:,1); b_tl = coef_tl(:,2);
        c_tl = coef_tl(:,3); d_tl = coef_tl(:,4);

                P = Pressure(I);

        A = b_ts + c_ts .* (P/8000);
        B = a_ts .* (P/8000) + d_ts;

        C = b_tl + c_tl .* (P/40e3);
        D = a_tl .* (P/40e3) + d_tl;

    cb2=cl2_1(i)*phi_1(i)+cs2_1(i)*(1-phi_1(i))+S_1(i);
    Ts=A(i)*cb2/Par_v/0.13+B(i);
    Tl=C(i)*cb2/Par_v/0.2+D(i);


C1= Tl;
B1= Ts-A1-C1;
% Parameters=[a0;b0;a1;b1; D1];
Cb=phi_1(i)*cl_1(i)+(1-phi_1(i))*cs_1(i);


    eps2=5e-2;
    solidus=(Cb/2*(1-tanh((T_1(i)-Ts)/eps2))-cs_1(i))/1e4