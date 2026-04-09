% A Newton's method assembler and solver
% total variable:
% um uf  length:N+1
% phi T cs cl: length N
% S cs2 cl2: length N
% HH, Mar 2026

%estimate the total degree of entries
% nnz=(3+1)*(N+1)+2*N... %momentum 
%     +2*(N+1)+N     ... %continuity
%     +4*(N+1)+4*N   ... %component transport
%     +2*(N+1)+3*N   ... %enthalpy transport
%     +4*N*2;           ... %chemical model
if Has_volatile==0
    Dof=(N+1)*2+N*4;
    nnz=10*(N+1)+4*(N+1)+13*N+8*N+8*N*2;
else
    Dof=(N+1)*2+N*7;
    nnz=12*(N+1)+4*(N+1)+13*N+8*N*2+8*N*2+14*N*2;
end

rows=zeros(nnz,1);
cols=zeros(nnz,1);
vals=zeros(nnz,1);

RHS=zeros(Dof,1);

% OLD_com=phi_old(1:N).*Cl_old+(1-phi_old).*Cs_old;
% OLD_ent=phi_old(1:N)*Lf+T_old*cp; 
Max_Newton_iter=200;
converged=false;

% Initial guess
if Has_volatile==0
    % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
    X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
else
    % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old; S_old; Cs2_old; Cl2_old];
    X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N); S_old; Cs2_old; Cl2_old];
end
X0=X;
X_pre=X;
if Has_volatile==1
    %preprocessing pressure index data
    Pressure=(nodez(end)-cellz)*g*rho_mean/1e5+1; %in bar
    Pressure=Pressure';
    index_pressure_nts=max(min(floor(Pressure/8000*N_Ts)+1,N_Ts),1);    
    index_pressure_ntl=max(min(floor(Pressure/40e3* N_Tl) + 1, N_Tl),1);

    Cb2_old=phi_old(1:N).*Cl2_old+(1-phi_old(1:N)).*Cs2_old+S_old;

    PT=[-4e-4,-5e-4,-6e-4,-13e-4, -15.5e-4,-17e-4,-16e-4,-5e-4, 0, 26e-4,5e-3, 5e-3];  % coefficient from (ref) for the temperature dependency of water saturation (Holtz et al,?) 
    PTx=[0,   0.12    0.2  0.3    0.5       1       2       3   4  5,    11, 20]; %in kbar
    dSdT = interp1(PTx, PT, Pressure/1000, 'linear', 'extrap');    
    
    P3=Pressure/10; %Pressure in Mpa
    b0_all=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100-8*dSdT;
end

Norm_pre=1e3;
Not_improve=0;
for iter=1:Max_Newton_iter
    %% Momontum assembly
    entry_count=1;
    if Has_volatile==2
        num_local=12;
    else
        num_local=10;
    end
        % update indes
    index_phi_nmus=max(min(floor(X((1:N)+(N+1)*2) * N_mus) + 1, N_mus),1);

    index_phi_nc=max(min(floor(X((1:N)+(N+1)*2) * N_C) + 1, N_C),1);
    index_cl_nc=max(min(floor(X((1:N)+(N+1)*2+N*3) * N_C) + 1, N_C),1);

    index_cs_nrhos=max(min(floor(X((1:N)+(N+1)*2+N*2) * N_rhos) + 1, N_rhos),1);   

    if Has_volatile==1
        index_cl2_nc=max(min(floor(X((1:N)+(N+1)*2+N*6) * N_C) + 1, N_C),1);
        % cb2=phi*cl2+(1-phi)*cs2+S;
        cb2=X((1:N)+(N+1)*2).*X((1:N)+(N+1)*2+N*6)+(1-X((1:N)+(N+1)*2)).*X((1:N)+(N+1)*2+N*5)+X((1:N)+(N+1)*2+N*4);
        index_cb2_nts=max(min(floor(cb2/0.13/Par_v * N_Ts) + 1, N_Ts),1);
        index_cb2_ntl=max(min(floor(cb2/0.2/Par_v * N_Ts) + 1, N_Ts),1);
    end

    for i=2:N
        %Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
        if Has_volatile==2
            variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6];
        else
            %     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
            variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3];
        end
        
        Variables=X(variable_indes);        
        %Solid coefficients
        Parameters=zeros(28,1);
        % index_phi_0=max(min(floor(X(i-1+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_0 index
        Parameters(1:2)=mus_coef_all(index_phi_nmus(i-1),:)  ; %mus_0 

        % index_phi_1=max(min(floor(X(i+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_1 index
        Parameters(3:4)=mus_coef_all(index_phi_nmus(i),:)  ; %mus_1 
        
        % Coupling coefficients
        % index_cl_1=min(floor(X(i+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_1 index
        % index_phi_1=max(min(floor(X(i+(N+1)*2) * N_C) + 1, N_C),1); %phi_1 index
        if Has_volatile==2
            index_cl2_1=min(floor(X(i+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
            Parameters(5:12)=C_coef_all(index_phi_nc(i),index_cl_nc(i),index_cl2_1,:);
        else
            Parameters(5:8)=C_coef_all(index_phi_nc(i),index_cl_nc(i),:);
        end
        
        % Density coefficients
        % index_cs_0=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_0 index
        Parameters(13:14)=rhos_coef_all(index_cs_nrhos(i-1),:);        
        % index_cs_1=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_1 index
        Parameters(15:16)=rhos_coef_all(index_cs_nrhos(i),:);  

        if Has_volatile==2
            % index_cl2_0=min(floor(X(i-1+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
            Parameters(17:20)=rhol_coef_all(index_cl_nc(i-1), index_cl2_nc(i-1),:);
            Parameters(21:24)=rhol_coef_all(index_cl_nc(i), index_cl2_nc(i),:);
        else
            % index_cl_0=min(floor(X(i-1+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_0 index
            Parameters(17:18)=rhol_coef_all(index_cl_nc(i-1),:);
            Parameters(21:22)=rhol_coef_all(index_cl_nc(i),:);
        end
        
        % dz
        Parameters(25:28)=[dz(i-1) dz(i) g mum_0];

        vals(entry_count:entry_count+num_local-1)=Jac_mom([Variables; Parameters]);
        RHS(i)=rhs_mom([Variables; Parameters]);

        rows(entry_count:entry_count+num_local-1)=i;
        cols(entry_count:entry_count+num_local-1)=variable_indes;
        entry_count=entry_count+num_local;
    end
    
    % Momentum node N, and boundary
    %velocity boundary conditions
    rows(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
    cols(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
    vals(entry_count+(0:3))=1;

    RHS([1 N+1 N+2 2*N+2])=X([1 N+1 N+2 2*N+2]);
    entry_count=entry_count+4;

    %% Continuity assembly
    num_local=4;
    for i=1:N+1
        % umi_1,ufi,phi0,phi1
        % in=[um(i), uf(i), phi(i-1), phi(i)];
        in=X([i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2]);
        vals(entry_count:entry_count+num_local-1)=Jac_con(in);
        RHS(i+N+1)=rhs_con(in);

        rows(entry_count:entry_count+num_local-1)=i+N+1;
        cols(entry_count:entry_count+num_local-1)=[i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2];

        entry_count=entry_count+num_local;
    end

    %% Major component transport assembly
    num_local=13;
    for i=2:N-1
        %umi_1,umi_2,ufi,ufi2,phi0,phi1,phi2,cs0,cs,cs2,cl0,cl,cl2,dt,dzi_1,Cb_old
        % in=[um(i), um(i+1), uf(i), uf(i+1), phi(i-1), phi(i), phi(i+1), Cs(i-1),Cs(i)            ,Cs(i+1)      ,Cl(i-1),Cl(i),Cl(i+1), dt, dz(i), Cb_old(i)];
        in=[X([i     ,i+1   ,i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N*2 ,i+(N+1)*2+N*2 ,i+1+(N+1)*2+N*2,i-1+(N+1)*2+N*3 ,i+(N+1)*2+N*3 ,i+1+(N+1)*2+N*3])', dt, dz(i), Cb_old(i)];

        vals(entry_count:entry_count+num_local-1)=Jac_ct(in');
        RHS(i+(N+1)*2)=rhs_ct(in');

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2;
        %umi_1 ,umi_2 ,ufi      ,ufi2     ,phi0       ,phi1     ,phi2       ,cs0             ,cs            ,cs2            ,cl0             ,cl            ,cl2
        cols(entry_count:entry_count+num_local-1)=[i     ,i+1   ,i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N*2 ,i+(N+1)*2+N*2 ,i+1+(N+1)*2+N*2,i-1+(N+1)*2+N*3 ,i+(N+1)*2+N*3 ,i+1+(N+1)*2+N*3];

        entry_count=entry_count+num_local;
    end

    %fix the first and last cell value
    rows(entry_count+(0:1))=[1 N]+(N+1)*2;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2)-phi_old([1 N]);
    entry_count=entry_count+2;

    %% Enthalpy transport assembly
    num_local=8;

    for i=2:N-1
        % ufi,ufi2,phi0,phi1,phi2,T0,T1,T2,dzi_0,dzi_1,dzi_2,dt,cp,Lf,kt0,H_old
        % in=[uf(i), uf(i+1), phi(i-1),  phi(i), phi(i+1), T(i-1), T(i), T(i+1), dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];
        in=[X([i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N,i+(N+1)*2+N,i+1+(N+1)*2+N  ])', dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];

        vals(entry_count:entry_count+num_local-1)=Jac_ent(in');
        RHS(i+(N+1)*2+N)=rhs_ent(in');

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N;
        %,ufi     ,ufi2     ,phi0       ,phi1     ,phi2       ,T0           ,T1         ,T2
        cols(entry_count:entry_count+num_local-1)=[i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N,i+(N+1)*2+N,i+1+(N+1)*2+N  ];
        entry_count=entry_count+num_local;
    end

    %fix the first and last cell value
    rows(entry_count+(0:1))=[1 N]+(N+1)*2+N;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2+N;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2+N)=X([1 N]+(N+1)*2+N)-T_old([1 N]);
    entry_count=entry_count+2;   

    %% Solidus and liquidus
    if Has_volatile==1
        num_local=7;
    else
        num_local=4;
    end
    for i=2:N-1
        if Has_volatile==1
            %[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1; ]
            column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3, i+(N+1)*2+N*4, i+(N+1)*2+N*5, i+(N+1)*2+N*6];
            in=X(column_index);
            
            a=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),1);
            b=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),2);
            c=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),3);
            d=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),4);

            A=b+c*Pressure(i)/8000;
            B=a*Pressure(i)/8000+d;

            a=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),1);
            b=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),2);
            c=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),3);
            d=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),4);

            C=b+c*Pressure(i)/40e3;
            D=a*Pressure(i)/40e3+d;
            Parameters=[A;B;C;D;Par_v];            
        else
            %[phi(i), T(i), Cs(i),  Cl(i)]
            column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3];
            in=X(column_index);
            Parameters=[Ts0(1);Tl0(1);0;0; 0];
        end

        vals(entry_count:entry_count+num_local-1)=Jac_solidus([in; Parameters;n_order]);
        RHS(i+(N+1)*2+N*2)=rhs_solidus([in; Parameters;n_order]);

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*2;
        % phi T cs cl

        cols(entry_count:entry_count+num_local-1)=column_index;
        entry_count=entry_count+num_local;


        
        vals(entry_count:entry_count+num_local-1)=Jac_liquidus([in; Parameters; A1]);
        RHS(i+(N+1)*2+N*3)=rhs_liquidus([in; Parameters; A1]);

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*3;
        cols(entry_count:entry_count+num_local-1)=column_index;
        entry_count=entry_count+num_local;
        
    end    

    %fix the first and last cell value
    rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2+N*2)=X([1 N]+(N+1)*2+N*2)-C_all_old([1 N]+N);
    entry_count=entry_count+2;

    
    rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2+N*3)=X([1 N]+(N+1)*2+N*3)-C_all_old([1 N]);
    entry_count=entry_count+2;
    %% Volatile component transport assembly
    if Has_volatile ==1
        num_local=16;
        for i=2:N-1
            %umi_1; umi_2; ufi; ufi2; phi_0; phi_1; phi_2; S_0; S_1; S_2; cs2_0; cs2_1; cs2_2; cl2_0; cl2_1; cl2_2
            in=[X([i ,i+1, i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2, i+(N+1)*2,i+1+(N+1)*2,...
                i-1+(N+1)*2+N*4, i+(N+1)*2+N*4, i+1+(N+1)*2+N*4,...
                i-1+(N+1)*2+N*5, i+(N+1)*2+N*5, i+1+(N+1)*2+N*5, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6, i+1+(N+1)*2+N*6]); dz(i-1);dz(i); dz(i+1); dt;  kf; Cb2_old(i)];

            vals(entry_count:entry_count+num_local-1)=Jac_ct2(in);
            RHS(i+(N+1)*2+N*4)=rhs_ct2(in);

            rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*4;
            cols(entry_count:entry_count+num_local-1)=[i ,i+1, i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2, i+(N+1)*2,i+1+(N+1)*2,...
                i-1+(N+1)*2+N*4, i+(N+1)*2+N*4, i+1+(N+1)*2+N*4,...
                i-1+(N+1)*2+N*5, i+(N+1)*2+N*5, i+1+(N+1)*2+N*5, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6, i+1+(N+1)*2+N*6];

            entry_count=entry_count+num_local;
        end

        %fix the first and last cell value
        rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*4;
        cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*4;
        vals(entry_count+(0:1))=1;
        RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*4)-S_old([1 N]);
        entry_count=entry_count+2;        

        %% Solid saturation equation
        num_local=4;
        for i=2:N-1
            %phi_1, cs_1;  cs2_1; cl2_1;
            column_index=[i+(N+1)*2, i+(N+1)*2+N*2, i+(N+1)*2+N*5, i+(N+1)*2+N*6 ];
            in=X(column_index);


            vals(entry_count:entry_count+num_local-1)=Jac_ssat([in;S_cap(1); S_cap(2); Par_v]);
            RHS(i+(N+1)*2+N*5)=rhs_ssat([in;S_cap(1); S_cap(2); Par_v]);
            rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*5;
            cols(entry_count:entry_count+num_local-1)=column_index;
            entry_count=entry_count+num_local;
        end
        %fix the first and last cell value
        rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*5;
        cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*5;
        vals(entry_count+(0:1))=1;
        RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*5)-Cs2_old([1 N]);
        entry_count=entry_count+2;
        %% Melt saturation equation
        num_local=4;
        % Sat=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100+(T-800)*dSdT/100;

        % a0=dSdT/100;
        for i=2:N-1
            %phi_1, T_1; S_1; cl2_1
            column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*4 ,i+(N+1)*2+N*6 ];
            in=X(column_index);

            a0=dSdT(i)/100;
            b0=b0_all(i);
            vals(entry_count:entry_count+num_local-1)=Jac_lsat([in; a0; b0]);
            RHS(i+(N+1)*2+N*6)=rhs_lsat([in;a0;b0]);
            rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*6;
            cols(entry_count:entry_count+num_local-1)=column_index;
            entry_count=entry_count+num_local;
        end

        %fix the first and last cell value
        rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*6;
        cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*6;
        vals(entry_count+(0:1))=1;
        RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*6)-Cl2_old([1 N]);
        entry_count=entry_count+2;


    end
%%
    norm_R = max(abs(RHS));
    if norm_R < Precision
        fprintf('Converged in %d iterations\n', iter);
        converged = true;
        break;
    end

    Matrix_A = sparse(rows(1:entry_count-1), cols(1:entry_count-1), vals(1:entry_count-1), Dof, Dof);
    
    du = Matrix_A \ (-RHS);
    % if any(isnan(du))
    %     dt=dt*0.5;   
    %     X=X0;
    %     disp(['Decrease dt, dt=' num2str(dt/Year)]) 
    %     continue
    % end
    % Force constant values on boundary
    du([[1 N]+(N+1)*2 [1 N]+(N+1)*2+N [1 N]+(N+1)*2+N*2 [1 N]+(N+1)*2+N*3])=0;
    
    alpha = 1.0;
    for ls = 1:10
        X = X_pre + alpha * du;
        % force 1>phi>0
        X(2*(N+1)+1:2*(N+1)+N)=max(X(2*(N+1)+1:2*(N+1)+N),-0.99e-2);
        % X(2*(N+1)+1:2*(N+1)+N)=min(X(2*(N+1)+1:2*(N+1)+N),1);

        % force 1350>T
        % X(2*(N+1)+N+1:2*(N+1)+N*2)=min(X(2*(N+1)+N+1:2*(N+1)+N*2),1350);
 
        RHS=zeros(Dof,1);
                index_phi_nmus=max(min(floor(X((1:N)+(N+1)*2) * N_mus) + 1, N_mus),1);

        index_phi_nc=max(min(floor(X((1:N)+(N+1)*2) * N_C) + 1, N_C),1);
        index_cl_nc=max(min(floor(X((1:N)+(N+1)*2+N*3) * N_C) + 1, N_C),1);

        index_cs_nrhos=max(min(floor(X((1:N)+(N+1)*2+N*2) * N_rhos) + 1, N_rhos),1);
        if Has_volatile==1
            index_cl2_nc=max(min(floor(X((1:N)+(N+1)*2+N*6) * N_C) + 1, N_C),1);
            % cb2=phi*cl2+(1-phi)*cs2+S;
            cb2=X((1:N)+(N+1)*2).*X((1:N)+(N+1)*2+N*6)+(1-X((1:N)+(N+1)*2)).*X((1:N)+(N+1)*2+N*5)+X((1:N)+(N+1)*2+N*4);
            index_cb2_nts=max(min(floor(cb2/0.13/Par_v * N_Ts) + 1, N_Ts),1);
            index_cb2_ntl=max(min(floor(cb2/0.2/Par_v * N_Ts) + 1, N_Ts),1);
        end
        for i=2:N
            %Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
            if Has_volatile==2
                variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6];
            else
                %     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
                variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3];
            end

            Variables=X(variable_indes);
            %parameters=[a0; b0; a1; b1; a2; b2; c2; d2; e2; f2; g2; h2; a3; b3; a4; b4; a5; b5; c5; d5; a6; b6; c6; d6;  dzi_0; dzi_1; g; mum_0];
            %Solid coefficients
            Parameters=zeros(28,1);
            % index_phi_0=max(min(floor(X(i-1+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_0 index
            Parameters(1:2)=mus_coef_all(index_phi_nmus(i-1),:)  ; %mus_0

            % index_phi_1=max(min(floor(X(i+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_1 index
            Parameters(3:4)=mus_coef_all(index_phi_nmus(i),:)  ; %mus_1

            % Coupling coefficients
            % index_cl_1=min(floor(X(i+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_1 index
            % index_phi_1=max(min(floor(X(i+(N+1)*2) * N_C) + 1, N_C),1); %phi_1 index
            if Has_volatile==2
                index_cl2_1=min(floor(X(i+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
                Parameters(5:12)=C_coef_all(index_phi_nc(i),index_cl_nc(i),index_cl2_1,:);
            else
                Parameters(5:8)=C_coef_all(index_phi_nc(i),index_cl_nc(i),:);
            end

            % Density coefficients
            % index_cs_0=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_0 index
            Parameters(13:14)=rhos_coef_all(index_cs_nrhos(i-1),:);
            % index_cs_1=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_1 index
            Parameters(15:16)=rhos_coef_all(index_cs_nrhos(i),:);

            if Has_volatile==2
                % index_cl2_0=min(floor(X(i-1+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
                Parameters(17:20)=rhol_coef_all(index_cl_nc(i-1), index_cl2_nc(i-1),:);
                Parameters(21:24)=rhol_coef_all(index_cl_nc(i), index_cl2_nc(i),:);
            else
                % index_cl_0=min(floor(X(i-1+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_0 index
                Parameters(17:18)=rhol_coef_all(index_cl_nc(i-1),:);
                Parameters(21:22)=rhol_coef_all(index_cl_nc(i),:);
            end

            % dz
            Parameters(25:28)=[dz(i-1) dz(i) g mum_0];


            RHS(i)=rhs_mom([Variables; Parameters]);
        end

        for i=1:N+1
            % umi_1,ufi,phi0,phi1
            % in=[um(i), uf(i), phi(i-1), phi(i)];
            in=X([i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2]);
            RHS(i+N+1)=rhs_con(in);
        end

        for i=2:N-1
            %umi_1,umi_2,ufi,ufi2,phi0,phi1,phi2,cs0,cs,cs2,cl0,cl,cl2,dt,dzi_1,Cb_old
            % in=[um(i), um(i+1), uf(i), uf(i+1), phi(i-1), phi(i), phi(i+1), Cs(i-1),Cs(i)            ,Cs(i+1)      ,Cl(i-1),Cl(i),Cl(i+1), dt, dz(i), Cb_old(i)];
            in=[X([i     ,i+1   ,i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N*2 ,i+(N+1)*2+N*2 ,i+1+(N+1)*2+N*2,i-1+(N+1)*2+N*3 ,i+(N+1)*2+N*3 ,i+1+(N+1)*2+N*3])', dt, dz(i), Cb_old(i)];

            RHS(i+(N+1)*2)=rhs_ct(in');
        end

        for i=2:N-1
            % ufi,ufi2,phi0,phi1,phi2,T0,T1,T2,dzi_0,dzi_1,dzi_2,dt,cp,Lf,kt0,H_old
            % in=[uf(i), uf(i+1), phi(i-1),  phi(i), phi(i+1), T(i-1), T(i), T(i+1), dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];
            in=[X([i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N,i+(N+1)*2+N,i+1+(N+1)*2+N  ])', dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];

            RHS(i+(N+1)*2+N)=rhs_ent(in');
        end
        
        for i=2:N-1
            if Has_volatile==1
                %[phi_1; T_1; cs_1; cl_1; S_1; cs2_1; cl2_1; ]
                column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3, i+(N+1)*2+N*4, i+(N+1)*2+N*5, i+(N+1)*2+N*6];
                in=X(column_index);

                a=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),1);
                b=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),2);
                c=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),3);
                d=Ts0_coefficient(index_pressure_nts(i), index_cb2_nts(i),4);

                A=b+c*Pressure(i)/8000;
                B=a*Pressure(i)/8000+d;

                a=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),1);
                b=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),2);
                c=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),3);
                d=Tl0_coefficient(index_pressure_ntl(i), index_cb2_ntl(i),4);

                C=b+c*Pressure(i)/40e3;
                D=a*Pressure(i)/40e3+d;
                Parameters=[A;B;C;D;Par_v];
            else
                %[phi(i), T(i), Cs(i),  Cl(i)]
                column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3];
                in=X(column_index);
                Parameters=[Ts0(1);Tl0(1);0;0; 0];
            end
            RHS(i+(N+1)*2+N*2)=rhs_solidus([in; Parameters; n_order]);
            RHS(i+(N+1)*2+N*3)=rhs_liquidus([in; Parameters; A1]);            
        end

        if Has_volatile==1
            for i=2:N-1
                %umi_1; umi_2; ufi; ufi2; phi_0; phi_1; phi_2; S_0; S_1; S_2; cs2_0; cs2_1; cs2_2; cl2_0; cl2_1; cl2_2
                in=[X([i ,i+1, i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2, i+(N+1)*2,i+1+(N+1)*2,...
                    i-1+(N+1)*2+N*4, i+(N+1)*2+N*4, i+1+(N+1)*2+N*4,...
                    i-1+(N+1)*2+N*5, i+(N+1)*2+N*5, i+1+(N+1)*2+N*5, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6, i+1+(N+1)*2+N*6]); dz(i-1);dz(i); dz(i+1); dt;  kf; Cb2_old(i)];

                RHS(i+(N+1)*2+N*4)=rhs_ct2(in);
            end
            for i=2:N-1
                %phi_1, cs_1;  cs2_1; cl2_1;
                column_index=[i+(N+1)*2, i+(N+1)*2+N*2, i+(N+1)*2+N*5, i+(N+1)*2+N*6 ];
                in=X(column_index);
                RHS(i+(N+1)*2+N*5)=rhs_ssat([in;S_cap(1); S_cap(2); Par_v]);

            end
            for i=2:N-1
                %phi_1, T_1; S_1; cl2_1
                column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*4 ,i+(N+1)*2+N*6 ];
                in=X(column_index);
                RHS(i+(N+1)*2+N*6)=rhs_lsat([in;a0;b0]);
            end
        end

        temp_norm=max(abs(RHS)); 
        if temp_norm < (1-alpha*1e-4)*norm_R           
            break;            
        end
        alpha = alpha * 0.5;
    end
    if temp_norm < Precision
        fprintf('Converged in %d iterations\n', iter);
        converged = true;
        break;
    end

    % we expect the Newton method to converge with certain speed, otherwise
    % the time step is too big, we reset everything with smaller time step
    if temp_norm<Norm_pre*0.99 %|| dt/dt0<1e-4
        Norm_pre=temp_norm;
    else
        Not_improve=Not_improve+1;
        if Not_improve>5           
            if Has_volatile==0
                % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
                X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
            else
                % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old; S_old; Cs2_old; Cl2_old];
                X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N); S_old; Cs2_old; Cl2_old];
            end
            dt=dt*0.5;   
            disp(['Decrease dt, dt=' num2str(dt/Year)])  
            Not_improve=0;
            Norm_pre=1e3;
        end
    end
    X_pre=X;
end
if iter<=5
    if dt<0.01*Max_dt
        dt=min(dt*1.5, Max_dt);
    elseif dt<0.05*Max_dt
        dt=min(dt*1.2, Max_dt);
    else
        dt=min(dt*1.1, Max_dt);
    end
elseif iter>=10
    dt=min(dt*0.9, Max_dt);
end
% um=X(1:N+1);
% uf=X(N+2:2*N+2);
u_all=[X(N+2:2*N+2); X(1:N+1)];


phi=[X((1:N)+2*N+2); 1-X((1:N)+2*N+2)];
T=X((1:N)+2*N+2+N);
% Cs=X((1:N)+2*N+2+N*2);
% Cl=X((1:N)+2*N+2+N*3);
C_all=[X((1:N)+2*N+2+N*3); X((1:N)+2*N+2+N*2)];

% not used for calculation
Cb=X((1:N)+2*N+2).*X((1:N)+2*N+2+N*3)+(1-X((1:N)+2*N+2)).*X((1:N)+2*N+2+N*2);
H=phi(1:N)*Lf+cp*T;

if Has_volatile==1
    S=X((1:N)+2*N+2+N*4);
    Cs2=X((1:N)+2*N+2+N*5);
    Cl2=X((1:N)+2*N+2+N*6);
    Cb2=Cl2.*phi(1:N)+Cs2.*phi(N+1:end)+S;

    Lsaturation=(2.859e-2*P3-1.495e-3*P3.^1.5+2.702e-5*P3.^2+0.257*P3.^0.5)/100+(max(T,500)-800).*dSdT/100;
    Ssaturation=S_cap(1)*Cs2+S_cap(2)*(1-Cs2);

    Ts0=zeros(N,1);
    Tl0=zeros(N,1);
    index_v_nts=max(min(floor(Cb2/13*100/Par_v*N_Ts)+1,N_Ts),1);
    % index_pressure_nts=max(min(floor(Pressure/8000*N_Ts)+1,N_Ts),1);
    for i=1:N
        coef=Ts0_coefficient(index_pressure_nts(i),index_v_nts(i),:);
        Ts0(i)=coef(1)*Pressure(i)/8000+coef(2)*Cb2(i)/Par_v/13*100+coef(3)*Pressure(i)*Cb2(i)/Par_v/8000/13*100+coef(4);
    end

    index_v_ntl=max(min(floor(Cb2/20*100/Par_v*N_Ts)+1,N_Ts),1);
    % index_pressure_ntl=max(min(floor(Pressure/40e3*N_Ts)+1,N_Ts),1);
    for i=1:N
        coef=Tl0_coefficient(index_pressure_ntl(i),index_v_ntl(i),:);
        Tl0(i)=coef(1)*Pressure(i)/40e3+coef(2)*Cb2(i)/Par_v/20*100+coef(3)*Pressure(i)*Cb2(i)/Par_v/400/20+coef(4);
    end
    Ts=Tl0-Cb.^(1/n_order).*(Tl0-Ts0);
    Tl=A1*Cb.^2+(Ts0-Tl0-A1).*Cb+Tl0;
else
    % Ts=(1-max(Cb,1e-12).^(1/n_order))*(Tl0-Ts0)+Ts0; 
    Ts=Tl0(1)-Cb.^(1/n_order).*(Tl0(1)-Ts0(1)); %local solidus
    Tl=A1*Cb.^2+B1*Cb+C1; %local liquidus
end





improve1=norm_R;
improve2=improve1;


%% variables that needs to be updated for Newton
rhom=zeros(N+1,1);
rhof=zeros(N+1,1);
rho_b=zeros(N,1);

% if min(Cl2)<0
%     aaa=1;
% end
    
% Total_cb2=sum(Cb2.*cellz);
% disp(['CB2 conservation:' num2str(Total_cb2/Total_cb20)])

if iter==Max_Newton_iter
    error('Newton solver not converged before max iteration');
end

