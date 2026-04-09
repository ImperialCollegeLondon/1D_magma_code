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
    nnz=12*(N+1)+4*(N+1)+13*N+8*N+8*N*2+8*N*2;
end






rows=zeros(nnz,1);
cols=zeros(nnz,1);
vals=zeros(nnz,1);

RHS=zeros(Dof,1);


% OLD_com=phi_old(1:N).*Cl_old+(1-phi_old).*Cs_old;
% OLD_ent=phi_old(1:N)*Lf+T_old*cp; 



Max_Newton_iter=100;
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

Norm_pre=1e3;
Not_improve=0;
for iter=1:Max_Newton_iter
    %% Momontum assembly
    entry_count=1;
    if Has_volatile==1
        num_local=12;
    else
        num_local=10;
    end

    for i=2:N
        %Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
        if Has_volatile==1
            variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6];
        else
            %     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
            variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3];
        end
        
        Variables=X(variable_indes);        
        %parameters=[a0; b0; a1; b1; a2; b2; c2; d2; e2; f2; g2; h2; a3; b3; a4; b4; a5; b5; c5; d5; a6; b6; c6; d6;  dzi_0; dzi_1; g; mum_0];
        %Solid coefficients
        Parameters=zeros(28,1);
        index_phi_0=max(min(floor(X(i-1+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_0 index
        Parameters(1:2)=mus_coef_all(index_phi_0,:)  ; %mus_0 

        index_phi_1=max(min(floor(X(i+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_1 index
        Parameters(3:4)=mus_coef_all(index_phi_1,:)  ; %mus_1 
        
        % Coupling coefficients
        % index1=min(floor(X(i+(N+1)*2) * N_C) + 1, N_C);%phi_1 index, reuse
        index_cl_1=min(floor(X(i+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_1 index
        index_phi_1=max(min(floor(X(i+(N+1)*2) * N_C) + 1, N_C),1); %phi_1 index
        if Has_volatile==0
            Parameters(5:8)=C_coef_all(index_phi_1,index_cl_1,:);
        else
            index_cl2_1=min(floor(X(i+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
            Parameters(5:12)=C_coef_all(index_phi_1,index_cl_1,index_cl2_1,:);
        end
        
        % Density coefficients
        index_cs_0=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_0 index
        Parameters(13:14)=rhos_coef_all(index_cs_0,:);        
        index_cs_1=min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos);%cs_1 index
        Parameters(15:16)=rhos_coef_all(index_cs_1,:);  

        if Has_volatile==0
            index_cl_0=min(floor(X(i-1+(N+1)*2+N*3) * N_C) + 1, N_C);%cl_0 index
            Parameters(17:18)=rhol_coef_all(index_cl_0,:);
            Parameters(21:22)=rhol_coef_all(index_cl_1,:);
        else
            index_cl2_0=min(floor(X(i-1+(N+1)*2+N*6) * N_C) + 1, N_C);%cl2_1 index
            Parameters(17:20)=rhol_coef_all(index_cl_0, index_cl2_0,:); 
            Parameters(21:24)=rhol_coef_all(index_cl_1, index_cl2_1,:); 
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

    % %% Solidus  and liquidus assembly
    % num_local=4;
    % for i=2:N-1
    % 
    %     if Has_volatile==1 % phi1,T1,cs,cl, S, cs2, cl2
    %         Variables=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3, i+(N+1)*2+N*4, i+(N+1)*2+N*5, i+(N+1)*2+N*6])];
    %     else % phi1,T1,cs,cl
    %         Variables=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])];
    %     end
    % 
    %     %Parameters=[a0; b0; a1; b1; ...
    %     %a2; b2; c2; d2; e2; f2;g2; h2; i2; j2; k2; l2; m2; n2; o2; p2...
    %     %D1];
    %     Parameters=zeros(21,1);
    %     if Has_volatile==1
    %         % when volatile is present Tl and Ts become functions of bulk water content (thus cb2=cl2*phi+(1-phi)*cs2)
    %         % presure dependency is constant for each time step.
    %         % cb2=cl2_1*phi_1+cs2_1*(1-phi_1)+S_1;
    %         cb2=X(i+(N+1)*2+N*6)*X(i+(N+1)*2)+X(i+(N+1)*2+N*5)*(1-X(i+(N+1)*2))+X(i+(N+1)*2+N*4);
    %         index_cb2=min(floor(cb2* N_vol) + 1, N_vol);
    %         Parameters(1:2)=Tl_coef_all(:,index_cb2);
    %         Parameters(3:4)=Ts_coef_all(:,index_cb2);
    % 
    %         Tl=Parameters(1)*cb2/D1+arameters(2);
    %         Ts=Parameters(3)*cb2/D1+arameters(4);
    %     else
    %         Parameters(1)=Tl0;
    %         Parameters(3)=Ts0;
    %         Tl=Tl0;
    %         Ts=Ts0;
    %     end
    %     T_scaled=min(max((X(i+(N+1)*2+N)-Ts)/(Tl-Ts),-0.1),1.1);
    % 
    %     %solidus
    %     index_phi=min(floor(X(i+(N+1)*2)    * N_sl) + 1, N_sl);
    %     index_T  =min(floor( (T_scaled+0.1)/1.2        * N_sl_T) + 1, N_sl_T);
    %     index_cs =min(floor(X(i+(N+1)*2+N*2)* N_sl) + 1, N_sl);
    %     index_cl =min(floor(X(i+(N+1)*2+N*3)* N_sl) + 1, N_sl);
    %     Parameters(5:20)=solidus_coef_all(index_phi,index_T,index_cs,index_cl, :);
    %     Parameters(21)=K1;
    % 
    %     vals(entry_count:entry_count+num_local-1)=Jac_solidus([Variables; Parameters]);
    %     RHS(i+(N+1)*2+N*2)=rhs_solidus([Variables; Parameters]);
    % 
    %     rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*2;
    %     % phi T cs cl
    %     cols(entry_count:entry_count+num_local-1)=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3];
    %     entry_count=entry_count+num_local;
    % 
    %     %liquidus
    %     % reuse the previous Parameters
    %     % Parameters=[a0; b0; a1; b1; a2; b2; D1];
    %     Parameters(5:20)=liquidus_coef_all(index_phi,index_T,index_cs,index_cl, :);
    %     vals(entry_count:entry_count+num_local-1)=Jac_liquidus([Variables; Parameters]);
    %     RHS(i+(N+1)*2+N*3)=rhs_liquidus([Variables; Parameters]);
    % end
    % 
    % %fix the first and last cell value
    % rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    % cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    % vals(entry_count+(0:1))=1;
    % RHS([1 N]+(N+1)*2+N*2)=X([1 N]+(N+1)*2+N*2)-C_all_old([1 N]+N);
    % entry_count=entry_count+2;
    % 
    % 
    % rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    % cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    % vals(entry_count+(0:1))=1;
    % RHS([1 N]+(N+1)*2+N*3)=X([1 N]+(N+1)*2+N*3)-C_all_old([1 N]);
    % entry_count=entry_count+2;
    %% solidus
    num_local=4;
    for i=2:N-1
        % phi1,T1,cs,cl,n,Ts,Tl
        % in=[phi(i), T(i), Cs(i),  Cl(i) n_order, Ts0, Tl0];
        in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', n_order, Ts0(1), Tl0(1)];

        vals(entry_count:entry_count+num_local-1)=Jac_solidus(in');
        RHS(i+(N+1)*2+N*2)=rhs_solidus(in');

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*2;
        % phi T cs cl
        cols(entry_count:entry_count+num_local-1)=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3];
        entry_count=entry_count+num_local;
    end
    
    %fix the first and last cell value
    rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*2;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2+N*2)=X([1 N]+(N+1)*2+N*2)-C_all_old([1 N]+N);
    entry_count=entry_count+2;
    %% liquidus
    num_local=4;
    for i=2:N-1
        % phi1,T1,cs,cl,a1,b1,c1
        % in=[phi(i), T(i),  Cs(i), Cl(i), a1, b1, c1];
        in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', A1, B1, C1];

        vals(entry_count:entry_count+num_local-1)=Jac_liquidus(in');
        RHS(i+(N+1)*2+N*3)=rhs_liquidus(in');

        rows(entry_count:entry_count+num_local-1)=i+(N+1)*2+N*3;
        % phi T cs cl
        cols(entry_count:entry_count+num_local-1)=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3];
        entry_count=entry_count+num_local;
    end
    
    %fix the first and last cell value
    rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*3;
    vals(entry_count+(0:1))=1;
    RHS([1 N]+(N+1)*2+N*3)=X([1 N]+(N+1)*2+N*3)-C_all_old([1 N]);
    entry_count=entry_count+2;
%%
    norm_R = norm(RHS);
    if norm_R < Precision
        fprintf('Converged in %d iterations\n', iter);
        converged = true;
        break;
    end

    Matrix_A = sparse(rows(1:entry_count-1), cols(1:entry_count-1), vals(1:entry_count-1), Dof, Dof);
    
    % if ~isreal(RHS)
    %     aaa=1; 
    % end
    du = Matrix_A \ (-RHS);
    % Force constant values on boundary
    du([[1 N]+(N+1)*2 [1 N]+(N+1)*2+N [1 N]+(N+1)*2+N*2 [1 N]+(N+1)*2+N*3])=0;
    
    alpha = 1.0;

    for ls = 1:20
        X = X0 + alpha * du;
        % force 1>phi>0
        X(2*(N+1)+1:2*(N+1)+N)=max(X(2*(N+1)+1:2*(N+1)+N),-0.99e-2);
        % X(2*(N+1)+1:2*(N+1)+N)=min(X(2*(N+1)+1:2*(N+1)+N),1);

        % force 1350>T
        % X(2*(N+1)+N+1:2*(N+1)+N*2)=min(X(2*(N+1)+N+1:2*(N+1)+N*2),1350);
 
        RHS=zeros(Dof,1);
        for i=2:N
            %Variables=[umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1];
            if Has_volatile==1
                variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3, i-1+(N+1)*2+N*6, i+(N+1)*2+N*6];
            else
                %     Variables=[umi_0; umi_1;  umi_2; ufi; ufi2;phi_0; phi_1; phi_2; T_0; T_1; T_2; cs_0; cs_1; cs_2; cl_0; cl_1; cl_2];
                variable_indes=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3];
            end

            Variables=X(variable_indes);
            %parameters=[a0; b0; a1; b1; a2; b2; c2; d2; e2; f2; g2; h2; a3; b3; a4; b4; a5; b5; c5; d5; a6; b6; c6; d6;  dzi_0; dzi_1; g; mum_0];
            %Solid coefficients
            Parameters=zeros(28,1);
            index_phi_0=max(min(floor(X(i-1+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_0 index
            Parameters(1:2)=mus_coef_all(index_phi_0,:)  ; %mus_0

            index_phi_1=max(min(floor(X(i+(N+1)*2) * N_mus) + 1, N_mus),1); %phi_1 index
            Parameters(3:4)=mus_coef_all(index_phi_1,:)  ; %mus_1

            % Coupling coefficients
            % index1=min(floor(X(i+(N+1)*2) * N_C) + 1, N_C);%phi_1 index, reuse
            index_cl_1=max(min(floor(X(i+(N+1)*2+N*3) * N_C) + 1, N_C),1);%cl_1 index
            index_phi_1=max(min(floor(X(i+(N+1)*2) * N_C) + 1, N_C),1); %phi_1 index
            if Has_volatile==0
                Parameters(5:8)=C_coef_all(index_phi_1,index_cl_1,:);
            else
                index_cl2_1=max(min(floor(X(i+(N+1)*2+N*6) * N_C) + 1, N_C),1);%cl2_1 index
                Parameters(5:12)=C_coef_all(index_phi_1,index_cl_1,index_cl2_1,:);
            end

            % Density coefficients
            index_cs_0=max(min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos),1);%cs_0 index
            Parameters(13:14)=rhos_coef_all(index_cs_0,:);
            index_cs_1=max(min(floor(X(i-1+(N+1)*2+N*2) * N_rhos) + 1, N_rhos),1);%cs_1 index
            Parameters(15:16)=rhos_coef_all(index_cs_1,:);

            if Has_volatile==0
                index_cl_0=max(min(floor(X(i-1+(N+1)*2+N*3) * N_rhol) + 1, N_rhol),1);%cl_0 index
                Parameters(17:18)=rhol_coef_all(index_cl_0,:);
                Parameters(21:22)=rhol_coef_all(index_cl_1,:);
            else
                index_cl2_0=max(min(floor(X(i-1+(N+1)*2+N*6) * N_rhol) + 1, N_rhol),1);%cl2_1 index
                Parameters(17:20)=rhol_coef_all(index_cl_0, index_cl2_0,:);
                Parameters(21:24)=rhol_coef_all(index_cl_1, index_cl2_1,:);
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
            % phi1,T1,cs,cl,n,Ts,Tl
            % in=[phi(i), T(i), Cs(i),  Cl(i) n_order, Ts0, Tl0];
            in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', n_order, Ts0(1), Tl0(1)];
            RHS(i+(N+1)*2+N*2)=rhs_solidus(in');
        end
        for i=2:N-1
            % phi1,T1,cs,cl,a1,b1,c1
            % in=[phi(i), T(i),  Cs(i), Cl(i), a1, b1, c1];
            in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', A1, B1, C1];
            RHS(i+(N+1)*2+N*3)=rhs_liquidus(in');
        end

        temp_norm=norm(RHS); 
        if temp_norm < norm_R
           
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
            % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
            X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
            dt=dt*0.5;   
            disp(['Decrease dt, dt=' num2str(dt/Year)])  
            Not_improve=0;
%             if dt/Year<1e-5
%                 aaa=1;
%             end
        end
    end
    % um=X(1:N+1);
    % uf=X(N+2:2*N+2);
    % phi=X((1:N)+2*N+2);
    % T=X((1:N)+2*N+2+N);
    % Cs=X((1:N)+2*N+2+N*2);
    % Cl=X((1:N)+2*N+2+N*3);
    X0=X;
end
if iter<6
    dt=min(dt*1.2, Max_dt);
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

Ts=(1-max(Cb,1e-12).^(1/n_order))*(Tl0(1)-Ts0(1))+Ts0(1); %local solidus
Tl=A1*Cb.^2+B1*Cb+C1; %local liquidus


improve1=Norm_pre;
improve2=improve1;


%% variables that needs to be updated for Newton
rhom=zeros(N+1,1);
rhof=zeros(N+1,1);
rho_b=zeros(N,1);

