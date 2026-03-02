% A Newton's method assembler
% total variable:
% um uf  length:N+1
% phi T cs cl length N
% HH, Mar 2026
Dof=(N+1)*2+N*4;

%estimate the total degree of entries
% nnz=(3+1)*(N+1)+2*N... %momentum 
%     +2*(N+1)+N     ... %continuity
%     +4*(N+1)+4*N   ... %component transport
%     +2*(N+1)+3*N   ... %enthalpy transport
%     +4*N*2;           ... %chemical model


nnz=10*(N+1)+4*(N+1)+13*N+8*N+4*N*2;

rows=zeros(nnz,1);
cols=zeros(nnz,1);
vals=zeros(nnz,1);

RHS=zeros(Dof,1);


% OLD_com=phi_old(1:N).*Cl_old+(1-phi_old).*Cs_old;
% OLD_ent=phi_old(1:N)*Lf+T_old*cp; 



Max_Newton_iter=100;
converged=false;

% Initial guess
% X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
X0=X;

Norm_pre=1e3;
Not_improve=0;
for iter=1:Max_Newton_iter
    %% momontum assembly
    entry_count=1;

    num_local=10;
    for i=2:N
        %umi_0,umi_1,umi_2,ufi,phi0,phi1,cs0,cs,cl0,cl,muf1,muf2,mum_0,c0,dzi_0,dzi_1,rhom_1,rhom_2,rhof_1,rhof_2,g
        % if Non_dimention==0
            in=[X([i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3])' mu_f1 mu_f2 mum_0 c0 dz(i-1) dz(i) rhom_1 rhom_2 rhof_1 rhof_2 g];
        % else
        %     in=[X([i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3])' mu_f1 mu_f2 1     1  dz(i-1) dz(i) rhom_1 rhom_2 rhof_1 rhof_2 g];
        % end

        vals(entry_count:entry_count+num_local-1)=Jac_mom(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15),in(16),in(17),in(18),in(19),in(20),in(21));
        RHS(i)=rhs_mom(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15),in(16),in(17),in(18),in(19),in(20),in(21));

        rows(entry_count:entry_count+num_local-1)=i;
        cols(entry_count:entry_count+num_local-1)=[i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3];
        entry_count=entry_count+num_local;
    end

    %velocity boundary conditions
    rows(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
    cols(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
    vals(entry_count+(0:3))=1;

    RHS([1 N+1 N+2 2*N+2])=X([1 N+1 N+2 2*N+2]);
    entry_count=entry_count+4;

    %% continuity assembly
    num_local=4;
    for i=1:N+1
        % umi_1,ufi,phi0,phi1
        % in=[um(i), uf(i), phi(i-1), phi(i)];
        in=X([i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2]);
        vals(entry_count:entry_count+num_local-1)=Jac_con(in(1),in(2),in(3),in(4));
        RHS(i+N+1)=rhs_con(in(1),in(2),in(3),in(4));

        rows(entry_count:entry_count+num_local-1)=i+N+1;
        cols(entry_count:entry_count+num_local-1)=[i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2];

        entry_count=entry_count+num_local;
    end

    %% component transport assembly
    num_local=13;
    for i=2:N-1
        %umi_1,umi_2,ufi,ufi2,phi0,phi1,phi2,cs0,cs,cs2,cl0,cl,cl2,dt,dzi_1,Cb_old
        % in=[um(i), um(i+1), uf(i), uf(i+1), phi(i-1), phi(i), phi(i+1), Cs(i-1),Cs(i)            ,Cs(i+1)      ,Cl(i-1),Cl(i),Cl(i+1), dt, dz(i), Cb_old(i)];
        in=[X([i     ,i+1   ,i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N*2 ,i+(N+1)*2+N*2 ,i+1+(N+1)*2+N*2,i-1+(N+1)*2+N*3 ,i+(N+1)*2+N*3 ,i+1+(N+1)*2+N*3])', dt, dz(i), Cb_old(i)];

        vals(entry_count:entry_count+num_local-1)=Jac_ct(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));
        RHS(i+(N+1)*2)=rhs_ct(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));

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

    %% enthalpy transport assembly
    num_local=8;

    for i=2:N-1
        % ufi,ufi2,phi0,phi1,phi2,T0,T1,T2,dzi_0,dzi_1,dzi_2,dt,cp,Lf,kt0,H_old
        % in=[uf(i), uf(i+1), phi(i-1),  phi(i), phi(i+1), T(i-1), T(i), T(i+1), dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];
        in=[X([i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N,i+(N+1)*2+N,i+1+(N+1)*2+N  ])', dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];

        vals(entry_count:entry_count+num_local-1)=Jac_ent(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));
        RHS(i+(N+1)*2+N)=rhs_ent(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));

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

    %% solidus
    num_local=4;
    for i=2:N-1
        % phi1,T1,cs,cl,n,Ts,Tl
        % in=[phi(i), T(i), Cs(i),  Cl(i) n_order, Ts0, Tl0];
        in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', n_order, Ts0, Tl0];

        vals(entry_count:entry_count+num_local-1)=Jac_solidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));
        RHS(i+(N+1)*2+N*2)=rhs_solidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));

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

        vals(entry_count:entry_count+num_local-1)=Jac_liquidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));
        RHS(i+(N+1)*2+N*3)=rhs_liquidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));

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
            in=[X([i-1, i, i+1, i+(N+1), i-1+(N+1)*2, i+(N+1)*2, i-1+(N+1)*2+N*2, i+(N+1)*2+N*2, i-1+(N+1)*2+N*3, i+(N+1)*2+N*3])' mu_f1 mu_f2 mum_0 c0 dz(i-1) dz(i) rhom_1 rhom_2 rhof_1 rhof_2 g];
            RHS(i)=rhs_mom(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15),in(16),in(17),in(18),in(19),in(20),in(21));
        end

        for i=1:N+1
            % umi_1,ufi,phi0,phi1
            % in=[um(i), uf(i), phi(i-1), phi(i)];
            in=X([i, i+(N+1),  i-1+(N+1)*2, i+(N+1)*2]);
            RHS(i+N+1)=rhs_con(in(1),in(2),in(3),in(4));
        end

        for i=2:N-1
            %umi_1,umi_2,ufi,ufi2,phi0,phi1,phi2,cs0,cs,cs2,cl0,cl,cl2,dt,dzi_1,Cb_old
            % in=[um(i), um(i+1), uf(i), uf(i+1), phi(i-1), phi(i), phi(i+1), Cs(i-1),Cs(i)            ,Cs(i+1)      ,Cl(i-1),Cl(i),Cl(i+1), dt, dz(i), Cb_old(i)];
            in=[X([i     ,i+1   ,i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N*2 ,i+(N+1)*2+N*2 ,i+1+(N+1)*2+N*2,i-1+(N+1)*2+N*3 ,i+(N+1)*2+N*3 ,i+1+(N+1)*2+N*3])', dt, dz(i), Cb_old(i)];

            RHS(i+(N+1)*2)=rhs_ct(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));
        end

        for i=2:N-1
            % ufi,ufi2,phi0,phi1,phi2,T0,T1,T2,dzi_0,dzi_1,dzi_2,dt,cp,Lf,kt0,H_old
            % in=[uf(i), uf(i+1), phi(i-1),  phi(i), phi(i+1), T(i-1), T(i), T(i+1), dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];
            in=[X([i+(N+1)  ,i+1+(N+1),i-1+(N+1)*2,i+(N+1)*2,i+1+(N+1)*2,i-1+(N+1)*2+N,i+(N+1)*2+N,i+1+(N+1)*2+N  ])', dz(i-1), dz(i), dz(i+1), dt, cp, Lf, kt0, H_old(i)];

            RHS(i+(N+1)*2+N)=rhs_ent(in(1),in(2),in(3),in(4),in(5),in(6),in(7),in(8),in(9),in(10),in(11),in(12),in(13),in(14),in(15), in(16));
        end

        for i=1:N
            % phi1,T1,cs,cl,n,Ts0,Tl0
            % in=[phi(i), T(i), Cs(i),  Cl(i) n_order, Ts, Tl];
            in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', n_order, Ts0, Tl0];

            RHS(i+(N+1)*2+N*2)=rhs_solidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));
        end

        for i=1:N
            % phi1,T1,cs,cl,a1,b1,c1
            % in=[phi(i), T(i),  Cs(i), Cl(i), a1, b1, c1];
            in=[X([i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*3])', A1, B1, C1];

            RHS(i+(N+1)*2+N*3)=rhs_liquidus(in(1),in(2),in(3),in(4),in(5),in(6),in(7));
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
            dt=dt*0.8;   
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

Ts=(1-max(Cb,1e-12).^(1/n_order))*(Tl0-Ts0)+Ts0; %local solidus
Tl=A1*Cb.^2+B1*Cb+C1; %local liquidus


improve1=temp_norm;
improve2=improve1;


%% variables that needs to be updated for Newton
rhom=zeros(N+1,1);
rhof=zeros(N+1,1);
rho_b=zeros(N,1);

