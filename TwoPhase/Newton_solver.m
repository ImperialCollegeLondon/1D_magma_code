% A Newton's method assembler and solver
% total variable:
% um uf  length:N+1
% phi T cs cl: length N
% S cs2 cl2: length N
% HH, Mar 2026






%%
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
Max_Newton_iter=300;
converged=false;

% Initial guess
if Has_volatile==0
    % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
    X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
else
    % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old; S_old; Cs2_old; Cl2_old];
    X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N); S_old; Cs2_old; Cl2_old];
end

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

% to_diffuse=ones(N,1)/1e2;



%% 

Kfreeze = 3;  
bad_count = zeros(N+1,1);
new_frozen=[];

Time_old=Time;
if Record_data2==1
    Time_gap=Record_period;
else
    Time_gap=Max_dtY*Year;
end
if Advance_time==1
    Time_intended=Time_intended+Time_gap;
end
while Time<Time_intended-Time_gap*1e-5
    % if dt_intended<0.5e-5*Year
    %     dt_intended=1*Year;
    % end
    if Advance_time==1
        dt=min(dt_intended, Time_intended-Time);
        dt=min(dt,Max_dtY*Year);
    else
        dt=0;
    end
    
    if Has_volatile==1
        System_top=find(X((1:N) + 2*(N+1))>1e-3,1,'last');    
        if ~isempty(System_top)
            temp=find(X((1:N) + (N+1)*2+N*4)>1e-3,1,'last');
            if ~isempty(temp)
                temp=max(temp,System_top);
            else
                temp=System_top;
            end
            System_top=find(cellz>cellz(temp)+max_dx/2,1,'first');
            System_bottom=find(X((1:N) + 2*(N+1))>1e-2,1,'first');
        else
            System_bottom=[];
        end
    end

    Norm_pre=1e3;
    Not_improve=0;
    X0=X;
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
        index_cl_nrhol=max(min(floor(X((1:N)+(N+1)*2+N*3) * N_rhol) + 1, N_rhol),1);
        if Has_volatile==1
            index_cl2_nc=max(min(floor(X((1:N)+(N+1)*2+N*6) * N_C) + 1, N_C),1);
            cb2=X((1:N)+(N+1)*2).*X((1:N)+(N+1)*2+N*6)+(1-X((1:N)+(N+1)*2)).*X((1:N)+(N+1)*2+N*5)+X((1:N)+(N+1)*2+N*4);
            index_cb2_nts=max(min(floor(cb2/0.13/Par_v * N_Ts) + 1, N_Ts),1);
            index_cb2_ntl=max(min(floor(cb2/0.2/Par_v * N_Tl) + 1, N_Tl),1);

            % index_cl2_nrhol=max(min(floor(X((1:N)+(N+1)*2+N*6)/0.1 * N_rhol) + 1, N_rhol),1);
        end
        I = (2:N)';
        nI = numel(I);

        i0 = I - 1;
        i1 = I;
        i2 = I + 1;
        um_0 = X(i0);
        um_1 = X(i1);
        um_2 = X(i2);
        uf_0 = X(i0 + (N+1));
        uf_1 = X(i1 + (N+1));
        phi_0 = X(i0 + 2*(N+1));
        phi_1 = X(i1 + 2*(N+1));
        cs_0 = X(i0 + 2*(N+1) + 2*N);
        cs_1 = X(i1 + 2*(N+1) + 2*N);

        cl_0 = X(i0 + 2*(N+1) + 3*N);
        cl_1 = X(i1 + 2*(N+1) + 3*N);

        if Has_volatile == 2
            cl2_0 = X(i0 + 2*(N+1) + 6*N);
            cl2_1 = X(i1 + 2*(N+1) + 6*N);
            cl2_2 = X(i2 + 2*(N+1) + 6*N);
        end

        if Has_volatile == 2
            %umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1
            Vars = [ ...
                um_0, um_1, um_2, ...
                uf_1, ...
                phi_0, phi_1, ...
                cs_0, cs_1,  ...
                cl_0, cl_1,  ...
                cl2_0, cl2_1,...
                ];
        else
            %umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1
            Vars = [ ...
                um_0, um_1, um_2, ...
                uf_1, ...
                phi_0, phi_1, ...
                cs_0, cs_1,  ...
                cl_0, cl_1 ...
                ];
        end

        Param = zeros(nI, 28);
        Param(:,1:2) = mus_coef_all(index_phi_nmus(i0),:);
        Param(:,3:4) = mus_coef_all(index_phi_nmus(i1),:);
        % coupling coefficients
        if Has_volatile == 2
            idx_cl2 = max(min(floor(X(i1 + (N+1)*2 + 6*N) * N_C) + 1, N_C),1);
            Param(:,5:12) = C_coef_all(index_phi_nc(i1), index_cl_nc(i1), idx_cl2, :);
        else
            lin = sub2ind([N_C, N_C], index_phi_nc(i1), index_cl_nc(i1));
            Param(:,5:8) = C_coef_all(lin, :);
            % Param(:,5:8) = C_coef_all(index_phi_nc(i1), index_cl_nc(i1), :);
        end
        % density solid
        Param(:,13:14) = rhos_coef_all(index_cs_nrhos(i0),:);
        Param(:,15:16) = rhos_coef_all(index_cs_nrhos(i1),:);

        % density liquid
        if Has_volatile == 2
            lin = sub2ind([N_rhol, N_rhol], index_cl_nrhol(i0), index_cl2_nrhol(i0));
            Param(:,17:20) = rhol_coef_all(lin, :);
            lin = sub2ind([N_rhol, N_rhol], index_cl_nrhol(i1), index_cl2_nrhol(i1));
            Param(:,21:24) = rhol_coef_all(lin, :);
        else
            Param(:,17:18) = rhol_coef_all(index_cl_nrhol(i0), :);
            Param(:,21:22) = rhol_coef_all(index_cl_nrhol(i1), :);
        end
       
        % mum_0_local=sqrt((mus_coef_all(index_phi_nmus(i0),1).*phi_0+mus_coef_all(index_phi_nmus(i0),2)+mus_coef_all(index_phi_nmus(i1),1).*phi_1+mus_coef_all(index_phi_nmus(i1),2)))./(dz(2:end)+dz(1+end-1))'.^2*2;
        mum_0_local=min(mus_coef_all(index_phi_nmus(i0),1).*phi_0+mus_coef_all(index_phi_nmus(i0),2))*ones(nI,1);
        % physical constants
        Param(:,25:28) = [dz(i0)', dz(i1)', g*ones(nI,1), mum_0_local]; %mum_0*ones(nI,1)
        RHS(I)=rhs_mom(Vars(:,1),Vars(:,2),Vars(:,3),Vars(:,4),Vars(:,5),Vars(:,6),Vars(:,7),Vars(:,8),Vars(:,9),Vars(:,10),...
            Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5),Param(:,6),Param(:,7),Param(:,8),Param(:,9),Param(:,10),...
            Param(:,11),Param(:,12),Param(:,13),Param(:,14),Param(:,15),Param(:,16),Param(:,17),Param(:,18),Param(:,19),Param(:,20),...
            Param(:,21),Param(:,22),Param(:,23),Param(:,24),Param(:,25),Param(:,26),Param(:,27),Param(:,28));

        J = Jac_mom(Vars(:,1),Vars(:,2),Vars(:,3),Vars(:,4),Vars(:,5),Vars(:,6),Vars(:,7),Vars(:,8),Vars(:,9),Vars(:,10),...
            Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5),Param(:,6),Param(:,7),Param(:,8),Param(:,9),Param(:,10),...
            Param(:,11),Param(:,12),Param(:,13),Param(:,14),Param(:,15),Param(:,16),Param(:,17),Param(:,18),Param(:,19),Param(:,20),...
            Param(:,21),Param(:,22),Param(:,23),Param(:,24),Param(:,25),Param(:,26),Param(:,27),Param(:,28));
        rows_block = repmat(I, 1, num_local);
        cols_block = [ ...
            i0, i1, i2, ...
            i1+(N+1), ...
            i0+2*(N+1), i1+2*(N+1),  ...
            i0+2*(N+1)+2*N, i1+2*(N+1)+2*N, ...
            i0+2*(N+1)+3*N, i1+2*(N+1)+3*N...
            ];
        idx = entry_count : entry_count + nI*num_local - 1;

        rows(idx) = rows_block(:);
        cols(idx) = cols_block(:);
        vals(idx) = J(:);

        entry_count = entry_count + nI*num_local;


        % Momentum node N, and boundary
        %velocity boundary conditions
        rows(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
        cols(entry_count+(0:3))=[1 N+1 N+2 2*N+2];
        vals(entry_count+(0:3))=1;

        RHS([1 N+1 N+2 2*N+2])=X([1 N+1 N+2 2*N+2]);
        entry_count=entry_count+4;

        %% Continuity assembly
        num_local=4;
        % umi_1,ufi,phi0,phi1
        

        r=(X((2:N-1) + 2*(N+1))-X((1:N-2) + 2*(N+1)))./(X((3:N) + 2*(N+1))-X((2:N-1) + 2*(N+1))+1e-15);
        flux_limiter=(r+abs(r))./(1+abs(r))/2*0+1;

        RHS(I+N+1)=rhs_con(um_1,uf_1,phi_0, phi_1, con_scaling,[flux_limiter;1]);
        J = Jac_con(um_1, uf_1, phi_0, phi_1, con_scaling, [flux_limiter;1]);

        
        rows_block = repmat(I + N + 1, 1, num_local);   % (N-1) × 4
        cols_block = [ ...
            I, ...
            I + (N+1), ...
            I-1 + (N+1)*2, ...
            I   + (N+1)*2 ...
            ];   % (N-1) × 4

        idx = entry_count : entry_count + num_local*(N-1) - 1;

        rows(idx) = rows_block(:);
        cols(idx) = cols_block(:);
        vals(idx) = J(:);

        entry_count = entry_count + num_local*(N-1);

        %% Major component transport assembly
        num_local=13;

        I = (2:N-1)';
        i0=I-1;
        i2=I+1;
        nI = length(I);
        um_1=X(I);
        um_2=X(I+1);
        uf_1=X(I+(N+1));
        uf_2=X(I+1+(N+1));
        phi_0=X(I-1+(N+1)*2);
        phi_1=X(I  +(N+1)*2);
        phi_2=X(I+1+(N+1)*2);
        cs_0=X(I-1+(N+1)*2+N*2);
        cs_1=X(I  +(N+1)*2+N*2);
        cs_2=X(I+1+(N+1)*2+N*2);
        cl_0=X(I-1+(N+1)*2+N*3);
        cl_1=X(I  +(N+1)*2+N*3);
        cl_2=X(I+1+(N+1)*2+N*3);

        % cb=phi_1.*cl_1+(1-phi_1).*cs_1;
        % enhance=(cb>0.995)*1e2;
        scale=ones(N,1);
        scale(Cb_old<2e-3)=1e3;
        k_stable1=k_stable_value*(phi_1>1e-2).*scale(I);
        k_stable2=[k_stable1(2:end);0];
        
        RHS(I+(N+1)*2)=rhs_ct(um_1,um_2, uf_1, uf_2, phi_0, phi_1, phi_2, cs_0,cs_1,cs_2, cl_0,cl_1,cl_2,dt, dz(i0)', dz(I)', dz(i2)',Cb_old(I), k_stable1, k_stable2, [1; flux_limiter(1:end-1)],flux_limiter);
        J = Jac_ct(um_1,um_2, uf_1, uf_2, phi_0, phi_1, phi_2, cs_0,cs_1,cs_2, cl_0,cl_1,cl_2, dt, dz(i0)', dz(I)', dz(i2)', Cb_old(I), k_stable1, k_stable2, [1; flux_limiter(1:end-1)],flux_limiter);


        rows_block = repmat(I+(N+1)*2, 1, num_local);
        cols_block = [ ...
            I, I+1, ...
            I+(N+1), I+1+(N+1), ...
            I-1+(N+1)*2, I+(N+1)*2, I+1+(N+1)*2, ...
            I-1+(N+1)*2+N*2, I+(N+1)*2+N*2, I+1+(N+1)*2+N*2, ...
            I-1+(N+1)*2+N*3, I+(N+1)*2+N*3, I+1+(N+1)*2+N*3 ...
            ];

        idx = entry_count : entry_count + num_local*(N-2) - 1;
        rows(idx) = rows_block(:);
        cols(idx) = cols_block(:);
        vals(idx) =J(:);
        entry_count = entry_count + num_local*(N-2);

        %fix the first and last cell value
        rows(entry_count+(0:1))=[1 N]+(N+1)*2;
        cols(entry_count+(0:1))=[1 N]+(N+1)*2;
        vals(entry_count+(0:1))=1;
        RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2)-phi_old([1 N]);
        entry_count=entry_count+2;

        %% Enthalpy transport assembly
        num_local=8;

        T_0=X(I-1+(N+1)*2+N);
        T_1=X(I  +(N+1)*2+N);
        T_2=X(I+1+(N+1)*2+N);
        RHS(I+(N+1)*2+N)=rhs_ent(uf_1, uf_2, phi_0, phi_1, phi_2, T_0, T_1, T_2, dz(i0)', dz(I)', dz(i2)', dt, cp, Lf, kt0, H_old(I),H_scaling);

        J = Jac_ent(uf_1, uf_2, phi_0, phi_1, phi_2, T_0, T_1, T_2, dz(i0)', dz(I)', dz(i2)', dt, cp, Lf, kt0, H_old(I),H_scaling);
        rows_block = repmat(I+(N+1)*2+N, 1, num_local);
        cols_block = [ ...
            I+(N+1), I+1+(N+1), ...
            I-1+(N+1)*2, I+(N+1)*2, I+1+(N+1)*2, ...
            I-1+(N+1)*2+N, I+(N+1)*2+N, I+1+(N+1)*2+N, ...
            ];

        idx = entry_count : entry_count + num_local*(N-2) - 1;
        rows(idx) = rows_block(:);
        cols(idx) = cols_block(:);
        vals(idx) =J(:);
        entry_count = entry_count + num_local*(N-2);

        %fix the first and last cell value
        rows(entry_count+(0:1))=[1 N]+(N+1)*2+N;
        cols(entry_count+(0:1))=[1 N]+(N+1)*2+N;
        vals(entry_count+(0:1))=1;
        RHS([1 N]+(N+1)*2+N)=X([1 N]+(N+1)*2+N)-T_old([1 N]);
        entry_count=entry_count+2;

        %% Solidus and liquidus
        if Has_volatile==1
            S_0=X(I-1+(N+1)*2+N*4);
            S_1=X(I  +(N+1)*2+N*4);
            S_2=X(I+1+(N+1)*2+N*4);

            cs2_0=X(I-1+(N+1)*2+N*5);
            cs2_1=X(I  +(N+1)*2+N*5);
            cs2_2=X(I+1+(N+1)*2+N*5);

            cl2_0=X(I-1+(N+1)*2+N*6);
            cl2_1=X(I  +(N+1)*2+N*6);
            cl2_2=X(I+1+(N+1)*2+N*6);
            num_local=7;
            lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts(I), index_cb2_nts(I));
            lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl(I), index_cb2_ntl(I));

            coef_ts = Ts0_coefficient(lin_ts, :);   % (nI × 4)
            coef_tl = Tl0_coefficient(lin_tl, :);   % (nI × 4)

            % unpack
            a_ts = coef_ts(:,1); b_ts = coef_ts(:,2);
            c_ts = coef_ts(:,3); d_ts = coef_ts(:,4);

            a_tl = coef_tl(:,1); b_tl = coef_tl(:,2);
            c_tl = coef_tl(:,3); d_tl = coef_tl(:,4);

            % -------------------------
            % 4. compute A,B,C,D
            % -------------------------
            P = Pressure(I);

            A = b_ts + c_ts .* (P/8000);
            B = a_ts .* (P/8000) + d_ts;

            C = b_tl + c_tl .* (P/40e3);
            D = a_tl .* (P/40e3) + d_tl;

            % -------------------------
            % 5. parameters
            % -------------------------
            Param = [A, B, C, D, Par_v*ones(nI,1)];
            RHS(I + (N+1)*2 + N*2)=rhs_solidus (phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
            RHS(I + (N+1)*2 + N*3)=rhs_liquidus(phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);

            J = Jac_solidus(phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
            rows_block = repmat(I+(N+1)*2+N*2, 1, num_local);
            cols_block = [I+(N+1)*2, I+(N+1)*2+N, I+(N+1)*2+N*2, I+(N+1)*2+N*3, I+(N+1)*2+N*4, I+(N+1)*2+N*5, I+(N+1)*2+N*6];
            idx = entry_count : entry_count + num_local*(N-2) - 1;
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);

            J = Jac_liquidus(phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);
            rows_block = repmat(I+(N+1)*2+N*3, 1, num_local);
            idx = entry_count : entry_count + num_local*(N-2) - 1;
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);
        else
            num_local=4;
            Param = repmat([Ts0(1), Tl0(1), 0, 0, 0], nI, 1);

            RHS(I + (N+1)*2 + N*2)=rhs_solidus (phi_1, T_1, cs_1, cl_1,  Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
            RHS(I + (N+1)*2 + N*3)=rhs_liquidus(phi_1, T_1, cs_1, cl_1,  Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);

            J = Jac_solidus(phi_1, T_1, cs_1, cl_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
            rows_block = repmat(I+(N+1)*2+N*2, 1, num_local);
            cols_block = [I+(N+1)*2, I+(N+1)*2+N, I+(N+1)*2+N*2, I+(N+1)*2+N*3];
            idx = entry_count : entry_count + num_local*(N-2) - 1;
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);

            J = Jac_liquidus(phi_1, T_1, cs_1, cl_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);
            rows_block = repmat(I+(N+1)*2+N*3, 1, num_local);
            idx = entry_count : entry_count + num_local*(N-2) - 1;
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);
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
        %%
        if Has_volatile ==1
            %% Volatile component transport assembly
            num_local=16;
            k2_stable1=kf*1e-1*(cellz(I)'<11e3);%*(phi_1>1e-3)*0
            k2_stable2=[k2_stable1(2:end); 0];

            if ~isempty(System_top)
                kf1=kf*ones(N-2,1);
                kf1(System_top:end)=0;
                kf1(1:System_bottom)=0;
                kf2=[kf1(2:end);0];

                % kf2=kf*(cb2(I)>1e-4);
                % kf2(System_top:end)=0;
                % kf2(1:System_bottom)=0;
                % kf1=[0; kf2(1:end-1)];
            else
                kf1=zeros(N-2,1);
                kf2=zeros(N-2,1);
            end


            RHS(I+(N+1)*2+N*4)=rhs_ct2(um_1, um_2, uf_1, uf_2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2, dz(i0)',dz(I)', dz(i2)', dt,  kf1, kf2, Cb2_old(I), k2_stable1, k2_stable2);

            cols_block = [ ...
                I, I+1, ...
                I+(N+1), I+1+(N+1), ...
                I-1+(N+1)*2, I+(N+1)*2, I+1+(N+1)*2, ...
                I-1+(N+1)*2+N*4, I+(N+1)*2+N*4, I+1+(N+1)*2+N*4, ...
                I-1+(N+1)*2+N*5, I+(N+1)*2+N*5, I+1+(N+1)*2+N*5, ...
                I-1+(N+1)*2+N*6, I+(N+1)*2+N*6, I+1+(N+1)*2+N*6];
            idx = entry_count : entry_count + num_local*(N-2) - 1;
            J=Jac_ct2(um_1, um_2, uf_1, uf_2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2, dz(i0)',dz(I)', dz(i2)', dt,  kf1, kf2, Cb2_old(I), k2_stable1, k2_stable2);
            rows_block = repmat(I+(N+1)*2+N*4, 1, num_local);
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);

            %fix the first and last cell value
            rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*4;
            cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*4;
            vals(entry_count+(0:1))=1;
            RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*4)-S_old([1 N]);
            entry_count=entry_count+2;
            %% solid saturation
            num_local=5;
            RHS(I+(N+1)*2+N*5)=rhs_ssat( T_1, cs_1, S_1, cs2_1, cl2_1, S_cap(1), S_cap(2), Par_v);
            cols_block = [I+(N+1)*2+N, I+(N+1)*2+N*2, I+(N+1)*2+N*4, I+(N+1)*2+N*5, I+(N+1)*2+N*6];

            idx = entry_count : entry_count + num_local*(N-2) - 1;
            J=Jac_ssat( T_1, cs_1, S_1, cs2_1, cl2_1, S_cap(1), S_cap(2), Par_v);
            J=repmat(J,nI,1); %for now, J is constant so can not automatically broadcasted

            rows_block = repmat(I+(N+1)*2+N*5, 1, num_local);
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);
            %% melt saturation
            num_local=3;
            RHS(I+(N+1)*2+N*6)=rhs_lsat(T_1, S_1,  cl2_1, dSdT(I)/100, b0_all(I));
            cols_block = [I+(N+1)*2+N, I+(N+1)*2+N*4,  I+(N+1)*2+N*6];

            idx = entry_count : entry_count + num_local*(N-2) - 1;
            J=Jac_lsat(T_1, S_1,  cl2_1, dSdT(I)/100, b0_all(I));
            rows_block = repmat(I+(N+1)*2+N*6, 1, num_local);
            rows(idx) = rows_block(:);
            cols(idx) = cols_block(:);
            vals(idx) =J(:);
            entry_count = entry_count + num_local*(N-2);

            rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*5;
            cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*5;
            vals(entry_count+(0:1))=1;
            RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*5)-Cs2_old([1 N]);
            entry_count=entry_count+2;

            rows(entry_count+(0:1))=[1 N]+(N+1)*2+N*6;
            cols(entry_count+(0:1))=[1 N]+(N+1)*2+N*6;
            vals(entry_count+(0:1))=1;
            RHS([1 N]+(N+1)*2)=X([1 N]+(N+1)*2+N*6)-Cl2_old([1 N]);
            entry_count=entry_count+2;
        end
        %%
        % RHS(new_frozen)=0;
        norm_R = max(abs(RHS));
        if norm_R < Precision
            con_scaling=min(1/max(abs(X(1:2*N+2)))/1e3,1e5);
            % if iter<10 && (dt>Min_dt || Advance_time==0)
            %     % con_scaling=con_scaling*2;
            %     % con_scaling=min(1e6,con_scaling);
            %     k_stable_value=k_stable_value*0.9;
            %     k_stable_value=max(k_stable_value,k_stable_value0);
            %     H_scaling=H_scaling*2;
            %     H_scaling=min(1e-2,H_scaling);
            %     if Has_volatile==1
            %         Precision=Precision/2;
            %         Precision=max(Precision,Precision0);
            %     end
            % else
            %     % con_scaling=con_scaling/10;
            %     % con_scaling=max(con_scaling,1e3);
            %     k_stable_value=k_stable_value*10;
            %     k_stable_value=min(k_stable_value,1e-9);
            %     H_scaling=H_scaling/10;
            %     H_scaling=max(1e-4,H_scaling);
            %     if Has_volatile==1
            %         Precision=Precision*5;
            %         Precision=min(Precision,Precision0*100);
            %     end
            % end
            prevStr=Display_monitor(Time/Year,Time_intended/Year, Time_gap/Year, dt/Year, iter, prevStr,conservation1,conservation2);
            converged = true;
            bad_count = zeros(N+1,1);
            new_frozen=[];
            break;
        end

        % bad_nodes=find(abs(RHS(1:N+1))>Precision); %detect_trouble_nodes(RHS,N,Precision);
        % 
        % bad_count=bad_count + ismember(1:N+1, bad_nodes)';
        % bad_count(~ismember(1:N+1,bad_nodes)) = max(bad_count(~ismember(1:N+1,bad_nodes)) - 1, 0);
        % new_frozen = find(bad_count >= Kfreeze);

        Matrix_A = sparse(rows(1:entry_count-1), cols(1:entry_count-1), vals(1:entry_count-1), Dof, Dof);
        
        % for k=1:length(new_frozen)
        %     Matrix_A(new_frozen(k),:)=0;
        %     Matrix_A(new_frozen(k),new_frozen(k))=1;            
        % end
        % RHS(new_frozen)=0;
        du = Matrix_A \ (-RHS);
        if ~isreal(du) || any(isnan(du))
            dt_intended=dt_intended*0.5;
            dt=dt/2;
            % if dt<1e-4*Year
            %     dt=1*Year;
            % end
            X=X0;
            Not_improve=0;
            Norm_pre=1e3;
            % disp(['Decrease dt, dt=' num2str(dt/Year)])
            continue
            % iter=Max_Newton_iter;
            % break
        end
        % Force constant values on boundary
        du([[1 N]+(N+1)*2 [1 N]+(N+1)*2+N [1 N]+(N+1)*2+N*2 [1 N]+(N+1)*2+N*3])=0;

        alpha = 1.0;
        for LS = 1:4
            X = X_pre + alpha * du;
            % force 1>phi>0
            X(2*(N+1)+1:2*(N+1)+N)=max(X(2*(N+1)+1:2*(N+1)+N),-0.99e-2);
            % X(2*(N+1)+1:2*(N+1)+N)=min(X(2*(N+1)+1:2*(N+1)+N),1);

            % force 1350>T
            % X(2*(N+1)+N+1:2*(N+1)+N*2)=min(X(2*(N+1)+N+1:2*(N+1)+N*2),1350);
            r=(X((2:N-1) + 2*(N+1))-X((1:N-2) + 2*(N+1)))./(X((3:N) + 2*(N+1))-X((2:N-1) + 2*(N+1))+1e-15);
            flux_limiter=(r+abs(r))./(1+abs(r))/2*0+1;

            RHS(:)=0;
            index_phi_nmus=max(min(floor(X((1:N)+(N+1)*2) * N_mus) + 1, N_mus),1);

            index_phi_nc=max(min(floor(X((1:N)+(N+1)*2) * N_C) + 1, N_C),1);
            index_cl_nc=max(min(floor(X((1:N)+(N+1)*2+N*3) * N_C) + 1, N_C),1);

            index_cs_nrhos=max(min(floor(X((1:N)+(N+1)*2+N*2) * N_rhos) + 1, N_rhos),1);
            index_cl_nrhol=max(min(floor(X((1:N)+(N+1)*2+N*3) * N_rhol) + 1, N_rhol),1);
            if Has_volatile==1
                index_cl2_nc=max(min(floor(X((1:N)+(N+1)*2+N*6) * N_C) + 1, N_C),1);
                cb2=X((1:N)+(N+1)*2).*X((1:N)+(N+1)*2+N*6)+(1-X((1:N)+(N+1)*2)).*X((1:N)+(N+1)*2+N*5)+X((1:N)+(N+1)*2+N*4);
                index_cb2_nts=max(min(floor(cb2/0.13/Par_v * N_Ts) + 1, N_Ts),1);
                index_cb2_ntl=max(min(floor(cb2/0.2/Par_v * N_Ts) + 1, N_Ts),1);

                % index_cl2_nrhol=max(min(floor(X((1:N)+(N+1)*2+N*6)/0.1 * N_rhol) + 1, N_rhol),1);
            end

            I = (2:N)';
            nI = numel(I);

            i0 = I - 1;
            i1 = I;
            i2 = I + 1;
            um_0 = X(i0);
            um_1 = X(i1);
            um_2 = X(i2);
            uf_0 = X(i0 + (N+1));
            uf_1 = X(i1 + (N+1));
            phi_0 = X(i0 + 2*(N+1));
            phi_1 = X(i1 + 2*(N+1));
            cs_0 = X(i0 + 2*(N+1) + 2*N);
            cs_1 = X(i1 + 2*(N+1) + 2*N);

            cl_0 = X(i0 + 2*(N+1) + 3*N);
            cl_1 = X(i1 + 2*(N+1) + 3*N);

            if Has_volatile == 2
                cl2_0 = X(i0 + 2*(N+1) + 6*N);
                cl2_1 = X(i1 + 2*(N+1) + 6*N);
                cl2_2 = X(i2 + 2*(N+1) + 6*N);
            end

            if Has_volatile == 2
                %umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1; cl2_0; cl2_1
                Vars = [ ...
                    um_0, um_1, um_2, ...
                    uf_1, ...
                    phi_0, phi_1, ...
                    cs_0, cs_1,  ...
                    cl_0, cl_1,  ...
                    cl2_0, cl2_1,...
                    ];
            else
                %umi_0; umi_1; umi_2; ufi; phi_0; phi_1; cs_0; cs_1; cl_0; cl_1
                Vars = [ ...
                    um_0, um_1, um_2, ...
                    uf_1, ...
                    phi_0, phi_1, ...
                    cs_0, cs_1,  ...
                    cl_0, cl_1 ...
                    ];
            end

            Param = zeros(nI, 28);
            Param(:,1:2) = mus_coef_all(index_phi_nmus(i0),:);
            Param(:,3:4) = mus_coef_all(index_phi_nmus(i1),:);
            % coupling coefficients
            if Has_volatile == 2
                idx_cl2 = max(min(floor(X(i1 + (N+1)*2 + 6*N) * N_C) + 1, N_C),1);
                Param(:,5:12) = C_coef_all(index_phi_nc(i1), index_cl_nc(i1), idx_cl2, :);
            else
                lin = sub2ind([N_C, N_C], index_phi_nc(i1), index_cl_nc(i1));
                Param(:,5:8) = C_coef_all(lin, :);
                % Param(:,5:8) = C_coef_all(index_phi_nc(i1), index_cl_nc(i1), :);
            end
            % density solid
            Param(:,13:14) = rhos_coef_all(index_cs_nrhos(i0),:);
            Param(:,15:16) = rhos_coef_all(index_cs_nrhos(i1),:);

            % density liquid
            if Has_volatile == 2
                lin = sub2ind([N_rhol, N_rhol], index_cl_nrhol(i0), index_cl2_nrhol(i0));
                Param(:,17:20) = rhol_coef_all(lin, :);
                lin = sub2ind([N_rhol, N_rhol], index_cl_nrhol(i1), index_cl2_nrhol(i1));
                Param(:,21:24) = rhol_coef_all(lin, :);
            else
                Param(:,17:18) = rhol_coef_all(index_cl_nrhol(i0), :);
                Param(:,21:22) = rhol_coef_all(index_cl_nrhol(i1), :);
            end
            % mum_0_local=sqrt((mus_coef_all(index_phi_nmus(i0),1).*phi_0+mus_coef_all(index_phi_nmus(i0),2)+mus_coef_all(index_phi_nmus(i1),1).*phi_1+mus_coef_all(index_phi_nmus(i1),2)))./(dz(2:end)+dz(1+end-1))'.^2*2;
            mum_0_local=min(mus_coef_all(index_phi_nmus(i0),1).*phi_0+mus_coef_all(index_phi_nmus(i0),2))*ones(nI,1);

            % physical constants
            Param(:,25:28) = [dz(i0)', dz(i1)', g*ones(nI,1), mum_0_local]; %mum_0*ones(nI,1)
            RHS(I)=rhs_mom(Vars(:,1),Vars(:,2),Vars(:,3),Vars(:,4),Vars(:,5),Vars(:,6),Vars(:,7),Vars(:,8),Vars(:,9),Vars(:,10),...
                Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5),Param(:,6),Param(:,7),Param(:,8),Param(:,9),Param(:,10),...
                Param(:,11),Param(:,12),Param(:,13),Param(:,14),Param(:,15),Param(:,16),Param(:,17),Param(:,18),Param(:,19),Param(:,20),...
                Param(:,21),Param(:,22),Param(:,23),Param(:,24),Param(:,25),Param(:,26),Param(:,27),Param(:,28));
            
            RHS(I+N+1)=rhs_con(um_1,uf_1,phi_0,phi_1, con_scaling,[flux_limiter; 1]);



            %
            I = (2:N-1)';
            i0=I-1;
            i2=I+1;
            nI = numel(I);
            um_1=X(I);
            um_2=X(I+1);
            uf_1=X(I+(N+1));
            uf_2=X(I+1+(N+1));
            phi_0=X(I-1+(N+1)*2);
            phi_1=X(I  +(N+1)*2);
            phi_2=X(I+1+(N+1)*2);
            cs_0=X(I-1+(N+1)*2+N*2);
            cs_1=X(I  +(N+1)*2+N*2);
            cs_2=X(I+1+(N+1)*2+N*2);
            cl_0=X(I-1+(N+1)*2+N*3);
            cl_1=X(I  +(N+1)*2+N*3);
            cl_2=X(I+1+(N+1)*2+N*3);

            scale=ones(N,1);
            scale(Cb_old<2e-3)=1e3;
            k_stable1=k_stable_value*(phi_1>1e-2).*scale(I);
            k_stable2=[k_stable1(2:end);0];
            RHS(I+(N+1)*2)=rhs_ct(um_1,um_2, uf_1, uf_2, phi_0, phi_1, phi_2, cs_0,cs_1,cs_2, cl_0,cl_1,cl_2, dt, dz(i0)', dz(I)', dz(i2)',Cb_old(I),k_stable1, k_stable2,[1; flux_limiter(1:end-1)],flux_limiter);

            T_0=X(I-1+(N+1)*2+N);
            T_1=X(I  +(N+1)*2+N);
            T_2=X(I+1+(N+1)*2+N);
            RHS(I+(N+1)*2+N)=rhs_ent(uf_1, uf_2, phi_0, phi_1, phi_2, T_0, T_1, T_2, dz(i0)', dz(I)', dz(i2)', dt, cp, Lf, kt0, H_old(I), H_scaling);

            %
            if Has_volatile==1
                S_0=X(I-1+(N+1)*2+N*4);
                S_1=X(I  +(N+1)*2+N*4);
                S_2=X(I+1+(N+1)*2+N*4);

                cs2_0=X(I-1+(N+1)*2+N*5);
                cs2_1=X(I  +(N+1)*2+N*5);
                cs2_2=X(I+1+(N+1)*2+N*5);

                cl2_0=X(I-1+(N+1)*2+N*6);
                cl2_1=X(I  +(N+1)*2+N*6);
                cl2_2=X(I+1+(N+1)*2+N*6);
                num_local=7;
                lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts(I), index_cb2_nts(I));
                lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl(I), index_cb2_ntl(I));

                coef_ts = Ts0_coefficient(lin_ts, :);   % (nI × 4)
                coef_tl = Tl0_coefficient(lin_tl, :);   % (nI × 4)

                % unpack
                a_ts = coef_ts(:,1); b_ts = coef_ts(:,2);
                c_ts = coef_ts(:,3); d_ts = coef_ts(:,4);

                a_tl = coef_tl(:,1); b_tl = coef_tl(:,2);
                c_tl = coef_tl(:,3); d_tl = coef_tl(:,4);

                % -------------------------
                % 4. compute A,B,C,D
                % -------------------------
                P = Pressure(I);

                A = b_ts + c_ts .* (P/8000);
                B = a_ts .* (P/8000) + d_ts;

                C = b_tl + c_tl .* (P/40e3);
                D = a_tl .* (P/40e3) + d_tl;

                % -------------------------
                % 5. parameters
                % -------------------------
                Param = [A, B, C, D, Par_v*ones(nI,1)];
                RHS(I + (N+1)*2 + N*2)=rhs_solidus (phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
                RHS(I + (N+1)*2 + N*3)=rhs_liquidus(phi_1, T_1, cs_1, cl_1, S_1, cs2_1, cl2_1, Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);
            else
                num_local=4;
                Param = repmat([Ts0(1), Tl0(1), 0, 0, 0], nI, 1);

                RHS(I + (N+1)*2 + N*2)=rhs_solidus (phi_1, T_1, cs_1, cl_1,  Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), n_order);
                RHS(I + (N+1)*2 + N*3)=rhs_liquidus(phi_1, T_1, cs_1, cl_1,  Param(:,1),Param(:,2),Param(:,3),Param(:,4),Param(:,5), A1);
            end


            if Has_volatile==1
                k2_stable1=kf*1e-1*(cellz(I)'<11e3);%*(phi_1>1e-3)*0
                k2_stable2=[k2_stable1(2:end); 0];
                if ~isempty(System_top)
                    kf1=kf*ones(N-2,1);                
                    kf1(System_top:end)=0;
                    kf1(1:System_bottom)=0;
                    kf2=[kf1(2:end);0];
                    % kf2=kf*(cb2(I)>1e-4);
                    % kf2(System_top:end)=0;
                    % kf2(1:System_bottom)=0;
                    % kf1=[0; kf2(1:end-1)];
                else
                    kf1=zeros(N-2,1);
                    kf2=zeros(N-2,1);
                end

                RHS(I+(N+1)*2+N*4)=rhs_ct2(um_1, um_2, uf_1, uf_2, phi_0, phi_1, phi_2, S_0, S_1, S_2, cs2_0, cs2_1, cs2_2, cl2_0, cl2_1, cl2_2, dz(i0)',dz(I)', dz(i2)', dt,  kf1, kf2, Cb2_old(I), k2_stable1, k2_stable2);

                RHS(I+(N+1)*2+N*5)=rhs_ssat( T_1, cs_1, S_1, cs2_1, cl2_1, S_cap(1), S_cap(2), Par_v);
                RHS(I+(N+1)*2+N*6)=rhs_lsat( T_1,       S_1,        cl2_1, dSdT(I)/100, b0_all(I));
                % for i=2:N-1
                %     column_index=[i+(N+1)*2, i+(N+1)*2+N, i+(N+1)*2+N*2, i+(N+1)*2+N*4,i+(N+1)*2+N*5, i+(N+1)*2+N*6 ];
                %     in=X(column_index);
                %     a0=dSdT(i)/100;
                %     b0=b0_all(i);
                %     RHS(i+(N+1)*2+N*5)=rhs_ssat([in;S_cap(1); S_cap(2); Par_v]);
                %     RHS(i+(N+1)*2+N*6)=rhs_lsat([in;a0;b0]);
                % end
            end

            temp_norm=max(abs(RHS));
            if temp_norm < (1-alpha*1e-4)*norm_R
                break;
            end
            alpha = alpha * 0.5;
        end
        
        if ~isreal(X) || any(isnan(X((1:N)  +(N+1)*2+N)))
            dt_intended=dt_intended*0.5;
            dt=dt/2;
            if dt<Min_dtY*Year
                dt=1*Year;
            end
            X=X0;
            Not_improve=0;
            Norm_pre=1e3;
            % disp(['Decrease dt, dt=' num2str(dt/Year)])
            continue
        end
        if temp_norm < Precision
            con_scaling=min(1/max(abs(X(1:2*N+2)))/1e3,1e5);
            % if iter<10 && (dt>Min_dt || Advance_time==0)
            %     % con_scaling=con_scaling*2;
            %     % con_scaling=min(1e6,con_scaling);
            %     k_stable_value=k_stable_value*0.9;
            %     k_stable_value=max(k_stable_value,k_stable_value0);
            %     H_scaling=H_scaling*2;
            %     H_scaling=min(1,H_scaling);
            %     if Has_volatile==1
            %         Precision=Precision/2;
            %         Precision=max(Precision,Precision0);
            %     end
            % else
            %     % con_scaling=con_scaling/10;
            %     % con_scaling=max(con_scaling,1e3);
            %     k_stable_value=k_stable_value*10;
            %     k_stable_value=min(k_stable_value,1e-9);
            %     H_scaling=H_scaling/10;
            %     H_scaling=max(1e-4,H_scaling);
            %     if Has_volatile==1
            %         Precision=Precision*5;
            %         Precision=min(Precision,Precision0*100);
            %     end
            % end

            prevStr=Display_monitor(Time/Year,Time_intended/Year, Time_gap/Year, dt/Year, iter, prevStr,conservation1, conservation2);

            converged = true;
            bad_count = zeros(N+1,1);
            new_frozen=[];
            break;
        end

        % we expect the Newton method to converge with certain speed, otherwise
        % the time step is too big, we reset everything with smaller time step
        if temp_norm<Norm_pre*0.99 %|| dt/dt0<1e-4
            Norm_pre=temp_norm;
        else
            Not_improve=Not_improve+1;
            if Not_improve>5 && Advance_time==1
                % if Has_volatile==0
                %     % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old];
                %     X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N)];
                % else
                %     % X=[um_old; uf_old; phi_old; T_old; Cs_old; Cl_old; S_old; Cs2_old; Cl2_old];
                %     X=[u_all_old(N+2:2*N+2); u_all_old(1:N+1); phi_old(1:N); T_old; C_all_old(N+1:2*N); C_all_old(1:N); S_old; Cs2_old; Cl2_old];
                % end
                X=X0;
                dt_intended=dt_intended*0.5;
                dt=dt/2;

                Not_improve=0;
                Norm_pre=1e3;
            end
        end
        X_pre=X;
    end
    Cb=X((1:N)+2*N+2).*X((1:N)+2*N+2+N*3)+(1-X((1:N)+2*N+2)).*X((1:N)+2*N+2+N*2);
    if Has_volatile==1
        Cb2=X((1:N)+2*N+2+N*6).*X((1:N)+2*N+2)+X((1:N)+2*N+2+N*5).*(1-X((1:N)+2*N+2))+X((1:N)+2*N+2+N*4);
    end

    if iter>=Max_iter
        error('Newton iteration fails')
    end
    Time=Time+dt;
    if Advance_time==0
        break
    elseif dt<=Min_dt*1e-1
        if Time-Last_adapted>Adaptive_step_time
            Adapt_mesh_master    
            Last_adapted=Time;
            dt_intend=Min_dt;
        end
    end

    
    if iter<=9
        if dt<0.01*Max_dt
            dt=min(dt*1.3, Max_dt);
        elseif dt<0.05*Max_dt
            dt=min(dt*1.2, Max_dt);
        else
            dt=min(dt*1.1, Max_dt);
        end
        dt_intended=max(dt,dt_intended);
    elseif iter>=10
        dt=min(dt*0.9, Max_dt);
        dt_intended=min(dt,dt_intended);
        dt_intended=max(Min_dt*1e-3,dt_intended);
    end
    if Time<Time_intended %update all old values
        Cb_old=Cb;
        if Has_volatile==1
            Cb2_old=Cb2; %X((1:N)+2*N+2+N*6).*X((1:N)+2*N+2)+X((1:N)+2*N+2+N*5).*(1-X((1:N)+2*N+2))+X((1:N)+2*N+2+N*4);
        end
        H_old=X((1:N)+2*N+2)*Lf+X((1:N)+2*N+2+N)*cp;      
    end
end
%%
dt=0; %MUST SET TO 0, so Time won't change outside of Newton solver.
%%

% um=X(1:N+1);
% uf=X(N+2:2*N+2);
u_all=[X(N+2:2*N+2); X(1:N+1)];


phi=[X((1:N)+2*N+2); 1-X((1:N)+2*N+2)];
T=X((1:N)+2*N+2+N);
Cs=X((1:N)+2*N+2+N*2);
Cl=X((1:N)+2*N+2+N*3);
C_all=[Cl; Cs];

% not used for calculation
% Cb=X((1:N)+2*N+2).*X((1:N)+2*N+2+N*3)+(1-X((1:N)+2*N+2)).*X((1:N)+2*N+2+N*2);
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
    lin_ts = sub2ind([N_Ts, N_Ts], index_pressure_nts, index_v_nts);   
    coef=Ts0_coefficient(lin_ts,:);
    for i=1:N
        % coef=Ts0_coefficient(index_pressure_nts(i),index_v_nts(i),:);        
        Ts0(i)=coef(i,1)*Pressure(i)/8000+coef(i,2)*Cb2(i)/Par_v/13*100+coef(i,3)*Pressure(i)*Cb2(i)/Par_v/8000/13*100+coef(i,4);
    end

    index_v_ntl=max(min(floor(Cb2/20*100/Par_v*N_Ts)+1,N_Ts),1);
    lin_tl = sub2ind([N_Tl, N_Tl], index_pressure_ntl, index_v_ntl);
    coef=Tl0_coefficient(lin_tl,:);
    for i=1:N
        % coef=Tl0_coefficient(index_pressure_ntl(i),index_v_ntl(i),:);        
        Tl0(i)=coef(i,1)*Pressure(i)/40e3+coef(i,2)*Cb2(i)/Par_v/20*100+coef(i,3)*Pressure(i)*Cb2(i)/Par_v/400/20+coef(i,4);
    end
    
    Tl=A1*Cb.^2+(Ts0-Tl0-A1).*Cb+Tl0;
    if is_eutectic==1
        Ts=Ts0;
    else         
        Ts=Tl0-max(Cb,0).^(1/n_order).*(Tl0-Ts0);    
    end
else
    
    Tl=A1*Cb.^2+B1*Cb+C1; %local liquidus
    if is_eutectic==1
        Ts=Ts0(1)*ones(N,1);
    else
        % Ts=(1-max(Cb,1e-12).^(1/n_order))*(Tl0-Ts0)+Ts0; 
        Ts=Tl0(1)-max(Cb,0).^(1/n_order).*(Tl0(1)-Ts0(1)); %local solidus        
    end
end




improve1=norm_R;
improve2=improve1;


%% variables that needs to be updated for Newton
rhom=zeros(N+1,1);
rhof=zeros(N+1,1);
rho_b=zeros(N,1);


phi_old=phi;
u_all_old=u_all;
H_old=H;
T_old=T;
C_all_old=C_all;
Cb_old=Cb;

if Has_volatile==1
    S_old=S;
    Cs2_old=Cs2;
    Cl2_old=Cl2;
end



function prevStr=Display_monitor(Time,Time_intended, Time_gap, dt, iter, prevStr, conservation1, conservation2)
    frac=(Time-Time_intended+Time_gap)/Time_gap;
    nfill = round(frac*40);
    bar = [repmat('=',1,nfill),repmat(' ',1,40-nfill)];
    str = sprintf([ ...
        'Time: %10.4f   Output time step: %6.1f\n' ...
        'dt:   %10.5f   Newton iter:  %6d\n' ...
        '[%s] %6.2f %%\n'...
        'Conservation: %3.8f  %3.8f'], ...
        Time, Time_gap, dt, iter, bar, 100*frac,conservation1,conservation2);
    fprintf(repmat('\b',1,length(prevStr)));
    fprintf('%s', str);
    prevStr = str;
end
% 
function bad_nodes=detect_trouble_nodes(RHS,N, tol)
    rnode=RHS(1:N+1);
    % Rc=RHS(N+2:2*N+2);
    % Rt=RHS(2*N+2+(1:N));
    % rnode=sqrt(Rm.^2+Rc.^2+Rt.^2);

    % rmed = median(abs(Rm));
    active = abs(rnode) > tol;    
    if max(active)<1
        bad_nodes=[];
        return;
    end
    r_active=RHS(active);
    r_sorted = sort(r_active);
    bulk = r_sorted(1:max(1, round(0.9*length(r_sorted))));
    r_bulk_max = max(bulk);

    alpha = 1e3;
    bad_nodes = find(active & rnode > alpha * r_bulk_max);
end