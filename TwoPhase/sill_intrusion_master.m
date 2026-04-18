%CAB 30/1/23
%Sill intrusion function file.
%% sill intruson master

%% Start and end nodes of sill depths

%Nodez for intrusion depth for (N+1)x1 arrays
Start_depthN = depthN;%-SillNodez; %N+1- depthN;

% Cellz for for intrusion depth for Nx1 arrays
Start_depthN_c =  depthCN;%-SillNodez;

%Nodez for intrusion depth for the second half of the 2(N+1)x1 arrays
Start_depth2N = N+1+depthN;%-SillNodez; %2*(N+1)-depthN;

%Cellz for intrusion depth for the second half of the 2N arrays
Start_depth2N_c = N+depthCN;%-SillNodez;%2*N-depthCN;

%% Update mesh parameters

% dz

dz=[dz(1:(Start_depthN_c)),dzF.*ones(1,SillNodez),dz(Start_depthN_c+1:end)];

% nodez
nodez=zeros(1,N+1+SillNodez);
for i=2:N+SillNodez+1
    nodez(i) = nodez(i-1)+dz(i-1);
end

%N
N=length(nodez)-1;

%cellz
cellz = (nodez(1:end-1)+nodez(2:end))/2; % cell center points. 


%% Intrusion in Nx1 arrays

%Bulk composition
Cb = [Cb(1:Start_depthN_c); Sillcomp.*ones(SillNodez,1);Cb(Start_depthN_c+1:end)];

%Enthalpy
H = [H(1:Start_depthN_c); SillH.*ones(SillNodez,1);H(Start_depthN_c+1:end)];

%Temperature
T = [T(1:Start_depthN_c); SillTemp.*ones(SillNodez,1);T(Start_depthN_c+1:end)];

if Has_volatile==1
    S=[S(1:Start_depthN_c); injection_S.*ones(SillNodez,1);S(Start_depthN_c+1:end)];
    Cl2=[Cl2(1:Start_depthN_c); injection_Cl2.*ones(SillNodez,1);Cl2(Start_depthN_c+1:end)];
    Cs2=[Cs2(1:Start_depthN_c); injection_Cs2.*ones(SillNodez,1);Cs2(Start_depthN_c+1:end)];
end

if FourMPD==1
    %Olivine
    ol = [ol(1:Start_depthN_c); zeros(SillNodez,1);ol(Start_depthN_c+1:end)];
    
    %Orthopyroxene
    opx = [opx(1:Start_depthN_c); 0.*ones(SillNodez,1);opx(Start_depthN_c+1:end)];
    
    %Clinopyroxene 
    cpx = [cpx(1:Start_depthN_c); 0.*ones(SillNodez,1);cpx(Start_depthN_c+1:end)];
    
    %Feldspar
    feld = [feld(1:Start_depthN_c); 0.*ones(SillNodez,1);feld(Start_depthN_c+1:end)];
    
    %Olivine mg
    ol_mg = [ol_mg(1:Start_depthN_c); zeros(SillNodez,1);ol_mg(Start_depthN_c+1:end)];
    
    %Pyroxene  mg 
    px_mg = [px_mg(1:Start_depthN_c); 0.*ones(SillNodez,1);px_mg(Start_depthN_c+1:end)];
    
    %Orthopyroxene mg
    opx_mg = [opx_mg(1:Start_depthN_c); 0.*ones(SillNodez,1);opx_mg(Start_depthN_c+1:end)];
    
    %Clinopyroxene  mg 
    cpx_mg = [cpx_mg(1:Start_depthN_c); 0.*ones(SillNodez,1);cpx_mg(Start_depthN_c+1:end)];
    
    %Feldspar mg 
    feld_mg = [feld_mg(1:Start_depthN_c); 0.*ones(SillNodez,1);feld_mg(Start_depthN_c+1:end)];
    
    %Solid Mn
    sol_mn = [sol_mn(1:Start_depthN_c); 0.*ones(SillNodez,1);sol_mn(Start_depthN_c+1:end)];
    
    %Olivine Mn
    ol_mn = [ol_mn(1:Start_depthN_c); 0.*ones(SillNodez,1);ol_mn(Start_depthN_c+1:end)];
end

%% Contributions
if To_cal_mc==1
   Con_T=[Con_T(1:Start_depthN_c); zeros(SillNodez,1); Con_T(Start_depthN_c+1:end)];  %contribution from temperature
   Con_H=[Con_H(1:Start_depthN_c); zeros(SillNodez,1); Con_H(Start_depthN_c+1:end)];  %contribution from entire enthalpy  
   Con_C=[Con_C(1:Start_depthN_c); zeros(SillNodez,1); Con_C(Start_depthN_c+1:end)];  %contribution from composition
   Con_C_dummy=[Con_C_dummy(1:Start_depthN_c); zeros(SillNodez,1); Con_C_dummy(Start_depthN_c+1:end)];  %contribution from composition dummy
   Con_M=[Con_M(1:Start_depthN_c); zeros(SillNodez,1); Con_M(Start_depthN_c+1:end)];  %contribution from compaction

   ACon_T=[ACon_T(1:Start_depthN_c); zeros(SillNodez,1); ACon_T(Start_depthN_c+1:end)]; %accumulated contribution from temperature
   ACon_M=[ACon_M(1:Start_depthN_c); zeros(SillNodez,1); ACon_M(Start_depthN_c+1:end)];   
   ACon_C=[ACon_C(1:Start_depthN_c); zeros(SillNodez,1); ACon_C(Start_depthN_c+1:end)]; 
   ACon_C_dummy=[ACon_C_dummy(1:Start_depthN_c); zeros(SillNodez,1); ACon_C_dummy(Start_depthN_c+1:end)]; 

   CD_M=[CD_M(1:Start_depthN_c); zeros(SillNodez,1); CD_M(Start_depthN_c+1:end)]; %CD caused by compaction
   CD_R=[CD_R(1:Start_depthN_c); zeros(SillNodez,1); CD_R(Start_depthN_c+1:end)];
   CD_R_dummy=[CD_R_dummy(1:Start_depthN_c); zeros(SillNodez,1); CD_R_dummy(Start_depthN_c+1:end)];
end
%% Intrusion in 2Nx1 arrays

%Melt Fraction / Solid Fraction
phi = [phi(1:Start_depthN_c); SillMF.*ones(SillNodez,1); phi(Start_depthN_c+1:Start_depth2N_c); (1-SillMF).*ones(SillNodez,1); phi(Start_depth2N_c+1:end)];

%C_all - Liquid comp / Solid comp
C_all = [C_all(1:Start_depthN_c); SillCf.*ones(SillNodez,1); C_all(Start_depthN_c+1:Start_depth2N_c); SillCm.*ones(SillNodez,1); C_all(Start_depth2N_c+1:end)];

%Cphi_all
Cphi_all = [C_all(1:N).*phi(1:N); C_all(N+1:2*N).*phi(N+1:2*N)];


%% Intrusion in 2(N+1)x1 arrays
u_all = [u_all(1:Start_depthN); zeros(SillNodez,1); u_all(Start_depthN+1:Start_depth2N); zeros(SillNodez,1); u_all(Start_depth2N+1:end)];

%% Update mass conservation
OG_composition=OG_composition+mean(Sillcomp)*SillNodez;
Conservative_pos(2)=Conservative_pos(2)+SillNodez;

%% Update base_crust
Base_crust=Base_crust+Sill_length/1000;
Initial_dN=Initial_dN+SillNodez;
Initial_dC=Initial_dC+SillNodez;
%depthTsCN=depthTsCN+SillNodez;

if HighTsF==1
    for i=1:1:N
        
        if i<=depthTsCN
            Ts(i)=TsD;
        else
         [Ts(i),Tl(i),~,~]=solid_state_Ts_Tl_local_master(HHJPet,SSPD,FourMPD, Cb(i),cp,Lf,DT,DC,lm1,lc1, crit_cb4,crit_cb5,s5T,...
                                                    liq_P2_C, liq_P3_C, liq_k1,liq_a1,liq_b1,liq_k2,liq_a2,liq_b2,liq_k3,liq_a3,liq_b3,A1,B1,C1,alpha, n_PD);
        end
        
    end
else
    for i=1:1:N
        [Ts(i),Tl(i),~,~]=solid_state_Ts_Tl_local_master(HHJPet,SSPD,FourMPD, Cb(i),cp,Lf,DT,DC,lm1,lc1, crit_cb4,crit_cb5,s5T,...
                                                    liq_P2_C, liq_P3_C, liq_k1,liq_a1,liq_b1,liq_k2,liq_a2,liq_b2,liq_k3,liq_a3,liq_b3,A1,B1,C1,alpha, n_PD);
    end
end


%%
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
dt=0;
Newton_solver2;

Just_intruded=1;
Last_intrusion_depth=nodez(depthN);
dt_intended=dt0Y*Year;