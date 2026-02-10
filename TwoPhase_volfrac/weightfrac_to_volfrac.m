function [phi,ol,opx,cpx,feld, rhof, rhom, rhool, rhoopx, rhocpx, rhofeld] = weightfrac_to_volfrac(T,phi_wt,ol_wt,opx_wt,cpx_wt, feld_wt,cl, ol_mg, opx_mg, cpx_mg, min_mg, max_mg, Ts_local,rhof_m, rhof_c, crit_mg_si_melt, rhof_2,rhool_m,rhool_c,...
                                                                                                                                    rhoopx_m, rhoopx_c, rhocpx_m, rhocpx_c, rhofeld_m, rhofeld_c,rhom_1)                                                                                                                
                                                                                                                                    
                                                                                                                               


%% This file needs to do the following
% 1. Calculate melt and mineral densities
% 2. Convert wt% to vol%
% 3. Calculate solid density.


%% Add these to input file
%rhof_m = 128.75; % rhof_m
%rhof_c = 2200; %rhof_c
%crit_mg_si_melt = 4;
%rhof_2=2715;

%rhool_m = -20.5;
%rhool_c = 4175;

%rhoopx_m = -21.4286;
%rhoopx_c = 3864.29;

%rhocpx_m = -12.5;
%rhocpx_c = 3562.5;

%rhofeld_m = 0.04819;
%rhofeld_c = 2621.45;

%% convert nd to dim

cl_dim = cl*(max_mg-min_mg)+min_mg;
ol_mg_dim = ol_mg*(max_mg-min_mg)+min_mg;
opx_mg_dim = opx_mg*(max_mg-min_mg)+min_mg;
cpx_mg_dim = cpx_mg*(max_mg-min_mg)+min_mg;





%% Calculating melt and mineral densities


    
% melt density 
if cl_dim>=crit_mg_si_melt
    rhof = rhof_2;
else
    rhof = rhof_m*cl_dim + rhof_c;
end

%ol density
rhool =  rhool_m*ol_mg_dim+rhool_c;

%opx density
rhoopx =  rhoopx_m*opx_mg_dim+rhoopx_c;

%cpx density 
rhocpx =  rhocpx_m*cpx_mg_dim+rhocpx_c;
    
if T>=Ts_local
    %feldspar density
    rhofeld =  rhofeld_m*T + rhofeld_c;
else
    rhofeld =  rhofeld_m*Ts_local + rhofeld_c;
end




%% Converting from wt% to vol%


phi_vol = phi_wt/rhof;
ol_vol = ol_wt/rhool;
opx_vol = opx_wt/rhoopx;
cpx_vol = cpx_wt/rhocpx;
feld_vol = feld_wt/rhofeld;

sum_vol = phi_vol + ol_vol + opx_vol + cpx_vol + feld_vol;

%% 
if sum_vol>0
    phi=phi_vol/sum_vol;
    ol=ol_vol/sum_vol;
    opx=opx_vol/sum_vol;
    cpx=cpx_vol/sum_vol;
    feld=feld_vol/sum_vol;
    sol = ol + opx + cpx + feld;
    test = ol+opx+cpx+feld+phi;
    %if abs(test-1)>1e-6
       % fprintf('%.15f\n',test)
       % error('oh no')
   % end
else
    phi=0;
    ol=0;
    opx=0;
    cpx=0;
    feld=0;
    sol=0;
    %sum_vol=1/3150;
end

%phi = phi_wt;
%ol = ol_wt;
%opx = opx_wt;
%cpx = cpx_wt;
%feld = feld_wt;
%sol = ol + opx + cpx + feld;


%% Calculate solid density

if sol>0
    rhom = ((ol/sol)*rhool + (opx/sol)*rhoopx + (cpx/sol)*rhocpx + (feld/sol)*rhofeld);
else
    rhom=rhom_1;
end


end