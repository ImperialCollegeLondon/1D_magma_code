%% HH 1/12/2022
% calculate the source term for the injected sill 
function [U_source,U_source_int,H_source,C_source,Index]=Cal_sources_master(injection_range,injection_rate,injection_t,cp,Lf,phi0,Cl0,Cs0,sink_location,nodez)
% This is to obtain all the source terms including : 
% 1. velocity source term and integral of velocity source
% over the distance. The integral represents the bulk velocity!
% 2. enthalpy and 
% 3. composition sources
% velocity source is on nodal position while enthalpy and composition
% source on cell position

    U_source=zeros(length(nodez),1);
    H_source=zeros(length(nodez)-1,1);
    C_source=zeros(length(nodez)-1,4);

    Index=zeros(size(injection_range,1)+1,2);
    S_amount=0;
    H_amount=0;
    C_amount=zeros(4,1);
    %building the source
    for i=1:size(injection_range,1)
        Index(i,1)=find(nodez<injection_range(i,1),1,'last');
        Index(i,2)=find(nodez>injection_range(i,2),1,'first');

        U_source(Index(i,1):Index(i,2))=injection_rate(i);
        H_source(Index(i,1):Index(i,2)-1)=injection_rate(i)*(cp*injection_t(i)+Lf*phi0(i));
        C_source(Index(i,1):Index(i,2)-1,1)=injection_rate(i)*phi0(i)*Cl0(i);
        C_source(Index(i,1):Index(i,2)-1,2)=injection_rate(i)*(1-phi0(i))*Cs0(i);
        C_source(Index(i,1):Index(i,2)-1,3)=injection_rate(i)*phi0(i)*(1-Cl0(i));
        C_source(Index(i,1):Index(i,2)-1,4)=injection_rate(i)*(1-phi0(i))*(1-Cl0(i));

        S_amount=S_amount+injection_rate(i)*(nodez(Index(i,2))-nodez(Index(i,1)));
        H_amount=H_amount+H_source(Index(i,1))*(nodez(Index(i,2))-nodez(Index(i,1)));
        C_amount(1)=C_amount(1)+C_source(Index(i,1),1)*(nodez(Index(i,2))-nodez(Index(i,1)));
        C_amount(2)=C_amount(2)+C_source(Index(i,1),2)*(nodez(Index(i,2))-nodez(Index(i,1)));
        C_amount(3)=C_amount(3)+C_source(Index(i,1),3)*(nodez(Index(i,2))-nodez(Index(i,1)));
        C_amount(4)=C_amount(4)+C_source(Index(i,1),4)*(nodez(Index(i,2))-nodez(Index(i,1)));
    end


   
    %building the sink
    Index(end,1)=find(nodez>sink_location(1),1,'first');
    Index(end,2)=find(nodez<sink_location(2),1,'last');
    U_source(Index(end,1):Index(end,2))=-S_amount/(nodez(Index(end,2))-nodez(Index(end,1)));
    H_source(Index(end,1):Index(end,2)-1)=-H_amount/(nodez(Index(end,2))-nodez(Index(end,1)));
    C_source(Index(end,1):Index(end,2)-1,1)=-C_amount(1)/(nodez(Index(end,2))-nodez(Index(end,1)));
    C_source(Index(end,1):Index(end,2)-1,2)=-C_amount(2)/(nodez(Index(end,2))-nodez(Index(end,1)));
    C_source(Index(end,1):Index(end,2)-1,3)=-C_amount(3)/(nodez(Index(end,2))-nodez(Index(end,1)));
    C_source(Index(end,1):Index(end,2)-1,4)=-C_amount(4)/(nodez(Index(end,2))-nodez(Index(end,1)));

    U_source_int=zeros(length(nodez),1);
    for i=2:length(nodez)
        U_source_int(i)=U_source_int(i-1)+(nodez(i)-nodez(i-1))*(U_source(i)+U_source(i-1))/2*(abs(U_source(i)*U_source(i-1))>0);
    end
end