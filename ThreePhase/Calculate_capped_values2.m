Mass_all=(sum(Mass_data,2));
Mass_solid=(sum(Mass_data(:,4:6),2));
% 
melt_fraction =sum(Mass_data(:,1:3),2)./Mass_all;
solid_fraction=(sum(Mass_data(:,4:6),2)-max(Mass_data(:,6)-Mass_solid.*S_cap,0))./Mass_all;
% indexS=[];


index=find(melt_fraction>Convective_cut);
if any(index)
    % indexL=find(index(1:end-1) & ~index(2:end)) + 1;
    % if index(end)
    %     indexL = [indexL; index(end)];
    % end

    % indexS=unique(index+(1:8));
    % indexS=unique(index+(-5:5));
    indexS=index;

    melt_fraction_capped=melt_fraction;
    solid_fraction_capped=solid_fraction;

    Mass_data_capped=Mass_data;
    for i=1:length(indexS)
        if melt_fraction_capped(indexS(i))>=Convective_cut
            melt_fraction_capped(indexS(i))=Convective_cut;
            solid_fraction_capped(indexS(i))=solid_fraction_capped(indexS(i))+melt_fraction(indexS(i))-Convective_cut;
            Ml=sum(Mass_data(:,1:3),2);
            fraction=Mass_data(:,1:3)./max(Ml,1e-5);
            transfer=(melt_fraction(indexS(i))-Convective_cut).*Ml(indexS(i)).*fraction(indexS(i),:);
            Mass_data_capped(indexS(i),1:3)=Mass_data_capped(indexS(i),1:3)-transfer;
            Mass_data_capped(indexS(i),4:6)=Mass_data_capped(indexS(i),4:6)+transfer;
        end
        % if  solid_fraction_capped(indexS(i))>=Convective_cut2
        %     solid_fraction_capped(indexS(i))=Convective_cut2;
        %     melt_fraction_capped(indexS(i))=melt_fraction_capped(indexS(i))+solid_fraction(indexS(i))-Convective_cut2;
        %     Ms=sum(Mass_data(:,4:6),2);
        %     fraction=Mass_data(:,4:6)./max(Ms,1e-5);
        %     transfer=Ms(indexS(i)).*(solid_fraction(indexS(i))-Convective_cut2).*fraction(indexS(i),:);
        %     Mass_data_capped(indexS(i),4:6)=Mass_data_capped(indexS(i),4:6)-transfer;
        %     Mass_data_capped(indexS(i),1:3)=Mass_data_capped(indexS(i),1:3)+transfer;
        % end
    end
else
    Mass_data_capped=Mass_data;
    melt_fraction_capped=melt_fraction;
    solid_fraction_capped=solid_fraction;
    indexS=[];
end


% Mass_data_capped=Mass_data;
% melt_fraction_capped=melt_fraction;
% solid_fraction_capped=solid_fraction;