Papers={'Blatter', 'Blatter2','Alonso','Melekhova','Nandedkar','Sisson'};

% Folder='/media/hh210/data/Matlab_Workplace/CV_1D_twophase/1D_threephase/Draw/Melting_data/';
% Folder='D:\1D_magma_code\1D_threephase\Draw\Melting_data\';
Folder='D:\Workspace\CV_1D_twophase\1D_threephase\Draw\Melting_data\';


Index=0;
results = struct();
for i=1:length(Papers)
    % File names
    modesFile  = [Folder Papers{i} '_modes.csv'];    % contains modal fractions per run
    oxidesFile = [Folder Papers{i} '_oxides.csv'];   % contains oxide compositions per phase

    % Read data
    modes  = readtable(modesFile);
    oxides = readtable(oxidesFile);

    % List of oxides to compute
    oxideList = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O','H2O'};

    % Initialize result structure
    

    % Loop over each Run
    uniqueRuns = unique(modes.Run);

    for r = 1:length(uniqueRuns)
        runID = uniqueRuns(r);
        % --- Extract mode and composition data for this run
        if iscell(modes.Run)
            runModes  = modes(strcmp(modes.Run, runID),:);
            runOxides = oxides(strcmp(oxides.Run, runID),:);
        else
            runModes  = modes(modes.Run == runID,:);
            runOxides = oxides(oxides.Run == runID,:);
        end
        % Extract numeric mode fractions (all columns after 3rd)
        modalCols = runModes(:, 4:end);
        modalFracs = table2array(modalCols);
        modalFracs(isnan(modalFracs)) = 0;
        % Normalize mode proportions to 1.0
        modalTotal = sum(modalFracs, 'all', 'omitnan');
        modalFracs = modalFracs / modalTotal;
        % Melt (glass) composition
        meltIdx = strcmpi(runOxides.Phase, 'Glass');
        meltComp = runOxides{meltIdx, oxideList};
        % Initialize solid composition accumulator
        solidComp = zeros(1, numel(oxideList));
        totalSolidWt = 0;
        % Loop over crystalline phases
        for p = 1:height(runOxides)
            phase = runOxides.Phase{p};
            % Skip the glass phase
            if strcmpi(phase, 'Glass')
                continue;
            end
            % Skip if phase not present in modal list
            if ~ismember(phase, runModes.Properties.VariableNames)
                continue;
            end
            modalVal = runModes{1, phase};
            if isempty(modalVal) || modalVal <= 0
                continue;
            end
            % Add weighted contribution
            thisComp = runOxides{p, oxideList};
            solidComp = solidComp + modalVal * thisComp;
            totalSolidWt = totalSolidWt + modalVal;
        end
        % Normalize solid composition
        solidComp = solidComp / totalSolidWt;
        % --- Store results
        results(Index+r).Run   = string(runID);
        results(Index+r).Melt  = array2table(meltComp, 'VariableNames', oxideList);
        results(Index+r).Solid = array2table(solidComp, 'VariableNames', oxideList);
    end
    Index=Index+length(uniqueRuns);
    % --- Display results summary
    % for r = 1:length(results)
    %     fprintf('Run %d:\n', results(r).Run);
    %     disp('  Melt composition (wt%%):');
    %     disp(results(r).Melt);
    %     disp('  Solid composition (wt%%):');
    %     disp(results(r).Solid);
    % end

end


%%
SiO2_melt=arrayfun(@(r) r.Melt.SiO2, results);
TiO2_melt=arrayfun(@(r) r.Melt.TiO2, results);
Al2O3_melt=arrayfun(@(r) r.Melt.Al2O3, results);
FeO_melt=arrayfun(@(r) r.Melt.FeO , results);
MnO_melt=arrayfun(@(r) r.Melt.MnO, results);
MgO_melt=arrayfun(@(r) r.Melt.MgO, results);
CaO_melt=arrayfun(@(r) r.Melt.CaO, results);
Na2O_melt=arrayfun(@(r) r.Melt.Na2O, results);
K2O_melt=arrayfun(@(r) r.Melt.K2O, results);


[SiO2_melt, idx] = sort(SiO2_melt);
TiO2_melt=TiO2_melt(idx);
Al2O3_melt=Al2O3_melt(idx);
FeO_melt=FeO_melt(idx);
MnO_melt=MnO_melt(idx);
MgO_melt=MgO_melt(idx);
CaO_melt=CaO_melt(idx);
Na2O_melt=Na2O_melt(idx);
K2O_melt=K2O_melt(idx);
%%
Fitting_order=3;
Oxides_coefficients=zeros(8,Fitting_order+1);


Font=15;
figure(1)
clf
set(gcf,'Units','normalized')
set(gcf,'color','w');
% set(gcf,'FontSize',20)

% TiO2 vs SiO2
subplot(3,3,1)
plot(SiO2_melt, TiO2_melt, 'o'); hold on
p = polyfit(SiO2_melt, TiO2_melt, Fitting_order);
Oxides_coefficients(1,:)=p;
xfit = linspace(min(SiO2_melt), max(SiO2_melt), 100);
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
set(gca,'FontSize',Font)
xlabel('SiO2 (wt%)'); ylabel('TiO2 (wt%)')
% title('TiO2 vs SiO2')

% Al2O3
subplot(3,3,2)
plot(SiO2_melt, Al2O3_melt, 'o'); hold on
p = polyfit(SiO2_melt, Al2O3_melt, Fitting_order);
Oxides_coefficients(2,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('Al2O3 (wt%)')
% title('Al2O3 vs SiO2')
set(gca,'FontSize',Font)

% FeO
subplot(3,3,3)
plot(SiO2_melt, FeO_melt, 'o'); hold on
p = polyfit(SiO2_melt, FeO_melt, Fitting_order);
Oxides_coefficients(3,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('FeO (wt%)')
% title('FeO vs SiO2')
set(gca,'FontSize',Font)

% MnO
subplot(3,3,4)
plot(SiO2_melt, MnO_melt, 'o'); hold on
p = polyfit(SiO2_melt, MnO_melt, Fitting_order);
Oxides_coefficients(4,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('MnO (wt%)')
% title('MnO vs SiO2')
set(gca,'FontSize',Font)

% MgO
subplot(3,3,5)
plot(SiO2_melt, MgO_melt, 'o'); hold on
p = polyfit(SiO2_melt, MgO_melt, Fitting_order);
Oxides_coefficients(5,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('MgO (wt%)')
% title('MgO vs SiO2')
set(gca,'FontSize',Font)

% CaO
subplot(3,3,6)
plot(SiO2_melt, CaO_melt, 'o'); hold on
p = polyfit(SiO2_melt, CaO_melt, Fitting_order);
Oxides_coefficients(6,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('CaO (wt%)')
% title('CaO vs SiO2')
set(gca,'FontSize',Font)

% Na2O
subplot(3,3,7)
plot(SiO2_melt, Na2O_melt, 'o'); hold on
p = polyfit(SiO2_melt, Na2O_melt, Fitting_order);
Oxides_coefficients(7,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('Na2O (wt%)')
% title('Na2O vs SiO2')
set(gca,'FontSize',Font)

% K2O
subplot(3,3,8)
plot(SiO2_melt, K2O_melt, 'o'); hold on
p = polyfit(SiO2_melt, K2O_melt, Fitting_order);
Oxides_coefficients(8,:)=p;
yfit = polyval(p, xfit);
plot(xfit, yfit, 'r-', 'LineWidth', 1.5); hold off
xlabel('SiO2 (wt%)'); ylabel('K2O (wt%)')
% title('K2O vs SiO2')
set(gca,'FontSize',Font)

% Empty panel
subplot(3,3,9)
axis off
% text(0.1,0.5,[num2str(Fitting_order) ' order least-squares fits'],'FontSize',12)
save('Oxides_coefficients.mat','Oxides_coefficients');
