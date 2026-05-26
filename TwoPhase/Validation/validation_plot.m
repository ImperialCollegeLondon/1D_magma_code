clc;
clear;

Step = 0;


% -------------------------------------------------
% Load initial data
% -------------------------------------------------
% Cases={'Non_2','N_mu2_13_2'};
% Labels={'Nonlinear','Newton'};
% Range = [0 5];

% Cases={'Non_linear_muM_13_20y_-4','Non_linear_muM_13_20y','Newton_muM_13_20y'};
% Labels={'Nonlinear-4','Nonlinear-6','Newton'};
% Range = [0 19];

% Cases={'Non_linear_muM_13_eu-4','Newton_muM_13_eu'};
% Labels={'Nonlinear-6','Newton'};
% Range = [0 70];
% dy=5;

Cases={'Non_linear_muC_eu','Non_linear_muC_eu_fine','Non_linear_muC_eu_fine2','Newton_muC_eu_stable0','Newton_muC_eu_fine'};
Labels={'Nonlinear-6','Nonlinear-8','Nonlinear-9-Cu0.1','Newton','Newton dt0.2'};
Range = [0 68];
dy=5;

y_range=[48 60];
x_range=[-15.15 -14.95];

data=cell(1,length(Cases));
for i=1:length(Cases)
    data{i} = readtable(['/media/hh210/data/Matlab_Workplace/1D_magma_code/TwoPhase/Validation/' Cases{i} '/output_' num2str(Step) '_CELLS.txt']);
end

Font = 20;
Linewidth=1.5;
% -------------------------------------------------
% Main Figure
% -------------------------------------------------
figure(11)
clf
set(gcf,'Color','w')

Windows=[1,2];  %1, melt fraction 2. SiO2 3.Temperature

Handle=zeros(1,length(Cases)*length(Windows)+1);
% ---------------- TOP PLOT ----------------
subplot(length(Windows),1,1)
hold on
box on
set(gca,'TickDir','out','FontSize',Font)

for i=1:length(Cases)
    Handle(i) = plot(data{i}.Depth_km_, data{i}.MeltFrac_, ...
        'linewidth',Linewidth);
end

Handle(end)=title(['Time=' num2str(Step*50) 'y'],'fontsize',Font);
ylabel('Melt fraction (-)','fontsize',Font)

legend(Labels,'fontsize',16)

xlim(x_range)
ylim([0, 1]);
% ---------------- BOTTOM PLOT ----------------
subplot(length(Windows),1,2)
hold on
box on
set(gca,'TickDir','out','FontSize',Font)

for i=1:length(Cases)
    Handle(i+length(Cases)) = plot(data{i}.Depth_km_, data{i}.CbSiO2, ...
        'linewidth',Linewidth);
end



ylabel('Bulk SiO_2 (%)','fontsize',Font)
ylim(y_range);
xlim(x_range)

if length(Windows)==3
subplot(length(Windows),1,3)
hold on
box on
set(gca,'TickDir','out','FontSize',Font)

for i=1:length(Cases)
    Handle(i+length(Cases)*2) = plot(data{i}.Depth_km_, data{i}.Temp__C_, ...
        'linewidth',Linewidth);
end
ylim([300 1360]);
xlim(x_range)
ylabel('temperature (^\circ C)','fontsize',Font)
end
xlabel('Depth (km)','fontsize',Font)
% -------------------------------------------------
% Slider UI
% -------------------------------------------------
uicontrol('Style','text', ...
          'Position',[150 80 120 20], ...
          'String','Select Step', ...
          'FontSize',12);

slider = uicontrol('Style','slider', ...
                   'Min',Range(1), ...
                   'Max',Range(2), ...
                   'Value',Step, ...
                   'SliderStep',[1/(Range(2)-Range(1)) 1/(Range(2)-Range(1))], ...
                   'Position',[50 50 300 20], ...
                   'Callback',@sliderCallback);

figLabel = uicontrol('Style','text', ...
                     'Position',[170 20 100 20], ...
                     'String',['Step = ' num2str(Step)], ...
                     'FontSize',12);

% -------------------------------------------------
% Store handles
% -------------------------------------------------
S.Handle   = Handle;
S.figLabel = figLabel;
S.Cases=Cases;
S.dy=dy;
S.Windows=Windows;

guidata(gcf,S);


% =================================================
% Callback Function
% =================================================
function sliderCallback(src,~)

    % Get stored data
    S = guidata(gcf);

    Handle   = S.Handle;
    figLabel = S.figLabel;
    Cases=S.Cases;
    dy=S.dy;
    Windows=S.Windows;

    % Current slider value
    Step = round(get(src,'Value'));

    % Snap slider to integer
    set(src,'Value',Step);

    % Update text label
    set(figLabel,'String',['Step = ' num2str(Step)])

    % -------------------------------------------------
    % Reload data
    % -------------------------------------------------
    data=cell(1,length(Cases));
    for i=1:length(Cases)
        data{i} = readtable(['/media/hh210/data/Matlab_Workplace/1D_magma_code/TwoPhase/Validation/' S.Cases{i} '/output_' num2str(Step) '_CELLS.txt']);
    end

    % -------------------------------------------------
    % Update plots
    % -------------------------------------------------
    for i=1:length(Cases)
    set(Handle(i), ...
        'XData',data{i}.Depth_km_, ...
        'YData',data{i}.MeltFrac_);
    end
    
    for i=1:length(Cases)
    set(Handle(i+length(Cases)), ...
        'XData',data{i}.Depth_km_, ...
        'YData',data{i}.CbSiO2);
    end
    
    if length(Windows)==3
        for i=1:length(Cases)
            set(Handle(i+length(Cases)*2), ...
                'XData',data{i}.Depth_km_, ...
                'YData',data{i}.Temp__C_);
        end
    end
    set(Handle(end),'string',['Time=' num2str(Step*dy) 'y'])
    drawnow;

end