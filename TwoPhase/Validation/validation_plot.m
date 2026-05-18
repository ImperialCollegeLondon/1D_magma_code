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

Cases={'Non_linear_muM_13_-6','Newton_muM_13'};
Labels={'Nonlinear-6','Newton'};
Range = [0 49];

data=cell(1,length(Cases));
for i=1:length(Cases)
    data{i} = readtable(['/media/hh210/data/Matlab_Workplace/1D_magma_code/TwoPhase/Validation/' Cases{i} '/output_' num2str(Step) '_CELLS.txt']);
end

Font = 20;

% -------------------------------------------------
% Main Figure
% -------------------------------------------------
figure(11)
clf
set(gcf,'Color','w')


Handle=zeros(1,length(Cases)*2+1);
% ---------------- TOP PLOT ----------------
subplot(2,1,1)
hold on
box on
set(gca,'TickDir','out','FontSize',Font)

for i=1:length(Cases)
    Handle(i) = plot(data{i}.Depth_km_, data{i}.MeltFrac_, ...
        'linewidth',3);
end

Handle(end)=title(['Time=' num2str(Step*50) 'y'],'fontsize',Font);
xlabel('Depth (km)','fontsize',Font)
ylabel('Melt fraction (-)','fontsize',Font)

legend(Labels,'fontsize',16)

xlim([-15.15 -14.95])

% ---------------- BOTTOM PLOT ----------------
subplot(2,1,2)
hold on
box on
set(gca,'TickDir','out','FontSize',Font)

for i=1:length(Cases)
    Handle(i+length(Cases)) = plot(data{i}.Depth_km_, data{i}.CbSiO2, ...
        'linewidth',3);
end


xlabel('Depth (km)','fontsize',Font)
ylabel('Bulk SiO_2 (%)','fontsize',Font)

legend({'Nonlinear','Newton'},'fontsize',16)

xlim([-15.2 -14.9])

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

    set(Handle(end),'string',['Time=' num2str(Step*50) 'y'])
    drawnow;

end