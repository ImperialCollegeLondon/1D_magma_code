%% ---- USER PARAMETERS ---------------------------------------------------
close all
% Depth axis: total depth range shown (km, positive = deeper)
depth_total_km = 50;      % total depth range to show (km)
depth_surface  = 0;       % surface depth label (km)
width=20; % half width of the figure in km
% ---- Reservoirs ----------------------------------------------------------
% Each row: [depth_center_km, half_width_km, half_height_km, melt_fraction]
% melt_fraction: 0 (solid) → 1 (pure melt); controls particle density
reservoirs = [
    15,  8, 3.5, 0.45;   % shallow reservoir  (depth, half-w, half-h, melt)
    35,  10, 2.5, 0.25;   % deep reservoir
];

% ---- Dikes (vertical conduits) ------------------------------------------
% Each row: [x_position , width_km, start,  end]
% x_position is fraction of figure width; 0 = center
dikes = [
   -2,  3, 45, 35;   % leftmost dike
    5,  2, 35, 15;   % central dike
];

% ---- Mantle layer --------------------------------------------------------
mantle_top_km    = 45;    % depth where red mantle begins (km)
mantle_color     = [0.85, 0.15, 0.10];

% ---- Background rock texture particles -----------------------------------
% These white/cream capsule-shapes represent crystalline matrix
% More particles = less melt (lower melt fraction in that region)
bg_particle_density = 180;   % particles per unit area in background
bg_melt_fraction    = 0.0;   % background is solid rock (no melt)

% ---- Colors --------------------------------------------------------------
col_bg_rock      = [0.96, 0.90, 0.75];   % sandy background
col_dike         = [0.90, 0.30, 0.10];   % orange-red dike
col_reservoir    = [0.95, 0.60, 0.45];   % light orange reservoir
col_particle_bg  = [1.00, 1.00, 1.00];   % white background particles
col_particle_res = [0.85, 0.45, 0.20];   % dark particles inside reservoir
col_arrow        = [0.85, 0.20, 0.08];   % upward flow arrows
col_scale        = [0.20, 0.20, 0.20];   % depth scale

% =========================================================================
%% ---- SETUP FIGURE ------------------------------------------------------
fig = figure('Color','white','Position',[100 50 900 900]);
ax  = axes('Position',[0.04 0.04 0.84 0.93]);
hold(ax,'on'); axis(ax,'off');

% Working coordinate system: x in [-1, 1], y in [0, depth_total_km]
% y increases downward (depth)
xlim(ax,[-width-1, width]);
ylim(ax,[0, depth_total_km]);
set(ax,'YDir','reverse');   % 0 at top, depth_total at bottom

% =========================================================================
%% ---- BACKGROUND ROCK ---------------------------------------------------
fill(ax,[-width width width -width],[0 0 depth_total_km depth_total_km], col_bg_rock, ...
     'EdgeColor','none');

% =========================================================================
%% ---- MANTLE LAYER (solid red bottom) -----------------------------------
fill(ax,[-width width width -width], ...
     [mantle_top_km mantle_top_km depth_total_km depth_total_km], ...
     mantle_color,'EdgeColor','none','FaceAlpha',0.92);

% Gradient effect on mantle top edge
for i = 1:20
    alpha_val = 0.05 + 0.04*i;
    y_top = mantle_top_km - 0.5 + 0.5*(i/20);
    fill(ax,[-width width width -width],[y_top y_top y_top+0.4 y_top+0.4], ...
         mantle_color,'EdgeColor','none','FaceAlpha',alpha_val);
end

% =========================================================================
%% ---- BACKGROUND TEXTURE PARTICLES (solid rock zones) ------------------
draw_particles(ax, -1, 1, 0, mantle_top_km, ...
               bg_particle_density, bg_melt_fraction, col_particle_bg, ...
               reservoirs, dikes, depth_total_km);

% =========================================================================
%% ---- DIKES (vertical conduits) -----------------------------------------
for d = 1:size(dikes,1)
    xc  = dikes(d,1);
    hw  = dikes(d,2);  % scale half-width to plot units
    % Dike goes full depth
    y_top    = dikes(d,4);
    y_bottom = dikes(d,3);
    
    
    fill(ax,[xc-hw/2, xc+hw/2, xc+hw/2, xc-hw/2, xc-hw/2 ],[y_bottom, y_bottom, y_top, y_top, y_bottom], col_dike, ...
         'EdgeColor','none','FaceAlpha',0.88);
    
    % Highlight / gradient inside dike
    % yy2 = linspace(y_top, y_bottom, 40);
    % xl2 = xc - hw*0.3 * ones(size(yy2));
    % xr2 = xc + hw*0.3 * ones(size(yy2));
    % fill(ax,[xl2, fliplr(xr2)],[yy2, fliplr(yy2)], [1 0.75 0.55], ...
    %      'EdgeColor','none','FaceAlpha',0.35);
end

% =========================================================================
%% ---- RESERVOIRS --------------------------------------------------------
for r = 1:size(reservoirs,1)
    dep_c  = reservoirs(r,1);
    hw_km  = reservoirs(r,2);
    hh_km  = reservoirs(r,3);
    melt_f = reservoirs(r,4);
    
    % Convert km to plot units
    hw_x = hw_km;   % x half-width in plot units
    hh_y = hh_km;               % y half-height in km
    
    % Draw ellipse
    theta = linspace(0, 2*pi, 200);
    xe = hw_x * cos(theta);
    ye = dep_c + hh_y * sin(theta);
    
    % Clip to depth range
    valid = ye >= 0 & ye <= depth_total_km;
    
    % Fill reservoir
    fill(ax, xe(valid), ye(valid), col_reservoir, ...
         'EdgeColor','none','FaceAlpha',0.85);
    
    % Inner gradient (lighter center)
    xe2 = hw_x*0.55 * cos(theta);
    ye2 = dep_c + hh_y*0.55 * sin(theta);
    fill(ax, xe2, ye2, [1 0.88 0.75], ...
         'EdgeColor','none','FaceAlpha',0.50);
    
    % Particles inside reservoir (fewer = more melt)
    draw_particles_ellipse(ax, 0, dep_c, hw_km, hh_km, melt_f);
end

%% ---- DEPTH SCALE BAR ---------------------------------------------------
x_scale = -width;
tick_x2 = -width-1;

% Main vertical line
line(ax,[x_scale x_scale],[0 depth_total_km],'Color',col_scale,'LineWidth',1.2);

% Depth tick marks and labels
tick_depths = 0:5:depth_total_km;
for td = tick_depths
    is_major = mod(td,10)==0;
    tx2 = tick_x2 + (is_major - 0.5)*0.04;
    line(ax,[x_scale tx2],[td td],'Color',col_scale, ...
         'LineWidth', 0.8 + is_major*0.6);
    if is_major
        text(ax, x_scale - 1.04, td, sprintf('%dkm', td), ...
             'HorizontalAlignment','right','FontSize',8, ...
             'Color',col_scale,'FontName','Helvetica');
    end
end

% =========================================================================
%% ---- TITLE & LABELS ----------------------------------------------------
text(ax, 0, -2.5, 'Magma Plumbing System', ...
     'HorizontalAlignment','center','FontSize',14,'FontWeight','bold', ...
     'Color',[0.3 0.1 0.05],'FontName','Helvetica');

% Reservoir labels
% for r = 1:size(reservoirs,1)
%     text(ax, 0.55, reservoirs(r,1), ...
%          sprintf('Reservoir %d\n(melt = %.0f%%)', r, reservoirs(r,4)*100), ...
%          'FontSize',7,'Color',[0.5 0.15 0.05],'FontName','Helvetica', ...
%          'HorizontalAlignment','left');
% end

text(ax, -0.3, (mantle_top_km + depth_total_km)/2, 'Mantle', ...
     'FontSize',9,'Color',[1 0.8 0.75],'FontWeight','bold', ...
     'FontName','Helvetica');

% =========================================================================
%% ---- FINALIZE ----------------------------------------------------------
axis(ax,'equal','off');
% xlim(ax,[-1.05, 1.05]);
ylim(ax,[-4, depth_total_km + 2]);
set(ax,'YDir','reverse');

% Save
exportgraphics(fig,'magma_plumbing.png','Resolution',200);
fprintf('Figure saved as magma_plumbing.png\n');



% =========================================================================
%% ---- HELPER: Draw background particles (capsule-shapes) ----------------
function draw_particles(ax, xmin, xmax, ymin, ymax, density, melt_f, ...
                        color, reservoirs, dikes, depth_total_km)
    % Particle count scales inversely with melt fraction
    n = round(density * (1 - melt_f));
    
    rng(42);  % reproducible
    for i = 1:n
        px = xmin + rand()*(xmax - xmin);
        py = ymin + rand()*(ymax - ymin);
        
        % Skip if inside any reservoir
        in_res = false;
        for r = 1:size(reservoirs,1)
            dep_c = reservoirs(r,1);
            hw_x  = reservoirs(r,2);
            hh_y  = reservoirs(r,3);
            if ((px/hw_x)^2 + ((py-dep_c)/hh_y)^2) < 1.1
                in_res = true; break;
            end
        end
        
        % Skip if inside any dike
        in_dike = false;
        for d = 1:size(dikes,1)
            xc = dikes(d,1);
            hw = dikes(d,2);
            if abs(px - xc) < hw * 1.1
                in_dike = true; break;
            end
        end
        
        % Skip if in mantle (bottom red area)
        mantle_top = reservoirs(end,1) + 10;  % approx
        
        if ~in_res && ~in_dike
            draw_capsule(ax, px, py, color, 0.65);
        end
    end
end


% =========================================================================
%% ---- HELPER: Draw particles inside an elliptical reservoir -------------
function draw_particles_ellipse(ax, xc, yc, hw_x, hh_y, melt_f)
    % More melt → fewer, darker particles
    n = round(120 * (1 - melt_f * 0.7));
    col_dark  = [0.65, 0.25, 0.12];   % dark mineral grains
    col_green = [0.30, 0.55, 0.20];   % green xenoliths
    col_white = [0.98, 0.97, 0.93];   % feldspar laths
    
    rng(99);
    for i = 1:n
        % Random point inside ellipse
        ang = rand()*2*pi;
        rad = sqrt(rand());
        px  = xc + hw_x * 0.90 * rad * cos(ang);
        py  = yc + hh_y * 0.85 * rad * sin(ang);
        
        % Random mineral type
        r_type = rand();
        if r_type < 0.55
            col = col_white;
            alpha = 0.80;
        elseif r_type < 0.80
            col = col_dark + rand()*[0.1 0.05 0.05];
            alpha = 0.85;
        else
            col = col_green;
            alpha = 0.75;
        end
        
        draw_capsule(ax, px, py, col, alpha);
    end
end


% =========================================================================
%% ---- HELPER: Draw a single oriented capsule (rounded rectangle) --------
function draw_capsule(ax, cx, cy, color, alpha_val)
    % Random orientation and size
    angle  = rand() * pi;
    length = 0.5 + rand() * 0.1;   % plot units
    width  = length * (0.25 + rand()*0.20);
    
    % Capsule approximated as rounded rectangle
    t = linspace(0, 2*pi, 24);
    % Ellipse vertices
    xe = (length/2) * cos(t);
    ye = (width/2)  * sin(t);
    
    % Rotate
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    pts = rot * [xe; ye];
    
    fill(ax, cx + pts(1,:), cy + pts(2,:), color, ...
         'EdgeColor','none','FaceAlpha',alpha_val);
end