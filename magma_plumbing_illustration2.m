%% ---- USER PARAMETERS ---------------------------------------------------
close all

n = 256;

% Build segments
n1 = round(0.2*n);  % dark red -> red
n2 = round(0.3*n);  % red -> orange
n3 = round(0.3*n);  % orange -> yellow
n4 = n - (n1+n2+n3); % yellow -> light yellow

% Segment 1: dark red -> red
r1 = linspace(0.5,1,n1)';
g1 = zeros(n1,1);
b1 = zeros(n1,1);

% Segment 2: red -> orange
r2 = ones(n2,1);
g2 = linspace(0,0.5,n2)';
b2 = zeros(n2,1);

% Segment 3: orange -> yellow
r3 = ones(n3,1);
g3 = linspace(0.5,1,n3)';
b3 = zeros(n3,1);

% Segment 4: yellow -> light yellow (add a bit of blue to soften)
r4 = ones(n4,1);
g4 = ones(n4,1);
b4 = linspace(0,0.4,n4)';

% Combine
cmap = [r1 g1 b1;
        r2 g2 b2;
        r3 g3 b3;
        r4 g4 b4];



% Depth axis: total depth range shown (km, positive = deeper)
depth_total_km = 50;      % total depth range to show (km)
depth_surface  = 0;       % surface depth label (km)
width=50; % 70 half width of the figure in km

% ---- Reservoirs ----------------------------------------------------------
% Each row: [depth_center_km, x_position, half_width_km, half_height_km, cb_min, cb_ax]
% melt_fraction: 0 (solid) → 1 (pure melt); controls particle density
reservoirs = [
    5,  -30,  5,  2,  0.0, 1;
    15,   0,  6, 1.5, 0.5, 0.8;   % shallow reservoir  (depth, half-w, half-h, melt)
    30,   -5, 10, 2.5, 0.0, 0.6;   % deep reservoir
    28,   28, 12,   4, 0.0,  0.7;   % deep reservoir
];

% reservoirs = [
%     5,  1,  5,  2,     0.0, 1;
%     20,   0,  6, 1.5,  0.2, 0.8;   
%     30,   -2, 7, 2, 0.0, 0.6;   
% ];

% reservoirs = [
%     20,   0,  6, 2.5,  0.5, 1;   
%     30,   -3, 4, 1.5, 0.0, 0.4;   
% ];

% reservoirs = [
%     5,  -0,  5,  2,  0.0, 1;
% ];
% ---- Dikes (vertical conduits) ------------------------------------------
% Each row: [x_position , width_km, start,  end, colortype]
% x_position is fraction of figure width; 0 = center
dikes = [
  -30,  1,   45,  6, 0;
   -5,  1.5, 45, 32, 0;   % leftmost dike
    0,  1.5, 31.5, 16, 0;   % central dike
   28,  2,   45, 31.5, 0;
];

% dikes = [
%     5,  1,   20,  5, 1;
%     0,  1.5,   45, 31.5, 0;
%     0,  0.01, 28.5, 21, 0;
% ];

% dikes = [
%     -3,  1,   45, 31.5, 0;
%     -1.5,  0.01, 29.5, 21, 0;
% ];
% 
% dikes = [
%   0,  1,   45,  6, 0;
% ];

% ---- Envelope parameters -------------------------------------------------
% Envelope = aureole/halo of partially-molten country rock around features
envelope_res_scale  = 1.8;   % how much larger than reservoir (multiplicative)
envelope_dike_pad   = 1.5;    % extra km padding on each side of dike (km)
envelope_color      = [0.98, 0.82, 0.62];  % warm cream — partially molten rock
envelope_melt_frac  = 0.9;   % envelope is partially molten → denser particles
envelope_particle_density = 150;  % particles per unit area in envelope

% ---- Mantle layer --------------------------------------------------------
mantle_top_km    = 45;    % depth where red mantle begins (km)
mantle_color     = [0.85, 0.15, 0.10];

% ---- Background rock texture particles -----------------------------------
bg_particle_density = 30;   % particles per unit area in background
bg_melt_fraction    = 0.0;   % background is solid rock (no melt)

% ---- Colors --------------------------------------------------------------
col_bg_rock      = [1, 1, 1];   % sandy background
col_dike         = [0.90, 0.30, 0.10];   % orange-red dike
col_dike2        = [1, 0.8, 0.0];   % orange-red dike
col_reservoir    = [0.95, 0.60, 0.45];   % light orange reservoir
col_particle_bg  = [1.00, 1.00, 1.00];   % white background particles
col_particle_res = [0.85, 0.45, 0.20];   % dark particles inside reservoir
col_arrow        = [0.85, 0.20, 0.08];   % upward flow arrows
col_scale        = [0.20, 0.20, 0.20];   % depth scale

% =========================================================================
%% ---- SETUP FIGURE ------------------------------------------------------
fig = figure('Color','white','Position',[100 50 1500 1500]);
ax  = axes('Position',[0.04 0.04 0.84 0.93]);
hold(ax,'on'); axis(ax,'off');
colormap(ax, cmap);
clim([0, 1])

xlim(ax,[-width-1, width]);
ylim(ax,[0, depth_total_km]);
set(ax,'YDir','reverse');   % 0 at top, depth_total at bottom

% =========================================================================
%% ---- BACKGROUND ROCK ---------------------------------------------------
fill(ax,[-width width width -width],[0 0 depth_total_km depth_total_km], col_bg_rock, ...
     'EdgeColor','none');

% =========================================================================
%% ---- MANTLE LAYER (solid red bottom) -----------------------------------
if depth_total_km>45
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
end
% =========================================================================
%% ---- BACKGROUND TEXTURE PARTICLES (solid rock zones) ------------------
% draw_particles(ax, -width, width, 0, mantle_top_km, ...
%                bg_particle_density, bg_melt_fraction, col_particle_bg, ...
%                reservoirs, dikes, depth_total_km, ...
%                envelope_res_scale, envelope_dike_pad);

% =========================================================================
%% ---- STEP 1: Draw dike envelopes FIRST (behind dikes) ------------------
for d = 1:size(dikes,1)
    xc       = dikes(d,1);
    hw       = dikes(d,2);
    y_bottom = dikes(d,3);
    y_top    = dikes(d,4);
    
    if hw<0.2
        scale=3.5;
    else
        scale=1;
    end
    env_hw   = hw + envelope_dike_pad*scale;   % padded half-width
    
    % Simple flat-ended rectangular envelope (no curved caps)
    fill(ax, [xc - env_hw,  xc + env_hw,  xc + env_hw,  xc - env_hw], ...
             [y_top,         y_top,         y_bottom,      y_bottom   ], ...
         envelope_color, 'EdgeColor','none','FaceAlpha',0.72);
    
    % Draw envelope particles
    draw_particles_dike_envelope(ax, xc, hw, y_top, y_bottom, ...
                                 env_hw, envelope_particle_density, ...
                                 envelope_melt_frac, col_particle_bg);
end

% =========================================================================
%% ---- STEP 2: Draw reservoir envelopes FIRST (behind reservoirs) --------
for r = 1:size(reservoirs,1)
    dep_c  = reservoirs(r,1);
    x_pos  = reservoirs(r,2);
    hw_km  = reservoirs(r,3);
    hh_km  = reservoirs(r,4);
    
    % Envelope ellipse — scaled up from reservoir
    env_hw = hw_km * envelope_res_scale;
    env_hh = hh_km * envelope_res_scale;
    
    theta = linspace(0, 2*pi, 200);
    xe = env_hw * cos(theta);
    ye = dep_c  + env_hh * sin(theta);
    
    fill(ax, x_pos+xe, ye, envelope_color, ...
         'EdgeColor','none','FaceAlpha',0.72);
    
    % Draw envelope particles (in ring between reservoir and envelope)
    draw_particles_res_envelope(ax, x_pos, dep_c, hw_km, hh_km, ...
                                env_hw, env_hh, envelope_particle_density, ...
                                envelope_melt_frac, col_particle_bg);
end

% =========================================================================
%% ---- STEP 3: Draw dikes ON TOP of envelopes ----------------------------
for d = 1:size(dikes,1)
    xc  = dikes(d,1);
    hw  = dikes(d,2);
    y_top    = dikes(d,4);
    y_bottom = dikes(d,3);
        
    if dikes(d,5)==0
        Color_code=col_dike;
    else
        Color_code=col_dike2;
    end
    if hw>0.2
        fill(ax,[xc-hw/2, xc+hw/2, xc+hw/2, xc-hw/2, xc-hw/2], ...
               [y_bottom, y_bottom, y_top,    y_top,    y_bottom], Color_code, ...
             'EdgeColor','none','FaceAlpha',0.88);
        
        draw_particles_dike(ax, xc, hw/2, y_top, y_bottom);
    end
end

% =========================================================================
%% ---- STEP 4: Draw reservoirs ON TOP of envelopes -----------------------
  

for r = 1:size(reservoirs,1)
    dep_c  = reservoirs(r,1);
    x_pos  = reservoirs(r,2);
    hw_km  = reservoirs(r,3);
    hh_km  = reservoirs(r,4);
    c_min  =reservoirs(r,5);
    c_max  =reservoirs(r,6);
    
    hw_x = hw_km;
    hh_y = hh_km;
    
    theta = linspace(0, 2*pi, 200);
    xe = hw_x * cos(theta);
    ye = dep_c + hh_y * sin(theta);
    
    valid = ye >= 0 & ye <= depth_total_km;
    
    xv = x_pos + xe(valid);
    yv = ye(valid);
    % Fill reservoir
    % fill(ax, xv, yv, col_reservoir, ...
    %      'EdgeColor','none','FaceAlpha',0.85);
    
    cv = ((-(yv-dep_c)/hh_y)/2+0.5)*(c_max-c_min)+c_min;   % color varies with vertical position

    patch(ax, xv, yv, cv, ...
          'FaceColor','interp', ...
          'EdgeColor','none', ...
          'FaceAlpha',0.85);
    

    
    % Optional: fix color scaling to your depth range
    % caxis(ax, [0 depth_total_km]);
        
    % Particles inside reservoir
    draw_particles_ellipse(ax, x_pos, dep_c, hw_km, hh_km, 0.99);
end
cb = colorbar(ax);

cb.Ticks = [0  0.5 1];              % positions (relative to caxis)
cb.TickLabels = {'45%','SiO_2 content','75%'};     % labels at ends
cb.FontSize=15;
%% ---- DEPTH SCALE BAR ---------------------------------------------------
x_scale = -width;
tick_x2 = -width-1;

line(ax,[x_scale x_scale],[0 depth_total_km],'Color',col_scale,'LineWidth',1.2);

tick_depths = 0:5:depth_total_km;
for td = tick_depths
    is_major = mod(td,10)==0;
    tx2 = tick_x2 + (is_major - 0.5)*0.04;
    line(ax,[x_scale tx2],[td td],'Color',col_scale, ...
         'LineWidth', 0.8 + is_major*0.6);
    if is_major
        text(ax, x_scale - 1.04, td, sprintf('%dkm', td), ...
             'HorizontalAlignment','right','FontSize',12, ...
             'Color',col_scale,'FontName','Helvetica');
    end
end

% =========================================================================
%% ---- TITLE & LABELS ----------------------------------------------------
% text(ax, 0, -2.5, 'Magma Plumbing System', ...
%      'HorizontalAlignment','center','FontSize',14,'FontWeight','bold', ...
%      'Color',[0.3 0.1 0.05],'FontName','Helvetica');

if depth_total_km>45
    text(ax, -0.3, (mantle_top_km + depth_total_km)/2, 'Mantle', ...
         'FontSize',12,'Color',[1 0.8 0.75],'FontWeight','bold', ...
         'FontName','Helvetica');
end
% =========================================================================
%% ---- FINALIZE ----------------------------------------------------------
axis(ax,'equal','off');
ylim(ax,[-4, depth_total_km + 2]);
set(ax,'YDir','reverse');

exportgraphics(fig,'magma_plumbing.png','Resolution',200);
fprintf('Figure saved as magma_plumbing.png\n');


% =========================================================================
%% ---- HELPER: Draw background particles ---------------------------------
% Skips regions inside envelopes (larger exclusion zones than before)
function draw_particles(ax, xmin, xmax, ymin, ymax, density, melt_f, ...
                        color, reservoirs, dikes, depth_total_km, ...
                        env_res_scale, env_dike_pad)
    n = round(density * (1 - melt_f));
    
    rng(42);
    for i = 1:n
        px = xmin + rand()*(xmax - xmin);
        py = ymin + rand()*(ymax - ymin);
        
        % Skip if inside any reservoir ENVELOPE (scaled ellipse)
        in_res = false;
        for r = 1:size(reservoirs,1)
            dep_c = reservoirs(r,1);
            hw_x  = reservoirs(r,2) * env_res_scale;
            hh_y  = reservoirs(r,3) * env_res_scale;
            if ((px/hw_x)^2 + ((py-dep_c)/hh_y)^2) < 1.05
                in_res = true; break;
            end
        end
        
        % Skip if inside any dike ENVELOPE (padded width)
        in_dike = false;
        for d = 1:size(dikes,1)
            xc      = dikes(d,1);
            env_hw  = dikes(d,2) + env_dike_pad;
            y_bot   = dikes(d,3);
            y_top   = dikes(d,4);
            if abs(px - xc) < env_hw && py >= y_top && py <= y_bot
                in_dike = true; break;
            end
        end
        
        if ~in_res && ~in_dike
            draw_capsule(ax, px, py, color, 0.65);
        end
    end
end


% =========================================================================
%% ---- HELPER: Particles in DIKE envelope (ring around dike) -------------
function draw_particles_dike_envelope(ax, xc, hw, y_top, y_bottom, ...
                                      env_hw, density, melt_f, color)
    % Sample in a box covering the envelope, reject if outside ring
    n = round(density * (1 - melt_f * 0.5) * (env_hw*2) * (y_bottom - y_top) / 100);
    
    rng(77);
    attempts = 0;
    drawn    = 0;
    while drawn < n && attempts < n * 20
        attempts = attempts + 1;
        px = (xc - env_hw) + rand() * (2*env_hw);
        py = y_top  + rand() * (y_bottom - y_top);
        
        % Must be inside envelope band but OUTSIDE the dike itself
        in_env  = abs(px - xc) <= env_hw;
        in_dike = abs(px - xc) <= hw/2;
        
        if in_env && ~in_dike
            draw_capsule(ax, px, py, color, 0.72);
            drawn = drawn + 1;
        end
    end
end


% =========================================================================
%% ---- HELPER: Particles in RESERVOIR envelope (ring around reservoir) ---
function draw_particles_res_envelope(ax, xc, yc, hw_res, hh_res, ...
                                     hw_env, hh_env, density, melt_f, color)
    % Sample in bounding box of envelope ellipse, keep ring region only
    n = round(density * (1 - melt_f * 0.5) * hw_env * hh_env * pi / 80);
    
    rng(55);
    attempts = 0;
    drawn    = 0;
    while drawn < n && attempts < n * 20
        attempts = attempts + 1;
        px = (xc - hw_env) + rand() * (2*hw_env);
        py = (yc - hh_env) + rand() * (2*hh_env);
        
        % Normalised ellipse coordinates
        norm_env = ((px-xc)/hw_env)^2 + ((py-yc)/hh_env)^2;
        norm_res = ((px-xc)/hw_res)^2 + ((py-yc)/hh_res)^2;
        
        % Inside envelope ellipse but outside reservoir ellipse
        if norm_env <= 1.0 && norm_res >= 1.05
            draw_capsule(ax, px, py, color, 0.72);
            drawn = drawn + 1;
        end
    end
end


% =========================================================================
%% ---- HELPER: Draw particles inside an elliptical reservoir -------------
function draw_particles_ellipse(ax, xc, yc, hw_x, hh_y, melt_f)
    n = round(120 * (1 - melt_f * 0.7))*0.7;
    col_dark  = [0.65, 0.25, 0.12];
    col_green = [0.30, 0.55, 0.20];
    col_white = [0.98, 0.97, 0.93];
    
    rng(99);
    for i = 1:n
        ang = rand()*2*pi;
        rad = sqrt(rand());
        px  = xc + hw_x * 0.90 * rad * cos(ang);
        py  = yc + hh_y * 0.85 * rad * sin(ang);
        
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
    angle  = rand() * pi;
    length = 0.5 + rand() * 0.1;
    width  = length * (0.25 + rand()*0.20);
    
    t = linspace(0, 2*pi, 24);
    xe = (length/2) * cos(t);
    ye = (width/2)  * sin(t);
    
    rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    pts = rot * [xe; ye];
    
    fill(ax, cx + pts(1,:), cy + pts(2,:), color, ...
         'EdgeColor','none','FaceAlpha',alpha_val);
end

% =========================================================================
%% ---- HELPER: Draw mineral particles inside a dike (same as reservoir) --
function draw_particles_dike(ax, xc, hw, y_top, y_bottom)
    % Particle count scales with dike area
    area = (2*hw) * (y_bottom - y_top);
    n = round(area * 4);   % density factor — tune as needed
    
    col_dark  = [0.65, 0.25, 0.12];
    col_green = [0.30, 0.55, 0.20];
    col_white = [0.98, 0.97, 0.93];
    
    rng(88);
    for i = 1:n
        px = (xc - hw) + rand() * (2*hw);
        py = y_top      + rand() * (y_bottom - y_top);
        
        r_type = rand();
        if r_type < 0.55
            col   = col_white;
            alpha = 0.80;
        elseif r_type < 0.80
            col   = col_dark + rand()*[0.1 0.05 0.05];
            alpha = 0.85;
        else
            col   = col_green;
            alpha = 0.75;
        end
        
        draw_capsule(ax, px, py, col, alpha);
    end
end