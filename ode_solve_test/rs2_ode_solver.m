%% ============================================================
%  two state variable RSF solver from Claude to check for apparent
% IC dependence of solution found using my Petc code
%% ============================================================
clear; clc;close all
function dydt = rs2_rhs(~, x, par)
expx = exp(x(1));
f1 =  (1.0 - expx) * par.k;
f2 = -expx * par.r * (par.b2 * x(1) + x(3));
f0 = expx * ((par.b1 - 1.0) * x(1) + x(2) - x(3)) + f1 -  f2;

dydt = [ f0; f1; f2 ];
end
% control parameter is kuse, the non dim stiffness
% vary from 1 (stable) to ~0.8515
%//par->knd = 0.93312950; // 2 orbit
%//par->knd = 0.86164850; // 4 orbit
%//par->knd = 0.85621575; // 8 orbit
%//par->knd = 0.85522170; // 16 orbit
%//par->knd = 0.85501330; // 32 orbit
%//par->knd = 0.85496885; // 64 orbit
%//par->knd = 0.85495935; // 128 orbit
%//par->knd = 0.85495730; // 256 orbit

%
% in weird IC range
%
kuse= 0.857651;
% close by orbit where there is no IC dependence
%kuse= 0.85765;
imode=0; % x-y-z perturbation of IC leads to ~5 orbit - the weird case
%imode=1; % x perturbation of IC leads to clean  4 orbit - the expected case found for surround knds

rtol=1e-10;atol=1e-12;% tolerances
%rtol=1e-12;atol=1e-13;% tolerances

% initial condition choice
eps = 1e-5;
% this does lead to stable behavior
%eps = 1e-3;
if imode==0
    IC         = [ eps; eps; eps];
else
    IC         = [eps; 0; 0];
end
% total integration time, will use last 1% for analysis
tEnd       = 1e5;tfrac=0.99
%tEnd = 1e6;tfrac=0.999

%
% parameters
%
par.b1=1.0;
par.b2=0.84;
par.r=0.048;
% absolute critical stiffnesses
par.kcr1 = par.b1 - 1.0;
par.kcr2=(par.kcr1+par.r*(2.*par.b1+(par.b2-1.)*(2.+par.r))+ ...
    sqrt(4.*par.r*par.r*(par.kcr1+par.b2)+ ...
    (par.kcr1+par.r*par.r*(par.b2-1.))^2))/(2.+2.*par.r);
% use dimensionalized
par.k = kuse * par.kcr2;


%% --- Time span & initial condition --------------------------
t0         = 0;
t_longtime  = tEnd*tfrac;   % only report x=0 crossings after this time
t_plot_start = tEnd*tfrac; % plot only the last fraction  of the time series

%% --- High-accuracy solver options ---------------------------
opts = odeset('RelTol',  rtol,  'AbsTol',  atol, 'MaxStep', 0.05,   ...
    'Events',  @(t,y) event_x_zero(t, y, t_longtime) );

%% --- Integrate ----------------------------------------------
fprintf('Integrating 2 RSF  system...\n');
tic
[t, Y, te, Ye, ie] = ode45(@(t,y) rs2_rhs(t, y, par), [t0, tEnd], IC, opts);
%[t, Y, te, Ye, ie] = ode89(@(t,y) rs2_rhs(t, y, par), [t0, tEnd], IC, opts);
%[t, Y, te, Ye, ie] = ode15s(@(t,y) rs2_rhs(t, y, par), [t0, tEnd], IC, opts);
fprintf('Done in %.2f s  (%d steps)\n\n', toc, length(t));

x = Y(:,1);  y_var = Y(:,2);  z = Y(:,3);

% Indices for plotting — only the tail of the simulation
i_plot  = t >= t_plot_start;
tp      = t(i_plot);
xp      = x(i_plot);
yp      = y_var(i_plot);
zp      = z(i_plot);

% Crossings within the plot window
te_plot = te(te >= t_plot_start);
Ye_plot = Ye(te >= t_plot_start, :);

%% --- Report crossings ---------------------------------------
fprintf('=== x=0 crossings after t=%.1f ===\n', t_longtime);
if isempty(te)
    fprintf('  None detected.\n');
else
    dirs = {'neg->pos', 'pos->neg'};
    fprintf('  %5s  %10s  %10s  %10s  %s\n','#','t','y','z','Dir');
    for k = 1:length(te)
        fprintf('  %5d  %10.5f  %10.5f  %10.5f  %s\n', ...
            k, te(k), Ye(k,2), Ye(k,3), dirs{ie(k)});
    end
    fprintf('\n  Total: %d crossings\n', length(te));
    if length(te) >= 2
        fprintf('  Mean interval: %.4f  |  Std: %.4f\n', ...
            mean(diff(te)), std(diff(te)));
    end
end

%% --- Figure 1 : 3-D Phase Portrait (fast surface coloring) --
figure('Name','3-D Phase Portrait','Color','w','Position',[100 100 860 680]);

% Fast colored 3-D line via the surface trick:
% reshape trajectory into 2-row matrices, use surf with EdgeColor interp
n  = length(t);
Xs = [xp,  xp ]';
Ys = [yp,  yp ]';
Zs = [zp,  zp ]';
Cs = [tp,  tp ]';

surf(Xs, Ys, Zs, Cs, ...
    'EdgeColor','interp', 'FaceColor','none', 'LineWidth', 0.5);

colormap(parula);
cb = colorbar;  cb.Label.String = 'Time';  cb.FontSize = 11;
clim([t_plot_start tEnd]);

hold on;
%plot3(xp(1),   yp(1),   zp(1),   'go','MarkerSize',10,'LineWidth',2,'DisplayName','Start');
%plot3(xp(end), yp(end), zp(end), 'rs','MarkerSize',10,'LineWidth',2,'DisplayName','End');
if ~isempty(te_plot)
    plot3(Ye_plot(:,1), Ye_plot(:,2), Ye_plot(:,3), 'k^', ...
        'MarkerSize',6,'MarkerFaceColor','k','DisplayName','x=0 crossing');
end
legend('Location','best','FontSize',10);
hold off;

xlabel('x','FontSize',14,'FontWeight','bold');
ylabel('y','FontSize',14,'FontWeight','bold');
zlabel('z','FontSize',14,'FontWeight','bold');
title('Phase Space','FontSize',15,'FontWeight','bold');
grid on; box on; view(30,25);
ax = gca; ax.FontSize = 12;
%% --- Figure 2 : y(t) ----------------------------------------
figure('Name','y(t)','Color','w','Position',[100 820 860 380]);

plot(tp, yp, 'b-', 'LineWidth', 0.7);
hold on;

yl = ylim;
patch([t_plot_start tEnd tEnd t_plot_start], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 1 0.88], 'EdgeColor','none','FaceAlpha',0.4);

if ~isempty(te_plot)
    plot(te_plot, Ye_plot(:,2), 'rv','MarkerSize',7,'LineWidth',1.5, ...
        'DisplayName','x=0 crossing');
end

xlabel('Time','FontSize',13,'FontWeight','bold');
ylabel('y(t)','FontSize',13,'FontWeight','bold');
title('y(t) — RSF 2SV System','FontSize',14,'FontWeight','bold');
legend('y(t)', sprintf('Plot window (t>%g)', t_plot_start), ...
    'x=0 crossing','Location','northwest','FontSize',10);
grid on; box on; xlim([t_plot_start tEnd]);
hold off;

%% --- Figure 3 : crossing interval histogram -----------------
if length(te) >= 5
    figure('Name','Crossing Intervals','Color','w','Position',[980 100 500 380]);
    histogram(diff(te), 30, 'FaceColor',[0.2 0.5 0.8],'EdgeColor','w');
    xlabel('Inter-crossing interval','FontSize',12,'FontWeight','bold');
    ylabel('Count','FontSize',12,'FontWeight','bold');
    title('x=0 crossing intervals (long times)', ...
        'FontSize',13,'FontWeight','bold');
    grid on; box on;
end

fprintf('\nAll plots rendered.\n');




function [value, isterminal, direction] = event_x_zero(t, y, t_long)
    if t < t_long
        value = [1; 1];
    else
        value = [y(1); y(1)];
    end
    isterminal = [0; 0];
    direction  = [1; -1];
end
