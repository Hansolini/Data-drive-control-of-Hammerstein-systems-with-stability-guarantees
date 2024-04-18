%% Some general settings
clc; close;
rng(2);

% Text setup
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');
set(0, 'defaultAxesLabelFontSizeMultiplier', 1.2);

% Figure setup
set(0, 'defaultFigureColor', 'w');
set(0, 'defaultAxesColor', 'w');
set(0, 'defaultFigurePosition', [100, 100, 1.1*624/2, 1.1*624*(2/4)/2]);
set(0, 'DefaultAxesClipping', 'on');
set(0, 'DefaultLineClipping', 'on');
set(0, 'defaultAxesLineWidth', 1.5);
set(0, 'defaultLineLineWidth', 2);

% Colors
NTNU_black  = [0,0,0];
NTNU_blue   = [0,80,158]/255;
NTNU_yellow = [247,208,25]/255;
NTNU_orange = [239,129,20]/255;
NTNU_brown  = [207,184,135]/255;
NTNU_purple = [176, 27, 129]/255;
NTNU_green  = [188, 208, 37]/255;


%% Loading data
load('data/convergence_10k_sat2_better_spacing.mat');

% Convergence plot
% Example bar plot
figure(1);
clf; grid on; hold on; box on;

dI = 1;
I = 2:dI:16;
h = bar(NUM_DATA(I), [sum(~IS_STABLE_LINEAR(:,I,:),3); sum(~IS_STABLE(:,I,:),3)]/(2*n_samples), 'stacked', 'LineWidth', 1.5, 'barwidth', 0.5/dI);

set(h(1), 'facecolor',NTNU_blue);
set(h(2), 'facecolor',NTNU_orange);

%plot(U*0.5)
xticks(NUM_DATA);
yticks(0:0.25:1);
yticklabels(num2str(100*(0:0.25:1)') + repmat("\%",5,1));
xlabel('\# Data');
ylabel('Unstable');

legend({'Linear','Hammerstein'},'location','best')
fancyLegend();

exportgraphics(gcf, 'figures/DMSD_convergence.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Plotting sample data
figure(3);
clf; 

grid on; hold on; box on;
plot(t,y_hammerstein,'color',NTNU_blue,'linewidth',1.5,'displayname','$y(t)$')
plot(t,u,'k--','linewidth',2,'displayname','$u(t)$')

xlabel('Time')
ylabel('Output')
legend('location','best')

fancyLegend();
exportgraphics(gcf, 'figures/DMSD_sample_data.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Plotting response
load('data\constrained_vrft.mat')

DISPLAYNAMES = {
    '$J_{mr}$';
    '$J_{est}$';
    '$J_{vrft}$';
    '$J_{vrft}$ + ideal con.';
    '$J_{vrft}$ + est. con.';
    '$J_{vrft}$ + est. Hamm. con.';
};

COLORS = {
    NTNU_blue;
    NTNU_blue;
    NTNU_blue;
    NTNU_orange;
    NTNU_green;
    NTNU_orange;
    NTNU_purple;
    NTNU_brown
};

STYLES = {
    '-';
    '-';
    '-';
    '-';
    '-';
    '-';
    '--';
    ':';
};

figure(4);
clf; grid on; hold on; box on;

figure(5);
clf; grid on; hold on; box on;
for i = 3:length(K)
    if i == 4
        continue;
    end
    % Simulate step response
    yl = step_closed_loop_hammerstein(G,K{i},@(x) x,t); % Using linear "nonlinearity" here
    yh = step_closed_loop_hammerstein(G,K{i},@(x) f(fi_hat(x)),t);
    
    % Plot response
    figure(4);
    plot(t,yl,STYLES{i}, 'Color', [COLORS{i},1], 'DisplayName', DISPLAYNAMES{i})
    
    figure(5);
    plot(t,yh,STYLES{i}, 'Color', [COLORS{i},1], 'DisplayName', DISPLAYNAMES{i})
end

figure(4);
plot(t,step(feedback(Kr*G,1),t),      '--', 'Color', [NTNU_black,1],  'DisplayName', '$y_{M_r}(t)$')
xlim([0,100])
ylim([-0.5,1.5])
xlabel('Time')
ylabel('Output')

fancyLegend();
exportgraphics(gcf, 'figures/DMSD_step_response_linear.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

figure(5)
plot(t,step(feedback(Kr*G,1),t),      '--', 'Color', [NTNU_black,1],  'DisplayName', '$y_{M_r}(t)$')
xlim([0,100])
ylim([-0.5,1.5])
xlabel('Time')
ylabel('Output')

fancyLegend();
exportgraphics(gcf, 'figures/DMSD_step_response_hammerstein.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
%%
figure(6)
clf; grid on; hold on; box on;
plot(t,step_closed_loop_hammerstein(G,K{3},@(x) x,t),STYLES{3}, 'Color', [COLORS{3},1], 'DisplayName', DISPLAYNAMES{3})

Kl = transpose(beta)*THETA_linear(:,3);
plot(t,step_closed_loop_hammerstein(G,Kl,@(x) x,t),'-', 'Color', [NTNU_green,1], 'DisplayName', '$J_{vrft}$  on linear data')

plot(t,step(feedback(Kr*G,1),t),      '--', 'Color', [NTNU_black,1],  'DisplayName', '$y_{M_r}(t)$')
xlim([0,100])
xlabel('Time')
ylabel('Output')

fancyLegend();
exportgraphics(gcf, 'figures/DMSD_step_response_linear_no_constraints.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');


figure(7)
clf; grid on; hold on; box on;
plot(t,step_closed_loop_hammerstein(G,K{3},@(x) f(fi_hat(x)),t),STYLES{3}, 'Color', [COLORS{3},1], 'DisplayName', DISPLAYNAMES{3})
plot(t,step_closed_loop_hammerstein(G,Kl,@(x) f(fi_hat(x)),t),'-', 'Color', [NTNU_green,1],  'DisplayName', '$J_{vrft}$  on linear data')

plot(t,step(feedback(Kr*G,1),t),      '--', 'Color', [NTNU_black,1],  'DisplayName', '$y_{M_r}(t)$')
xlim([0,100])
ylim([-0.5,1.5])

xlabel('Time')
ylabel('Output')

fancyLegend();
exportgraphics(gcf, 'figures/DMSD_step_response_hammerstein_no_constraints.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Functions
function fancyLegend()
    l = legend('location','southeast','color','none','EdgeColor', 'none');
    pos = l.Position;
    
    % Adjust legend's parent to the figure
    l.Parent = gcf;
    
    % Update position to be in terms of the figure
    l.Units = 'normalized';
    pos = l.Position;
    
    ax = gca;
    ax_pos = ax.Position;
    
    % Compute the absolute position and size of the legend in axis units
    x = ax.XLim(1) + diff(ax.XLim) * ((pos(1) - ax_pos(1)) / ax_pos(3));
    y = ax.YLim(1) + diff(ax.YLim) * ((pos(2) - ax_pos(2)) / ax_pos(4));
    w = diff(ax.XLim) * (pos(3) / ax_pos(3));
    h = diff(ax.YLim) * (pos(4) / ax_pos(4));

    xs = x + 0.05*w/2;
    ys = y - 0.075*h;

    rectangle('Position',[xs,ys,w,h],'FaceColor',[0.5,0.5,0.5],'Curvature',0.25,'EdgeColor','none','LineWidth',1.5)
    rectangle('Position',[x,y,w,h],'FaceColor','w','Curvature',0.25,'LineWidth',1.5)
end