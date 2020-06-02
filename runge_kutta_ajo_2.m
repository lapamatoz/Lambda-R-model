%% Cosmology Calculations
% The script has a function for each plot. This way toggling the plots on/off is
% easy. The script is done with Matlab version 2020a, which supports
% functions inside scripts, thus with older versions each function
% must be separated to Matlab function files. The script part is first
% and all the functions are declared afterwards.
%% Calculation parameters
% Here we declare global variables, First is Hubble constant, with
% correction factor, in unit (10 Gyr)^-1.
global H_T; H_T = 0.6715 * 1.0226;
%%
% Precision value: smaller gives more precise result. This single value
% affects various tolerance values in different algorithms.
global a0; a0 = 10^-12;
%%
% Maximum t-value for differential equation solver, in giga years and
% starting from the big bang.
% In some cases it is wanted that the solver is terminated after point
% $a(T) = 1$ is reached.
% This is the limit, if we don't want to stop the solver at this point.
global t_end; t_end = 69;
%%
% Parameters for the models. $\Omega^{\Lambda}$ for the $\Lambda R$-model
% is solved from equation
% $\Omega^{B} + \Omega^{\Lambda} + \Omega^{\Lambda R}_T = 1$. Function
% for |solve_omegaL| is provided later.
global F_omegaB;  F_omegaB  = 0.049;
global F_omegaD;  F_omegaD  = 0.268;
global F_omegaL;  F_omegaL  = 0.683;
global LR_omegaB; LR_omegaB = F_omegaB;
global LR_omegaL; LR_omegaL = solve_omegaL(LR_omegaB);

% Prints parameters at which Ω^{ΛR}_T is maximised.
global omega_LR_opt
global omega_B_opt
global omega_L_opt

[omega_LR_opt,omega_B_opt,omega_L_opt] = solve_optimum();
%disp(legend_text('Optimal parameters', omega_LR_opt, NaN, omega_B_opt, omega_L_opt, NaN, NaN))

%% Calling the Plotting Functions
% Here you can select which Figures you want to plot. Comment out the ones
% that you don't want.

%[a_res, t_res, omega_LR] = LR_model(1,0, 1,true);
%-t_res(1) - 10* 2/(3*H_T)

% Kuva 1 / Fig 1
final_plot_1()

% Kuva 2 / Fig 3
final_plot_50()

% Kuva 3 / Fig 4
final_plot_3()

% Kuva 4 / Fig 5-7
final_plot_4()

% Kuva 5 / Fig 8-16
final_plot_many_B()

% Kuva 6 / 19 yhtälö 25
[omegaB, omegaLR, omegaL, T, omegaLR_lin, omegaL_lin] = final_plot_RB();
print_table_two_col([omegaB,omegaLR,omegaL,T])
print_table_two_col([omegaB, omegaLR_lin, omegaL_lin])

%% Differential equations
% Here we define the differential equations used.
% In both models, we use a $b$-function to solve for
% $\Omega^{\Lambda R}_t$ with
% 
% $$\Omega^{\Lambda R}_t = \frac{H \sqrt{\Omega^{\Lambda}}}{a} b.$$
% 
% First the Friedmann
% model.
% 
% $$\dot{a} = H_T \sqrt{\frac{\Omega^{B+D}}{a} + a^2 \Omega^{\Lambda}}$$
% 
% $$\dot{b} = a e^{-H_T \sqrt{\Omega^{\Lambda}} t}$$
% 
% Matlab uses single vector to represent $a$ and $b$, so
% in the code, |ab(1)| is $a$ and |ab(2)| is $b$.

function res = diff_eq_F(t,ab,BD,L)
    global H_T
    res = [H_T * sqrt(BD/ab(1) + L*ab(1)^2)
           ab(1) * exp(-H_T * sqrt(L) * t)];
end

%%
% Next, $\Lambda R$-model,
% 
% $$\dot{a} = H_T \sqrt{\frac{\Omega^{B+D}}{a} + \frac{\alpha H_T \sqrt{\Omega^{\Lambda}}}{a^3} b + a^2 \Omega^{\Lambda}}$$
% 
% $$\dot{b} = a e^{-H_T \sqrt{\Omega^{\Lambda}} t}$$
% 
function res = diff_eq_E(t,ab,BD,L,alpha)
    global H_T
    res = [H_T * sqrt(BD/ab(1) + H_T * alpha * sqrt(L)/ab(1)^2 * ab(2) + L*ab(1)^2)
           ab(1) * exp(-H_T * sqrt(L) * t)];
end

%% Differential equation solvers
% We use Matlab's ode45 -algorithm. We solve the equations with time unit
% 10 Gyrs, and the unit conversion to Gyrs is done after that.
% As above, |ab(:,1)| is $a$, |ab(:,2)| is $b$, and we solve
% $\Omega^{\Lambda R}_t$ with
% 
% $$\Omega^{\Lambda R}_t = \frac{H \sqrt{\Omega^{\Lambda}}}{a} b.$$

function [a_res, t_res, omega_LR] = F_model(omegaBD,omegaL, terminate_T)
    % Load global variables
    global H_T t_end a0;
    
    % Some options regarding accuracy and the point a(T) == 1
    options = odeset('Events',@(t,y)eventFunction(t,y,terminate_T),...
                     'MaxStep',a0^(1/6),...
                     'InitialStep',a0/8,...
                     'RelTol',a0^0.6);
                 
    % Run the solver
    [t,ab,T,~,~] = ode45(@(t,ab)diff_eq_F(t,ab,omegaBD,omegaL),...
                         [0, t_end/10],... % integration area
                         [a0,0],... % initial values for a and b
                         options); % load options set above
    
    % Alert the user if solver didn't reach a(T) == 1
    if isempty(T)
        disp(['t_end is too small, so that a(T) = 1 is not reached.'...
            'Increase the value in parameter listing.'...
            'You can also try to increase accuracy.'])
    end
    
    % Set T = 0 and do a unit conversion
    t_res = 10*(t - T);
    
    % Final values
    a_res = ab(:,1);
    omega_LR = H_T .* sqrt(omegaL) .* ab(:,2) ./ a_res;
end

%%
% Very similar to Friedmann algrorithm. See comments above.
function [a_res, t_res, omega_LR] = LR_model(omegaB,omegaL, alpha, terminate_T)
    global H_T t_end a0;
    options = odeset('Events',@(t,ab)eventFunction(t,ab,terminate_T),...
                     'MaxStep',a0^(1/6),...
                     'InitialStep',a0/8,...
                     'RelTol',a0^0.6);
                 
    % For some reason the solver gets stuck when omegaB == 0.
    % That's why we start with a small initial b, when omegaB == 0.
    % When omegaB ~= 0 our initial b is 0.
    
    [t,ab,T,~,~] = ode45(@(t,ab)diff_eq_E(t,ab,omegaB,omegaL,alpha),...
                         [0, t_end/10],... % integration area
                         [a0; (omegaB == 0)*a0/100],... % initial values for a and b
                         options); % load options set above
    
    if isempty(T)
        disp(['t_end is too small, so that a(T) = 1 is not reached.'...
              'Increase the value in parameter listing.'...
              'You can also try to increase accuracy.'])
    end
    t_res = 10*(t - T);
    a_res = ab(:,1);
    omega_LR = H_T .* sqrt(omegaL) .* ab(:,2) ./ a_res;
end

%%
% An event function for ode45-algorithm. With this function we record
% the event $a(T) = 1$, and possibly stop the integration when that point
% is reached.

function [position,isterminal,direction] = eventFunction(t,y,terminate_T)
    position = y(1) - 1; % The value that we want to be zero
    isterminal = terminate_T;  % Halt integration at the event
    direction = 0;   % The event can be approached from either direction
end

%%
% Simple function that returns $\Omega^{\Lambda R}_T$. This is useful
% with solvers furher down.
function omegaLR_T = LR_model_omegaLR_T(omegaB,omegaL,alpha)
    % 'true' below means, that we stop the integration at a(T) == 1.
    [~, ~, omega_LR] = LR_model(omegaB, omegaL, alpha, true);
    omegaLR_T = omega_LR(end);
end

%% Other Solver Functions
% Solves $\Omega^{\Lambda}$ from the flatness algorithm
% $\Omega^B + \Omega^{\Lambda} + \Omega^{\Lambda R}_T = 1$.
function omegaL = solve_omegaL(omegaB)
    omegaL_0 = (1-omegaB)^1.6; % merely an initial guess
    global a0
    options = optimset('TolFun',a0);
    omegaL = fzero(@(omegaL) LR_model_omegaLR_T(omegaB, omegaL, 1) + omegaB + omegaL - 1, omegaL_0, options);
end

%%
% Same as above, except $a$ is approximated with $a = H_T t$.
% No integration is needed in this algorithm.
function omegaL = solve_omegaL_linear(omegaB)
    if omegaB == 1
        omegaL = 0;
    else
        global a0
        options = optimset('TolFun',a0);
        omegaL = fzero(@(omegaL)real(1 - omegaB - (1-exp(-sqrt(omegaL))*(1+sqrt(omegaL)))/sqrt(omegaL) - omegaL), 0.7*(1-omegaB), options);
    end
end

%%
% Solves $\alpha$ from $\Omega^B + \Omega^{\Lambda} + \alpha \Omega^{\Lambda R}_T  + \Omega^D = 1$
function alpha = solve_alpha(omegaB, omegaD, omegaL)
    global a0
    options = optimset('TolFun',a0);
    alpha = fzero(@(alpha)omegaB + omegaD + omegaL + alpha * LR_model_omegaLR_T(omegaB + omegaD, omegaL, alpha) - 1, 0.5, options);
end

%%
% Solver for optimum $\max_{\Omega^B}{\Omega^{\Lambda R}_T}$. 
% $\Omega^{\Lambda}$ is solved from the flatness equation,
% $\Omega^B + \Omega^{\Lambda} + \Omega^{\Lambda R}_T = 1$.
function [omegaLR, omegaB, omegaL] = solve_optimum()
    global a0
    options = optimoptions(@fminunc, 'OptimalityTolerance', a0,'Display','off');
    [omegaB, omegaLR] = fminunc(@(omegaB)-LR_model_omegaLR_T(omegaB, solve_omegaL(omegaB),1), 0.0458, options);
    omegaLR = -omegaLR;
    omegaL = 1-omegaB-omegaLR;
end

%% Tools for plotting
% This section contains some general functions for plotting.

% Saves current plot to current directory. This function adds grid.
% |leg| is legend object.
function save_plot(name, leg)
    grid on
    if ~isa(leg,'double')
        leg.ItemTokenSize = [15,9]; % default [30,18]
    end
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 3.5*screenposition(3:4)],...
        'PaperSize',[3.5*screenposition(3:4)]);
    %print -dpdf -painters num2str(name)
    print(name,'-dpdf','-fillpage')
end

% Function for creating legend texts from numerical data
function res = legend_text(name,R,T,B,D,alpha,L)
    res = [name ': '];
    if ~isnan(R)
        res = [res 'Ω^{Λ{\itR}}_{\itT} = ' print_num(R,4) ', '];
    end
    if ~isnan(T)
        res = [res '{\itT} = ' print_num(T,5) ', '];
    end
    if ~isnan(B) && isnan(D)
        res = [res 'Ω^{\itB} = ' print_num(B,4) ', '];
    end
    if ~isnan(B) && ~isnan(D)
        BD = 0;
        if ~isnan(B)
            BD = B;
        end
        if ~isnan(D)
            BD = BD + D;
        end
        res = [res 'Ω^{{\itB}+{\itD}} = ' print_num(BD,4) ', '];
    end
    if ~isnan(alpha)
        res = [res '{\itα} = ' print_num(alpha,4) ', '];
    end
    if ~isnan(L)
        res = [res 'Ω^Λ = ' print_num(L,4) ', '];
    end
    res = res(1:end-2);
end

% Returns current date in string format
function today = print_today()
    formatOut = '(dd.mm.yyyy)';
    today = datestr(datetime('now'),formatOut);
end

% Prints a two-column table in LaTeX-format
function print_table_two_col(matrix)
    et = char(ones(length(matrix(:,1)),1) * ' & ');
    res = char(ones(length(matrix(:,1)),1) * ' ');
    for p = 1:(length(matrix(1,:))-1)
        res = [res, print_num(matrix(:,p),4), et];
    end
    res = [res, print_num(matrix(:,end),4)];
    lines = ceil(length(matrix(:,1))/2);
    % If we have uneven number of rows, we need to add a blank row to the end
    if mod(length(matrix(:,1)),2) == 1
        blanks = char(ones(1,length(res(1,:))) * ' ');
    else
        blanks = [];
    end
    line = char(ones(lines,1) * ' \\');
    et = char(ones(lines,1) * ' &');
    disp([res(1:lines,:), et, [res(lines+1:end,:);blanks], line])
end

% Define colors
function c = plotColors(str)
    if str == "blue"
        c = [50,110,230]/255;
    elseif str == "red"
        c = [250,42,30]/255;
    elseif str == "green"
        c = [24,185,50]/255;
    else
        c = [0,0,0];
    end
end

% Define line thickness in the plots
function w = plotLineWidth()
    w = 1;
end

% Draws vertical line x = 0 to a plot
function draw_y_axis()
    hold on
    plot([0 0], ylim, '-k', 'LineWidth',0.5,'HandleVisibility','off')
    hold off
end

% Number rounding function
function res = print_num(x, n)
    res = num2str(round(x,n,'significant'));
end

%% Figure 1
% Plot showing Friedmann model and ΛR-model until a(T) = 1.
function final_plot_1()
    % Create new figure frame
    figure; hold on;
    
    % Load global variables
    global F_omegaB F_omegaD F_omegaL
    global LR_omegaB LR_omegaL
    
    % We only need to plot functions until a(T) = 1,
    % so we stop the differential equation solver at that point
    terminate_T = true;
    
    % ΛR-model alpha = 1
    alpha = 1;
    
    % Run the Friedmann model
    [a0, t0, omega_LR0] = F_model(F_omegaB + F_omegaD, F_omegaL, terminate_T);
    p0 = plot(t0, a0, '--', 'LineWidth', plotLineWidth(), 'Color', plotColors('blue'));
    
    % Run the ΛR-model
    [a1, t1, omega_LR1] = LR_model(LR_omegaB, LR_omegaL, alpha, terminate_T);
    p1 = plot(t1, a1, 'LineWidth', plotLineWidth(), 'Color', plotColors('black'));
    
    % Set plot legend
    leg = legend([p1 p0],...
        {legend_text('Λ{\itR}-model', omega_LR1(end), -t1(1), LR_omegaB, NaN,      NaN, LR_omegaL),...
         legend_text('F-model',  omega_LR0(end), -t0(1), F_omegaB,  F_omegaD, NaN, F_omegaL)},...
        'Location',...
        'northwest');
    
    % Set title etc., and save the figure
    title(['Kuva 1 ', print_today()])
    xlabel('Time {\itt} in Gyr'); ylabel('Scale factor {\ita}({\itt})')
    axis([-18 0.5 0 1.1])
    draw_y_axis()
    save_plot('kuva_1', leg)
end

%% Figure 2
function final_plot_50()
    % Create new figure frame
    figure; hold on;
    
    % Load global variables
    global H_T
    global F_omegaB F_omegaD F_omegaL
    global LR_omegaB LR_omegaL
    
    % Create t-axis points for f and g
    t = linspace(0,51,200);
    
    % Plot f
    pf = plot(t, exp(t/10 * H_T), '--', 'LineWidth', plotLineWidth(), 'Color', plotColors('green'));
    % Plot g
    pg = plot(t, exp(t/10 * H_T * sqrt(LR_omegaL)), 'LineWidth', plotLineWidth(), 'Color', plotColors('green'));
    
    % ΛR-model alpha = 1
    alpha = 1;
    
    % F-model
    [a0, t0, ~] = F_model(F_omegaB + F_omegaD, F_omegaL, false);
    p0 = plot(t0, a0, '--', 'LineWidth', plotLineWidth(), 'Color', plotColors('blue'));
    
    % ΛR-model
    [a1, t1, ~] = LR_model(LR_omegaB, LR_omegaL, alpha, false);
    p1 = plot(t1, a1, 'LineWidth', plotLineWidth(), 'Color', plotColors('black'));
    
    % Add legend to the plot
    leg = legend([p1 p0 pf pg],...
        {legend_text('Λ{\itR}-model',NaN,-t1(1),LR_omegaB, NaN,NaN, LR_omegaL),...
         legend_text('F-model', NaN,-t0(1),F_omegaB,F_omegaD,NaN,F_omegaL),...
         'Function {\itf}',...
         'Function {\itg}'},...
        'Location',...
        'northwest');
    
    % Add title etc. and save the plot
    title(['Kuva 2 ', print_today()])
    xlabel('Time {\itt} in Gyr'); ylabel('Scale factor {\ita}({\itt})')
    draw_y_axis()
    axis([-20 50 0 20])
    save_plot('kuva_2', leg)
end

%% Figure 3
function final_plot_3()
    figure; hold on;
    global F_omegaB F_omegaD F_omegaL
    global LR_omegaB LR_omegaL
    
    % Go pass the point a(T) == 1
    terminate_T = false;
    
    % Friedmann
    [~,t0,omega_LR0] = F_model(F_omegaB + F_omegaD, F_omegaL, terminate_T);
    p0 = plot(t0, omega_LR0, '--', 'Color', plotColors('blue'), 'LineWidth', plotLineWidth());
    
    % ΛR-model alpha = 1
    alpha = 1;
    
    % ΛR-model
    [~, t1, omega_LR1] = LR_model(LR_omegaB, LR_omegaL, alpha, terminate_T);
    p1 = plot(t1, omega_LR1, 'Color', plotColors('black'), 'LineWidth', plotLineWidth());
    
    leg = legend([p1, p0],...
            {legend_text('Λ{\itR}-model',NaN, -t1(1), LR_omegaB, NaN,      NaN, LR_omegaL),...
             legend_text('F-model', NaN, -t0(1), F_omegaB,  F_omegaD, NaN, F_omegaL)},...
             'Location',...
             'southeast');
    
    % Separation of the omegaLR_T texts
    textOffset = 0.003;
    
    % Solve for the omegaLR_T and draw the values on the plot
    [~,~,omegaLR] = LR_model(LR_omegaB, LR_omegaL, 1, true);
    draw_x_value('', omegaLR(end), -12, 0, plotColors('black'),-textOffset)
    [~,~,omegaLR] = F_model(F_omegaB + F_omegaD, F_omegaL, true);
    draw_x_value('', omegaLR(end), -12, 0, plotColors('blue'), textOffset)
    
    % Set title etc. and save the figure
    title(['Kuva 3 ', print_today()])
    xlabel('Time {\itt} in Gyr'); ylabel('Ω^{Λ{\itR}}({\itt})')
    axis([-20 50 0 0.3])
    draw_y_axis()
    save_plot('kuva_3', leg)
end

function draw_x_value(name,value,xmin,xIntersect,color,offset)
    hold on
    plot([xmin,xIntersect],[value value], 'Color', color, 'LineWidth',0.5,'HandleVisibility','off')
    if name == ""
        t = text(xmin,value+offset,[print_num(value,4), ' '],'HorizontalAlignment','right','Rotation',0,'FontSize',9);
    else
        t = text(xmin,value+offset,[name, ' = ', print_num(value,4), ' '],'HorizontalAlignment','right','Rotation',0,'FontSize',9);
    end
    t.Color = color;
    hold off
end

%% Figure 4
function final_plot_4()
    % a(T) kuvaaja 50 Gyr asti
    figure; hold on;
    global F_omegaB
    global LR_omegaB LR_omegaL
    global Dopt
    
    % F-model
    [a_res, t_res0, omega_LR0] = F_model(F_omegaB + Dopt, LR_omegaL, true);
    p0 = plot(t_res0, a_res,'-','LineWidth',plotLineWidth(),'Color',plotColors('blue'));
    
    % ΛR-model
    D1 = 0.2 - LR_omegaB;
    alpha1 = solve_alpha(LR_omegaB, D1, LR_omegaL);
    [a_res, t_res1, omega_LR1] = LR_model(LR_omegaB+D1, LR_omegaL, alpha1,true);
    p1 = plot(t_res1, a_res, ':','LineWidth',plotLineWidth(),'Color',plotColors('black'));
    
    % ΛR-model
    D2 = 0.27 - LR_omegaB;
    alpha2 = solve_alpha(LR_omegaB, D2, LR_omegaL);
    [a_res, t_res2, omega_LR2] = LR_model(LR_omegaB+D2, LR_omegaL, alpha2,true);
    p2 = plot(t_res2, a_res, '-.','LineWidth',plotLineWidth(),'Color',plotColors('black'));
    
    % ΛR-model
    D3 = 0.3 - LR_omegaB;
    alpha3 = solve_alpha(LR_omegaB, D3, LR_omegaL);
    [a_res, t_res3, omega_LR3] = LR_model(LR_omegaB+D3, LR_omegaL, alpha3,true);
    p3 = plot(t_res3, a_res, '--','LineWidth',plotLineWidth(),'Color',plotColors('black'));
    
    leg = legend([p1 p2 p3 p0],...
        {legend_text('Λ{\itR}-model', omega_LR1(end),-t_res1(1),LR_omegaB, D1,alpha1, LR_omegaL),...
         legend_text('Λ{\itR}-model', omega_LR2(end),-t_res2(1),LR_omegaB, D2,alpha2, LR_omegaL),...
         legend_text('Λ{\itR}-model', omega_LR3(end),-t_res3(1),LR_omegaB, D3,alpha3, LR_omegaL),...
         legend_text('F-model',  omega_LR0(end),-t_res0(1),F_omegaB, Dopt, NaN, LR_omegaL)},...
        'Location',...
        'northwest');
    
    title(['Kuva 4 ', print_today()])
    xlabel('Time {\itt} in Gyr'); ylabel('Scale factor {\ita}({\itt})')
    axis([-18 0.5 0 1.1])
    draw_y_axis()
    save_plot('kuva_4', leg)
end

%% Figure 5
function final_plot_many_B()
    figure;
    global Bopt
    hold on
    t = 51;
    drawNewB(0,plotColors('red'),':',true);
    drawNewB(0.02,plotColors('red'),'--',false); % -10
    
    drawNewB(0.1,plotColors('blue'),'--',true);
    drawNewB(0.2,plotColors('blue'),'-.',true);
    drawNewB(0.4,plotColors('blue'),':',true);
    drawNewB(0.8,plotColors('blue'),':',true);

    drawNewB(Bopt,plotColors('black'),'-',true);
    xlabel('Time {\itt} in Gyr'); ylabel('Ω^{Λ{\itR}}({\itt})')
    title(['Kuva 5, ', print_today()])
    axis([-20 50 0 0.3])
    draw_y_axis()
    save_plot('kuva_5', NaN);
end

function p1 = drawNewB(B,color,mark,textB)
    textRotation = -25;
    L = solve_omegaL(B);
    [~,t_res,omega_LR] = LR_model(B, L, 1, false);
    %(t_res < 50)
    p1 = plot(t_res,omega_LR,mark,'LineWidth',plotLineWidth());
    p1.Color = color;
    R_omegatext = omega_LR(t_res < 50);
    R_omegatext = R_omegatext(end);
    if textB
        txt = [' Ω^{\itB} = ', print_num(B,4)];
        t = text(50.5,R_omegatext,txt,'HorizontalAlignment','left','Rotation',textRotation,'FontSize',9);
        t.Color = color;
    end
end

%% Figure 6
% The function outputs data which can be visualised later.
function [omegaB_table, omegaLR_table, omegaL_table, T_table, omegaLR_table_lin, omegaL_table_lin] = final_plot_RB()
    % New figure frame
    figure; hold on;
    
    % Load global variables
    global Bopt LR_omegaB
    
    % Create a vector for values of Ω^B, which we want to put in a table
    omegaB_table = [0:0.01:0.04, Bopt, LR_omegaB, 0.05, 0.06:0.01:0.1, 0.12:0.02:0.9, 0.91:0.01:1].';
    
    % Add some points for the plot. We'll merge this vector with the
    % values for the table
    omegaB_plot = [0:0.00125:0.04, 0.8:0.05:1].';
    omegaB_plot = sort(unique([omegaB_table; omegaB_plot]));
    
    % Solve for Ω^L with each Ω^B
    % We also calculate another set of Omega^L using the linear model
    omegaL_plot     = arrayfun(@(omegaB) solve_omegaL(omegaB),        omegaB_plot);
    omegaL_plot_lin = arrayfun(@(omegaB) solve_omegaL_linear(omegaB), omegaB_plot);
    
    % Calculate Ω^ΛR using flatness equation
    omegaLR_plot     = 1 - omegaB_plot - omegaL_plot;
    omegaLR_plot_lin = 1 - omegaB_plot - omegaL_plot_lin;
    
    % Plotting
    p2 = plot(omegaB_plot,omegaLR_plot_lin, '--', 'Color', plotColors('green'), 'LineWidth', plotLineWidth());
    p1 = plot(omegaB_plot,omegaLR_plot,     '-',  'Color', plotColors('black'), 'LineWidth', plotLineWidth());
    
    % Find indices of Bopt and LR_omegaB
    [~, ind] =    min(abs(LR_omegaB - omegaB_plot));
    [~, indOpt] = min(abs(Bopt      - omegaB_plot));
    
    % Draw red x at Bopt and LR_omegaB
    plot(omegaB_plot([ind,indOpt]), omegaLR_plot_lin([ind,indOpt]), '+', 'Color', plotColors('red'), 'LineWidth', plotLineWidth());
    plot(omegaB_plot([ind,indOpt]), omegaLR_plot([ind,indOpt]),     '+', 'Color', plotColors('red'), 'LineWidth', plotLineWidth());
    
    leg = legend([p1 p2],...
        {'Λ{\itR}-model',...
         'Λ{\itR}-model with linear approximation for {\ita}({\itt})'},...
        'Location',...
        'northeast');
    
    % Set title etc. and save the figure
    xlabel('Ω^{\itB}'); ylabel('Ω^{Λ{\itR}}_{\itT}')
    title(['Kuva 6, ', print_today()])
    save_plot('kuva_6', leg)
    
    %daspect([1,1,1])
    %axis([0 0.1 0.23 0.28])
    %title(['Kuva 7, ', print_today()])
    %save_plot('kuva_7', NaN);
    
    % Let's create the table
    omegaL_table = zeros(length(omegaB_table),1);
    T_table = zeros(length(omegaB_table),1); % this is for tables, not plot
    
    for p = 1:length(omegaB_table)
        [~, ind] = min(abs(omegaB_table(p) - omegaB_plot));
        omegaL_table(p) = omegaL_plot(ind);
    end
    omegaLR_table = 1 - omegaB_table - omegaL_table;
    for p = 1:length(omegaB_table)
        [~, t_res, ~] = LR_model(omegaB_table(p),omegaL_table(p), 1, true);
        T_table(p) = -t_res(1);
    end
    
    omegaL_table_lin = zeros(length(omegaB_table),1);
    T_table_lin = zeros(length(omegaB_table),1); % this is for tables, not plot
    for p = 1:length(omegaB_table)
        [~, ind] = min(abs(omegaB_table(p) - omegaB_plot));
        omegaL_table_lin(p) = omegaL_plot_lin(ind);
    end
    omegaLR_table_lin = 1 - omegaB_table - omegaL_table_lin;
end
