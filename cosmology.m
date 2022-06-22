%% Cosmology Calculations
% Lauri Jokinen, 2022.
% 
% The script is made with Matlab version 2020a, with Live Script enabled. The 
% code is much more readable with Live Script swiched on.
% 
% The script has a function for each plot, and the functions are called in section 
% 'Calling the plotting Functions'. This way toggling the plots on/off is easy.
% 
% The functions have many similar commands, so some of the code further down 
% is not commented as well as the code above it.

close all    % Close all open figures
format long  % With this Matlab won't round results in console
%% Calculation Parameters and Global Variables
% Here we declare "global" variables in struct G. This struct is passed to nearly 
% all functions.
% 
% First is Hubble constant. With the correction factor, the constant is in unit 
% (Gyr)^-1.

G.save = true;
G.correction = 3.15576 / (1.495978707 * 6.48 / pi) / 10;
G.H = 0.6726 * G.correction;
%% 
% Precision value: smaller gives more precise result. This single value affects 
% the precision of various algorithms. 

G.precision = 1e-6; % 1e-12 on aika jees
%% 
% Set options for Matlab's equation system solver |fsolve| and optimizer |fminunc|.

G.fsolveOptions  = optimoptions('fsolve', ...
            'OptimalityTolerance',G.precision, ...
            'FunctionTolerance', G.precision, ...
            'StepTolerance',G.precision, ...
            'MaxFunctionEvaluations', ceil(sqrt(1/G.precision)), ...
            'MaxIterations', ceil(sqrt(1/G.precision)), ...
            'Display', 'off');

G.fminuncOptions = optimoptions('fminunc', ...
            'OptimalityTolerance', G.precision, ...
            'StepTolerance', G.precision, ...
            'MaxFunctionEvaluations', ceil(sqrt(1/G.precision)), ...
            'MaxIterations', ceil(sqrt(1/G.precision)), ...
            'Display', 'off');
%% 
% Maximum time point for differential equation solver, in unit Gyr and starting 
% from the big bang. In some cases it is wanted that the solver is terminated 
% after point $a(T) = 1$ is reached. This is the limit, if we don't want to stop 
% the solver at that point.

G.t_end = 80;
%% 
% General parameters for the models. 

G.omegaB    = 0.049;
G.FH_omegaD = 0.268; G.FH_omegaL = 0.683;
G.F_omegaD  = 0.151; G.F_omegaL  = 0.8;
%% 
% $\Omega^{\Lambda}$ for the $\Lambda \text{R}$-model is solved from equation 
% $\Omega^{B} + \Omega^{\Lambda} + \Omega^{\Lambda R}_T = 1$. Function |flatness_solve_omegaL| 
% is defined further down.

G.LR_omegaL = flatness_solve_omegaL(G, G.omegaB, 1);
%% 
% Record the age of the universe, using the benchmark Friedmann-model.

[~, tH] = F_model(G, G.omegaB + G.FH_omegaD, G.FH_omegaL, true);
G.benchmarkAge = -tH(1);
%% 
% Parameters at which $\Omega^{\Lambda \text{R}}_T$ is maximised.

[G.omegaLR_opt, G.omegaB_opt, G.omegaL_opt] = optimize_omegaB(G);
%% 
% Parameter at which the $\Lambda \text{R}$-model returns the correct age of 
% the universe, and satistfies the equation, $\Omega^D + \alpha \Omega^{\Lambda 
% R}_T = \frac{\Omega^{\Lambda R}_{T,opt}}{\Omega^B_{opt}} \Omega^B$.

[G.omegaD_opt, G.alpha_opt] = optimize_omegaD_alpha(G);
%% Code for plotting
% Various graphs can be plotted here.

% Tests the precision of the differential equation

%[a_res, t_res, ~, y_res] = LR_model(G, 1, 0); % a0 = 1e-12 => e = 1e-10.
%disp(-t_res(1) - 2/(3*G.H));

[a_res, t_res] = LR_model(G, 0.04, 0.6); % a0 = 1e-12 => e = 1e-10.
plot(t_res(1:10), a_res(1:10),'o')
hold on
plot(t_res(1:10), a_res(1:10))

%% Differential Equation Solvers
% The ODE system is,
% 
% $\dot{A} = H_T \sqrt{\frac{1}{At^{2}} \Omega^{B+D} +  \frac{ \alpha H_T \sqrt{\Omega^{\Lambda}}}{A^2t}  
% B + A^2 \Omega^{\Lambda}} - \frac{2}{3 t}A$,
% 
% $\dot{B} = \frac{1}{t} \left[ A \cdot e^{-t \cdot H_T \sqrt{\Omega^{\Lambda}}} 
% - \frac{5}{3} B \right]$.
% 
% 

function [a, t, y] = LR_model(G, omegaBD, omegaL, alpha, terminate_T)
    if nargin==4
        terminate_T = true;
    elseif nargin==3
        terminate_T = true;
        alpha = 1;
    end

    precision = G.precision;
    
    if omegaBD ~= 0
        a0 = G.precision*10000;
        t0 = a0^(3/2) / (3/2 * G.H * sqrt(omegaBD));
        b0 = 3/5 * a0 * t0;
    else
        t0 = G.precision*10000;
        a0 = G.H * (alpha * sqrt(omegaL) / 2)^(1/3) * t0;
        b0 = G.H * (alpha * sqrt(omegaL) / 2)^(1/3) * t0^2 / 2;
    end
        
    L = G.H*sqrt(omegaL);

    % The differential equations in an array form
    ay_dot = @(t,ay)[G.H * sqrt(omegaBD/ay(1) + alpha*L*ay(2)/ay(1)^2 + ay(1)^2*omegaL);
                     ay(1) * exp(-t * L)];


    opts = odeset('Events',@(t,ay)ode_events(t, ay, ay_dot, terminate_T),'RelTol',1e-1*precision,'AbsTol',1e-2*precision); % Create events function
    
    % ay_res = [a_res; y_res]
    [t, ay, t_events, ay_events] = ode23s(ay_dot, [t0, G.t_end], [a0;b0], opts);
    
    % Adds initial points to result vectors
    ay = [[0,0]; ay];
    t = [0; t];
    
    % Add 'events' to the result vectors
    for p = 1:length(t_events)
        ind = sum(t < t_events(p));
        t = [t(1:ind);   t_events(p);   t((ind+1):end)];
        ay = [ay(1:ind,:); ay_events(p,:); ay((ind+1):end,:)];
    end

    % Remove identical points
    p = 1;
    while p < length(t)
        if t(p) == t(p+1) || ay(p,1) == ay(p+1,1) || ay(p,2) == ay(p+1,2)
            t(p+1) = [];
            ay(p+1,:) = [];
        else
            p = p+1;
        end
    end
    
    a = ay(:,1); y = ay(:,2);

    % Shift t, such that a(t=0) = 1.
    [~, ind_T] = min(abs(a - 1)); % find index of a(t) = 1
    ind_T = ind_T(1); % in case there are many such events
    t = t - t(ind_T); % shift values of t
end

function [value,isterminal,direction] = ode_events(t, y, y_dot, terminate_T)
    y_dot_val = y_dot(t,y);
    value = [y(1) - 1; y(1)*y_dot_val(2) - y(2)*y_dot_val(1)];
    isterminal = [terminate_T; 0];
    direction = [0; 0];
end
%% 
% The Friedmann model's equations are as above, but $\alpha = 0$,

function [a_res, t_res, y_res] = F_model(G, omegaBD, omegaL, terminate_T)
    if nargin==3
        terminate_T = true;        
    end
    [a_res, t_res, y_res] = LR_model(G, omegaBD, omegaL, 0, terminate_T);
end
%% 
% Simple function that returns $\Omega^{\Lambda \text{R}}_T$. This function 
% is useful with solvers furher down.

function omegaLR_T = LR_model_omegaLR_T(G, omegaB, omegaL, alpha)
    if nargin==3
        alpha = 1;
    end
    [a, ~, y] = LR_model(G, omegaB, omegaL, alpha);
    omegaLR_T = G.H * sqrt(omegaL) * y(end) ./ a(end);
end
%% Other Equation Solver Functions
% Solves $\Omega^{\Lambda}$ from the flatness algorithm $\Omega^B + \Omega^{\Lambda} 
% + \alpha \Omega^{\Lambda \text{R}}_T = 1$. Evaluation of $ \Omega^{\Lambda \text{R}}_T$ 
% requires integration.

function omegaL = flatness_solve_omegaL(G, omegaB, alpha)
    if omegaB==1
        omegaL = 0;
        return
    end
    if nargin==2
        alpha = 1;
    end
    omegaL = fsolve(@(omegaL) alpha*LR_model_omegaLR_T(G, omegaB, abs(omegaL), alpha) + omegaB + abs(omegaL) - 1, ...
                    (1-omegaB)^1.6,... % initial guess
                    G.fsolveOptions);
    omegaL = abs(omegaL);
end
%% 
% Solver for optimum $\max_{\Omega^B}{\Omega^{\Lambda \text{R}}_T}$.  $\Omega^{\Lambda}$ 
% is solved from the flatness equation, $\Omega^B + \Omega^{\Lambda} + \Omega^{\Lambda 
% \text{R}}_T = 1$. In the function we solve $\min_{\Omega^B}{(-\Omega^{\Lambda 
% \text{R}}_T)}$

function [omegaLR, omegaB, omegaL] = optimize_omegaB(G)
    [omegaB, omegaLRneg] = fminunc(@(omegaB)-LR_model_omegaLR_T(G, omegaB, flatness_solve_omegaL(G,omegaB)), ...
                                0.0458, ... % initial guess
                                G.fminuncOptions);
    omegaLR = -omegaLRneg;
    omegaL = 1 - omegaB - omegaLR;
end
%% 
% Solver for $\Omega^D$ and $\alpha$, such that age of the universe matches 
% with the Friedmann model, and $\Omega^D + \alpha \Omega^{\Lambda \text{R}}_T 
% = \frac{\Omega^{\Lambda \text{R}}_{T,opt}}{\Omega^B_{opt}} \Omega^B$.

function [omegaD_opt, alpha_opt] = optimize_omegaD_alpha(G)
    x = fsolve(@(x)targetValTime(G, x(1), x(2)),[0.26,0.083], G.fsolveOptions);
    omegaD_opt = x(1); alpha_opt = x(2);
end

function values = targetValTime(G, omegaD, alpha)
    omegaL = 1 - G.omegaLR_opt/G.omegaB_opt * G.omegaB - G.omegaB;
    [a, t, y] = LR_model(G, G.omegaB + omegaD, omegaL, alpha);
    omegaLR_T = G.H * sqrt(omegaL) * y(end) ./ a(end);
    values = [(omegaD + alpha*omegaLR_T) / G.omegaB - G.omegaLR_opt/G.omegaB_opt,...
                 t(1) + G.benchmarkAge];
 end