%%% DECLARE GLOBALLY USED VARIABLES

% They are not defined as global variables, but a struct that is passed
% to functions.

% Hubble constant, and unit conversion factor
G.convertion_factor = 0.001 * 3.15576 / (1.495978707 * 6.48 / pi);
G.Hubble_constant = 67.26 * G.convertion_factor;

% Integration termination time
% (put a low value for speed; high value for robustness)
G.ode_t_termination = 30;

% Tune this to vary the precision of all numerical algorithms
G.precision = 1e-3;

% Options for numerical algorithms
G.fsolveOptions  = optimoptions('fsolve', ...
            'OptimalityTolerance',G.precision, ...
            'FunctionTolerance', G.precision, ...
            'StepTolerance', G.precision, ...
            'MaxFunctionEvaluations', ceil(sqrt(1/G.precision)), ...
            'MaxIterations', ceil(sqrt(1/G.precision)), ...
            'Display', 'off');
G.fminuncOptions = optimoptions('fminunc', ...
            'OptimalityTolerance', G.precision, ...
            'StepTolerance', G.precision, ...
            'MaxFunctionEvaluations', ceil(sqrt(1/G.precision)), ...
            'MaxIterations', ceil(sqrt(1/G.precision)), ...
            'Display', 'off');




%%% Example


[a, t, Omega_LR, b] = LR_model(G,0.049,0.6881,1);
plot(t, a)
title('Example \LambdaR-model')
xlabel('Time in Gyrs')
ylabel('a(t)')




%%% Run algorithms to find various solutions, optimums etc.

% These functions defined below


G.Omega_B = 0.049; % this is used in the functions below

% Solve flatness equation with specific parameters
G.LR_Omega_L = flatness_solve_Omega_L(G, G.Omega_B, 1);

% Evaluate the age of the Universe, with benchmark model
[~, t_benchmark] = LR_model(G, 0.049 + 0.268, 0.683, 0);
G.benchmarkAge = -t_benchmark(1);

disp('Evaluating optimal Omega^B... this may take a while')
[G.Omega_LR_opt, G.Omega_B_opt, G.Omega_L_opt] = optimal_Omega_B(G);
disp('Done!')

[G.Omega_D_opt, G.alpha_opt] = optimal_Omega_D_and_alpha(G);







%%% ODE SYSTEM

% Main function
function [a, t, Omega_LR, b, T_index] = LR_model(G, Omega_BD, Omega_L, alpha)
    % G is the global struct
    % T_index has the property a(T_index) == 1 (or closest to it)

    a0 = 1e-16;
    b0 = 0;
    t0 = 0;

    y_dot = @(t,y)odes(t, y, G.Hubble_constant, Omega_BD, Omega_L, alpha);

    opts = odeset('Events',@(t,y)ode_events(t, y, y_dot), ...
        'RelTol',1e-1*G.precision, ...
        'AbsTol',1e-2*G.precision); % Create events function
    
    [t, y, t_events, y_events] = ode23s(y_dot, [t0, G.ode_t_termination], [a0;b0], opts);
    
    % Add the 'events' to the result vectors, one by one
    for p = 1:length(t_events)
        ind = sum(t < t_events(p));
        t = [t(1:ind);   t_events(p);   t((ind+1):end)];
        y = [y(1:ind,:); y_events(p,:); y((ind+1):end,:)];
    end
    
    a = y(:,1);
    b = y(:,2);
    
    % find the event, a(T) = 1, and shift time axis, s.t., a(t=0) = 1
    [~, T_index] = min(abs(a - 1));
    T_index = T_index(1);
    t = t - t(T_index);
    
    % evaluate Omega_Lambda_R(t)
    Omega_LR = G.Hubble_constant * sqrt(Omega_L) * b ./ a;
end

% ODE functions
function y_dot = odes(t, y, H, omegaBD, omegaL, alpha) % differential equation
    y_dot = [H * sqrt(omegaBD/y(1) + alpha*H*sqrt(omegaL)*y(2)/(y(1))^2 + (y(1))^2*omegaL);
             y(1) * exp(-t * H * sqrt(omegaL))];
end

% Declare ODE events
function [value,isterminal,direction] = ode_events(t, y, y_dot)
    y_dot_value = y_dot(t,y);
    value = [y(1) - 1;                                    % event a(t) = 1
             y(1)*y_dot_value(2) - y(2)*y_dot_value(1)];  % event max Omega^{Lambda R}_T
    % ignore further definitions
    isterminal = [0; 0];
    direction = [0; 0];
end

% This is used in the functions below. It returns a single value, omegaLR_T
function omegaLR_T = LR_model_Omega_LR_T(G, omegaB, omegaL, alpha)
    [~, ~, omega_LR, ~, ind_T] = LR_model(G, omegaB, omegaL, alpha);
    omegaLR_T = omega_LR(ind_T);
end



%%% EQUATION SOLVERS


% Solves omega^Lambda from the flatness equation
function omegaL = flatness_solve_Omega_L(G, omegaB, alpha)
    if omegaB==1
        omegaL = 0;
        return
    end
    omegaL = fsolve(@(omegaL) alpha*LR_model_Omega_LR_T(G, omegaB, abs(omegaL), alpha) + omegaB + abs(omegaL) - 1, ...
                    (1-omegaB)^1.6,... % initial guess
                    G.fsolveOptions);
    omegaL = abs(omegaL); % this is to suppress complex result
end

% Solves optimal Omega^B
function [Omega_LR, Omega_B, Omega_L] = optimal_Omega_B(G)
    [Omega_B, Omega_LR] = fminunc(@(Omega_B)-LR_model_Omega_LR_T(G, Omega_B, flatness_solve_Omega_L(G,Omega_B,1), 1), ...
                                   0.0458, ... % initial guess
                                   G.fminuncOptions);
    Omega_LR = -Omega_LR; % change sign because of minimization instead of max.
    Omega_L = 1 - Omega_B - Omega_LR;
end

% Solves optimal Omega_D and alpha
function [Omega_D_opt, alpha_opt] = optimal_Omega_D_and_alpha(G)
    x = fsolve(@(x)objective_function_optimal_Omega_D_and_alpha(G, x(1), x(2)),[0.26,0.083], G.fsolveOptions);
    Omega_D_opt = x(1);
    alpha_opt = x(2);
end

% This is used in the above solver
function res = objective_function_optimal_Omega_D_and_alpha(G, Omega_D, alpha)
    Omega_L = 1 - G.Omega_LR_opt/G.Omega_B_opt * G.Omega_B - G.Omega_B;
    [~, t, Omega_LR, ~, T_index] = LR_model(G, G.Omega_B + Omega_D, Omega_L, alpha);
    res = [(Omega_D + alpha*Omega_LR(T_index)) / G.Omega_B - G.Omega_LR_opt/G.Omega_B_opt;
           t(1) + G.benchmarkAge];
 end
