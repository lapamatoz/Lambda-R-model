function [a_res,b_res,t_res] = runge_kutta_vpa(dadt,dbdt,final_t,steps,a0,b0,t0)
%RUNGE_KUTTA Solver for a system of differential equations
%   Solves dadt and dbdt up to a point when a(T) == 1.
%   final_t should be more than T, because
%   steps are calculated from t0 to final_t, but the solver stops when
%   above condition is true.
%   The final point is iterated with a secant method.
%   The method uses quadraticly increasing step size. 

% Initialize the solution arrays
%a = vpa(zeros(steps,1));  a(1) = vpa(a0);
%b = vpa(zeros(steps,1));  b(1) = vpa(b0);
%t = vpa(zeros(steps,1));  t(1) = vpa(t0);

a = vpa(zeros(2,1));  a(1) = vpa(a0);
b = vpa(zeros(2,1));  b(1) = vpa(b0);
t = vpa(zeros(2,1));  t(1) = vpa(t0);

% k is a factor that we use to fix the number of steps
k = (final_t - t0) * (6*steps - 6) / (steps * (2*steps - 1));
h = (1 / (steps-1))^2 * k;

for p = 1:steps-1
    % We stop if we pass the point where a(T)==1.
    if a(2) > 1
        last_index = 2;
        break;
    end
    % New step size
    h = (p / (steps-1))^2 * k;
    if a(2) ~= 0
        a(1) = a(2); b(1) = b(2); t(1) = t(2);
    end
    % Evaluate new point and add it to the arrays
    [a(2), b(2), t(2)] = next_point(dadt,dbdt,a(1),b(1),t(1),h);
end
% Now we have iterated RK4 one iteration over the point a(T)==1.
% Let's iterate the last point, so condition a(T)==1 is more excact.

% Initial values for secant method
h_array = [h,h/2];

a_array = [a(last_index), 0];
b_array = [b(last_index), 0];
t_array = [t(last_index), 0];

for q = 0:100
    % Calculate new point with RK4
    [a_array(end), b_array(end), t_array(end)] = next_point(dadt,dbdt,a(last_index-1),b(last_index-1),t(last_index-1),h_array(end));
    % Evaluate new step size with previous iterations
    h_new = h_array(end) - (a_array(end)-1) * (h_array(end-1) - h_array(end)) / (a_array(end-1) - a_array(end));
    % Conditions for breaking the loop
    if h_array(end) == h_new || abs(a_array(end)-1)<1e-28
        break;
    else
        % Add the evaluated h
        a_array = [a_array(2:end), 0];
        b_array = [b_array(2:end), 0];
        t_array = [t_array(2:end), 0];
        h_array = [h_array(2:end), h_new];
    end
end

% Remove zeros from the arrays and add the iterated point
%a = [a(1:last_index-1); a_array(end)]; a_res = a(end);
%b = [b(1:last_index-1); b_array(end)]; b_res = b(end); 
%t = [t(1:last_index-1); t_array(end)]; t_res = t(end); 
a_res = a_array(end);
b_res = b_array(end);
t_res = t_array(end);
    
    % This function evaluates the next point with 4. order Runge-Kutta
    function [a_1,b_1,t_1] = next_point(dadt,dbdt,a,b,t,h)
    	k_0 = h * dadt(t, a, b);
    	l_0 = h * dbdt(t, a, b);

        k_1 = h * dadt(t + h/2, a + k_0/2, b + l_0/2);
        l_1 = h * dbdt(t + h/2, a + k_0/2, b + l_0/2);

        k_2 = h * dadt(t + h/2, a + k_1/2, b + l_1/2);
        l_2 = h * dbdt(t + h/2, a + k_1/2, b + l_1/2);

        k_3 = h * dadt(t + h, a + k_2, b + l_2);
        l_3 = h * dbdt(t + h, a + k_2, b + l_2);

        a_1 = a + 1/6 * (k_0 + 2*k_1 + 2*k_2 + k_3);
        b_1 = b + 1/6 * (l_0 + 2*l_1 + 2*l_2 + l_3);
        t_1 = t + h;
    end

end

