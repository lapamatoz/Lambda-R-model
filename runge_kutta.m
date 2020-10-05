function [a, b, t_converted] = runge_kutta(a_dot, b_dot, ...
%% RUNGE_KUTTA Differential Equation Solver
% Solves a two equation system of first order differential equations with a 
% fourth order Runge-Kutta algorithm. The algorithm records an event where $a(T) 
% = 1$, which is solved with a secant method.
% 
% The algorithm is designed to solve differential equation, where derivative 
% $\dot{a}$ is close to beign singular at initial time, and almost a constant 
% further in time, like Friedmann differential equation. With this kind of problem, 
% the algorithm produces more precise results that Matlab's own |ode45| -algorithm. 
% This is because of cubicly increasing step size. The algorithm also solves the 
% event $a(T)=1$ more precisely than Matlab's event-feature can produce.
% 
% Algorithm integrates from $t_1=0$, but instead of returning $\bf{t}$, a converted 
% array $10({\bf t}-T)$ is returned.
% 
% Input variables are
%% 
% * |a_dot|, function $\dot{a}$
% * |b_dot|, function $\dot{b}$
% * initial value for $a$, |a_0| (initial values for $b$ and $t$ are set to 
% 0)
% * final value for $t$, |t_n|
% * initial step size $h_1$, |initial_step|
% * a choise to stop at point $a(T)=1$, |terminate_T| (true / false)
%% 
% The algorithm has three phases: integration until $a_p>1$, solving for the 
% point $a(T)=1$, and finally integrating from $T$ to $t_n$. The last phase might 
% not be run, depending on |terminate_T|.
                                 a_0, t_n, initial_step, terminate_T, findMax)
%% 
% Let's calculate the number of steps; Matlab can work with arrays that increase 
% in size at every iteration, but the code will run faster if we initialize the 
% arrays to a specific size and don't increase their lengths later.
% 
% Step sizes are $h_p = h_1 p^3$, where $p = 2,3,\ldots,n$, and $h_1$ is given. 
% We want to set the final value for time, $t_n < \sum_{p=1}^n h_p = h_1 \sum_{p=1}^n 
% p^3$, and solve for $n$. Sum of cubes is $1^3+2^3+\cdots+n^3 = \left( \frac{n(n+1)}{2} 
% \right)^2$, so $t_n < h_1 \left( \frac{n(n+1)}{2} \right)^2$ and finally we 
% can solve the number of steps $n > \frac{\sqrt{1+8\sqrt{t_n / h_1}}-1}{2}$. 
% We'll round this expression up, and because we stop along the way to calculate 
% point $a(T)=1$, we also add one to the expression.
% Solve for step count, n
n = sqrt(t_n / initial_step);
n = sqrt(1 + 8*n) - 1;
n = ceil(n / 2) + 1;
% Initialise the solution arrays
a = zeros(n,1);  a(1) = a_0;
b = zeros(n,1);
t = zeros(n,1);
% Iteration step index
p = 1;
%% Integration from $t=0$ until $a_p > 1$
% At the end of this file, we have function |next_point|, which calculates next 
% point with given condition and step size.
while p < n
    
    if a(p) > 1
        % break out of the while-loop
        break;
    end
    
    % Calculate step size
    h = p^3 * initial_step;
    
    % Evaluate new point with RK4 and add it to the arrays
    [a(p+1), b(p+1), t(p+1)] = next_point(a_dot, b_dot, a(p), b(p), t(p), h);
    
    % Index increment
    p = p+1;
    
end
%% Solving for $a(T)=1$
% Now we have $a_{p-1}<1$ and $a_p>1$. Now it's time to solve $a(T)=1$ with 
% secant method. We start from previous iteration condition, $p-1$, and try to 
% find a step size such that $a$ of the following RK4 iteration is close to 1. 
% Let's denote these iteration variables with primes, so we won't mix them with 
% the previous integration steps.
% 
% Our goal is to have $a' - 1=0$ by changing $h'$.
% 
% Let's create vectors containing last two iterations, $[h'_{k-1} ~ h'_k]$ and 
% $[a'_{k-1} ~ a'_k]$. Our initial values come from iteration indices $p-1$ and 
% $p$. The initial values are,
% 
% $h'_1 = t_{p-1} - t_{p-1}=0$,         $h'_2 = t_{p} - t_{p-1} = h_{p-1}$,
% 
% $a'_1 = a_{p-1}$,                         $a'_2=a_p$.
% 
% With Friedmann differential equation typically only 6-7 steps is needed, due 
% to its almost linear shape at $a(T)=1$.
% 
% Function for the secant method is defined below.
[a(p), b(p), t(p)] = secantMethod(@(t,a,b)a-1, t(p-1), h, a(p-1), a(p), b(p-1), b(p), a_dot, b_dot);
% Save t(p) for the conversion later on
T = t(p);
%% Integrating from $T$ to $t_n$
% If |terminate_T| is true, we slice the arrays and don't iterate anymore. If 
% |terminate_T| is false, we continue integrating the same way as before.
if terminate_T
    % Shorten the arrays
    a = a(1:p);
    b = b(1:p);
    t = t(1:p);
else
    while p < n
        % Calculate step size
        h = p^3 * initial_step;
        
        % Evaluate new point and add it to the arrays
        [a(p+1), b(p+1), t(p+1)] = next_point(a_dot, b_dot, a(p), b(p), t(p), h);
        
        % Increment p
        p = p+1;
    end
end
%% Finding maximum of b / a
% Finds maximum of $b/a \iff \frac{d (b/a)}{dt} = 0$. By chain rule, we get
% 
% $$0=\frac{d(b/a)}{dt} = \frac{a \dot{b} - b \dot{a}}{a^2}$$
% 
% and multiplying both sides with $a^2 \neq 0$, we get $a \dot{b} - b \dot{a} 
% = 0$.
if findMax
    equation = @(t,a,b)b_dot(t,a,b)*a - a_dot(t,a,b)*b;
    % Let's find the indexes, q-1 and q, between which the wanted point lies
    q = 1;
    while equation(t(q), a(q), b(q)) > 0 && q < length(a)
        q = q+1;
    end
    
    [a(q), b(q), t(q)] = secantMethod(equation, ...
                            t(q-1), (q-1)^3 * initial_step, a(q-1), a(q), b(q-1), b(q), a_dot, b_dot);
end
%% 
% Finally, we do a unit conversion. The algorithm is finished after this line.
t_converted = t - T;
%% Evaluating a Runge-Kutta step
% Here is the function that evaluates next point, with RK4. The code is straight 
% forward and can be found here
% 
% England, Roland. "Error estimates for Runge-Kutta type solutions to systems 
% of ordinary differential equations." _The Computer Journal_ 12.2 (1969): 166-170.
function [a_res, b_res, t_res] = next_point(a_dot, b_dot, a, b, t, h)
    k_0 = h * a_dot(t, a, b);
    l_0 = h * b_dot(t, a, b);
    k_1 = h * a_dot(t + h/2, a + k_0/2, b + l_0/2);
    l_1 = h * b_dot(t + h/2, a + k_0/2, b + l_0/2);
    k_2 = h * a_dot(t + h/2, a + k_1/2, b + l_1/2);
    l_2 = h * b_dot(t + h/2, a + k_1/2, b + l_1/2);
    k_3 = h * a_dot(t + h, a + k_2, b + l_2);
    l_3 = h * b_dot(t + h, a + k_2, b + l_2);
    a_res = a + 1/6 * (k_0 + 2*k_1 + 2*k_2 + k_3);
    b_res = b + 1/6 * (l_0 + 2*l_1 + 2*l_2 + l_3);
    
    t_res = t + h;
end
%% Secant Method
% We want to solve $f\left(t_0 +h,a\left(t_0 +h\right),b\left(t_0 +h\right)\right)=0$, 
% with respect to $h$.
% 
% Initial values to be provided are,
% 
% $a_0 =a\left(t_0 \right)$,      $a_1 =a\left(t_0 +h_0 \right)$,
% 
% $b_0 =b\left(t_0 \right)$,      $b_1 =b\left(t_0 +h_0 \right)$,
% 
% along with $t_0$ and $h_0$.
function [a,b,t] = secantMethod(fun, t0, h0, a0, a1, b0, b1, a_dot, b_dot)
%% 
% Let's declare two dimensional vectors for $h$ and $f$.
    h_k = [0, ...
           h0];
    
    f_k = [fun(t0, a0, b0), ...
           fun(t0+h0, a1, b1)];
%% 
% Set $a$, $b$ and $t$ to something, in case the secant equation is singular 
% right from the start.
    a = a0; b = b0; t = t0;
    
    % 100 step limit for the secant method
    for w = 1:100
%% 
% Let's calculate, if the secant equation produces singular value or $h_{k+2} 
% = h_{k+1}$. The secant equation is
% 
% $h_{k+2} = h_{k+1} - f(h_{k+1}) \frac{h_{k+1} - h_{k}}{f(h_{k+1}) - f(h_{k})}$.
% 
% Condition $h_{k+2} = h_{k+1} \iff f(h_{k+1}) =0$ or $h_{k+1} - h_{k} =0 $.
% 
% The secant equation is singular, if $f(h_{k+1}) - f(h_{k}) = 0$.
        if f_k(2) == 0 || h_k(1) - h_k(2) == 0 || f_k(1) - f_k(2) == 0
            break;
%% 
% If the equation is ok, we continue:
        else
            % Evaluate new step size with previous iteration with secant equation
            h_new = h_k(2) - ...
                f_k(2) * (h_k(1)-h_k(2)) / (f_k(1)-f_k(2));
            
            % Update the older array values
            f_k(1) = f_k(2);
            h_k(1) = h_k(2);
            
            % Calculate new point with RK4
            [a, b, t] =...
                    next_point(a_dot, b_dot, a0, b0, t0, h_new);
            
            % Update the newer aray values
            h_k(2) = h_new;
            f_k(2) = fun(t, a, b);
        end
    end
end
end