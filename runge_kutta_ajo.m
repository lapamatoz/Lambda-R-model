% Kuva 1 / Fig 1
%final_plot_1()

% Kuva 2 / Fig 3
%final_plot_50()

% Kuva 3 / Fig 4
%final_plot_3()

% Kuva 4 / Fig 5-7
final_plot_4()

% Kuva 5 / Fig 8-16
%final_plot_many_B()

% Kuva 6 / 19 yhtälö 25
%[B_array, R_array, L_array, T_array] = final_plot_RB();

%et = char(ones(length(B_array),1) * ' & ');
%[print_num(B_array,4), et, print_num(R_array,4), et, print_num(L_array,4), et, print_num(T_array,5)]

% (Kuva 7 / Landscape)
%[B_array, R_array, L_array, T_array] = final_plot_RB_cut();


%[a_res, b_res, t_res, R_res] = R_T_all(1,0, 20, 1);
%-t_res(1)- 9.70860731629968

function c = f_color()
    c = [10,60,248]/255;
end

function c = e_color()
    c = [0,0,0];
end

function c = red_color()
    c = [255,8,18]/255;
end

function w = plot_lineWidth()
    w = .8;
end

function [a_res, b_res, t_res, R_res] = friedmann(B,L, t_end)
    [H, ~, a0] = parameters();
    options = odeset('Events',@(t,y)eventFcn(t,y),'MaxStep',0.03,'InitialStep',a0/8,'RelTol',1e-10);
    [t,y,T,ab_T,~] = ode45(@(t,y)diff_eq_F(t,y,B,L),[0, t_end/10], [a0,0], options);
    a_res = y(:,1); b_res = y(:,2); t_res = t;
    t_res = 10*(t_res - T);
    left = sum(t_res < 0);
    t_res = [t_res(1:left); 0; t_res(left+1:end)];
    a_res = [a_res(1:left); ab_T(1); a_res(left+1:end)];
    b_res = [b_res(1:left); ab_T(2); b_res(left+1:end)];
    R_res = H .* sqrt(L) .* b_res ./ a_res;
end

function [a_res, b_res, t_res, R_res] = R_T_all(B,L, t_end, alpha)
    [H, ~, a0] = parameters();
    options = odeset('Events',@(t,y)eventFcn(t,y),'MaxStep',0.03,'InitialStep',a0/8,'RelTol',1e-10);
    [t,y,T,ab_T,~] = ode45(@(t,y)diff_eq_E(t,y,B,L,alpha),[0, t_end/10],[a0; (B == 0)*a0/100], options);
    a_res = y(:,1); b_res = y(:,2); t_res = t;
    if length(T) == 0
        't_end on liian pieni. a(T) = 1 ei sisälly ratkaisuun.'
    end
    t_res = 10*(t_res - T);
    left = sum(t_res < 0);
    t_res = [t_res(1:left); 0; t_res(left+1:end)];
    a_res = [a_res(1:left); ab_T(1); a_res(left+1:end)];
    b_res = [b_res(1:left); ab_T(2); b_res(left+1:end)];
    R_res = H .* sqrt(L) .* b_res ./ a_res;
end

function res = diff_eq_E(t,y,B,L,alpha)
    [H, ~, ~] = parameters();
    res = [H * sqrt(B/y(1) + H * alpha * sqrt(L)/y(1)^2 * y(2) + L*y(1)^2)
           y(1) * exp(-H * sqrt(L) * t)];
end

function res = diff_eq_F(t,y,B,L)
    [H, ~, ~] = parameters();
    res = [H * sqrt(B/y(1) + L*y(1)^2)
           y(1) * exp(-H * sqrt(L) * t)];
end

function [position,isterminal,direction] = eventFcn(t,y)
    position = y(1) - 1; % The value that we want to be zero
    isterminal = 0;  % Halt integration
    direction = 0;   % The zero can be approached from either direction
end

function final_plot_4()
    % a(T) kuvaaja 50 Gyr asti
    figure; hold on;
    
    % F-model
    B0 = 0.317; L0 = 0.683;
    [a_res, ~, t_res0, ~] = friedmann(B0,L0, 17);
    p0 = plot(t_res0, a_res,'--','LineWidth',plot_lineWidth(),'Color',f_color());
    
    % ΛR-model
    alpha1 = 0.414;
    L = 0.6881; B1 = 0.049; D1 = 0.151;
    [a_res, ~, t_res1, ~] = R_T_all(B1+D1, L, 17, alpha1);
    p1 = plot(t_res1, a_res, '-.','LineWidth',plot_lineWidth(),'Color',e_color());
    
    % ΛR-model
    alpha2 = 0.1548;
    L = 0.6881; B2 = 0.049; D2 = 0.221;
    [a_res, ~, t_res2, ~] = R_T_all(B2+D2, L, 17, alpha2);
    p2 = plot(t_res2, a_res, '--','LineWidth',plot_lineWidth(),'Color',e_color());
    
    % ΛR-model
    alpha3 = 0.04398;
    L = 0.6881; B3 = 0.049; D3 = 0.251;
    [a_res, ~, t_res3, ~] = R_T_all(B3+D3, L, 17, alpha3);
    p3 = plot(t_res3, a_res, '-','LineWidth',plot_lineWidth(),'Color',e_color());
    
    legend([p1 p2 p3 p0],...
        {['ΛR-model: T = ', print_num(-t_res1(1),5), ', Ω^B = ', print_num(B1,4), ', Ω^Λ = ', print_num(L,4), ', Ω^D = ', print_num(D1,4), ', α = ', print_num(alpha1,4)],...
        ['ΛR-model: T = ', print_num(-t_res2(1),5), ', Ω^B = ', print_num(B2,4), ', Ω^Λ = ', print_num(L,4), ', Ω^D = ', print_num(D2,4), ', α = ', print_num(alpha2,4)],...
        ['ΛR-model: T = ', print_num(-t_res3(1),5), ', Ω^B = ', print_num(B3,4), ', Ω^Λ = ', print_num(L,4), ', Ω^D = ', print_num(D3,4), ', α = ', print_num(alpha3,4)],...
        ['F-model: T = ', print_num(-t_res0(1),5), ', Ω^B = ', print_num(B0,4), ', Ω^Λ = ', print_num(L0,4)]},...
        'Location',...
        'northwest')
    
    title(['Kuva 4 ', print_today()])
    xlabel('Time t in Gyr'); ylabel('Scale factor a(t)')
    axis([-18 0 0 1])
    draw_y_axis()
    save_plot('kuva_4')
end

function final_plot_3()
    figure; hold on;
    
    % Friedmann
    B0 = 0.317; L0 = 0.683;
    [~,~,T0,R] = friedmann(B0,L0, 70);
    p2 = plot(T0,R,'--','Color',f_color(),'LineWidth',plot_lineWidth());
    
    % ΛR-model
    B = 0.049;
    L = eval_L(B,20);
    [~,~,T,R] = R_T_all(B,L, 70,1);
    p1 = plot(T,R,'Color',e_color(),'LineWidth',plot_lineWidth());
    
    legend([p1, p2],...
    {['ΛR-model: T = ', print_num(-T(1),5), ', Ω^B = ', print_num(B,4), ', Ω^Λ = ', print_num(L,4)],...
    ['F-model: T = ', print_num(-T0(1),5), ', Ω^B = ', print_num(B0,4), ', Ω^Λ = ', print_num(L0,4)]},...
    'Location',...
    'northeast')
    axis([-20 50 0 0.3])
    xlabel('Time t in Gyr'); ylabel('Ω^{ΛR}(t)')
    draw_y_axis()
    draw_x_value('',0.2629,-12,0,e_color())
    draw_x_value('',0.2696,-12,0,f_color())
    draw_x_value('',0.2722,-5,4.8,e_color())
    draw_x_value('',0.2925,-5,6.78,f_color())
    
    draw_y_value('',6.78, 0.17, 0.2925, f_color())
    draw_y_value('',4.80, 0.16, 0.2722, e_color())
    title(['Kuva 3 ', print_today()])
    save_plot('kuva_3')
end

function draw_x_value(name,value,xmin,xIntersect,color)
    hold on
    plot([xmin,xIntersect],[value value], 'Color', color, 'LineWidth',0.5,'HandleVisibility','off')
    if name == ""
        t = text(xmin,value+0.002,[print_num(value,4), ' '],'HorizontalAlignment','right','Rotation',0,'FontSize',9);
    else
        t = text(xmin,value+0.002,[name, ' = ', print_num(value,4), ' '],'HorizontalAlignment','right','Rotation',0,'FontSize',9);
    end
    t.Color = color;
    hold off
end

function draw_y_value(name,value,ymin,yIntersect,color)
    hold on
    plot([value value],[ymin,yIntersect], 'Color', color, 'LineWidth',0.5,'HandleVisibility','off')
    if name == ""
        t = text(value,ymin,[' ', print_num(value,4)],'HorizontalAlignment','left','Rotation',0,'FontSize',9);
    else
        t = text(value,ymin,[name, ' = ', print_num(value,4), ' '],'HorizontalAlignment','left','Rotation',0,'FontSize',9);
    end
    t.Color = color;
    hold off
end

function [R,B,L] = solve_optimum()
    [~, steps, a0] = parameters();
    % Optimipisteen ratkaisu
    options = optimoptions(@fminunc, 'OptimalityTolerance', 1e-3);
    [B, R] = fminunc(@(B)real(-R_T(B, eval_L(B, 1,1,0), a0, steps, 1)), 0.045812, options);
    R = -R;
    L = 1-B-R;
end

function [R,L] = solve_B_zero()
    [~, steps, a0] = parameters();
    B = 0;
    L = fzero(@(L)R_T(B,L, 40,1) + B + L - 1,0.7);
    R = R_T(B,L,a0,steps,1);
end

function [B_array, R_array, L_array, T_array] = final_plot_RB()
    B_array = [0:0.01:0.04, 0.045812890174276, 0.049, 0.05:0.01:0.1, 0.12:0.02:0.9, 0.91:0.01:1].';
    L_array = arrayfun(@(B)eval_L(B, 19), B_array);
    R_array = 1 - B_array - L_array;
    T_array = zeros(length(B_array),1); % uncomment if T for tables
    for p = 1:length(B_array)
        [~, ~, t_res, ~] = R_T_all(B_array(p),L_array(p), 19, 1);
        T_array(p) = -t_res(1);
    end
    figure
    hold on
    plot(B_array,R_array, '+', 'Color', e_color(), 'LineWidth', plot_lineWidth())
    plot(B_array(6:7),R_array(6:7), '+', 'Color', red_color(), 'LineWidth', plot_lineWidth())
    xlabel('Ω^B'); ylabel('Ω^{ΛR}_T')
    title(['Kuva 6, ', print_today()])
    save_plot('kuva_6')
end

function [B_array, R_array, L_array, T_array] = final_plot_RB_cut()
    B_array = [0:0.01:0.04, 0.045812890174276, 0.049, 0.05:0.01:0.1].';%, 0.12:0.02:0.9, 0.91:0.01:1].';
    L_array = arrayfun(@(B)eval_L(B, 19), B_array);
    R_array = 1 - B_array - L_array;
    T_array = zeros(length(B_array),1); % uncomment if T for tables
    for p = 1:length(B_array)
        [~, ~, t_res, ~] = R_T_all(B_array(p),L_array(p), 19, 1);
        T_array(p) = -t_res(1);
    end
    figure
    hold on
    plot(B_array,R_array, '+', 'Color', e_color(), 'LineWidth', plot_lineWidth())
    plot(B_array(6:7),R_array(6:7), '+', 'Color', red_color(), 'LineWidth', plot_lineWidth())
    xlabel('Ω^B'); ylabel('Ω^{ΛR}_T')
    title(['Kuva 7, ', print_today()])
    daspect([1 1 1])
    save_plot('kuva_7')
end

function final_plot_50()
    % a(T) kuvaaja 50 Gyr asti
    figure; hold on;
    
    % f and g
    [H, ~, ~] = parameters();
    t = linspace(0,51/10,100);
    p4 = plot(t*10, exp(t*H), '--', 'Color', red_color(), 'LineWidth', plot_lineWidth()/2); % f
    t = linspace(0,51/10,100);
    p3 = plot(t*10, exp(H * sqrt(0.6881)* t), 'Color', red_color(), 'LineWidth', plot_lineWidth()/2);
    
    % F-model
    B2 = 0.317; L2 = 0.683;
    [a_res, ~, t_res2, ~] = friedmann(B2,L2, 65);
    p2 = plot(t_res2, a_res,'--','LineWidth',plot_lineWidth(),'Color',f_color());
    
    % ΛR-model
    L = 0.6881; B = 0.049;
    [a_res, ~, t_res1, ~] = R_T_all(B,L, 68, 1);
    p1 = plot(t_res1, a_res,'LineWidth',plot_lineWidth(),'Color',e_color());
    
    legend([p1 p2],...
        {['ΛR-model: T = ', print_num(-t_res1(1),5), ', Ω^B = ', print_num(B,4), ', Ω^Λ = ', print_num(L,4)],...
        ['F-model: T = ', print_num(-t_res2(1),5), ', Ω^B = ', print_num(B2,4), ', Ω^Λ = ', print_num(L2,4)]},...
        'Location',...
        'northwest')
    
    title(['Kuva 2 ', print_today()])
    xlabel('Time t in Gyr'); ylabel('Scale factor a(t)')
    draw_y_axis()
    axis([-20 50 0 20])
    save_plot('kuva_2')
end

function final_plot_1()
    figure; hold on;
    
    % F-model
    B2 = 0.317; L2 = 0.683;
    [a_res, ~, t_res2, ~] = friedmann(B2,L2, 15);
    p2 = plot(t_res2, a_res,'--','LineWidth',plot_lineWidth(),'Color',f_color());
    
    % ΛR-model
    L = 0.6881; B = 0.049;
    [a_res, ~, t_res1, ~] = R_T_all(B,L, 18, 1);
    p1 = plot(t_res1, a_res,'LineWidth',plot_lineWidth(),'Color',e_color());
    
    legend([p1 p2],...
        {['ΛR-model: T = ', print_num(-t_res1(1),5), ', Ω^B = ', print_num(B,4), ', Ω^Λ = ', print_num(L,4)],...
        ['F-model: T = ', print_num(-t_res2(1),5), ', Ω^B = ', print_num(B2,4), ', Ω^Λ = ', print_num(L2,4)]},...
        'Location',...
        'northwest')
    
    title(['Kuva 1 ', print_today()])
    xlabel('Time t in Gyr'); ylabel('Scale factor a(t)')
    draw_y_axis()
    axis([-18 0 0 1])
    save_plot('kuva_1')
end

function draw_y_axis()
    hold on
    plot([0 0], ylim, '-k', 'LineWidth',0.5,'HandleVisibility','off')
    hold off
end

function res = print_num(x, n)
    res = num2str(round(x,n,'significant'));
end

function final_plot_many_B()
    %[~, steps, a0] = parameters();
    figure;
    hold on
    t = 51;
    drawNewB(0,red_color(),':',t+18.3,true);
    drawNewB(0.02,red_color(),'--',t+17.4-0.6,false); % -10
    
    drawNewB(0.1,f_color(),'--',t+15.5,true);
    drawNewB(0.2,f_color(),'-.',t+14.06,true);
    drawNewB(0.4,f_color(),':',t+12.27,true);
    drawNewB(0.8,f_color(),':',t+10.31,true);

    drawNewB(0.046,e_color(),'-',t+16.64,true);
    xlabel('Time t in Gyr'); ylabel('Ω^{ΛR}(t)')
    title(['Kuva 5, ', print_today()])
    axis([-20 60 0 0.3])
    draw_y_axis()
    save_plot('kuva_5')
end

function p1 = drawNewB(B,color,mark,t_end,textB)
    %[~, steps, a0] = parameters();
    textRotation = -32;
    L = eval_L(B,19);
    %[~, ~, t_res, ~] = R_T_all(B,L, 1,1);
    %T_offset = t_res(end);
    [~,~,t_res,R_res] = R_T_all(B, L, t_end, 1);
    %ind = plot_indexes(t_res);
    p1 = plot(t_res,R_res,mark,'LineWidth',plot_lineWidth());
    p1.Color = color;
    if textB
        txt = [' Ω^B = ', num2str(B)];
        t = text(t_res(end),R_res(end),txt,'HorizontalAlignment','left','Rotation',textRotation);
        t.Color = color;
    end
end

function L = eval_L(B, t_end)
    %options = optimset('Display', 'iter');
    L_0 = (1-B)^1.6; % just a guess
    L = fzero(@(L) R_T(B, L, t_end) + B + L - 1, L_0);
end

function res = R_T(B,L, t_end)
    [H, steps, a0] = parameters();
    [~, b_res, t_res, ~] = R_T_all(B,L, t_end, 1);
    [~,tind] = min(abs(t_res));
    res = H * sqrt(L) * b_res(tind);
end

function res = R_T2(B,L, a_end)
    [H, steps, a0] = parameters();
    a_dot = @(t,a,b)H*(B/a + H * sqrt(L)/a^2 * b + L*a^2)^(1/2);
    b_dot = @(t,a,b)a*exp(-H*sqrt(L)*t);
    b0 = 0; t0 = 0; t_final = 10;
    [~,b_res,~] = runge_kutta(a_dot, b_dot, t_final, steps, a0, b0, t0, a_end);
    res = H * sqrt(L) * b_res(end);
end

function [a_res, b_res, t_res, R_res] = friedmann2(B,L, a_end)
    [H, steps, a0] = parameters();
    a_dot = @(t,a,b)H * sqrt(B/a + L*a^2);
    b_dot = @(t,a,b)0;
    b0 = 0; t0 = 0; t_final = 9;
    [a_res,b_res,t_res] = runge_kutta(a_dot, b_dot, t_final, steps, a0, b0, t0, a_end);
    R_res = H .* sqrt(L) .* b_res ./ a_res;
end

function [a_res, b_res, t_res, R_res] = R_T_all2(B,L, a_end, alpha)
    [H, steps, a0] = parameters();
    a_dot = @(t,a,b)H * sqrt(B/a + H * alpha * sqrt(L)/a^2 * b + L*a^2);
    b_dot = @(t,a,b)a * exp(-H*sqrt(L)*t);
    b0 = 0; t0 = 0; t_final = 9;
    [a_res,b_res,t_res] = runge_kutta(a_dot, b_dot, t_final, steps, a0, b0, t0, a_end);
    R_res = H .* sqrt(L) .* b_res ./ a_res;
end


function save_plot(name)
    grid on
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 4*screenposition(3:4)],...
        'PaperSize',[4*screenposition(3:4)]);
    %print -dpdf -painters num2str(name)
    print(name,'-dpdf','-fillpage')
end

function today = print_today()
    formatOut = '(dd.mm.yyyy)';
    today = datestr(datetime('now'),formatOut);
end

function res = plot_indexes(x)
    inc = floor(length(x) / 200);
    inc = max(1,inc);
    res = 1:inc:length(x);
    if res(end) ~= length(x)
        res = [res, length(x)];
    end
end

function [H, steps, a0] = parameters()
    H = 0.6866759;
    steps = 800;
    a0 = 10^-12;
end
