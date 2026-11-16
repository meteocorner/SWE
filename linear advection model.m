% ==================================================================================================================================================
% MAIN scripts of the exercise for linear advection model
% ==================================================================================================================================================
% GIVEN parameter settings:
params.idim = 99;
params.del_x = 10000; % m
params.del_t = 300; % s
params.U_const = 10; % m/s
% and the initial height field: 
% h(x) = h_0*exp(-(x-x_0)^2/(2*sigma^2)), 
h_0 = 10; % 10 m
sigma = 100*1000; % 100 km

%% exercise 1. initial height distribution
% initialize the grids on x-axis, 
% for convenience, grid is centered at x_0 = 0
half_len = ((params.idim - 1) / 2) * params.del_x; % scalar
x_grid = linspace(-half_len, half_len, params.idim); % vector, linspace(point interval pair-begin, point interval pair-end, number of points)
h_init = h_0 * exp(-x_grid.^2 / (2 * sigma^2)); % vector operation

figure
plot(x_grid,h_init);
title("initial height distribution");
xlabel("x (m)");
ylabel("h (m)");
ylim([0 10]);

save_figure("ex1")

%% exercise 2. initial height tendency
h_tend = height_tendency(h_init, params.del_x, params.U_const); % function height_tendency(h_t, del_x, U_const) is defined at the end

figure
fig = plot(x_grid(2:end-1),h_tend);
title("initial height tendency");
xlabel("x (m)");
ylabel("∂h/∂t (m/s)");

save_figure('ex2')

%% exercise 3.1 run the model for 10 time steps applying constant b.c.; 
% return a plot of h at the initial and final times (in the same plot).
params.ts = 10; % total number of time steps
h_tx = linear_advection_model(h_init, params, 'constant');

figure
fig_0 = plot(x_grid,h_tx(1,:));
hold on
fig_ts = plot(x_grid,h_tx(end,:));
title("height distribution");
xlabel("x (m)");
ylabel("h (m)");
ylim([0 10]);
legend([fig_0,fig_ts],"initial","final")

save_figure('ex3-1')

%% exercise 3.2 run the model for 100 time steps applying constant b.c.; 
% return a plot of h at the initial and final times (in the same plot).
params.ts = 100;
h_tx = linear_advection_model(h_init, params, 'constant');

figure
fig_0 = plot(x_grid,h_tx(1,:));
hold on
fig_ts = plot(x_grid,h_tx(end,:));
title("height distribution");
xlabel("x (m)");
ylabel("h (m)");
ylim([0 10]);
legend([fig_0,fig_ts],"initial","final")

save_figure('ex3-2')

%% exercise 4.1 run the model for 100 time steps applying periodic b.c.; 
% return a plot of h at the initial and final times (in the same plot).
params.ts = 100;
h_tx = linear_advection_model(h_init, params, 'periodic');

figure
fig_0 = plot(x_grid,h_tx(1,:));
hold on
fig_ts = plot(x_grid,h_tx(end,:));
title("height distribution");
xlabel("x (m)");
ylabel("h (m)");
ylim([0 10]);
legend([fig_0,fig_ts],"initial","final")

save_figure('ex4-1')

%% Exercise 4.2 make a long run so that the perturbation propagates many times around the grid. 
% plot distribution of h at the initial time and after [320, 640, 960,…, 2880] time steps in the same plot.
params.ts = 2880;
h_tx = linear_advection_model(h_init, params, 'periodic');

figure
for time_step = 0:320:2880
    fig_t = plot(x_grid,h_tx(time_step+1,:));
    hold on
end
title("height distribution");
xlabel("x (m)");
ylabel("h (m)");
legend("initial","320","640","960","1280","1600","1920","2240","2560","2880");

save_figure('ex4-2')

%% =================================================================================================================================================
% FUNCTIONS: ∂h/∂t + U * ∂h/∂x = 0, U_const is constant
% ==================================================================================================================================================
% height tendency at time t, ∂h/∂t = -U * ∂h/∂x
function h_tend = height_tendency(h_t, del_x, U_const)
    % h_tend : height tendency at time t
    % h_t : height distribution at t time step
    % del_x : grid spacing on x-axis
    % U_const : constant velocity on x-axis

    % height gradient on x at time t
    h_grad_x = (h_t(3:end) - h_t(1:end-2)) / (2 * del_x); 
    % discretization at central difference, vector operation
    % ∂h/∂x ≈ (h(x + Δx) - h(x - Δx)) / (2 * Δx)

    % height tendency at time t
    h_tend = -U_const * h_grad_x;
    % ∂h/∂t = -U * ∂h/∂x
end

% time extrapolation
% forward Euler: h(t+1) = h(t) + del_t * h_tend, 2:end-1
function h_next = forward_euler(h_t, h_tend, del_t)
    % h_next : height distribution at t+1 time step
    % h_t : height distribution at t time step
    % h_tend : height tendency at time t
    % del_t : time step

    h_next = h_t(2:end-1) + del_t * h_tend;
end

% leap-frog: h(t+1) = h(t-1) + 2 * del_t * h_tend, 2:end-1
function h_next = leap_frog(h_prev, h_tend, del_t)
    % h_next : height distribution at t+1 time step
    % h_prev : height distribution at t-1 time step
    % h_tend : height tendency at time t
    % del_t : time step

    h_next = h_prev(2:end-1) + 2 * del_t * h_tend;
end

% boundary conditions
function h_t = apply_boundary_condition(h_init, h_t, bc_type)
    % h_t : height distribution at t time step, boundary values are updated in this function
    % h_init : height distribution at initial time step
    % bc_type : boundary condition types, string, 'constant'||'periodic'
    % constant b.c. : the height values at the boundary are fixed 
    % periodic b.c. : the height values at the boundary are updated by the penultimate points

    switch lower(bc_type)
        % convert 'text' to lower case

        case 'constant'
            % constant b.c.
            h_t(1)   = h_init(1);     
            h_t(end) = h_init(end);   
            
        case 'periodic' 
            % periodic b.c.
            % h(1) = h(end-1)
            % h(end) = h(2)
            h_t(1)   = h_t(end-1); 
            h_t(end) = h_t(2);     
    end
end


% linear advection model
function h_tx = linear_advection_model(h_init, params, bc_type)
    % h_tx : h(time step, x-axis), array, size (ts+1) x idim
    % h_init : height distribution at initial time step
    % params : parameters, struct
    % params.del_x : grid spacing
    % params.del_t : time step 
    % params.U_const : constant velocity on x-axis
    % params.ts : total number of time steps
    % bc_type : boundary condition types, string, 'constant'||'periodic'

    % 1. preallocate
    h_tx = zeros(params.ts + 1, params.idim); % 1-based index in MATLAB
    h_tx(1, :) = h_init;

    % 2. first step (t = 0 -> 1): forward Euler
    % time extrapolation
    h_tend = height_tendency(h_tx(1, :), params.del_x, params.U_const);
    h_tx(2, 2:end-1) = forward_euler(h_tx(1, :), h_tend, params.del_t); 
    % update boundary
    h_tx(2, :) = apply_boundary_condition(h_init, h_tx(2, :), bc_type);

    % 3. steps 2 and later: leap_frog
    for t = 2:params.ts
        % time extrapolation
        h_tend = height_tendency(h_tx(t, :), params.del_x, params.U_const);
        h_tx(t+1, 2:end-1) = leap_frog(h_tx(t-1, :), h_tend, params.del_t);
        % update boundary
        h_tx(t+1, :) = apply_boundary_condition(h_init, h_tx(t+1, :), bc_type);
    end
end

% Save current figure 
function save_figure(file_name)
    fig = gcf; % get current figure
    set(fig, 'Units', 'inches', 'Position', [1, 1, 10, 5]);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 10, 5]);
    
    print(fig, file_name, '-dpng', '-r192'); % ., ., format, resolution
end