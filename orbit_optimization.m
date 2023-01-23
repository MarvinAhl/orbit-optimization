% These parameters are not real but from Kerbal Space Program
const.gravity = 6.67430e-11;  % N m^2 kg^-2
const.sun_mass = 1.7565459e28;  % kg
const.earth_mass = 5.2915158e22;  % kg
const.earth_sem_maj_ax = 13599840256;  % m
const.earth_period = 9203545;  % s

t0 = 0;
tmax = const.earth_period;  % Max time before objective must be completed
t_act = 100;  % Min time between maneuvres in s

% Set initial conditions:
%r0 = [const.earth_sem_maj_ax; -700000; 0; 2000; 9285; 0];
r0 = [const.earth_sem_maj_ax*0.95; 0; 0; 0; 9500; 0];
%target_r0 = [0; -const.earth_sem_maj_ax*0.8; 0; 10000; 0; 0];
target_r0 = [const.earth_sem_maj_ax*0.8/sqrt(2); -const.earth_sem_maj_ax*0.8/sqrt(2); 0; 10000/sqrt(2); 10000/sqrt(2); 0];

% Solve target path
target_trajectory = simulate_target(t0, tmax, target_r0, const);

% ~~~ Optimization ~~~

problem = optimproblem('ObjectiveSense', 'min');

phases = 10;  % Number of shooting phases (minimum 2)
% These values are scaled, always remember to unscale them before use
ts_sc = optimvar('ts_sc', phases, LowerBound=sc_t(t0+t_act, t0, tmax), UpperBound=sc_t(tmax, t0, tmax));  % Phase start times
rs_sc = optimvar('rs_sc', 6, phases-1);  % Every phases initial states. Column contains: x, y, z, vx, vy, vz

% Simulate shooting phases with simulate and fcn2optimexpr
sim_rs_sc = fcn2optimexpr(@simulate_sc, ts_sc, rs_sc, t0, r0, tmax, const, 'OutputSize', [6, phases], 'ReuseEvaluation', true, 'Analysis', 'off');

% Use this code instead of the line above if you want to split up the simulate function instead of having it vectorized
% sim_rs_sc = optimexpr(6, phases);
% for i = 1 : phases
%     sim_rs_sc(:, i) = fcn2optimexpr(@simulate_single_sc, i, ts_sc, rs_sc, t0, r0, tmax, const, 'OutputSize', [6, 1], 'ReuseEvaluation', true, 'Analysis', 'off');
% end

% Final target state vector
target_rf_sc = fcn2optimexpr(@target_sc, ts_sc(end), target_trajectory, t0, tmax, const, 'OutputSize', [6, 1]);

% Velocities objective
obj_vel = sum( sqrt(sum((rs_sc(4:6, :) - sim_rs_sc(4:6, 1:end-1)).^2)) ) + ...  % Intermediate Delta v burns
    sqrt(sum((target_rf_sc(4:6) - sim_rs_sc(4:6, end)).^2));  % Final Delta v burn
obj_pos = sqrt(sum((target_rf_sc(1:3) - sim_rs_sc(1:3, end)).^2));  % How closely target and vehicle meet

% Add weighted objectives to optimization problem
problem.Objective = obj_vel;  % + obj_pos;

% Time constraints
t_act_sc = sc_t(t_act, t0, tmax);
t_constr = ts_sc(2:end) >= ts_sc(1:end-1) + t_act_sc;
problem.Constraints.t_constr = t_constr;

% Simulation constraints (x(i)' == x(i+1))
r_constr = sqrt(sum((rs_sc(1:3, :) - sim_rs_sc(1:3, 1:end-1)).^2)) == zeros(1, phases-1);
problem.Constraints.r_constr = r_constr;

% Final target position constraint
rf_constr = sqrt(sum((target_rf_sc(1:3) - sim_rs_sc(1:3, end)).^2)) == 0;
problem.Constraints.rf_constr = rf_constr;

% Initial values: Solve trajectory with initial conditions and use intermediate points
t_guess = (t0+(tmax-t0)/phases : (tmax-t0)/phases : tmax)';
%t_guess(1) = 200;  % Time isn't altered by the solver very much so I have to set the first time to something reasonable manually

% r Guess 1
% r_guess_traj = ode89(@(t, r) force(t, r, const), [t0, tmax], r0);
% r_guess = deval(t_guess(1:end-1), r_guess_traj);

% r Guess 2
rf = target(tmax, target_trajectory, t0, tmax);
r_guess = initial_r(r0, rf, phases);

x0.rs_sc = sc_r(r_guess, const);
x0.ts_sc = sc_t(t_guess, t0, tmax);

% Finite Difference Step is set a little bigger than default because of simulation errors
% Problem is already scaled manually but there is no harm in more scaling
prob_options = optimoptions(problem, 'UseParallel', true, 'Display', 'iter', 'MaxFunctionEvaluations', 50000, 'FiniteDifferenceStepSize', 1e-7, 'ScaleProblem', true);

%show(problem)
tic
[sol, fval, eflag, output] = solve(problem, x0, 'Options', prob_options);
toc

sol_ts = usc_t(sol.ts_sc, t0, tmax);
sol_rs = usc_r(sol.rs_sc, const);

% Debug
%disp(num2str(sol_ts - t_guess, t0, tmax));

% ~~~ Plot ~~~
plot3(0, 0, 0, 'r*');  % Sun

title('Trajectory');
xlabel('x / m');
ylabel('y / m');
zlabel('z / m');
axis equal
axis([-1 1 -1 1 -1 1] * const.earth_sem_maj_ax*1.5);
hold on

steps = 1000;  % Only for displaying earth and target
tf = sol_ts(end);
t = t0 : (tf-t0)/steps : tf;
earth_pos_arr = earth_pos(t, const);
plot3(earth_pos_arr(1, :), earth_pos_arr(2, :), earth_pos_arr(3, :), 'b--');  % Earth

% Vehicle Plot
[~, eval_rs, burn_points] = evaluate(sol_ts, sol_rs, t0, r0, target(tf, target_trajectory, t0, tmax), const);
plot3(eval_rs(1, :), eval_rs(2, :), eval_rs(3, :), 'k-');
plot3(burn_points(1, :), burn_points(2, :), burn_points(3, :), 'r.');

% Target Plot
%idx = find(target_trajectory.x <= tmax);
%plot3(target_trajectory.y(1, 1:idx(end)), target_trajectory.y(2, 1:idx(end)), target_trajectory.y(3, 1:idx(end)), 'g-');
target_traj = target(t, target_trajectory, t0, tmax);
plot3(target_traj(1, :), target_traj(2, :), target_traj(3, :), 'g-');

hold off

% Analytically compute earth's position on a circular orbit
function position = earth_pos(t, const)
    position = zeros(3, length(t));
    angles = 2*pi/const.earth_period * t;
    position(1, :) = const.earth_sem_maj_ax .* cos(angles);
    position(2, :) = const.earth_sem_maj_ax .* sin(angles);
end

% Computes spiral from initial to final point as initial guess
function r_init = initial_r(r0, rf, phases)
    r_init = zeros(6, phases-1);

    % Get position as equally spaced spiral points around the sun towards target
    [r0_az, r0_el, r0_r] = cart2sph(r0(1), r0(2), r0(3));
    [rf_az, rf_el, rf_r] = cart2sph(rf(1), rf(2), rf(3));
    if rf_az <= r0_az
        rf_az = rf_az + 2*pi;
    end
    azs = linspace(r0_az, rf_az, phases+1);
    els = linspace(r0_el, rf_el, phases+1);
    rs = linspace(r0_r, rf_r, phases+1);
    [xs, ys, zs] = sph2cart(azs(2:end-1), els(2:end-1), rs(2:end-1));
    r_init(1:3, :) = [xs; ys; zs];

    % Velocity will just point from one point to the next
    v0 = norm(r0(4:6));
    comb_rs = [r_init(1:3, :), rf(1:3)];
    pointing = comb_rs(:, 2:end) - comb_rs(:, 1:end-1);
    r_init(4:6, :) = pointing ./ vecnorm(pointing) * v0;
end

% The ode's equations of motion
function dr = force(t, r, const)  % dr/dt = force(r), r = [x; y; z; vx; vy; vz]
    dr = zeros(6, 1);
    dr(1:3) = r(4:6);
    vec_obj_sun = -r(1:3);
    dr(4:6) = dr(4:6) + const.gravity*const.sun_mass / norm(vec_obj_sun)^3 * vec_obj_sun;  % Sun gravity
    vec_obj_earth = earth_pos(t, const) - r(1:3);
    dr(4:6) = dr(4:6) + const.gravity*const.earth_mass / norm(vec_obj_earth)^3 * vec_obj_earth;  % Earth gravity
end

% Does the multiple-shooting simulation
function sim_rs_sc = simulate_sc(ts_sc, rs_sc, t0, r0, tmax, const)
    % Solve ODE in times given by ts from initial values rs to sim_rs and
    % return that. Try-Catch in case of ode errors
    ts = usc_t(ts_sc, t0, tmax);
    rs = usc_r(rs_sc, const);
    phases = length(ts);
    start_ts = [t0; ts(1:end-1)];
    end_ts = ts;
    init_rs = [r0, rs];
    sim_rs = zeros(6, phases);
    parfor i = 1 : phases  % for instead of parfor may be faster depending on the circumstances
        try
            traj = ode89(@(t, r) force(t, r, const), [start_ts(i), end_ts(i)], init_rs(:, i));
            sim_rs(:, i) = deval(end_ts(i), traj);
        catch
            warning(['Integration failed at phase ', num2str(i), '. Setting solution to NaN']);
            sim_rs(:, i) = [NaN; NaN; NaN; NaN; NaN; NaN];
        end
    end
    sim_rs_sc = sc_r(sim_rs, const);
end

% Like simulate_sc but only simulates one phase
function sim_r_sc = simulate_single_sc(phase, ts_sc, rs_sc, t0, r0, tmax, const)
    ts = usc_t(ts_sc, t0, tmax);
    rs = usc_r(rs_sc, const);
    start_ts = [t0; ts(1:end-1)];
    end_ts = ts;
    init_rs = [r0, rs];
    try
        traj = ode89(@(t, r) force(t, r, const), [start_ts(phase), end_ts(phase)], init_rs(:, phase));
        sim_r = deval(end_ts(phase), traj);
    catch
        warning('Integration failed. Setting solution to NaN');
        sim_r = [NaN; NaN; NaN; NaN; NaN; NaN];
    end
    sim_r_sc = sc_r(sim_r, const);
end

% Computes the trajectory for final plot and prints useful information
function [eval_ts, eval_rs, burn_points] = evaluate(ts, rs, t0, r0, rf, const)
    phases = length(ts);
    start_ts = [t0; ts(1:end-1)];
    end_ts = ts;
    init_rs = [r0, rs];
    next_rs = [rs, rf];
    next_vs = next_rs(4:6, :);

    eval_ts = [];
    eval_rs = [];
    burn_points = zeros(6, phases);
    parfor i = 1 : phases
        try
            traj = ode89(@(t, r) force(t, r, const), [start_ts(i), end_ts(i)], init_rs(:, i));
            eval_ts = [eval_ts; traj.x'];
            eval_rs = [eval_rs, traj.y];
            burn_points(:, i) = deval(end_ts(i), traj);
        catch
            warning('Integration failed at final evaluation');
        end
    end

    % Print required dv text
    total_dv = 0;
    for i = 1 : phases
       dv = chop(next_vs(:, i) - burn_points(4:6, i));
       total_dv = total_dv + norm(dv);
       disp(['t = ', num2str(round(end_ts(i))), ' s;', 9, 'dv = ', num2str(norm(dv)), ' m/s', 9, '(x=', num2str(dv(1)), '; y=', num2str(dv(2)), '; z=', num2str(dv(3)), ')']); 
    end
    disp(['Total dv = ', num2str(total_dv), ' m/s']);
    
    dr = burn_points(1:3, end) - rf(1:3, end);
    disp(['Final deviation: ', num2str(norm(dr)), ' m (x=', num2str(dr(1)), '; y=', num2str(dr(2)), '; z=', num2str(dr(3)), ')']);
end

% Compute target trajectory in advance
function target_trajectory = simulate_target(t0, tmax, r0, const)
    % Simulate target trajectory
    target_trajectory = ode89(@(t, r) force(t, r, const), [t0, tmax], r0);
end

function target_r_sc = target_sc(t_sc, target_trajectory, t0, tmax, const)
    % Return Target position at time t given the precomputed target
    % simulation
    t = usc_t(t_sc, t0, tmax);
    t_ = max(min(t, tmax), t0);
    target_r = deval(t_, target_trajectory);
    target_r_sc = sc_r(target_r, const);
end

function target_r = target(t, target_trajectory, t0, tmax)
    % Return Target position at time t given the precomputed target
    % simulation
    t_ = max(min(t, tmax), t0);
    target_r = deval(t_, target_trajectory);
end

% Scale r to have values in the order of 1
function scaled_r = sc_r(r, const)
    scaled_r = r;
    % Just some arbitrary scaling factors to have all variables in about the
    % range from 0 to 1
    scaled_r(1:3, :) = scaled_r(1:3, :) / const.earth_sem_maj_ax;
    scaled_r(4:6, :) = scaled_r(4:6, :) / 15000;
end

% Unscale r back to normal values
function r = usc_r(scaled_r, const)
    r = scaled_r;
    r(1:3, :) = r(1:3, :) * const.earth_sem_maj_ax;
    r(4:6, :) = r(4:6, :) * 15000;
end

function scaled_t = sc_t(t, t0, tmax)
    scaled_t = (t - t0) / (tmax - t0);
end

function t = usc_t(scaled_t, t0, tmax)
    t = scaled_t * (tmax - t0) + t0;
end

function b = chop(a, eps)
    if (nargin < 2)
        eps = 1e-5;
    end

    b = a;
    b(abs(b) < eps) = 0;
end