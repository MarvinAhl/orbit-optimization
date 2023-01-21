% These parameters are not real but from Kerbal Space Program
const.gravity = 6.67430e-11;  % N m^2 kg^-2
const.sun_mass = 1.7565459e28;  % kg
const.earth_mass = 5.2915158e22;  % kg
const.earth_sem_maj_ax = 13599840256;  % m
const.earth_period = 9203545;  % s

t0 = 0;
tmax = const.earth_period;  % Max time before objective must be completed
t_act = 100;  % Min time between maneuvres in s

% Todo:
%r0 = [const.earth_sem_maj_ax; -700000; 0; 2000; 9285; 0];
r0 = [const.earth_sem_maj_ax*0.95; 0; 0; 0; 9500; 0];
target_r0 = [0; -const.earth_sem_maj_ax*0.8; 0; 10000; 0; 0];

% Solve target path
target_trajectory = simulate_target(t0, tmax, target_r0, const);

% ~~~ Optimization ~~~

problem = optimproblem('ObjectiveSense', 'min');

phases = 5;  % Number of shooting phases
ts = optimvar('ts', phases, LowerBound=t0+t_act, UpperBound=tmax);  % Phase start times
rs = optimvar('rs', 6, phases-1);  % Every phases initial states. Column contains: x, y, z, vx, vy, vz

% Simulate shooting phases with simulate and fcn2optimexpr
sim_rs = fcn2optimexpr(@simulate, ts, rs, t0, r0, const, 'OutputSize', [6, phases], 'ReuseEvaluation', true, 'Analysis', 'off');
% Final target state vector
target_rf = fcn2optimexpr(@target, ts(end), target_trajectory, t0, tmax, 'OutputSize', [6, 1]);

% Velocities objective
obj_vel = sum( sqrt(sum((rs(4:6, :) - sim_rs(4:6, 1:end-1)).^2)) ) + ...  % Intermediate Delta v burns
    sqrt(sum((target_rf(4:6) - sim_rs(4:6, end)).^2));  % Final Delta v burn
obj_pos = sqrt(sum((target_rf(1:3) - sim_rs(1:3, end)).^2));  % How closely target and vehicle meet

% Add weighted objectives to optimization problem
problem.Objective = obj_vel + obj_pos;

% Time constraints
t_constr = ts(2:end) >= ts(1:end-1) + t_act;
problem.Constraints.t_constr = t_constr;

% Simulation constraints (x(i)' == x(i+1))
r_constr = rs(1:3, :) == sim_rs(1:3, 1:end-1);
problem.Constraints.r_constr = r_constr;

% Initial values: Solve trajectory with initial conditions and use intermediate points
r_guess_traj = ode89(@(t, r) force(t, r, const), [t0, tmax], r0);
t_guess = t0+tmax/phases : tmax/phases : tmax;
r_guess = deval(t_guess(1:end-1), r_guess_traj);

x0.rs = r_guess;
x0.ts = t_guess;

%solver = 'ga';

prob_options = optimoptions(problem, 'ConstraintTolerance', 10, 'UseParallel', true, 'Display', 'iter');
%prob_options = optimset('Display', 'iter');

%show(problem)
[sol, fval, eflag, output] = solve(problem, x0, 'Options', prob_options);

% ~~~ Plot ~~~
plot3(0, 0, 0, 'r.');  % Sun

axis equal
axis([-1 1 -1 1 -1 1] * const.earth_sem_maj_ax*1.5);
hold on

steps = 1000;  % Only for displaying the earth
t = t0:(tf-t0)/steps:tf;
earth_pos_arr = earth_pos(t, const);
plot3(earth_pos_arr(1, :), earth_pos_arr(2, :), earth_pos_arr(3, :), 'b-');  % Earth

% Todo: Vehicle Plot
% Todo: Target Plot

hold off


function position = earth_pos(t, const)
    position = zeros(3, length(t));
    angles = 2*pi/const.earth_period * t;
    position(1, :) = const.earth_sem_maj_ax .* cos(angles);
    position(2, :) = const.earth_sem_maj_ax .* sin(angles);
end

function dr = force(t, r, const)  % dr/dt = force(r), r = [x; y; z; vx; vy; vz]
    dr = zeros(6, 1);
    dr(1:3) = r(4:6);
    vec_obj_sun = -r(1:3);
    dr(4:6) = dr(4:6) + const.gravity*const.sun_mass / norm(vec_obj_sun)^3 * vec_obj_sun;  % Sun gravity
    vec_obj_earth = earth_pos(t, const) - r(1:3);
    dr(4:6) = dr(4:6) + const.gravity*const.earth_mass / norm(vec_obj_earth)^3 * vec_obj_earth;  % Earth gravity
end

function sim_rs = simulate(ts, rs, t0, r0, const)
    % Solve ODE in times given by ts from initial values rs to sim_rs and
    % return that. Try-Catch in case of ode errors
    phases = length(ts);
    start_ts = [t0; ts(1:end-1)];
    end_ts = ts;
    init_rs = [r0, rs];
    sim_rs = zeros(6, phases);
    parfor i = 1 : phases
        try
            traj = ode89(@(t, r) force(t, r, const), [start_ts(i), end_ts(i)], init_rs(:, i));
            sim_rs(:, i) = deval(end_ts(i), traj);
        catch
            warning(['Integration failed at phase ', num2str(i), ' failed. Setting solution to NaN']);
            sim_rs(:, i) = [NaN; NaN; NaN; NaN; NaN; NaN];
        end
    end
end

function target_trajectory = simulate_target(t0, tmax, r0, const)
    % Simulate target trajectory
    target_trajectory = ode89(@(t, r) force(t, r, const), [t0, tmax], r0);
end

function target_r = target(t, target_trajectory, t0, tmax)
    % Return Target position at time t given the precomputed target
    % simulation
    t_ = max(min(t, tmax), t0);
    target_r = deval(t_, target_trajectory);
end