% These parameters are not real but from Kerbal Space Program
const.gravity = 6.67430e-11;  % N m^2 kg^-2
const.sun_mass = 1.7565459e28;  % kg
const.earth_mass = 5.2915158e22;  % kg
const.earth_sem_maj_ax = 13599840256;  % m
const.earth_period = 9203545;  % s

t0 = 0;
tmax = const.earth_period;  % Max time before objective must be completed
t_act = 100;  % Min time between maneuvres

% Todo: Solve target path

problem = optimproblem();

phases = 20;  % Number of shooting phases
ts = optimvar('ts', phases, LowerBound=t0+t_act, UpperBound=tmax);  % Phase start times
rs = optimvar('rs', 6, phases);  % Every phases initial states. Column contains: x, y, z, vx, vy, vz

% Todo: Simulate shooting phases with simulate and fcn2optimexpr
% Todo: Objective

% Time constraints
t_constr = optimconstr(phases-1);
for i = 1 : phases-1
    t_constr(i) = ts(i+1) >= ts(i) + t_act;
end
problem.Constraints.t_constr = t_constr;

% Todo: Collocation constraints (x(i)' == x(i+1))

% ~~~ Plot ~~~
plot3(0, 0, 0, 'r.');  % Sun

axis equal
axis([-const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5 ...
   -const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5 ...
   -const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5]);
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

function new_r = simulate(t, r, const)
    % Solve ODE in times given by to from initial values r to new_r and
    % return that. Try-Catch for ode errors
    new_r = NaN;
end

function target_trajectory = simulate_target(t0, tmax, r0, const)
    % Simulate target trajectory
    target_trajectory = ode89(@(t, r) force(t, r, const), [t0, tmax], r0);
end

function target_r = target(t, target_trajectory)
    % Return Target position at time t given the precomputed target
    % simulation
    target_r = deval(t, target_trajectory);
end