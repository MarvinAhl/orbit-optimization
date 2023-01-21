% These parameters are not real but from Kerbal Space Program
const.gravity = 6.67430e-11;  % N m^2 kg^-2
const.sun_mass = 1.7565459e28;  % kg
const.earth_mass = 5.2915158e22;  % kg
const.earth_sem_maj_ax = 13599840256;  % m
const.earth_period = 9203545;  % s

t0 = 0;
tmax = const.earth_period;  % Max time before objective must be completed

r0 = [const.earth_sem_maj_ax; -700000; 0; 2000; 9285; 0];
r0s = [r0 r0 r0];
t0s = [0; const.earth_period/3; const.earth_period*2/3];

sim_rs = simulate(t0s(2:3), r0s(:, 2:3), t0s(1), r0s(:, 1), const);

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
    start_ts = [t0; ts(1:phases-1)];
    end_ts = ts;
    init_rs = [r0, rs(:, 1:phases-1)];
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