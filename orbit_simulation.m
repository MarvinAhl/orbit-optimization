% These parameters are not real but from Kerbal Space Program
const.gravity = 6.67430e-11;  % N m^2 kg^-2
const.sun_mass = 1.7565459e28;  % kg
const.earth_mass = 5.2915158e22;  % kg
const.earth_sem_maj_ax = 13599840256;  % m
const.earth_period = 9203545;  % s

t0 = 0;
tf = const.earth_period / 2;
steps = 1000;  % Only for displaying the earth

t = t0:(tf-t0)/steps:tf;
earth_pos_arr = earth_pos(t, const);

%r0 = [0; -const.earth_sem_maj_ax*0.8; 0; 10000; 0; 0];
r0 = [const.earth_sem_maj_ax*0.95; 0; 0; 0; 9500; 0];

[t_ode, vehicle_r_t] = ode89(@(t, r) force(t, r, const), [t0, tf], r0);
vehicle_r = vehicle_r_t';

plot3(0, 0, 0, 'r.');  % Sun
axis equal
axis([-const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5 ...
   -const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5 ...
   -const.earth_sem_maj_ax*1.5 const.earth_sem_maj_ax*1.5]);
hold on
plot3(earth_pos_arr(1, :), earth_pos_arr(2, :), earth_pos_arr(3, :), 'b-');  % Earth
plot3(vehicle_r(1, :), vehicle_r(2, :), vehicle_r(3, :), 'k-');  % Vehicle
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