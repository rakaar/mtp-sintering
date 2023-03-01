n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:t_max;

% particle params
m = 1;
k=1;
damping_frac = 0.8;
v0 = 1e-6;
radius = 1;

% initial vals
distance = zeros(length(t_durn),1);
velocity = zeros(length(t_durn),1);
force_e = zeros(length(t_durn),1);
force_d = zeros(length(t_durn),1);

% initialise
distance(1) = 2*radius*rand;
velocity(1) = rand;
force_e(1) = k*distance(1);
force_d(1) = damping_frac*sqrt(2*m*k)*velocity(1);

for t=2:length(t_durn)
    force_e(t) = k*distance(t-1);
    force_d(t) = damping_frac*sqrt(2*m*k)*velocity(t-1);

    velocity(t) = velocity(t-1) + (  (force_e(t) + force_d(t))/m   )*dt;
    distance(t) = distance(t-1) + velocity(t)*dt;
    
end
disp('done')
%%
close all
figure
    subplot(2,2,1)
    plot(distance)
    title('distance')

    subplot(2,2,2)
    plot(velocity)
    title('velocity')

    subplot(2,2,3)
    plot(force_e)
    title('force e')

    subplot(2,2,4)
    plot(force_d)
    title('force d')
    
