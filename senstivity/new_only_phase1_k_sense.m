clear; close all;
n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:1;

% params
R = 0.1;
m = 1;

m = 1;
k=1;
damping_frac = 0.8;

ord = 2;
ord1=2;
% initial vals
x = zeros(n_particles,length(t_durn));
y = zeros(n_particles,length(t_durn));

vx = zeros(n_particles,length(t_durn));
vy = zeros(n_particles,length(t_durn));
vrn_vec = zeros(1,length(t_durn));


force_e_x = zeros(1,length(t_durn));
force_e_y = zeros(1,length(t_durn));

force_d_x = zeros(1,length(t_durn));
force_d_y = zeros(1,length(t_durn));

force_total_x = zeros(1, length(t_durn));
force_total_y = zeros(1, length(t_durn));
% intialise
x(1,1) = rand;
y(1,1) = rand;

target_distance = 0.95 * 2 * R;
theta = 2 * pi * rand();

x(2,1) = x(1,1) + target_distance * cos(theta);
y(2,1) = y(1,1) + target_distance * sin(theta);


while x(2,1) < 0 || x(2,1) > 1 || y(2,1) < 0 || y(2,1) > 1
    theta = 2 * pi * rand();
    x(2,1) = x(1,1) + target_distance * cos(theta);
    y(2,1) = y(1,1) + target_distance * sin(theta);
end

vx(1,1) = rand;
vx(2,1) = rand;

vy(1,1) = rand;
vy(2,1) = rand;

% diff_coef = 3.832e-29;
% atomic_vol = 1.18e-29;
diff_coef = 3.832e-19;
atomic_vol = 1.18e-1;
surface_energy = 1.72;
dihedral_angle = 146*(pi/180); % radians
density = 8920;
kT = 1.38e-23*1e3;

a = zeros(1,length(t_durn));
    t=2;

    a(t-1) = 0.7*rand;
    d = sqrt( (x(1,t-1) - x(2,t-1))^2  +  (y(1,t-1) - y(2,t-1))^2 );
    force_e = k*(d^ord);
    dx = ( x(1,t-1) - x(2,t-1) );
    dy = ( y(1,t-1) - y(2,t-1) );
     force_e_x(t-1) = force_e*(abs(dx)/d);
    force_e_y(t-1) = force_e*(abs(dy)/d);

    v1x = vx(1,t-1);
    v1y = vy(1,t-1);

    v2x = vx(2,t-1);
    v2y = vy(2,t-1);

    unit_vec = [-dx, -dy]./d;
    v1n = v1x*(unit_vec(1)) + v1y*(unit_vec(2));
    v2n = v2x*(unit_vec(1)) + v2y*((unit_vec(2)));
    vrn = -v1n + v2n;
    vrn_vec(t-1) = vrn;
    force_d = damping_frac*(d^ord1)*sqrt(2*m*k)*vrn;
    force_d_x(t-1) = force_d*(abs(dx)/d);
    force_d_y(t-1) = force_d*(abs(dy)/d);

param_range = linspace(0.1,1.5,10);
all_force_diff_param = [];

for param=param_range
        k =  param;
        for t=2:length(t_durn)
            d = sqrt( (x(1,t-1) - x(2,t-1))^2  +  (y(1,t-1) - y(2,t-1))^2 );
            force_e = k*(d^ord);
            dx = ( x(1,t-1) - x(2,t-1) );
            dy = ( y(1,t-1) - y(2,t-1) );
            force_e_x(t) = force_e*(abs(dx)/d);
            force_e_y(t) = force_e*(abs(dy)/d);
            
            
            v1x = vx(1,t-1);
            v1y = vy(1,t-1);
            
            v2x = vx(2,t-1);
            v2y = vy(2,t-1);
            
            unit_vec = [-dx, -dy]./d;
            v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
            v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
            vrn = -v1n + v2n;
            vrn_vec(t) = vrn;
            force_d = damping_frac*(d^ord1)*sqrt(2*m*k)*vrn;
            force_d_x(t) = force_d*(abs(dx)/d);
            force_d_y(t) = force_d*(abs(dy)/d);
            
            force_total_x(t) = force_e_x(t) + force_d_x(t);
            force_total_y(t) = force_e_y(t) + force_d_y(t);

            for n=1:n_particles
                vx(n,t) =  vx(n,t-1)+ (  (force_e_x(t) + force_d_x(t))/m   )*dt;
                vy(n,t) =  vy(n,t-1) + (  (force_e_y(t) + force_d_y(t))/m   )*dt;
            end
            
            for n=1:n_particles
                x(n,t) = x(n,t-1) + vx(n,t)*dt;
                y(n,t) = y(n,t-1) + vy(n,t)*dt;
            end
            a(t) = sqrt(a(t-1)^2 -  2*R*vrn*dt);
        end
        all_force_diff_param = [all_force_diff_param; sqrt( (force_total_x).^2 + (force_total_y).^2  )];
end % param range


%%
figure
    plot(sqrt(force_e_x.^2 + force_e_y.^2))
    title('Force e')
    xlabel('Time')
    ylabel('Force')


figure
     plot(sqrt(force_d_x.^2 + force_d_y.^2))
    title('Force d')
    xlabel('Time')
    ylabel('Force')

  %%
  figure
    plot( sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 ) )
    title('Distance between particles')
    xlabel('Time')
    ylabel('Distance')
    %%
    figure
        plot(vrn_vec)
        title('Relative Normal Velocity')
        xlabel('time')
        ylabel('velocity')

        %%
        figure
            plot(all_force_diff_param)
            title('all force')
  %%
  figure
    imagesc(all_force_diff_param)
    colorbar
    ylabel('Parameter Range');
yLabels = round(param_range,2);

% Update the Y-axis ticks and labels
ax = gca;
ax.YTick = 1:numel(yLabels);
ax.YTickLabel = yLabels;
xlabel('Time')
title('For different values of k')
