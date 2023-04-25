clear;close all
n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:1;

% params
R = 0.1;
m = 75;

% initial vals
x = zeros(n_particles,length(t_durn));
y = zeros(n_particles,length(t_durn));

vx = zeros(n_particles,length(t_durn));
vy = zeros(n_particles,length(t_durn));
vrn_vec = zeros(1,length(t_durn));

force_v_x = zeros(1,length(t_durn));
force_v_y = zeros(1,length(t_durn));

force_s_x = zeros(1,length(t_durn));
force_s_y = zeros(1,length(t_durn));

all_force_diff_param = [];

a_vec = [];
% intialise
x(1,1) = rand;
y(1,1) = rand;

target_distance = 0.05 * 2 * R;
theta = 2 * pi * rand();

x(2,1) = x(1,1) + target_distance * cos(theta);
y(2,1) = y(1,1) + target_distance * sin(theta);


while x(2,1) < 0 || x(2,1) > 1 || y(2,1) < 0 || y(2,1) > 1
    theta = 2 * pi * rand();
    x(2,1) = x(1,1) + target_distance * cos(theta);
    y(2,1) = y(1,1) + target_distance * sin(theta);
end

rand_v_in_x_dir = rand;
if x(1,1) < x(2,1)
    vx(1,1) = rand_v_in_x_dir;
    vx(2,1) = -rand_v_in_x_dir;
else
    vx(1,1) = -rand_v_in_x_dir;
    vx(2,1) = rand_v_in_x_dir;
end

rand_v_in_y_dir = rand;
if y(1,1) < y(2,1)
    vy(1,1) = rand_v_in_y_dir;
    vy(2,1) = -rand_v_in_y_dir;
else
    vy(1,1) = -rand_v_in_y_dir;
    vy(2,1) = rand_v_in_y_dir;
end


% diff_coef = 3.832e-29;
% atomic_vol = 1.18e-29;
diff_coef = 3.832e-13;
atomic_vol = 1.18e-13;
surface_energy = 1.72;
dihedral_angle = 146*(pi/180); % radians
density = 8920;
kT = 1.38e-23*1e3;

a = zeros(1,length(t_durn));
    t=2;

    a(t-1) = 0.03*R;
    a_vec = [a_vec a(t-1)];
    d = sqrt( (x(1,t-1) - x(2,t-1))^2  +  (y(1,t-1) - y(2,t-1))^2 );
    force_v = (pi*(a(t-1)^4))/(8*((diff_coef*atomic_vol)/(kT)));
    dx = ( x(1,t-1) - x(2,t-1) );
    dy = ( y(1,t-1) - y(2,t-1) );
    force_v_x(t-1) = force_v*(abs(dx)/d);
    force_v_y(t-1) = force_v*(abs(dy)/d);

    v1x = vx(1,t-1);
    v1y = vy(1,t-1);

    v2x = vx(2,t-1);
    v2y = vy(2,t-1);

    unit_vec = [-dx, -dy]./d;
    v1n = v1x*(unit_vec(1)) + v1y*(unit_vec(2));
    v2n = v2x*(unit_vec(1)) + v2y*((unit_vec(2)));
    vrn = -v1n + v2n;
    vrn_vec(t-1) = vrn;
    force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
    force_s_x(t-1) = force_s*(abs(dx)/d);
    force_s_y(t-1) = force_s*(abs(dy)/d);
  
%   param_range = linspace(1,2,10);
  param_range = linspace(144*(pi/180),148*(pi/180),10);

  for param=param_range
      dihedral_angle = param;
      for t=2:length(t_durn)
          d = sqrt( (x(1,t-1) - x(2,t-1))^2  +  (y(1,t-1) - y(2,t-1))^2 );
          force_v = (pi*(a(t-1)^4))/(8*((diff_coef*atomic_vol)/(kT)));
          dx = ( x(1,t-1) - x(2,t-1) );
          dy = ( y(1,t-1) - y(2,t-1) );
          force_v_x(t) = force_v*(abs(dx)/d);
          force_v_y(t) = force_v*(abs(dy)/d);


          v1x = vx(1,t-1);
          v1y = vy(1,t-1);

          v2x = vx(2,t-1);
          v2y = vy(2,t-1);

          unit_vec = [-dx, -dy]./d;
          v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
          v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
          vrn = -v1n + v2n;
          vrn_vec(t) = vrn;
          force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
          force_s_x(t) = force_s*(abs(dx)/d);
          force_s_y(t) = force_s*(abs(dy)/d);

          for n=1:n_particles
              vx(n,t) =  vx(n,t-1) +  (  (force_v_x(t) + force_s_x(t))/m   )*dt;
              vy(n,t) =  vy(n,t-1) + (  (force_v_y(t) + force_s_y(t))/m   )*dt;
          end

          for n=1:n_particles
              x(n,t) = x(n,t-1) + vx(n,t)*dt;
              y(n,t) = y(n,t-1) + vy(n,t)*dt;
          end
          %     if a(t-1)^2 -  2*R*vrn*dt < 0
          %         disp(['f=', num2str(a(t-1)^2),' ' ,num2str(2*R*vrn*dt)])
          %         disp(a(t-1)^2 -  2*R*vrn*dt)
          %     end
          a(t) = sqrt(abs(a(t-1)^2 -  2*R*vrn*dt));
          a_vec = [a_vec a(t)];
          %     if a(t-1)^2 -  2*R*vrn*dt < 0
          %         disp('return')
          %         return
          %     end
      end % t
      all_force_diff_param = [all_force_diff_param; sqrt((force_v_x + force_s_x).^2 + (force_v_y + force_s_y).^2)];
  end % param


%%
figure
    plot(sqrt(force_v_x.^2 + force_v_y.^2))
    title('force v')

figure
     plot(sqrt(force_s_x.^2 + force_s_y.^2))
    title('force s')

  %%
  figure
    plot( sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 ) )
    title('dist')
    %%
    figure
        plot(vrn_vec)
        title('vrn')
        %%
figure
    plot(a_vec)
    title('neck radius')
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
title('For different values of dihedral angle')
