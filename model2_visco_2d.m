    n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:1;

% params
R = 10;
m = 1;

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

% intialise
x(1,1) = rand;
x(2,1) = rand;

y(1,1) = rand;
y(2,1) = rand;

d = sqrt( (x(1,1) - x(2,1))^2  +  (y(1,1) - y(2,1))^2 );
if d >= 2*R
    while d > 2*R
        x(1,1) = rand;
        x(2,1) = rand;
        
        y(1,1) = rand;
        y(2,1) = rand;
        
        d = sqrt( (x(1,1) - x(2,1))^2  +  (y(1,1) - y(2,1))^2 );
    end
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
    vrn = v1n - v2n;
    vrn_vec(t-1) = vrn;
    force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
    force_s_x(t-1) = force_s*(abs(dx)/d);
    force_s_y(t-1) = force_s*(abs(dy)/d);
    
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
    vrn = v1n - v2n;
    vrn_vec(t) = vrn;
    force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
    force_s_x(t) = force_s*(abs(dx)/d);
    force_s_y(t) = force_s*(abs(dy)/d);

    for n=1:n_particles
        vx(n,t) =  vx(n,t-1)+ (  (force_v_x(t) + force_s_x(t))/m   )*dt;
        vy(n,t) =  vy(n,t-1) + (  (force_v_y(t) + force_s_y(t))/m   )*dt;
    end
    
    for n=1:n_particles
        x(n,t) = x(n,t-1) + vx(n,t)*dt;
        y(n,t) = y(n,t-1) + vy(n,t)*dt;
    end
    a(t) = sqrt(a(t-1)^2 -  2*R*vrn*dt);
end

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
  