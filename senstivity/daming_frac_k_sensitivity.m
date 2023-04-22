clear all
% close all
n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:t_max;


% particle params
m = 1;
k=0.75;
damping_frac = 0.9;
R = 10;

%
diff_coef = 3.832e-10;
atomic_vol = 1.18e-7;
surface_energy = 2.5*0.01;
dihedral_angle = 146*(pi/180); % radians
density = 8920;
kT = 1.38e-23*1e3;

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


force_v_x = zeros(1,length(t_durn));
force_v_y = zeros(1,length(t_durn));

force_s_x = zeros(1,length(t_durn));
force_s_y = zeros(1,length(t_durn));

force_total_x = zeros(1,length(t_durn));
force_total_y = zeros(1,length(t_durn));

a = zeros(1,length(t_durn));


forces_all_order = zeros(length(-0.5:0.5:3), length(t_durn));
e_forces_all_order = zeros(length(-0.5:0.5:3), length(t_durn));
d_forces_all_order = zeros(length(-0.5:0.5:3), length(t_durn));
v_forces_all_order = zeros(length(-0.5:0.5:3), length(t_durn));
s_forces_all_order = zeros(length(-0.5:0.5:3), length(t_durn));


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
% initialize
t=2; ord = 1;
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

unit_vec = [dx, dy]./d;
v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
vrn = v1n - v2n;
vrn_vec(t-1) = vrn;

%     force_d = damping_frac*sqrt(2*m*k)*vrn;
% according to other paper, so that penetration power can be kept at both
% terms
force_d = damping_frac*(d^ord)*sqrt(2*m*k)*vrn;
force_d_x(t-1) = force_d*(abs(dx)/d);
force_d_y(t-1) = force_d*(abs(dy)/d);

force_total_x(t-1) = force_e_x(t-1) + force_d_x(t-1);
force_total_y(t-1) = force_e_y(t-1) + force_d_y(t-1);


ord = 1;
% sensitivity
param_range = linspace(0.1,0.9,10);
param_range2 = linspace(0.1,0.9,10);

force_at_single_t = zeros(length(param_range), length(param_range2));

for param=1:length(param_range)
    for param2=1:length(param_range2)
        k = param_range(param);
        damping_frac = param_range2(param2);
        t_regime_change = 0;
        for t=2:length(t_durn)

            d = sqrt( (x(1,t-1) - x(2,t-1))^2  +  (y(1,t-1) - y(2,t-1))^2 );

            if d < R/5
                force_e = k*(d^ord);
                dx = ( x(1,t-1) - x(2,t-1) );
                dy = ( y(1,t-1) - y(2,t-1) );
                force_e_x(t) = force_e*(abs(dx)/d);
                force_e_y(t) = force_e*(abs(dy)/d);


                v1x = vx(1,t-1);
                v1y = vy(1,t-1);

                v2x = vx(2,t-1);
                v2y = vy(2,t-1);

                unit_vec = [dx, dy]./d;
                v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
                v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
                vrn = v1n - v2n;
                vrn_vec(t) = vrn;
                %     force_d = damping_frac*sqrt(2*m*k)*vrn;
                % according to other paper, so that penetration power can be kept at both
                % terms
                force_d = damping_frac*(d^ord)*sqrt(2*m*k)*vrn;
                force_d_x(t) = force_d*(abs(dx)/d);
                force_d_y(t) = force_d*(abs(dy)/d);

                force_total_x(t) = force_e_x(t) + force_d_x(t);
                force_total_y(t) = force_e_y(t) + force_d_y(t);

                e_forces_all_order(ord*2 + 2,:) = sqrt( (force_e_x).^2 + (force_e_y).^2  );
                d_forces_all_order(ord*2 + 2,:) = sqrt( (force_d_x).^2 + (force_d_y).^2  );
                t_regime_change = t;
            else % d > r/2 else
                if a(t-1) == 0
                    a(t-1) = 1.5*rand;
                end

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
                a(t) = sqrt(a(t-1)^2 -  2*R*vrn*dt);

                force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
                force_s_x(t) = force_s*(abs(dx)/d);
                force_s_y(t) = force_s*(abs(dy)/d);


                force_total_x(t) = force_v_x(t) + force_s_x(t);
                force_total_y(t) = force_v_y(t) + force_s_y(t);


                v_forces_all_order(ord*2 + 2,:) = sqrt( (force_v_x).^2 + (force_v_y).^2  );
                s_forces_all_order(ord*2 + 2,:) = sqrt( (force_s_x).^2 + (force_s_y).^2  );


            end % d > r/2

            for n=1:n_particles
                vx(n,t) =  vx(n,t-1)+ (  (-1^n)*(force_total_x(t))/m   )*dt;
                vy(n,t) =  vy(n,t-1) + (  (-1^n)*(force_total_y(t))/m   )*dt;
            end

            for n=1:n_particles
                x(n,t) = x(n,t-1) + vx(n,t)*dt;
                y(n,t) = y(n,t-1) + vy(n,t)*dt;
            end
        end % t

        force_at_single_t(param, param2) = sqrt(force_total_x(t_regime_change)^2 + force_total_y(t_regime_change)^2);
    end
   
    
 
end % end of param



figure
    imagesc(force_at_single_t)
    title('c k change')
    xlabel('k', 'FontSize',20)
    ylabel('c', 'FontSize',20)




