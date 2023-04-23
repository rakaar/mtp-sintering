clear; close all
tic
n_particles = 9;
[rows, cols] = ind2sub([sqrt(n_particles) sqrt(n_particles)], linspace(1,n_particles, n_particles));
t_max = 10;
dt = 0.01;
t_durn = 0:dt:1;


% params
m = 1;
k=0.2;
damping_frac = 0.2;

R = 10;
m = 1;

% initial vals
x = zeros(n_particles,length(t_durn));
y = zeros(n_particles,length(t_durn));

vx = zeros(n_particles,length(t_durn));
vy = zeros(n_particles,length(t_durn));

force_total_x = zeros(n_particles, n_particles, length(t_durn));
force_total_y = zeros(n_particles, n_particles, length(t_durn));

vrn_vec = zeros(n_particles, n_particles,length(t_durn));

force_v_x = zeros(n_particles, n_particles,length(t_durn));
force_v_y = zeros(n_particles, n_particles,length(t_durn));

force_s_x = zeros(n_particles, n_particles,length(t_durn));
force_s_y = zeros(n_particles, n_particles,length(t_durn));

force_d_x = zeros(n_particles, n_particles,length(t_durn));
force_d_y = zeros(n_particles, n_particles,length(t_durn));

force_e_x = zeros(n_particles, n_particles,length(t_durn));
force_e_y = zeros(n_particles, n_particles,length(t_durn));

d_over_time = zeros(n_particles, n_particles, length(t_durn));

sep = 2 * R;
counter = 1;
for i = 1:sqrt(n_particles)
    for j = 1:sqrt(n_particles)
        x(counter,1) = (i - 1) * sep + R;
        y(counter,1) = (j - 1) * sep + R;
        counter = counter + 1;
    end
end


for n=1:n_particles
    vx(n,1) = rand;
    vx(n,1) = rand;

    vy(n,1) = rand;
    vy(n,1) = rand;
end % n

% diff_coef = 3.832e-29;
% atomic_vol = 1.18e-29;
diff_coef = 3.832e-10;
atomic_vol = 1.18e-1;
surface_energy = 1.72;
dihedral_angle = 146*(pi/180); % radians
density = 8920;
kT = 1.38e-23*1e3;

a = zeros(n_particles, n_particles ,length(t_durn));
t=2;


t_second_phase = -1;
ord = 1;
for n1=1:n_particles
    for n2=1:n_particles
        % no if needed for 4
        if n1 == n2
            continue
        end
        n1_row = rows(n1);
        n1_col = cols(n1);

        n2_row = rows(n2);
        n2_col = cols(n2);

        if ~( (abs(n1_row - n2_row) == 1 && abs(n1_col - n2_col) == 0)  || (abs(n1_row - n2_row) == 0 && abs(n1_col - n2_col) == 1)  || (abs(n1_row - n2_row) == 1 && abs(n1_col - n2_col) == 1))
            continue
        end
        d = sqrt( (x(n1,t-1) - x(n2,t-1))^2  +  (y(n1,t-1) - y(n2,t-1))^2 );
        d_over_time(n1,n2,t-1) = d;
        force_e = k*(d^ord);
        dx = -( x(n1,t-1) - x(n2,t-1) );
        dy = -( y(n1,t-1) - y(n2,t-1) );
        force_e_x(n1,n2,t-1) = force_e*(abs(dx)/d);
        force_e_y(n1,n2,t-1) = force_e*(abs(dy)/d);

        v1x = vx(n1,t-1);
        v1y = vy(n1,t-1);

        v2x = vx(n2,t-1);
        v2y = vy(n2,t-1);

        unit_vec = [dx, dy]./d;
        v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
        v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
        vrn = -(v1n - v2n);
        vrn_vec(n1,n2,t-1) = vrn;
        %     force_d = damping_frac*sqrt(2*m*k)*vrn;
        % according to other paper, so that penetration power can be kept at both
        % terms
        force_d = damping_frac*(d^ord)*sqrt(2*m*k)*vrn;
        force_d_x(n1,n2,t-1) = force_d*(abs(dx)/d);
        force_d_y(n1,n2,t-1) = force_d*(abs(dy)/d);

        force_total_x(n1,n2,t-1) = force_e_x(n1,n2,t-1) + force_d_x(n1,n2,t-1);
        force_total_y(n1,n2,t-1) = force_e_y(n1,n2,t-1) + force_d_y(n1,n2,t-1);

    end
end

 
for t=2:length(t_durn)
    for n1=1:n_particles
        for n2=1:n_particles
            % no if needed for 4
            if n1 == n2 || n1 < n2
                continue
            end
            n1_row = rows(n1);
            n1_col = cols(n1);

            n2_row = rows(n2);
            n2_col = cols(n2);

            if ~( (abs(n1_row - n2_row) == 1 && abs(n1_col - n2_col) == 0)  || (abs(n1_row - n2_row) == 0 && abs(n1_col - n2_col) == 1)  || (abs(n1_row - n2_row) == 1 && abs(n1_col - n2_col) == 1))
                continue
            end
             d = sqrt( (x(n1,t-1) - x(n2,t-1))^2  +  (y(n1,t-1) - y(n2,t-1))^2 );
             d_over_time(n1,n2,t-1) = d;
            % ---------- both phases with if statement
            if d > 0.05*2*R
                force_e = k*(d^ord);
                dx = -( x(n1,t-1) - x(n2,t-1) );
                dy = -( y(n1,t-1) - y(n2,t-1) );
                force_e_x(n1,n2,t) = force_e*(abs(dx)/d);
                force_e_y(n1,n2,t) = force_e*(abs(dy)/d);

                v1x = vx(n1,t-1);
                v1y = vy(n1,t-1);

                v2x = vx(n2,t-1);
                v2y = vy(n2,t-1);

                unit_vec = [dx, dy]./d;
                v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
                v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
                vrn = -(v1n - v2n);
                vrn_vec(n1,n2,t) = vrn;
                %     force_d = damping_frac*sqrt(2*m*k)*vrn;
                % according to other paper, so that penetration power can be kept at both
                % terms
                force_d = damping_frac*(d^ord)*sqrt(2*m*k)*vrn;
                force_d_x(n1,n2,t) = force_d*(abs(dx)/d);
                force_d_y(n1,n2,t) = force_d*(abs(dy)/d);

                force_total_x(n1,n2,t) = force_e_x(n1,n2,t) + force_d_x(n1,n2,t);
                force_total_y(n1,n2,t) = force_e_y(n1,n2,t) + force_d_y(n1,n2,t);

            else % d > r/2 else
                if a(t-1) == 0
                    a(t-1) = 3*rand;
                end
                if t_second_phase == -1
                    t_second_phase = t;
                end

                force_v = (pi*(a(t-1)^4))/(8*((diff_coef*atomic_vol)/(kT)));
                dx = -( x(n1,t-1) - x(n2,t-1) );
                dy = -( y(n1,t-1) - y(n2,t-1) );
                force_v_x(n1,n2,t) = force_v*(abs(dx)/d);
                force_v_y(n1,n2,t) = force_v*(abs(dy)/d);


                v1x = vx(n1,t-1);
                v1y = vy(n1,t-1);

                v2x = vx(n2,t-1);
                v2y = vy(n2,t-1);

                unit_vec = [dx, dy]./d;
                v1n = v1x*unit_vec(1) + v1y*unit_vec(2);
                v2n = v2x*unit_vec(1) + v2y*unit_vec(2);
                vrn = -(v1n - v2n);
                vrn_vec(n1,n2,t) = vrn;
                a(n1,n2,t) = sqrt(a(n1,n2,t-1)^2 -  2*R*vrn*dt);

                force_s = pi*(surface_energy)*(4*R*( 1 - 0.5*cos(dihedral_angle/2) ) + a(t-1)*sin(dihedral_angle/2) );
                force_s_x(n1,n2,t) = force_s*(abs(dx)/d);
                force_s_y(n1,n2,t) = force_s*(abs(dy)/d);


                force_total_x(n1,n2,t) = force_v_x(n1,n2,t) + force_s_x(n1,n2,t);
                force_total_y(n1,n2,t) = force_v_y(n1,n2,t) + force_s_y(n1,n2,t);
            end % d > r/2

            % ---------- end of both if

        end
    end

    for n=1:n_particles
        vx(n,t) =  vx(n,t-1)+ (  sum(force_total_x(n,:,t)./m) )*dt;
        vy(n,t) =  vy(n,t-1) + (  sum(force_total_y(n,:,t)./m) )*dt;
    end

    for n=1:n_particles
        x(n,t) = x(n,t-1) + vx(n,t)*dt;
        y(n,t) = y(n,t-1) + vy(n,t)*dt;
    end

end

toc
disp(['t 2nd phase starts at = ', num2str(t_second_phase)])
%%
for t=1:length(t_durn)
    imagesc( (squeeze(force_total_x(:,:,t)).^2 + squeeze(force_total_y(:,:,t)).^2).^0.5 )
%     caxis([0 ]);
    colorbar()
    title(num2str(t))
    pause(0.50)  
end

%%
for t=1:length(t_durn)
     imagesc( squeeze(d_over_time(:,:,t)) );
     title(num2str(t))
     colorbar()
     pause(0.5)
end