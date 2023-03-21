clear all
close all
n_particles = 2;
t_max = 10;
dt = 0.01;
t_durn = 0:dt:t_max;


% particle params
m = 1;
k=1;
damping_frac = 0.8;
R = 1;

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
figure
hold on
 for ord=-0.5:0.5:3
 
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

    for n=1:n_particles
        vx(n,t) =  vx(n,t-1)+ (  (-1^n)*(force_e_x(t) + force_d_x(t))/m   )*dt;
        vy(n,t) =  vy(n,t-1) + (  (-1^n)*(force_e_y(t) + force_d_y(t))/m   )*dt;
    end
    
    for n=1:n_particles
        x(n,t) = x(n,t-1) + vx(n,t)*dt;
        y(n,t) = y(n,t-1) + vy(n,t)*dt;
    end
 end % t
 plot(sqrt( (force_e_x + force_d_x).^2 + (force_e_y + force_d_y).^2  ))
 
 end% end of ord
    plot( sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 ), 'LineWidth',3 )
    plot(vrn_vec,'LineWidth',3)

    legend([strsplit(num2str(-0.5:0.5:3)), 'dist', 'vrn'])
    ylim([0, 5])
%     xlim([0. 5])
    title('forces')
    hold off

    

 