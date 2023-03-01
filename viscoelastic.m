%%
N = 5;  
R = 1; 
k = 1e4; 
eta = 0.1; % Viscosity coeff
T = 100; 
dt = 0.01; 

x = rand(N,1); 
y = rand(N,1); 
vx = zeros(N,1); 
vy = zeros(N,1); 

% spring connections
d = pdist2([x,y],[x,y]); % Pairwise distances between particles
d(d==0) = inf; % self-interactions

I = find(d < 2*R); 
[I1,I2] = ind2sub([N,N],I); 
L = d(I); % dist btn centers for particles in contact 
dx = (x(I1) - x(I2)) ./ L; % x-components of connection unit vectors
dy = (y(I1) - y(I2)) ./ L; % y-components of connection unit vectors
fx = zeros(N,1); % Zero x-forces
fy = zeros(N,1); % Zero y-forces

%%
% Main simulation loop
for t = 0:dt:T
    % Calculate spring forces
    F = k * (L - 2*R);
    Fx = F .* dx;
    Fy = F .* dy;
    
    % Calculate viscoelastic forces
    Fv = eta * (vx(I1)-vx(I2)) + eta * (vy(I1)-vy(I2));
    Fvx = Fv .* dx;
    Fvy = Fv .* dy;
    
    % Update particle velocities
    fx(:) = 0; % Reset forces
    fy(:) = 0;
    for i = 1:length(I)
        fx(I1(i)) = fx(I1(i)) + Fx(i) + Fvx(i);
        fy(I1(i)) = fy(I1(i)) + Fy(i) + Fvy(i);
        fx(I2(i)) = fx(I2(i)) - Fx(i) - Fvx(i);
        fy(I2(i)) = fy(I2(i)) - Fy(i) - Fvy(i);
    end
    vx = vx + fx*dt;
    vy = vy + fy*dt;
    
    % Update particle positions
    x = x + vx*dt;
    y = y + vy*dt;
    
    % Update spring connections
    d = pdist2([x,y],[x,y]); % Pairwise distances between particles
    d(d==0) = inf; % Exclude self-interactions
    I = find(d < 2*R); % Indices of connected pairs
    [I1,I2] = ind2sub([N,N],I); % Convert linear indices to subscripts
    L = d(I); % Lengths of connections
    dx = (x(I1) - x(I2)) ./ L; % x-components of connection unit vectors
    dy = (y(I1) - y(I2)) ./ L; % y-components of connection unit vectors
    
end 