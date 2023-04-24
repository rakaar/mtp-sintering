% Define R
R = 0.1; % You can replace this with any desired value within the range

% Calculate the target distance as > 90% of 2 times R, but < 2 times R
target_distance = 0.95 * 2 * R;

% Initialize the first point (x1, y1) within the 0 to 1 range
x1 = rand();
y1 = rand();

% Calculate the angle between the two points in radians


% Calculate the second point (x2, y2) using the target distance and angle
x2 = x1 + target_distance * cos(theta);
y2 = y1 + target_distance * sin(theta);

% Check if the points are within the 0 to 1 range, if not, try again
while x2 < 0 || x2 > 1 || y2 < 0 || y2 > 1
    theta = 2 * pi * rand();
    x2 = x1 + target_distance * cos(theta);
    y2 = y1 + target_distance * sin(theta);
end

% Print the two points
disp(['Point 1: (' num2str(x1) ', ' num2str(y1) ')'])
disp(['Point 2: (' num2str(x2) ', ' num2str(y2) ')'])
