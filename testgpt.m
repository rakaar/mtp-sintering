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

%%
% Define your data as 2D matrices or vectors
x = linspace(0, 2*pi, 100);
y1 = sin(x);
y2 = cos(x);
y3 = sin(2*x);
y4 = cos(2*x);
y5 = sin(3*x);
y6 = cos(3*x);
y7 = sin(4*x);
y8 = cos(4*x);

% Combine data into a single matrix
Y = [y1; y2; y3; y4; y5; y6; y7; y8];

% Define custom colors (you can choose any RGB values you like)
colors = [
    0, 0, 1;  % blue
    0, 1, 0;  % green
    1, 0, 0;  % red
    0, 1, 1;  % cyan
    1, 0, 1;  % magenta
    0.5, 0.5, 0.5;  % gray
    1, 0.5, 0;  % orange
    0.5, 0, 0.5;  % purple
];

% Create a new figure
figure;

% Plot the lines using a loop
for i = 1:8
    plot(x, Y(i, :), 'Color', colors(i, :), 'LineWidth', 1.5);
    hold on;
end

% Add a legend or other plot elements if desired
legend('Line 1', 'Line 2', 'Line 3', 'Line 4', 'Line 5', 'Line 6', 'Line 7', 'Line 8');
xlabel('x-axis');
ylabel('y-axis');
title('Eight Lines with Different Colors');

% Hold off to stop adding more plots to the figure
hold off;
