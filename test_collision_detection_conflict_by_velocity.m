% test_collision_detection_conflict_by_velocity

cpslist = zeros(100, 2, 2);
vcpslist = zeros(100, 2, 2);

% construct cpslist with overlapped trajectories but same velocities

cpslist(:, 1, 1) = linspace(15.1, 30.2, 100);
cpslist(:, 2, 1) = linspace(10.9, 10.9, 100);
vcpslist(:, 1, 1) = 1.57 * cpslist(:, 1, 1) + 0.2 * cpslist(:, 2, 1) + 3.1; 
vcpslist(:, 2, 1) = -0.4 * cpslist(:, 1, 1) + 1.2 * cpslist(:, 2, 1) + 2.7; 

cpslist(:, 1, 2) = linspace(22.3, 36.9, 100);
cpslist(:, 2, 2) = linspace(10.9, 10.9, 100);
vcpslist(:, 1, 2) = 1.57 * cpslist(:, 1, 2) + 0.2 * cpslist(:, 2, 2) + 3.1; 
vcpslist(:, 2, 2) = -0.4 * cpslist(:, 1, 2) + 1.2 * cpslist(:, 2, 2) + 2.7; 



g.nb_T = 100;
g.nb_cps = 2;
g.r = 0.25;
g.dv = 2; % difference of velocity

clist = collision_detection_conflicted_velocity(g, cpslist, vcpslist)

vcpslist = vcpslist * 0.001;

figure; clf;
hold on;
plot(cpslist(:, 1, 1), cpslist(:, 2, 1), 'b.');
quiver(cpslist(:, 1, 1), cpslist(:, 2, 1), vcpslist(:, 1, 1), vcpslist(:, 2, 1), 0, 'b');

plot(cpslist(:, 1, 2), cpslist(:, 2, 2), 'r.');
quiver(cpslist(:, 1, 2), cpslist(:, 2, 2), vcpslist(:, 1, 2), vcpslist(:, 2, 2), 0, 'r');

hold off;

axis equal;


% construct cpslist with overlapped trajectories but different velocities

clist = collision_detection(g, cpslist)