x_data = csvread("x_data.csv");
y_data = csvread("y_data.csv");
z_data = csvread("z_data.csv");
xx = x_data(:,1);
xy = x_data(:,2);
yx = y_data(:,1);
yy = y_data(:,2);
zx = z_data(:,1);
zy = z_data(:,2);
figure(1);
plot(xx,xy,'-.r',yx,yy,'b*',zx,zy,'g:','LineWidth',3)
grid on;
xlabel('Time(s)');
ylabel('Position(m)');
title('Optimal trajectory for X Y and Z')
legend('X-position','Y-position','Z-position','Location','southeast')
%%
% xy = csvread("xy_traj.csv");
% x = xy(:,1);
% y = xy(:,2);
% N = length(x);
% z = linspace(2, 2.5, N)';
% 
% start_point = [x(1), y(1), z(1)];
% end_point = [x(end), y(end), z(end)];
% plot3(x, y, z,'o');;
yy_new = interp1(yx,yy,xx);
plot3(xy,yy_new,zy,'-.b','LineWidth',3)
hold on;
plot3(xy(1),yy_new(1),zy(1), 'ro','LineWidth',4)
plot3(xy(end),yy_new(end),zy(end), 'ro','LineWidth',4)
xlabel('X-posiiton(m)');
ylabel('Y-posiiton(m)');
zlabel('Z-posiiton(m)');
grid on;
title('Trajectory in 3D space')
text(xy(1),yy_new(1) + 1,zy(1), 'Initial Position', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'red');
text(xy(end),yy_new(end) - 0.5,zy(end), 'Terminal Position', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'red');

view(3)