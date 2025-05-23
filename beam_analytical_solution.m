% Analytical solution for cantilever beam with moment and shear force
% This script generates analytical solution plots for a cantilever beam
% with moment and shear force boundary conditions

clear all;
close all;
clc;

% Create output directory if it doesn't exist
if ~exist('./output', 'dir')
    mkdir('./output');
end

%% Parameters from pre_process.m
% Material properties
E = 2.1e11;  % Young's modulus (Pa)
nu = 0.3;    % Poisson's ratio
G = E/(2*(1+nu));  % Shear modulus

% Beam dimensions
L = 0.1;      % Length (m)
h = 0.02;       % Height (m) - from -1 to 1
t = 0.01;       % Thickness (m)

% 选择载荷类型: 1=剪力, 2=力矩
load_type = 2;  % 修改此值以切换载荷类型

% Boundary conditions (from pre_process.m)
if load_type == 1
    % 仅施加剪力
    P = 1e6;     % Shear force (Pa) - from pre_process.m
    M = 0;       % No moment
    load_name = 'Shear Force';
    fprintf('仅施加剪力 P = %g Pa\n', P);
else
    % 仅施加力矩
    P = 0;       % No shear force
    M = -20;      % Moment (N·m) - reduced for better visualization
    load_name = 'Moment';
    fprintf('仅施加力矩 M = %g N·m\n', M);
end

% Moment of inertia
I = t*h^3/12;  % Second moment of area (m^4)

%% Create mesh for visualization
nx = 100;  % Number of points in x direction
ny = 20;   % Number of points in y direction

x = linspace(0, L, nx);
y = linspace(-h/2, h/2, ny);
[X, Y] = meshgrid(x, y);

%% Analytical solutions
% Displacement field for a cantilever beam with end moment and shear force
% v = -M/(E*I)*(x^2 + nu*y^2) - P/(6*E*I)*(3*L*x^2 - x^3) - nu*P/(2*E*I)*y^2 - P*h^2/(8*I*G)*x

% Calculate displacement components
U1 = zeros(size(X));  % u displacement (x-direction)
U2 = zeros(size(X));  % v displacement (y-direction)

for i = 1:nx
    for j = 1:ny
        % Current coordinates
        xi = X(j,i);
        yi = Y(j,i);

        % Displacement in x-direction
        % For moment M: u = M/(EI)·xy
        % For force P: u = P/(2EI)·(2Lx-x²)y - νP/(6EI)·y³ + P/(6IG)·y³
        U1(j,i) = M/(E*I)*xi*yi + P/(2*E*I)*(2*L*xi-xi^2)*yi - nu*P/(6*E*I)*yi^3 + P/(6*I*G)*yi^3;

        % Displacement in y-direction (main component)
        % For moment M: v = -M/(2EI)·(x² + νy²)
        % For force P: v = -P/(6EI)·(3Lx²-x³) - νP/(2EI)·(L-x)y² - Ph²/(8IG)·x
        U2(j,i) = -M/(2*E*I)*(xi^2 + nu*yi^2) - P/(6*E*I)*(3*L*xi^2 - xi^3) - nu*P/(2*E*I)*(L-xi)*yi^2 - P*h^2/(8*I*G)*xi;
    end
end

% Calculate total displacement magnitude
U = sqrt(U1.^2 + U2.^2);

%% Calculate stress components
S11 = zeros(size(X));  % σxx
S22 = zeros(size(X));  % σyy
S12 = zeros(size(X));  % τxy

for i = 1:nx
    for j = 1:ny
        % Current coordinates
        xi = X(j,i);
        yi = Y(j,i);

        % Normal stress in x-direction
        % For moment M: σₓ = M/I·y
        % For force P: σₓ = P/I·(L-x)y
        S11(j,i) = M*yi/I + P*(L-xi)*yi/I;

        % Normal stress in y-direction (zero for beam theory)
        % For both moment M and force P: σᵧ = 0
        S22(j,i) = 0;

        % Shear stress
        % For moment M: τₓᵧ = 0
        % For force P: τₓᵧ = -P/(2I)·[(h/2)² - y²]
        S12(j,i) = -P/(2*I)*((h/2)^2 - yi^2);
    end
end

% Calculate principal stresses
S1 = zeros(size(X));
S2 = zeros(size(X));
Mises = zeros(size(X));

for i = 1:nx
    for j = 1:ny
        % Average stress
        Savg = (S11(j,i) + S22(j,i))/2;

        % Radius of Mohr's circle
        R = sqrt((S11(j,i) - S22(j,i))^2/4 + S12(j,i)^2);

        % Principal stresses
        S1(j,i) = Savg + R;
        S2(j,i) = Savg - R;

        % Von Mises stress
        Mises(j,i) = sqrt(S1(j,i)^2 + S2(j,i)^2 - S1(j,i)*S2(j,i));
    end
end

%% Plot displacement field
figure(1);
set(gcf, 'Position', [100, 100, 1000, 800]);

% Plot U1 (x-displacement)
subplot(2,2,1);
contourf(X, Y, U1, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('U1 (x-displacement)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum displacement
[maxU1, maxU1Idx] = max(abs(U1(:)));
[maxU1_y, maxU1_x] = ind2sub(size(U1), maxU1Idx);
hold on;
plot(X(maxU1_y, maxU1_x), Y(maxU1_y, maxU1_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxU1_y, maxU1_x), Y(maxU1_y, maxU1_x), [' max = ', num2str(U1(maxU1_y, maxU1_x), '%.3e')], 'FontSize', 8);
fprintf('U1最大值: %g (位置: x=%g, y=%g)\n', U1(maxU1_y, maxU1_x), X(maxU1_y, maxU1_x), Y(maxU1_y, maxU1_x));

% Plot U2 (y-displacement)
subplot(2,2,2);
contourf(X, Y, U2, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('U2 (y-displacement)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum displacement
[maxU2, maxU2Idx] = max(abs(U2(:)));
[maxU2_y, maxU2_x] = ind2sub(size(U2), maxU2Idx);
hold on;
plot(X(maxU2_y, maxU2_x), Y(maxU2_y, maxU2_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxU2_y, maxU2_x), Y(maxU2_y, maxU2_x), [' max = ', num2str(U2(maxU2_y, maxU2_x), '%.3e')], 'FontSize', 8);
fprintf('U2最大值: %g (位置: x=%g, y=%g)\n', U2(maxU2_y, maxU2_x), X(maxU2_y, maxU2_x), Y(maxU2_y, maxU2_x));

% Plot total displacement
subplot(2,2,3);
contourf(X, Y, U, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('U (total displacement)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum displacement
[maxU, maxUIdx] = max(U(:));
[maxU_y, maxU_x] = ind2sub(size(U), maxUIdx);
hold on;
plot(X(maxU_y, maxU_x), Y(maxU_y, maxU_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxU_y, maxU_x), Y(maxU_y, maxU_x), [' max = ', num2str(maxU, '%.3e')], 'FontSize', 8);
fprintf('总位移U最大值: %g (位置: x=%g, y=%g)\n', maxU, X(maxU_y, maxU_x), Y(maxU_y, maxU_x));

% Plot deformed shape
subplot(2,2,4);
% Scale factor for visualization - more reasonable value
scale = 500;  % Adjust this value to make deformation visible but not distorted

% Create more points along the beam length for smoother visualization
x_beam = linspace(0, L, 100);
y_top = ones(size(x_beam)) * h/2;
y_bottom = ones(size(x_beam)) * (-h/2);

% Calculate analytical displacement for these points
u1_top = zeros(size(x_beam));
u2_top = zeros(size(x_beam));
u1_bottom = zeros(size(x_beam));
u2_bottom = zeros(size(x_beam));

% For pure moment loading (M), the beam deforms into a circular arc
% The radius of curvature is R = EI/M
if load_type == 2  % Pure moment case
    % Calculate radius of curvature
    R = E*I/M;

    % Calculate the angle subtended by the beam
    theta_max = L/R;

    % Generate points along the circular arc
    theta = linspace(0, theta_max, 100);

    % Calculate coordinates of the neutral axis (circular arc)
    x_neutral = R*sin(theta);
    y_neutral = R*(1-cos(theta));

    % Calculate coordinates of top and bottom edges
    % Rotate the h/2 distance around the neutral axis
    x_top_def = x_neutral - (h/2)*sin(theta);
    y_top_def = y_neutral + (h/2)*cos(theta);

    x_bottom_def = x_neutral + (h/2)*sin(theta);
    y_bottom_def = y_neutral - (h/2)*cos(theta);
else  % Shear force case
    for i = 1:length(x_beam)
        xi = x_beam(i);
        yi_top = h/2;
        yi_bottom = -h/2;

        % Displacement for top edge
        % For force P: u = P/(2EI)·(2Lx-x²)y - νP/(6EI)·y³ + P/(6IG)·y³
        u1_top(i) = P/(2*E*I)*(2*L*xi-xi^2)*yi_top - nu*P/(6*E*I)*yi_top^3 + P/(6*I*G)*yi_top^3;

        % For force P: v = -P/(6EI)·(3Lx²-x³) - νP/(2EI)·(L-x)y² - Ph²/(8IG)·x
        u2_top(i) = -P/(6*E*I)*(3*L*xi^2 - xi^3) - nu*P/(2*E*I)*(L-xi)*yi_top^2 - P*h^2/(8*I*G)*xi;

        % Displacement for bottom edge
        u1_bottom(i) = P/(2*E*I)*(2*L*xi-xi^2)*yi_bottom - nu*P/(6*E*I)*yi_bottom^3 + P/(6*I*G)*yi_bottom^3;
        u2_bottom(i) = -P/(6*E*I)*(3*L*xi^2 - xi^3) - nu*P/(2*E*I)*(L-xi)*yi_bottom^2 - P*h^2/(8*I*G)*xi;
    end

    % Calculate deformed coordinates
    x_top_def = x_beam + scale * u1_top;
    y_top_def = y_top + scale * u2_top;
    x_bottom_def = x_beam + scale * u1_bottom;
    y_bottom_def = y_bottom + scale * u2_bottom;
end

% Plot original shape
plot([0 L L 0 0], [-h/2 -h/2 h/2 h/2 -h/2], 'k--', 'LineWidth', 1.5);
hold on;

% Plot deformed shape - top and bottom edges
plot(x_top_def, y_top_def, 'b-', 'LineWidth', 1.5);
plot(x_bottom_def, y_bottom_def, 'b-', 'LineWidth', 1.5);

% Connect the ends with vertical lines
plot([x_top_def(1), x_bottom_def(1)], [y_top_def(1), y_bottom_def(1)], 'b-', 'LineWidth', 1.5); % Left end (fixed)
plot([x_top_def(end), x_bottom_def(end)], [y_top_def(end), y_bottom_def(end)], 'b-', 'LineWidth', 1.5); % Right end (free)

% Add some intermediate vertical lines for better visualization
for i = [25, 50, 75]
    plot([x_top_def(i), x_bottom_def(i)], [y_top_def(i), y_bottom_def(i)], 'b-', 'LineWidth', 0.5);
end

title('Deformed shape (scale factor = 500)');
xlabel('x (m)');
ylabel('y (m)');

% Set reasonable axis limits to focus on the beam
if load_type == 2  % Pure moment case
    % For moment loading, adjust limits to show the circular arc
    max_x = max(max(x_top_def), max(x_bottom_def));
    max_y = max(max(y_top_def), max(y_bottom_def));
    min_x = min(min(x_top_def), min(x_bottom_def));
    min_y = min(min(y_top_def), min(y_bottom_def));

    % Add some padding
    xlim([min_x-0.01, max_x+0.01]);
    ylim([min_y-0.01, max_y+0.01]);
else
    xlim([-0.5, L+0.5]);
    ylim([-h, h]);
end

axis equal;
grid on;
legend('Original', 'Deformed', 'Location', 'best');

% Add text to indicate fixed end and load
text(-0.3, 0, 'Fixed', 'FontSize', 10, 'HorizontalAlignment', 'right');
text(L+0.3, 0, sprintf('P = %g Pa\nM = %g N·m', P, M), 'FontSize', 10, 'HorizontalAlignment', 'left');

% Save displacement figure with load type in filename
if load_type == 1
    saveas(gcf, './output/analytical_displacement_shear.png');
else
    saveas(gcf, './output/analytical_displacement_moment.png');
end

%% Create a separate figure for clearer deformation visualization
figure(3);
set(gcf, 'Position', [100, 100, 800, 400]);

% Plot original shape
plot([0 L L 0 0], [-h/2 -h/2 h/2 h/2 -h/2], 'k--', 'LineWidth', 2);
hold on;

% Plot deformed shape - top and bottom edges
plot(x_top_def, y_top_def, 'b-', 'LineWidth', 2);
plot(x_bottom_def, y_bottom_def, 'b-', 'LineWidth', 2);

% Connect the ends with vertical lines
plot([x_top_def(1), x_bottom_def(1)], [y_top_def(1), y_bottom_def(1)], 'b-', 'LineWidth', 2); % Left end (fixed)
plot([x_top_def(end), x_bottom_def(end)], [y_top_def(end), y_bottom_def(end)], 'b-', 'LineWidth', 2); % Right end (free)

% Add some intermediate vertical lines for better visualization
if load_type == 2  % Pure moment case
    % For circular arc, use evenly spaced angles
    for i = [10, 25, 40, 55, 70, 85]
        plot([x_top_def(i), x_bottom_def(i)], [y_top_def(i), y_bottom_def(i)], 'b-', 'LineWidth', 0.8);
    end
else  % Shear force case
    for i = [10, 25, 40, 55, 70, 85]
        plot([x_top_def(i), x_bottom_def(i)], [y_top_def(i), y_bottom_def(i)], 'b-', 'LineWidth', 0.8);
    end
end

% Add arrows to indicate forces/moments at the right end based on load type
if load_type == 1
    % Draw shear force arrow
    arrow_length = 0.5;
    quiver(L, 0, 0, -arrow_length, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(L+0.1, -arrow_length/2, 'P', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
else
    % Draw moment curved arrow
    t = linspace(pi/2, 3*pi/2, 20);
    r = 0.3;
    xc = L + r*cos(t);
    yc = r*sin(t);
    plot(xc, yc, 'r-', 'LineWidth', 2);

    % Add arrowhead for moment
    arrow_angle = 3*pi/2;
    ah_x = [L+r*cos(arrow_angle), L+r*cos(arrow_angle)+0.1, L+r*cos(arrow_angle)-0.1];
    ah_y = [r*sin(arrow_angle), r*sin(arrow_angle)+0.1, r*sin(arrow_angle)+0.1];
    fill(ah_x, ah_y, 'r');

    % Add M label
    text(L+0.4, 0, 'M', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
end

% Add fixed end symbol (multiple short lines)
for i = -h/2:0.1:h/2
    line([-0.1, 0], [i, i], 'Color', 'k', 'LineWidth', 1);
end

title(sprintf('Cantilever Beam Deformation - %s Only (scale factor = 500)', load_name));
xlabel('x (m)');
ylabel('y (m)');

% Set reasonable axis limits to focus on the beam
if load_type == 2  % Pure moment case
    % For moment loading, adjust limits to show the circular arc
    max_x = max(max(x_top_def), max(x_bottom_def));
    max_y = max(max(y_top_def), max(y_bottom_def));
    min_x = min(min(x_top_def), min(x_bottom_def));
    min_y = min(min(y_top_def), min(y_bottom_def));

    % Add some padding
    xlim([min_x-0.01, max_x+0.01]);
    ylim([min_y-0.01, max_y+0.01]);
else
    xlim([-1, L+1]);
    ylim([-h, h]);
end

grid on;
legend('Original Shape', 'Deformed Shape', 'Location', 'best');

% Display relevant parameters based on load type
if load_type == 1
    text(L/2, -h*0.8, sprintf('Parameters: E = %.1e Pa, P = %.1e Pa', E, P), 'HorizontalAlignment', 'center');
else
    text(L/2, -h*0.8, sprintf('Parameters: E = %.1e Pa, M = %.1f N·m', E, M), 'HorizontalAlignment', 'center');
end

% Save the dedicated deformation figure with load type in filename
if load_type == 1
    saveas(gcf, './output/beam_deformation_shear.png');
else
    saveas(gcf, './output/beam_deformation_moment.png');
end

%% Plot stress field
figure(2);
set(gcf, 'Position', [100, 100, 1000, 800]);

% Plot S11 (normal stress in x-direction)
subplot(2,2,1);
contourf(X, Y, S11, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('S11 (σ_x_x)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum stress
[maxS11, maxS11Idx] = max(abs(S11(:)));
[maxS11_y, maxS11_x] = ind2sub(size(S11), maxS11Idx);
hold on;
plot(X(maxS11_y, maxS11_x), Y(maxS11_y, maxS11_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxS11_y, maxS11_x), Y(maxS11_y, maxS11_x), [' max = ', num2str(S11(maxS11_y, maxS11_x), '%.3e')], 'FontSize', 8);
fprintf('S11最大值: %g (位置: x=%g, y=%g)\n', S11(maxS11_y, maxS11_x), X(maxS11_y, maxS11_x), Y(maxS11_y, maxS11_x));

% Plot S12 (shear stress)
subplot(2,2,2);
contourf(X, Y, S12, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('S12 (τ_x_y)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum stress
[maxS12, maxS12Idx] = max(abs(S12(:)));
[maxS12_y, maxS12_x] = ind2sub(size(S12), maxS12Idx);
hold on;
plot(X(maxS12_y, maxS12_x), Y(maxS12_y, maxS12_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxS12_y, maxS12_x), Y(maxS12_y, maxS12_x), [' max = ', num2str(S12(maxS12_y, maxS12_x), '%.3e')], 'FontSize', 8);
fprintf('S12最大值: %g (位置: x=%g, y=%g)\n', S12(maxS12_y, maxS12_x), X(maxS12_y, maxS12_x), Y(maxS12_y, maxS12_x));

% Plot S1 (maximum principal stress)
subplot(2,2,3);
contourf(X, Y, S1, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('S1 (Maximum Principal Stress)');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum stress
[maxS1, maxS1Idx] = max(abs(S1(:)));
[maxS1_y, maxS1_x] = ind2sub(size(S1), maxS1Idx);
hold on;
plot(X(maxS1_y, maxS1_x), Y(maxS1_y, maxS1_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxS1_y, maxS1_x), Y(maxS1_y, maxS1_x), [' max = ', num2str(S1(maxS1_y, maxS1_x), '%.3e')], 'FontSize', 8);
fprintf('S1最大值: %g (位置: x=%g, y=%g)\n', S1(maxS1_y, maxS1_x), X(maxS1_y, maxS1_x), Y(maxS1_y, maxS1_x));

% Plot von Mises stress
subplot(2,2,4);
contourf(X, Y, Mises, 50, 'LineStyle', 'none');
colormap(jet(50));
colorbar;
title('von Mises Stress');
xlabel('x (m)');
ylabel('y (m)');
axis equal tight;

% Find and mark maximum stress
[maxMises, maxMisesIdx] = max(Mises(:));
[maxMises_y, maxMises_x] = ind2sub(size(Mises), maxMisesIdx);
hold on;
plot(X(maxMises_y, maxMises_x), Y(maxMises_y, maxMises_x), 'ko', 'MarkerFaceColor', 'k');
text(X(maxMises_y, maxMises_x), Y(maxMises_y, maxMises_x), [' max = ', num2str(maxMises, '%.3e')], 'FontSize', 8);
fprintf('Mises应力最大值: %g (位置: x=%g, y=%g)\n', maxMises, X(maxMises_y, maxMises_x), Y(maxMises_y, maxMises_x));

% Save stress figure with load type in filename
if load_type == 1
    saveas(gcf, './output/analytical_stress_shear.png');
else
    saveas(gcf, './output/analytical_stress_moment.png');
end

if load_type == 1
    fprintf('剪力载荷下的悬臂梁解析解云图已生成并保存在output文件夹中。\n');
else
    fprintf('力矩载荷下的悬臂梁解析解云图已生成并保存在output文件夹中。\n');
end
