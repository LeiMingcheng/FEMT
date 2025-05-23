% post_processing program for FEMT
%
% input: PROGRAM.CTR, PROGRAM.OUT
%
% Last Edited by  GAO HEXUAN at 2016-2-25
% Modified to support 9-node elements

%% initial
clear all;
close all;
clc;

% 设置异常值过滤的阈值因子（均值±threshold_factor*标准差之外的值被视为异常值）
threshold_factor = 6;

% 是否显示网格信息
show_mesh_info = true;

% 创建output目录（如果不存在）
if ~exist('./output', 'dir')
    mkdir('./output');
end

%% open files
fprintf('open files\n');

project_name = input('Input the program''s name:','s');  % 保存干净的项目名
tmpline = project_name;  % 用于后续处理
ctr = fopen(strcat(tmpline,'.ctr'), 'r');
out = fopen(strcat(tmpline,'.out'), 'r');


%% read coordinates
fprintf('read coordinates\n');
while ~strncmp(tmpline, ' ELEMENTS',9)
    tmpline = fgetl(ctr);
    if strncmp(tmpline, ' COORDINATES',12)
        [~,Nnodes] = strtok(tmpline);
        tmpline = fgetl(ctr);
        for i = 1 : str2num(Nnodes)
            P(i, :) = str2num(fgetl(ctr));
        end
    end
end

[~,Nelements] = strtok(tmpline);

%% read elements
fprintf('read elements\n');
while ~strncmp(deblank(tmpline), '    ELEMENT_MATERIAL',20)
    tmpline = fgetl(ctr);
    if strcmp(tmpline, '    ELEMENT_NODES')
        for i = 1 : str2num(Nelements)
            T(i, :) = str2num(fgetl(ctr));
        end
    end
end
fclose(ctr);

%% read displacement
fprintf('read displacement\n');
while isempty(strfind(tmpline, '*** NODAL STRESS ***'))
    tmpline = fgetl(out);
    if ~isempty(strfind(tmpline, '*** DISPLACEMENT ***'))
            tmpline = fgetl(out);
            tmpline = fgetl(out);
        for i = 1 : str2num(Nnodes)
            disp(i, :) = str2num(fgetl(out));
        end
    end
end
frewind(out);

%% read stress
fprintf('read stress\n');
while isempty(strfind(tmpline, 'PROGRAM STARTED'))
    tmpline = fgetl(out);
    if ~isempty(strfind(tmpline, '*** NODAL STRESS ***'))
            tmpline = fgetl(out);
            tmpline = fgetl(out);
        for i = 1 : str2num(Nnodes)
            stress(i, :) = str2num(fgetl(out));
        end
    end
end
fclose(out);

%% check mesh

%   viewmsh(P(:,2:3),T(:,2:end));

% 显示网格信息
if show_mesh_info
    [m, u] = size(T(:,2:end));
    fprintf('\n网格信息:\n');
    fprintf('节点数: %d\n', size(P, 1));
    fprintf('单元数: %d\n', m);

    % 确定元素类型
    if u == 3
        fprintf('元素类型: TRIANGLE3 (3节点三角形)\n');
    elseif u == 4
        fprintf('元素类型: RECTANGLE4 (4节点四边形)\n');
    elseif u == 6
        fprintf('元素类型: TRIANGLE6 (6节点三角形)\n');
    elseif u == 8
        fprintf('元素类型: RECTANGLE8 (8节点四边形)\n');
    elseif u == 9
        fprintf('元素类型: RECTANGLE9 (9节点四边形)\n');
    else
        fprintf('元素类型: 未知 (%d节点)\n', u);
    end
    fprintf('\n');
end


%% draw displacement
fprintf('calculate displacement\n');
scale_factor = max(max(abs(P(:,2:end))))/max(max(abs(disp(:,2:end))))*0.1;
post_P(:,2:3) = scale_factor*disp(:,2:3) + P(:,2:3);
post_P(:,1) = P(:,1);

U1 = disp(:,2);
U2 = disp(:,3);
U  = sqrt(disp(:,2).^2+disp(:,3).^2);
r = sqrt(P(:,2).^2 + P(:,3).^2);
n(:,1:2) = [P(:,2)./r,P(:,3)./r];
Un = n(:,1).*disp(:,2)+n(:,2).*disp(:,3);
t(:,1:2) = [-P(:,3)./r,P(:,2)./r];
Ut = t(:,1).*disp(:,2)+t(:,2).*disp(:,3);

x = P(:, 2);
y = P(:, 3);
x1= post_P(:,2);
y1= post_P(:,3);
tt = T(:,2:end);

[m, u] = size(tt);
if u == 6
    tt = [tt(:, 1), tt(:, 4), tt(:, 2), tt(:, 5), tt(:, 3), tt(:, 6)];
elseif u == 8
    tt = [tt(:, 1), tt(:, 5), tt(:, 2), tt(:, 6), tt(:, 3), tt(:, 7), tt(:, 4), tt(:, 8)];
elseif u == 9
    % 对于9节点元素，我们需要重新排列节点以便正确绘制
    % 9节点元素的节点顺序：4个角点，4个中点，1个中心点
    % 为了绘图，我们只使用8个外部节点，忽略中心节点
    tt = [tt(:, 1), tt(:, 5), tt(:, 2), tt(:, 6), tt(:, 3), tt(:, 7), tt(:, 4), tt(:, 8)];
end

tt = tt';


fprintf('draw displacement\n');

% Create figure with specific size to ensure enough space for text
figure(4)
set(gcf, 'Position', [100, 100, 900, 700])  % Wider figure
subplot(2,2,4)
patch(x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
patch(-x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(-x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;

patch(x1(tt), y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
patch(-x1(tt), y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(x1(tt), -y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(-x1(tt), -y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;

title('U-displacement', 'FontSize', 12, 'FontWeight', 'bold')
legend('before','after','Location','NorthEast')
axis image;
%colorbar;

% Find maximum displacement and print it
[Umax_disp, UI_disp] = max(U);
fprintf('【总位移U-displacement】最大值: %g (节点 %d)\n', Umax_disp, UI_disp);

% Add text annotation for maximum displacement at the right edge of the plot
ax = gca;
x_pos = max(x1(:));  % Right edge of the plot
y_pos = max(y1(:));  % Top of the plot
max_disp_text = sprintf('max = %.4e', Umax_disp);
text(x_pos, y_pos, max_disp_text, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

subplot(2,2,1)
drawc_filter(x,y,tt,U1,threshold_factor,'U1位移');
title('U1', 'FontSize', 12, 'FontWeight', 'bold');

subplot(2,2,2)
drawc_filter(x,y,tt,U2,threshold_factor,'U2位移');
title('U2', 'FontSize', 12, 'FontWeight', 'bold');

subplot(2,2,3)
drawc_filter(x,y,tt,U,threshold_factor,'总位移U');
title('U', 'FontSize', 12, 'FontWeight', 'bold');

% 保存位移云图
saveas(gcf, ['./output/' project_name '_displacement.png']);

%% draw stress
fprintf('calculate stress\n');
S11 = stress(:,2);
S22 = stress(:,3);
S12 = stress(:,4);
S1  = stress(:,5);
S2  = stress(:,6);

for i = 1:length(P(:,2))
    SN = [n(i,1),n(i,2);t(i,1),t(i,2)]*[stress(i,2),stress(i,4);stress(i,4),stress(i,3)]*[n(i,1),n(i,2);t(i,1),t(i,2)]';
    SNN(i) = SN(1,1);
    STT(i) = SN(2,2);
    SNT(i) = SN(1,2);
end

Mises  = sqrt(3/2*((S1-(S1+S2)/3).^2+(S2-(S1+S2)/3).^2));


fprintf('draw stress\n');

% Create figure with specific size to ensure enough space for text
figure(5)
set(gcf, 'Position', [100, 100, 900, 700])  % Wider figure
subplot(2,2,1)
drawc_filter(x,y,tt,S11,threshold_factor,'S11应力');
title('S11', 'FontSize', 12, 'FontWeight', 'bold');

subplot(2,2,2)
drawc_filter(x,y,tt,S22,threshold_factor,'S22应力');
title('S22', 'FontSize', 12, 'FontWeight', 'bold');

subplot(2,2,3)
drawc_filter(x,y,tt,S12,threshold_factor,'S12应力');
title('S12', 'FontSize', 12, 'FontWeight', 'bold');


subplot(2,2,4)
% Identify and filter outliers in Mises stress
U_mean = mean(Mises, 'omitnan');
U_std = std(Mises, 'omitnan');
U_var = var(Mises, 'omitnan');
upper_threshold = U_mean + threshold_factor * U_std;
lower_threshold = U_mean - threshold_factor * U_std;

% Print mean and variance
fprintf('【Mises应力】统计信息：均值 = %g，方差 = %g\n', U_mean, U_var);

% Identify outliers
outlier_mask = Mises > upper_threshold | Mises < lower_threshold;
outlier_nodes = find(outlier_mask);
outlier_values = Mises(outlier_mask);

% Create a copy of Mises with outliers replaced by mean value
Mises_filtered = Mises;
Mises_filtered(outlier_mask) = U_mean;

% Display information about outliers
if ~isempty(outlier_nodes)
    fprintf('【Mises应力】发现 %d 个异常节点：\n', length(outlier_nodes));
    for i = 1:length(outlier_nodes)
        fprintf('节点 %d: 值 = %g\n', outlier_nodes(i), outlier_values(i));
    end
else
    fprintf('【Mises应力】未发现异常节点。\n');
end

% Draw the filtered plot
patch(x1(tt), y1(tt), Mises_filtered(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
patch(-x1(tt), y1(tt), Mises_filtered(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;

title('Mises', 'FontSize', 12, 'FontWeight', 'bold');
%colorbar;
axis image;
box off;

% Find maximum absolute value (excluding outliers)
Mises_no_outliers = Mises;
Mises_no_outliers(outlier_mask) = 0;
[Umax, UI] = max(abs(Mises_no_outliers));

% Add text annotation for maximum value at the right edge of the plot
ax = gca;
x_pos = max(x1(:));  % Right edge of the plot
y_pos = max(y1(:));  % Top of the plot
max_text = sprintf('max = %.4e', sign(Mises(UI))*Umax);
text(x_pos, y_pos, max_text, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

% Print the maximum value to the console
fprintf('【Mises应力】最大值: %g (节点 %d)\n', sign(Mises(UI))*Umax, UI);

% 保存应力云图
saveas(gcf, ['./output/' project_name '_stress.png']);

