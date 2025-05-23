function [outlier_nodes, outlier_values] = drawc_filter(x, y, tt, U, threshold_factor, plot_title)
    % drawc_filter - Draw contour plot with outlier filtering
    %
    % Inputs:
    %   x, y - Coordinates of nodes
    %   tt - Element connectivity matrix (transposed)
    %   U - Values at nodes
    %   threshold_factor - Factor to determine outliers
    %                     Values outside mean ± threshold_factor*std are considered outliers
    %   plot_title - Title of the plot (for printing statistics)
    %
    % Outputs:
    %   outlier_nodes - Node indices of outliers
    %   outlier_values - Values at outlier nodes

    if nargin < 5
        threshold_factor = 3; % Default threshold factor
    end

    if nargin < 6
        plot_title = 'Unknown';
    end

    % Find outliers using mean and standard deviation
    U_mean = mean(U, 'omitnan');
    U_std = std(U, 'omitnan');
    U_var = var(U, 'omitnan');
    upper_threshold = U_mean + threshold_factor * U_std;
    lower_threshold = U_mean - threshold_factor * U_std;

    % Print mean and variance
    fprintf('【%s】统计信息：均值 = %g，方差 = %g\n', plot_title, U_mean, U_var);

    % Identify outliers
    outlier_mask = U > upper_threshold | U < lower_threshold;
    outlier_nodes = find(outlier_mask);
    outlier_values = U(outlier_mask);

    % Create a copy of U with outliers replaced by mean value
    U_filtered = U;
    U_filtered(outlier_mask) = U_mean;

    % Display information about outliers
    if ~isempty(outlier_nodes)
        fprintf('【%s】发现 %d 个异常节点：\n', plot_title, length(outlier_nodes));
        for i = 1:length(outlier_nodes)
            fprintf('节点 %d: 值 = %g\n', outlier_nodes(i), outlier_values(i));
        end
    else
        fprintf('【%s】未发现异常节点。\n', plot_title);
    end

    % Draw the filtered plot (using the original tt but with filtered values)
    patch(x(tt), y(tt), U_filtered(tt), 'EdgeColor', 'none');
    hold on;

    colormap(jet(50));
    %colorbar;
    axis image;
    box off;

    % Find maximum absolute value (excluding outliers)
    U_no_outliers = U;
    U_no_outliers(outlier_mask) = 0;
    [Umax, UI] = max(abs(U_no_outliers));

    % Add text annotation for maximum value at the right edge of the plot
    ax = gca;
    x_pos = max(x(:));  % Right edge of the plot
    y_pos = max(y(:));  % Top of the plot
    max_text = sprintf('max = %.4e', sign(U(UI))*Umax);
    text(x_pos, y_pos, max_text, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);

    % Print the maximum value to the console
    fprintf('【%s】最大值: %g (节点 %d)\n', plot_title, sign(U(UI))*Umax, UI);
end
