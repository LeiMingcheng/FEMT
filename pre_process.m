% pre_process program for FEMT
%
% input: PROGRAM.CTR
% output: PROGRAM.CTR
%
% Last Edited by  GAO HEXUAN at 2016-2-24

%% initial
clear all;
close all;
clc;

%% 统一几何和材料参数设置
% -------------------------------------------------------------------------
% 材料参数
E_modulus = 2.1E11;           % 弹性模量 (Pa)，钢材为210 GPa
Poisson_ratio = 0.3;          % 泊松比

% 几何参数
element_thickness = 0.01;      % 平面单元厚度 (m)

%% 边界条件设置
% -------------------------------------------------------------------------
% 边界条件类型设置
% 1=拉力/压力(原程序默认), 2=剪力, 3=力矩, 4=分布力P
boundary_type = 4;  % 修改此值以改变边界条件类型

% 拉力/压力边界设置 (boundary_type = 1)
Force_Line_Numbers1 = [2,4,6,8];  % 设置第一组力边界的Line numbers数组
Pa_Tensile1 = 1e6;            % 设置第一组力的大小 (Pa)
Force_Line_Numbers2 = [1,3];  % 设置第二组力边界的Line numbers数组
Pa_Tensile2 = 1e6;            % 设置第二组力的大小 (Pa)

% 剪力边界设置 (boundary_type = 2)
Shear_Line_Number = 2;        % 设置应用剪力的边界线编号
Shear_Stress = 1e6;           % 设置剪应力大小 (Pa)

% 力矩边界设置 (boundary_type = 3)
Moment_Line_Number = 2;       % 设置应用力矩的边界线编号
Moment_Value = 20.0;          % 设置力矩大小 (N·m)
Beam_Height = 0.02;           % 设置梁的高度 (m)，用于计算力矩产生的应力

% 分布力P边界设置 (boundary_type = 4)
P_Force_Line_Number = 2;      % 设置应用分布力P的边界线编号
P_Force_Value = -3e6;        % 设置力P的大小 (N)

% -------------------------------------------------------------------------
% USER INPUT FOR BOUNDARIES - END


%% open files
fprintf('open files\n');
project = input('Input the program''s name:','s');
inp = fopen(strcat(project,'.inp'), 'r');
ctr = fopen(strcat(project,'.ctr'), 'w');
bnd = fopen(strcat(project,'.bnd'), 'w');

%% read sub_items
fprintf('Read sub_items\n');
tic;
numl = 0;
line_num = 0;
Nset_num = 0;
elem_sta = 0;
elem_num = 0;
ei=0;
tmpline = fgetl(inp);numl = numl+1;

while tmpline ~= -1
    % read coordinates
    if strcmp(tmpline, '*NODE')
    fprintf('\tRead coordinates ');
        node_sta = numl+1;
        i=0;
        tmpline = fgetl(inp);numl = numl+1;
        while ~strncmp(tmpline, '*',1)
            i=i+1;
            P(i, :) = str2num(tmpline);
            tmpline = fgetl(inp);numl = numl+1;
        end
        node_end = numl-1;
     fprintf(strcat(num2str(node_end-node_sta+1),'\n'));

    % read line set
    elseif strncmp(tmpline(max(strfind(tmpline,'='))+1:end-1),'Line',4) % 修改: 使用strncmp检查前4个字符是否为'Line'
        line_num = line_num+1;
        line_sta(line_num) = numl+1;
        i=0;
        tmpline = fgetl(inp);numl = numl+1;
        while ~strncmp(tmpline, '*',1)
            i=i+1;
            Line(i, :,line_num) = str2num(tmpline);
            tmpline = fgetl(inp);numl = numl+1;
        end
        line_end(line_num) = numl-1;

     fprintf('\t\t %d elements from x%5.2f y%5.2f to  x%5.2f y%5.2f\n'...
            ,line_end(line_num)-line_sta(line_num)+1 ...
            ,[P(Line(1, 2,line_num),2),P(Line(1, 2,line_num),3)]...
            ,[P(Line(line_end(line_num)-line_sta(line_num)+1, end,line_num),2),P(Line(line_end(line_num)-line_sta(line_num)+1, end,line_num),3)]);


    % read element
    elseif startsWith(tmpline, '*ELEMENT') % Check if it's an element definition line
        idx_elset_equals = strfind(tmpline, 'ELSET=');
        elset_name_actual = ''; % Initialize to ensure ELSET is not found
        if ~isempty(idx_elset_equals)
            % Start from after "ELSET=" and take up to the trailing comma, or to the end if no comma
            name_part_after_elset = strtrim(tmpline(idx_elset_equals(1)+length('ELSET='):end));
            comma_position = strfind(name_part_after_elset, ',');
            if ~isempty(comma_position)
                elset_name_actual = strtrim(name_part_after_elset(1:comma_position(1)-1));
            else
                elset_name_actual = name_part_after_elset; % ELSET name is the last part of the line
            end
        end

        if ~isempty(elset_name_actual) && startsWith(elset_name_actual, 'Surface')
            % Currently, only process ELSETs starting with "Surface" (e.g., Surface, Surface9, Surface11, MySurface)
            % This part can be expanded later
            fprintf('	Read elements from ELSET %s ', elset_name_actual); % Print the actual ELSET name being read

            % --- User's original element reading logic ---
            elem_sta = numl+1; % Points to the line after *ELEMENT ... ELSET=...
            tmpline = fgetl(inp);numl = numl+1; % Read the first element data line

            elements_read_in_this_block = 0; % Counter for elements read in this block
            while ~strncmp(tmpline, '*',1) % User's original inner loop condition
                ei=ei+1; % Global element counter, correct logic
                T(ei, :) = str2num(tmpline);
                elements_read_in_this_block = elements_read_in_this_block + 1;
                tmpline = fgetl(inp);numl = numl+1;
            end
            temp = elem_num; % Store the element count before this block
            % Use user's original counting method (based on line numbers), this method ensures correctness for this block
            elem_num = elem_num + (numl-1-elem_sta+1);
            fprintf('%d to %d\n',temp+1,elem_num); % Print the range of elements read in this block
        end
    elseif strcmp(strtok(tmpline,','),'*NSET')
        % --- MODIFIED NSET READING LOGIC ---
        idx_nset_equals = strfind(tmpline, 'NSET=');
        is_target_nset = false;
        actual_nset_name_for_print = '';

        if ~isempty(idx_nset_equals)
            name_part_nset = tmpline(idx_nset_equals(1)+length('NSET='):end);
            comma_pos_nset = strfind(name_part_nset, ',');
            if ~isempty(comma_pos_nset)
                actual_nset_name = strtrim(name_part_nset(1:comma_pos_nset(1)-1));
            else
                actual_nset_name = strtrim(name_part_nset);
            end
            actual_nset_name_for_print = actual_nset_name; % For printing

            % Check for both PhysicalPoint and PhysicalLine
            if startsWith(actual_nset_name, 'PhysicalPoint') || startsWith(actual_nset_name, 'PhysicalLine')
                is_target_nset = true;
            end
        end

        if is_target_nset
            Nset_num = Nset_num+1; % Increment NSET counter

            if Nset_num == 1 % First Physical Point/Line NSET
                fprintf('	Read boundary set(U=0) for %s ', actual_nset_name_for_print);
            elseif Nset_num == 2 % Second Physical Point/Line NSET
                fprintf('	Read boundary set(V=0) for %s ', actual_nset_name_for_print);
            else % Other Physical Point/Line NSETs
                fprintf('	Read NSET %s (processed as set index %d) ', actual_nset_name_for_print, Nset_num);
            end

            j = 0; % Node counter for the current NSET
            tmpline = fgetl(inp);numl = numl+1; % Read the first line of nodes under NSET
            while ~strncmp(tmpline, '*',1)
                % Nodes under NSET in .inp file are comma-separated on one line
                current_nodes_on_line = str2num(tmpline);
                for k_node = 1:length(current_nodes_on_line)
                   j = j+1;
                   Nset(Nset_num,j) = current_nodes_on_line(k_node); % Use incremented Nset_num as row index
                end
                tmpline = fgetl(inp);numl = numl+1;
            end
            Nset_size(Nset_num) = j; % Store the number of nodes in this NSET (identified by Nset_num)
            fprintf('%d nodes\n',j);

        % Keep user's original elseif strcmp(tmpline(20:23),'Surf') to handle other NSET types
        % For example, in hw1.inp, this prevents errors.
        elseif strcmp(tmpline(20:23),'Surf')
            break;
        end

    else
    tmpline = fgetl(inp);numl = numl+1;
    end
end
toc;

p = P(:,2:3);
t = T(:,2:end);
L = Line(:,2:end,:); % L will contain all line sets read from the .inp file
clear P T Line;


%% semiband optimazition
fprintf('Semiband optimazition\n');
K = zeros(length(p), length(p));

semiBandwidth = 0;
semiBandwidthOpt = 0;

[m, n] = size(t);

for i = 1 : m
    tmp = max(t(i, :)) - min(t(i, :));
    if tmp > semiBandwidth
        semiBandwidth = tmp;
    end
end


semiBandwidth = (semiBandwidth + 1) * 2;

% 与原始pre_process_ori.m完全一致：为K矩阵构建时添加一列
t = [t, t(:, 1)];
for i = 1 : length(p)
    K(i, i) = 1;
end
for i = 1 : m
    for j = 1 : n
        K(t(i, j), t(i, j + 1)) = 1;
        K(t(i, j + 1), t(i, j)) = 1;
    end
end

% Use Cuthill-McKee optimazition
figure(4)
subplot(1,2,1)
spy(K);
title('Original Adjacency Matrix before Optimization')
axis equal;
axis([0 length(p) 0 length(p)])

subplot(1,2,2)
q = symrcm(K);
spy(K(q,q));
title('Optimized Adjacency Matrix (Cuthill-McKee)')
axis equal;
axis([0 length(p) 0 length(p)])
clear K;

p = p(q, :);
q1 = q * 0;
for i = 1 : length(q)
    q1(q(i)) = i;
end

% 严格按照pre_process_ori.m处理t矩阵
[m, n] = size(t);
t(:, n) = [];
for i = 1 : m
    for j = 1 : n - 1
        t(i, j) = q1(t(i, j));
    end
end

[m, n_nodes_per_line_segment, o_num_line_sets] = size(L);
for k = 1 : o_num_line_sets % Iterate over all line sets
    for i = 1 : m % Iterate over segments in a line set
        for j = 1 : n_nodes_per_line_segment % Iterate over nodes in a segment
            if L(i, j, k) ~= 0
                L(i, j, k) = q1(L(i, j, k));
            else
                L(i, j, k) = 0; % Should already be 0 if not used
            end
        end
    end
end

[m,n] = size(Nset);
for i = 1 : m
    for j = 1 : n
        if Nset(i,j) ~= 0
            Nset(i,j) = q1(Nset(i, j));
        end
    end
end

% 计算优化后的半带宽 - 完全按照ori版本
[m, n] = size(t);
for i = 1 : m
    tmp = max(t(i, :)) - min(t(i, :));
    if tmp > semiBandwidthOpt
        semiBandwidthOpt = tmp;
    end
end

semiBandwidthOpt = (semiBandwidthOpt + 1) * 2;


%% Set up force boundaries based on the selected boundary type
fprintf('Setting up force boundaries with boundary_type = %d\n', boundary_type);

% Initialize variables for different boundary types
ForceSets1 = [];
ForceSets2 = [];

% Only process tensile force boundaries if boundary_type is 1
if boundary_type == 1
    fprintf('Processing tensile/compressive force boundaries\n');
    fprintf('Group 1 line numbers: %s with force magnitude: %5.2f\n', mat2str(Force_Line_Numbers1), Pa_Tensile1);
    fprintf('Group 2 line numbers: %s with force magnitude: %5.2f\n', mat2str(Force_Line_Numbers2), Pa_Tensile2);

    % Initialize empty cell arrays to store force sets for each group
    ForceSets1 = cell(length(Force_Line_Numbers1), 1);
    ForceSets2 = cell(length(Force_Line_Numbers2), 1);

    % Process each line number in the first Force_Line_Numbers array
    for line_idx = 1:length(Force_Line_Numbers1)
        line_num = Force_Line_Numbers1(line_idx);

        fprintf('Force %5.2f applied on No.%d line (Group 1)\n', Pa_Tensile1, line_num);

        % Extract the line data for the current line number
        current_line_set = L(:,:,line_num);

        % Add column for element index
        current_force_set = [zeros(size(current_line_set,1),1), current_line_set];

        % Replace 0 with NaN for easier handling, except for first column
        current_force_set(current_force_set == 0) = NaN;

        % Find elements included in the current force set
        for i=1:elem_num
            switch (length(t(1,:))) % Number of nodes per element
                case {4,3} % Linear quadrilateral or triangle (2 nodes per edge)
                    for j = 1:size(current_force_set,1)
                        if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) % Ensure nodes are valid
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ~isempty(find(t(i,:) == current_force_set(j,3), 1))
                                current_force_set(j,1) = i;
                            end
                        end
                    end
                case {6,8,9} % Quadratic triangle or quadrilateral (3 nodes per edge for force application)
                    for j = 1:size(current_force_set,1)
                         if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) && ~isnan(current_force_set(j,4)) % Ensure nodes are valid
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,3), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,4), 1))
                                    current_force_set(j,1) = i;
                            end
                        end
                    end
            end
        end

        % Remove rows where no element was found or nodes were NaN
        current_force_set = current_force_set(~isnan(current_force_set(:,1)) & ~isnan(current_force_set(:,2)),:);

        % Store the processed force set in the cell array
        ForceSets1{line_idx} = current_force_set;

        fprintf('  Processed %d force elements for line %d (Group 1)\n', size(current_force_set,1), line_num);
    end

    % Process each line number in the second Force_Line_Numbers array
    for line_idx = 1:length(Force_Line_Numbers2)
        line_num = Force_Line_Numbers2(line_idx);

        fprintf('Force %5.2f applied on No.%d line (Group 2)\n', Pa_Tensile2, line_num);

        % Extract the line data for the current line number
        current_line_set = L(:,:,line_num);

        % Add column for element index
        current_force_set = [zeros(size(current_line_set,1),1), current_line_set];

        % Replace 0 with NaN for easier handling, except for first column
        current_force_set(current_force_set == 0) = NaN;

        % Find elements included in the current force set
        for i=1:elem_num
            switch (length(t(1,:))) % Number of nodes per element
                case {4,3} % Linear quadrilateral or triangle (2 nodes per edge)
                    for j = 1:size(current_force_set,1)
                        if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) % Ensure nodes are valid
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ~isempty(find(t(i,:) == current_force_set(j,3), 1))
                                current_force_set(j,1) = i;
                            end
                        end
                    end
                case {6,8,9} % Quadratic triangle or quadrilateral (3 nodes per edge for force application)
                    for j = 1:size(current_force_set,1)
                         if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) && ~isnan(current_force_set(j,4)) % Ensure nodes are valid
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,3), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,4), 1))
                                    current_force_set(j,1) = i;
                            end
                        end
                    end
            end
        end

        % Remove rows where no element was found or nodes were NaN
        current_force_set = current_force_set(~isnan(current_force_set(:,1)) & ~isnan(current_force_set(:,2)),:);

        % Store the processed force set in the cell array
        ForceSets2{line_idx} = current_force_set;

        fprintf('  Processed %d force elements for line %d (Group 2)\n', size(current_force_set,1), line_num);
    end
elseif boundary_type == 2
    fprintf('Processing shear force boundary on line %d with magnitude: %10.3E Pa\n', Shear_Line_Number, Shear_Stress);
    % 实际处理在FORCE BOUNDARY SECTION中完成
elseif boundary_type == 3
    fprintf('Processing moment boundary on line %d with magnitude: %5.2f N·m\n', Moment_Line_Number, Moment_Value);
    % 实际处理在FORCE BOUNDARY SECTION中完成
elseif boundary_type == 4
    fprintf('Processing distributed force P boundary on line %d with magnitude: %5.2f N\n', P_Force_Line_Number, P_Force_Value);
    % 实际处理在FORCE BOUNDARY SECTION中完成
end


%% write to ctr and bnd
fprintf('Write to ctr\n');

% Debug information for t matrix
fprintf('Debug: t matrix size = [%d, %d]\n', size(t,1), size(t,2));
fprintf('Debug: length of first row = %d\n', length(t(1,:)));
fprintf('Debug: length(t(1,:)) = %d, which should determine element type\n', length(t(1,:)));
if length(t(1,:)) == 4
    fprintf('Debug: Should use RECTANGLE4\n');
elseif length(t(1,:)) == 8
    fprintf('Debug: Should use RECTANGLE8\n');
elseif length(t(1,:)) == 9
    fprintf('Debug: Should use RECTANGLE9\n');
elseif length(t(1,:)) == 3
    fprintf('Debug: Should use TRIANGLE3\n');
elseif length(t(1,:)) == 6
    fprintf('Debug: Should use TRIANGLE6\n');
else
    fprintf('Debug: WARNING - No matching element type for length(t(1,:)) = %d\n', length(t(1,:)));
    fprintf('Debug: First few values of t(1,:): ');
    if ~isempty(t)
        fprintf('%d ', t(1,1:min(10,size(t,2))));
    else
        fprintf('t is empty!');
    end
    fprintf('\n');
end

% BASIC DATA
fprintf(ctr,    ' BASIC_DATA MECHANICAL\r\n');
fprintf(ctr,    '\r\n');

% NODE COORDINATES
fprintf(ctr,    ' COORDINATES %d\r\n',node_end - node_sta + 1);
fprintf(ctr,    ' !   NODE             X                   Y                     Z\r\n');
for i = 1:node_end-node_sta+1
    fprintf(ctr,'%8d    %18.15f    %18.15f    %18.15f\r\n', i,p(i,1),p(i,2),0 );
end
fprintf(ctr, '\r\n');

% ELEMENT
fprintf(ctr,    ' ELEMENTS  %d\r\n',elem_num);
fprintf(ctr,    '    ELEMENT_TYPE\r\n');
switch(length(t(1,:)))
    case 4
        ELEMTYPE = 'RECTANGLE4';
        fprintf(ctr,  '       1 TO  %d  TYPE %s\r\n',elem_num,ELEMTYPE);
        fprintf(ctr,  '    END ELEMENT_TYPE\r\n');
        fprintf(ctr, '\r\n');
        fprintf(ctr,  '    ELEMENT_NODES\r\n');
        for i = 1:elem_num
            fprintf(ctr,'%8d    %6d    %6d    %6d    %6d\r\n', i,t(i,1),t(i,2),t(i,3),t(i,4));
        end
        fprintf(ctr, '\r\n');

    case 8
        ELEMTYPE = 'RECTANGLE8';
        fprintf(ctr,  '       1 TO  %d  TYPE %s\r\n',elem_num,ELEMTYPE);
        fprintf(ctr,  '    END ELEMENT_TYPE\r\n');
        fprintf(ctr, '\r\n');
        fprintf(ctr,  '    ELEMENT_NODES\r\n');
        for i = 1:elem_num
            fprintf(ctr,'%8d    %6d    %6d    %6d    %6d    %6d    %6d    %6d    %6d\r\n',...
                i,t(i,1),t(i,2),t(i,3),t(i,4),t(i,5),t(i,6),t(i,7),t(i,8));
        end
        fprintf(ctr, '\r\n');

    case 9
        ELEMTYPE = 'RECTANGLE9';
        fprintf(ctr,  '       1 TO  %d  TYPE %s\r\n',elem_num,ELEMTYPE);
        fprintf(ctr,  '    END ELEMENT_TYPE\r\n');
        fprintf(ctr, '\r\n');
        fprintf(ctr,  '    ELEMENT_NODES\r\n');
        for i = 1:elem_num
            fprintf(ctr,'%8d    %6d    %6d    %6d    %6d    %6d    %6d    %6d    %6d    %6d\r\n',...
                i,t(i,1),t(i,2),t(i,3),t(i,4),t(i,5),t(i,6),t(i,7),t(i,8),t(i,9));
        end
        fprintf(ctr, '\r\n');
    case 3
        ELEMTYPE = 'TRIANGLE3';
        fprintf(ctr,  '       1 TO  %d  TYPE %s\r\n',elem_num,ELEMTYPE);
        fprintf(ctr,  '    END ELEMENT_TYPE\r\n');
        fprintf(ctr, '\r\n');
        fprintf(ctr,  '    ELEMENT_NODES\r\n');
        for i = 1:elem_num
            fprintf(ctr,'%8d    %6d    %6d    %6d\r\n', i,t(i,1),t(i,2),t(i,3));
        end
        fprintf(ctr, '\r\n');
    case 6
        ELEMTYPE = 'TRIANGLE6';
        fprintf(ctr,  '       1 TO  %d  TYPE %s\r\n',elem_num,ELEMTYPE);
        fprintf(ctr,  '    END ELEMENT_TYPE\r\n');
        fprintf(ctr, '\r\n');
        fprintf(ctr,  '    ELEMENT_NODES\r\n');
        for i = 1:elem_num
            fprintf(ctr,'%8d    %6d    %6d    %6d    %6d    %6d    %6d\r\n',...
                i,t(i,1),t(i,2),t(i,3),t(i,4),t(i,5),t(i,6));
        end
        fprintf(ctr, '\r\n');
end

%-- ELEMENT MATERIAL --
fprintf(ctr,'    ELEMENT_MATERIAL\r\n');
fprintf(ctr,'      1 TO    %d    MATERIAL 1\r\n', elem_num);
fprintf(ctr,'    END ELEMENT_MATERIAL\r\n');
fprintf(ctr, '\r\n');


%-- ELEMENT GEOMETRY --
fprintf(ctr,'    ELEMENT_GEOMETRY\r\n');
fprintf(ctr,'       1 TO    %d    GEOMETRY 1\r\n', elem_num);
fprintf(ctr,'    END ELEMENT_GEOMETRY\r\n');
fprintf(ctr,' END ELEMENTS\r\n');
fprintf(ctr, '\r\n');


%-- MATERIALS --
fprintf(ctr,' MATERIALS 1\r\n');
fprintf(ctr,'    MATERIAL 1 TYPE ISOTROPIC\r\n');
fprintf(ctr,'       E  %g\r\n', E_modulus); % 使用统一定义的弹性模量
fprintf(ctr,'       v  %g\r\n', Poisson_ratio); % 使用统一定义的泊松比
fprintf(ctr,'    END MATERIAL\r\n');
fprintf(ctr,' END MATERIALS\r\n');
fprintf(ctr, '\r\n');


%-- GEOMETRIES --
fprintf(ctr,' GEOMETRIES 1\r\n');
fprintf(ctr,'    GEOMETRY 1 TYPE PLANE\r\n');
fprintf(ctr,'       THICKNESS %g\r\n', element_thickness);
fprintf(ctr,'    END GEOMETRY\r\n');
fprintf(ctr,' END GEOMETRIES\r\n');
fprintf(ctr,'  \r\n');
fprintf(ctr,' END BASIC_DATA\r\n');
fprintf(ctr,'  \r\n');


%-- SOLUTION --
fprintf(ctr,' SOLUTION\r\n');
fprintf(ctr,'    STIFFNESS PLANE_STRESS\r\n');
fprintf(ctr,'    BOUNDARY\r\n');
fprintf(ctr,'    FORCE\r\n');
fprintf(ctr,'    SOLVE MECHANICAL STATIC DIR\r\n');
fprintf(ctr,'    STRESS NODAL PLANE_STRESS\r\n');
fprintf(ctr,'    OUTPUT DISPLACEMENT NODAL_STRESS\r\n');
fprintf(ctr,' END SOLUTION\r\n');

fclose(ctr);

%-- BOUNDARY CONDITIONS --
fprintf('Write to bnd\n');
fprintf(bnd,  '! BOUNDARY CONDITIONS\r\n');
fprintf(bnd,  '  STEP 1\r\n');

% DISPLACEMENT CONSTRAINTS
fprintf(bnd,  '    DISPLACEMENT\r\n');

% Line boundary (Default from original code, ensure Nset_size is valid)
if Nset_num >= 2 % Check if at least two NSETs (for U and V) are defined
    % Apply necessary constraints to ensure the stiffness matrix is positive definite
    % For 1/4 symmetric model:
    % 1. X=0 boundary constrains U direction
    % 2. Y=0 boundary constrains V direction

    % Calculate total number of constraints
    total_constraints = Nset_size(1) + Nset_size(2);
    fprintf(bnd, '       GIVEN_DISP    %d\r\n', total_constraints);

    % Apply U constraints to X=0 boundary (first NSET)
    for i = 1:Nset_size(1)
        fprintf(bnd, '        %4d   U  0.0\r\n', Nset(1,i));
    end

    % Apply V constraints to Y=0 boundary (second NSET)
    for i = 1:Nset_size(2)
        fprintf(bnd, '        %4d   V  0.0\r\n', Nset(2,i));
    end
else % Fallback or error if not enough NSETs
    fprintf('Warning: Not enough NSETs defined for default displacement constraints.\n');
    % fprintf(bnd,  '       GIVEN_DISP    0\r\n'); % Or handle as error
end

fprintf(bnd,  '       END GIVEN_DISP\r\n');
fprintf(bnd,  '    END DISPLACEMENT\r\n');
fprintf(bnd, '\r\n');
fprintf(bnd,  '    FORCE\r\n');


% FORCE BOUNDARY SECTION - START
% -------------------------------------------------------------------------

% 根据边界条件类型选择不同的处理方式
if boundary_type == 1
    % 原始拉力/压力边界条件处理 (原程序默认)
    % Calculate total number of force elements across all force sets from both groups
    total_line_forces = 0;

    % Count forces from Group 1
    for i = 1:length(ForceSets1)
        total_line_forces = total_line_forces + size(ForceSets1{i}, 1);
    end

    % Count forces from Group 2
    for i = 1:length(ForceSets2)
        total_line_forces = total_line_forces + size(ForceSets2{i}, 1);
    end

    fprintf(bnd,'       LINE_FORCE  %d\r\n', total_line_forces);

    % Write forces for Group 1 force boundaries
    % Count total elements in Group 1
    total_elements_group1 = 0;
    for i = 1:length(ForceSets1)
        total_elements_group1 = total_elements_group1 + size(ForceSets1{i}, 1);
    end
    fprintf('Writing %d force elements for Group 1 (magnitude: %5.2f)\n', total_elements_group1, Pa_Tensile1);
    for set_idx = 1:length(ForceSets1)
        current_force_set = ForceSets1{set_idx};

        for i = 1:size(current_force_set, 1)
            if length(t(1,:)) == 4 || length(t(1,:)) == 3
                fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_N     %10.3E    %10.3E\r\n',...
                current_force_set(i,1), current_force_set(i,2), current_force_set(i,3), Pa_Tensile1, Pa_Tensile1);
            elseif length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
                fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_N     %10.3E    %10.3E    %10.3E\r\n',...
                current_force_set(i,1), current_force_set(i,2), current_force_set(i,3), current_force_set(i,4), Pa_Tensile1, Pa_Tensile1, Pa_Tensile1);
            end
        end
    end

    % Write forces for Group 2 force boundaries
    % Count total elements in Group 2
    total_elements_group2 = 0;
    for i = 1:length(ForceSets2)
        total_elements_group2 = total_elements_group2 + size(ForceSets2{i}, 1);
    end
    fprintf('Writing %d force elements for Group 2 (magnitude: %5.2f)\n', total_elements_group2, Pa_Tensile2);
    for set_idx = 1:length(ForceSets2)
        current_force_set = ForceSets2{set_idx};

        for i = 1:size(current_force_set, 1)
            if length(t(1,:)) == 4 || length(t(1,:)) == 3
                fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_N     %10.3E    %10.3E\r\n',...
                current_force_set(i,1), current_force_set(i,2), current_force_set(i,3), Pa_Tensile2, Pa_Tensile2);
            elseif length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
                fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_N     %10.3E    %10.3E    %10.3E\r\n',...
                current_force_set(i,1), current_force_set(i,2), current_force_set(i,3), current_force_set(i,4), Pa_Tensile2, Pa_Tensile2, Pa_Tensile2);
            end
        end
    end

elseif boundary_type == 2
    % 剪力边界条件处理
    % 获取指定边界线上的元素
    shear_force_set = [];

    % 从指定的Line中获取边界元素
    if Shear_Line_Number <= size(L, 3)
        current_line_set = L(:,:,Shear_Line_Number);

        % 添加元素索引列
        current_force_set = [zeros(size(current_line_set,1),1), current_line_set];

        % 替换0为NaN以便处理
        current_force_set(current_force_set == 0) = NaN;

        % 查找包含在当前力集中的元素
        for i=1:elem_num
            switch (length(t(1,:))) % 每个元素的节点数
                case {4,3} % 线性四边形或三角形（每边2个节点）
                    for j = 1:size(current_force_set,1)
                        if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ~isempty(find(t(i,:) == current_force_set(j,3), 1))
                                current_force_set(j,1) = i;
                            end
                        end
                    end
                case {6,8,9} % 二次三角形或四边形（每边3个节点用于力的应用）
                    for j = 1:size(current_force_set,1)
                         if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) && ~isnan(current_force_set(j,4)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,3), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,4), 1))
                                    current_force_set(j,1) = i;
                            end
                        end
                    end
            end
        end

        % 移除未找到元素或节点为NaN的行
        shear_force_set = current_force_set(~isnan(current_force_set(:,1)) & ~isnan(current_force_set(:,2)),:);
    end

    % 写入剪力边界条件
    fprintf(bnd,'       LINE_FORCE  %d\r\n', size(shear_force_set, 1));

    % 对每个边界元素应用剪应力
    fprintf('Writing %d shear stress elements (magnitude: %10.3E Pa)\n', size(shear_force_set, 1), Shear_Stress);

    % 使用统一定义的几何参数（厚度）

    for i = 1:size(shear_force_set, 1)
        if length(t(1,:)) == 4 || length(t(1,:)) == 3
            % 对于线性元素（每边2个节点）
            fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_T     %10.3E    %10.3E\r\n',...
            shear_force_set(i,1), shear_force_set(i,2), shear_force_set(i,3), -Shear_Stress, -Shear_Stress);
        elseif length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
            % 对于高阶元素（每边3个节点）
            fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_T     %10.3E    %10.3E    %10.3E\r\n',...
            shear_force_set(i,1), shear_force_set(i,2), shear_force_set(i,3), shear_force_set(i,4), -Shear_Stress, -Shear_Stress, -Shear_Stress);
        end
    end

elseif boundary_type == 3
    % 力矩边界条件处理
    % 获取指定边界线上的元素
    moment_force_set = [];

    % 从指定的Line中获取边界元素
    if Moment_Line_Number <= size(L, 3)
        current_line_set = L(:,:,Moment_Line_Number);

        % 添加元素索引列
        current_force_set = [zeros(size(current_line_set,1),1), current_line_set];

        % 替换0为NaN以便处理
        current_force_set(current_force_set == 0) = NaN;

        % 查找包含在当前力集中的元素
        for i=1:elem_num
            switch (length(t(1,:))) % 每个元素的节点数
                case {4,3} % 线性四边形或三角形（每边2个节点）
                    for j = 1:size(current_force_set,1)
                        if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ~isempty(find(t(i,:) == current_force_set(j,3), 1))
                                current_force_set(j,1) = i;
                            end
                        end
                    end
                case {6,8,9} % 二次三角形或四边形（每边3个节点用于力的应用）
                    for j = 1:size(current_force_set,1)
                         if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) && ~isnan(current_force_set(j,4)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,3), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,4), 1))
                                    current_force_set(j,1) = i;
                            end
                        end
                    end
            end
        end

        % 移除未找到元素或节点为NaN的行
        moment_force_set = current_force_set(~isnan(current_force_set(:,1)) & ~isnan(current_force_set(:,2)),:);
    end

    % 写入力矩边界条件
    fprintf(bnd,'       LINE_FORCE  %d\r\n', size(moment_force_set, 1));



    % 计算边界的高度范围
    y_coords = [];
    for i = 1:size(moment_force_set, 1)
        y_coords = [y_coords; p(moment_force_set(i,2), 2); p(moment_force_set(i,3), 2)];
        if length(t(1,:)) >= 6 && ~isnan(moment_force_set(i,4))
            y_coords = [y_coords; p(moment_force_set(i,4), 2)];
        end
    end

    y_min = min(y_coords);
    y_max = max(y_coords);
    height = y_max - y_min;
    center_y = (y_max + y_min) / 2;

    % 使用统一定义的几何参数（厚度）

    % 计算力矩所需的法向应力分布
    fprintf('Writing %d moment stress elements (magnitude: %5.2f N·m)\n', size(moment_force_set, 1), Moment_Value);

    for i = 1:size(moment_force_set, 1)
        % 获取节点的y坐标
        node1_y = p(moment_force_set(i,2), 2);
        node2_y = p(moment_force_set(i,3), 2);

        % 计算法向应力（F_N）以产生力矩
        % 经典梁理论：σ = M*y/I，其中I = b*h³/12
        % 对于矩形截面，σ = 12*M*y/(b*h³)，其中b是厚度，h是梁高
        % 应力的大小与到中心线的距离成正比

        % 确保height不为零，避免除以零错误
        if height < 1e-10
            height = 1e-10;  % 设置一个小的非零值
        end

        % 使用实际计算出的梁高度而不是预设值，确保计算更准确
        % 使用统一定义的厚度参数和我们之前记住的正确公式
        I = element_thickness * height^3 / 12; % 矩形截面的惯性矩

        % 确保I不为零，避免除以零错误
        if abs(I) < 1e-15
            I = 1e-15;  % 设置一个小的非零值
        end

        stress1 = Moment_Value * (node1_y - center_y) / I;
        stress2 = Moment_Value * (node2_y - center_y) / I;

        if length(t(1,:)) == 4 || length(t(1,:)) == 3
            fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_N     %10.3E    %10.3E\r\n',...
            moment_force_set(i,1), moment_force_set(i,2), moment_force_set(i,3), stress1, stress2);
        elseif length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
            % 对于高阶元素，需要计算中间节点的应力
            node3_y = p(moment_force_set(i,4), 2);

            % 使用与前面相同的惯性矩计算公式
            % 确保计算稳定性
            if abs(I) < 1e-15
                I = 1e-15;  % 设置一个小的非零值
            end

            stress3 = Moment_Value * (node3_y - center_y) / I;

            fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_N     %10.3E    %10.3E    %10.3E\r\n',...
            moment_force_set(i,1), moment_force_set(i,2), moment_force_set(i,3), moment_force_set(i,4), stress1, stress2, stress3);
        end
    end
elseif boundary_type == 4
    % 分布力P边界条件处理
    % 获取指定边界线上的元素
    p_force_set = [];

    % 从指定的Line中获取边界元素
    if P_Force_Line_Number <= size(L, 3)
        current_line_set = L(:,:,P_Force_Line_Number);

        % 添加元素索引列
        current_force_set = [zeros(size(current_line_set,1),1), current_line_set];

        % 替换0为NaN以便处理
        current_force_set(current_force_set == 0) = NaN;

        % 查找包含在当前力集中的元素
        for i=1:elem_num
            switch (length(t(1,:))) % 每个元素的节点数
                case {4,3} % 线性四边形或三角形（每边2个节点）
                    for j = 1:size(current_force_set,1)
                        if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ~isempty(find(t(i,:) == current_force_set(j,3), 1))
                                current_force_set(j,1) = i;
                            end
                        end
                    end
                case {6,8,9} % 二次三角形或四边形（每边3个节点用于力的应用）
                    for j = 1:size(current_force_set,1)
                         if ~isnan(current_force_set(j,2)) && ~isnan(current_force_set(j,3)) && ~isnan(current_force_set(j,4)) % 确保节点有效
                            if ~isempty(find(t(i,:) == current_force_set(j,2), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,3), 1)) && ...
                               ~isempty(find(t(i,:) == current_force_set(j,4), 1))
                                    current_force_set(j,1) = i;
                            end
                        end
                    end
            end
        end

        % 移除未找到元素或节点为NaN的行
        p_force_set = current_force_set(~isnan(current_force_set(:,1)) & ~isnan(current_force_set(:,2)),:);
    end

    % 写入分布力P边界条件
    fprintf(bnd,'       LINE_FORCE  %d\r\n', size(p_force_set, 1));

    % 计算边界的高度范围
    y_coords = [];
    for i = 1:size(p_force_set, 1)
        y_coords = [y_coords; p(p_force_set(i,2), 2); p(p_force_set(i,3), 2)];
        if length(t(1,:)) >= 6 && ~isnan(p_force_set(i,4))
            y_coords = [y_coords; p(p_force_set(i,4), 2)];
        end
    end

    y_min = min(y_coords);
    y_max = max(y_coords);
    height = y_max - y_min;

    % 计算分布力P所需的法向应力分布
    % 根据公式 p_y = -0.75P(1-y²)，其中y是归一化坐标
    fprintf('Writing %d distributed force P elements (magnitude: %5.2f N)\n', size(p_force_set, 1), P_Force_Value);

    for i = 1:size(p_force_set, 1)
        % 获取节点的y坐标
        node1_y = p(p_force_set(i,2), 2);
        node2_y = p(p_force_set(i,3), 2);

        % 将y坐标归一化到[-1, 1]范围
        % 确保height不为零，避免除以零错误
        if height < 1e-10
            height = 1e-10;  % 设置一个小的非零值
        end

        % 修正归一化计算
        norm_y1 = 2 * ((node1_y - y_min) / height) - 1;
        norm_y2 = 2 * ((node2_y - y_min) / height) - 1;

        % 确保归一化值在[-1, 1]范围内
        norm_y1 = max(-1, min(1, norm_y1));
        norm_y2 = max(-1, min(1, norm_y2));

        % 计算分布力
        % 修正公式: p_y = -0.75P(1-y²)
        % 使用更稳定的计算方式
        force1 = -0.75 * P_Force_Value * (1 - norm_y1*norm_y1);
        force2 = -0.75 * P_Force_Value * (1 - norm_y2*norm_y2);

        if length(t(1,:)) == 4 || length(t(1,:)) == 3
            fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_T     %10.3E    %10.3E\r\n',...
            p_force_set(i,1), p_force_set(i,2), p_force_set(i,3), force1, force2);
        elseif length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
            % 对于高阶元素，需要计算中间节点的力
            node3_y = p(p_force_set(i,4), 2);

            % 使用与前面相同的归一化方法
            norm_y3 = 2 * ((node3_y - y_min) / height) - 1;

            % 确保归一化值在[-1, 1]范围内
            norm_y3 = max(-1, min(1, norm_y3));

            % 使用相同的公式计算力
            force3 = -0.75 * P_Force_Value * (1 - norm_y3*norm_y3);

            fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_T     %10.3E    %10.3E    %10.3E\r\n',...
            p_force_set(i,1), p_force_set(i,2), p_force_set(i,3), p_force_set(i,4), force1, force2, force3);
        end
    end
end
% -------------------------------------------------------------------------
% FORCE BOUNDARY SECTION - END

% 为所有边界条件类型添加结束标记
fprintf(bnd,  '       END LINE_FORCE\r\n');
fprintf(bnd,  '    END FORCE\r\n');
fprintf(bnd,  '  END STEP\r\n');
fclose(bnd);

fprintf('Pre-processing complete. Output files: %s.ctr, %s.bnd\n', project, project);