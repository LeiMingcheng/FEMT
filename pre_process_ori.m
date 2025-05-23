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

Pa = -1e5;  % distributed load  (left hand direction) line start----line end------positive force
Fnum = 3;   % Load line number (loads on nodes need to be set by User in '%% write to bnd')


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

    % read line set    chifan
    elseif strcmp(tmpline(max(strfind(tmpline,'='))+1:end-1),'Line')
    fprintf('\tRead line %d\n',line_num+1);
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
    elseif strcmp(tmpline(max(strfind(tmpline,'='))+1:end-1),'Surface')
    fprintf('\tRead elements ');
        elem_sta = numl+1;
        tmpline = fgetl(inp);numl = numl+1;
        while ~strncmp(tmpline, '*',1)
            ei=ei+1;
            T(ei, :) = str2num(tmpline);
            tmpline = fgetl(inp);numl = numl+1;
        end
        temp = elem_num;
        elem_num = elem_num + numl-1-elem_sta+1;
    fprintf('%d to %d\n',temp+1,elem_num);

    elseif strcmp(strtok(tmpline,','),'*NSET')

        if strcmp(tmpline(20:23),'Line')
            if Nset_num+1 == 1
                fprintf('\tRead boundary set(U=0) ');
            elseif Nset_num+1 == 2
                fprintf('\tRead boundary set(V=0) ');
            else
                fprintf('\tRead %d set  ',Nset_num-1);
            end
            Nset_num = Nset_num+1;
            j = 0;
            tmpline = fgetl(inp);numl = numl+1;
            while ~strncmp(tmpline, '*',1)
                temp = str2num(tmpline);
                while ~isempty(temp)
                   j = j+1;
                   Nset(Nset_num,j) = temp(1);
                   temp = temp(2:end);
                end
                tmpline = fgetl(inp);numl = numl+1;
            end
        Nset_size(Nset_num) = j;
        fprintf('%d nodes\n',j);

%         Line(:,1,line_num) = [1:line_end(line_num)-line_sta(line_num)+1]';
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
L = Line(:,2:end,:);
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
title('������Ż�ǰ��������նȾ���')
axis equal;
axis([0 length(p) 0 length(p)])

subplot(1,2,2)
q = symrcm(K); %���� K �ĶԳ� Cuthill-McKee ����q��ʹ K �ķ� 0 Ԫ�ؼ��������Խ��߸���
spy(K(q,q));
title('������Ż��󡪡�����նȾ���')
axis equal;
axis([0 length(p) 0 length(p)])
clear K;
% use new permutation on P & T
% permutate T
p = p(q, :);
% sequence of q
q1 = q * 0;
for i = 1 : length(q)
    q1(q(i)) = i;
end
% permutate P
[m, n] = size(t);
t(:, n) = [];
for i = 1 : m
    for j = 1 : n - 1
        t(i, j) = q1(t(i, j));
    end
end
% permutate L
[m, n, o] = size(L);
for k = 1 : o
    for i = 1 : m
        for j = 1 : n
            if L(i, j, k) ~= 0
                L(i, j, k) = q1(L(i, j, k));
            else
                L(i, j, k) = 0;
            end
        end
    end
end
% permutate Nset
[m,n] = size(Nset);
for i = 1 : m
    for j = 1 : n
        if Nset(i,j) ~= 0
            Nset(i,j) = q1(Nset(i, j));
        end
    end
end

% new semiband width
[m, n] = size(t);
for i = 1 : m
    tmp = max(t(i, :)) - min(t(i, :));
    if tmp > semiBandwidthOpt
        semiBandwidthOpt = tmp;
    end
end

semiBandwidthOpt = (semiBandwidthOpt + 1) * 2;

%% view mesh
%    viewmsh(p,t);

%% set the Fset

fprintf('Force %5.2f applied on No.%d line (Left hand direction)\n',Pa,Fnum);
Fset = L(:,:,Fnum);
Fset = [Fset(:,1),Fset];
Fset(Fset == 0) = NaN;

% find element included in Fset
for i=1:elem_num
switch (length(t(1,:)))
    case {4,3}
        for j = 1:length(Fset)
            if ~isempty(find(t(i,:) == Fset(j,2), 1)) && ~isempty(find(t(i,:) == Fset(j,3), 1))
                    Fset(j,1) = i;
            end
        end
    case {6,8,9}
        for j = 1:length(Fset)
            if ~isempty(find(t(i,:) == Fset(j,2), 1)) && ...
               ~isempty(find(t(i,:) == Fset(j,3), 1)) && ...
               ~isempty(find(t(i,:) == Fset(j,4), 1))
                    Fset(j,1) = i;
            end
        end
end
end

%% write to ctr and bnd
fprintf('Write to ctr\n');

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
fprintf(ctr,'       E  2.1E9\r\n');
fprintf(ctr,'       v  0.3\r\n');
fprintf(ctr,'    END MATERIAL\r\n');
fprintf(ctr,' END MATERIALS\r\n');
fprintf(ctr, '\r\n');


%-- GEOMETRIES --
fprintf(ctr,' GEOMETRIES 1\r\n');
fprintf(ctr,'    GEOMETRY 1 TYPE PLANE\r\n');
fprintf(ctr,'       THICKNESS 1\r\n');
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

%         % Line boundary (Default)
        fprintf(bnd,  '       GIVEN_DISP    %d\r\n',Nset_size(1) + Nset_size(2));

        for i = 1:Nset_size(1)
            fprintf(bnd,'        %4d   U  0.0\r\n', Nset(1,i));
        end
        for i = 1:Nset_size(2)
            fprintf(bnd,'        %4d   V  0.0\r\n', Nset(2,i));
        end

        % Node boundary (User defined)
%         fprintf(bnd,  '       GIVEN_DISP    %d\r\n',Nset_size(1) + 1);   % total number of given displacement
%
%         for i = 1:Nset_size(1)
%             fprintf(bnd,'        %4d   U  0.0\r\n', Nset(1,i));
%         end
%         for i = 1:1
%             fprintf(bnd,'        %4d   V  0.0\r\n', Nset(2,end));
%         end


fprintf(bnd,  '       END GIVEN_DISP\r\n');
fprintf(bnd,  '    END DISPLACEMENT\r\n');
fprintf(bnd, '\r\n');
fprintf(bnd,  '    FORCE\r\n');


% LINE FORCE
fprintf(bnd,'       LINE_FORCE  %d\r\n', line_end(Fnum)-line_sta(Fnum)+1);
for i = 1: line_end(Fnum)-line_sta(Fnum)+1
    if length(t(1,:)) == 4 || length(t(1,:)) == 3
        fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_N     %4.1f    %4.1f\r\n',...
        Fset(i,1), Fset(i,2), Fset(i,3), Pa, Pa);
    elseif  length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
        fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_N     %4.1f    %4.1f    %4.1f\r\n',...
        Fset(i,1), Fset(i,2), Fset(i,3), Fset(i,4) ,Pa , Pa, Pa);
    end
end

fprintf(bnd,  '       END LINE_FORCE\r\n');
fprintf(bnd,  '    END FORCE\r\n');
fprintf(bnd,  '  END STEP\r\n');
fclose(bnd);
