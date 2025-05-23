# 9节点网格支持实现总结

为了在FEMT有限元分析系统中添加对9节点网格的支持，我们对整个系统的三个主要部分（计算程序、前处理程序和后处理程序）进行了一系列修改。以下是所有修改的详细总结：

## 一、计算程序修改

### 1. 基础数据定义修改
- **BASIC_DATA.F90**：
  - 将`NUM_ELE_TYPE`常量从4增加到5，以包含新的元素类型
  ```fortran
  INTEGER(INT_KIND), PARAMETER:: NUM_ELE_TYPE = 5  ! THE TOTAL NUMBER OF ELEMENT TYPES
  ```

### 2. 元素库初始化修改
- **INIT_ELE_LIB.F90**：
  - 添加了RECTANGLE9元素类型的定义
  ```fortran
  CASE (5)
     ELEMENT_LIB(I)%ELEMENT_TYPE = 'RECTANGLE9'
     ELEMENT_LIB(I)%NUM_NODE     = 9
     ELEMENT_LIB(I)%NUM_INTP     = 9
     ELEMENT_LIB(I)%DISP_DOF     = 2
     ELEMENT_LIB(I)%NUM_STIFF    = 171
  ```

### 3. 积分点初始化修改
- **INIT_INTEGRATION.F90**：
  - 修改了积分点初始化，使RECTANGLE9使用与RECTANGLE8相同的积分点
  ```fortran
  CASE ('RECTANGLE8', 'RECTANGLE9')
     N_GAUSS = 0
     DO I = 1, 3          ! KESI
        DO J = 1, 3       ! ETA
           N_GAUSS = N_GAUSS + 1
           ELEMENT_LIB%INTE_COORD(N_GAUSS)%COORD(1) = GAUSS_COORD3(I)
           ELEMENT_LIB%INTE_COORD(N_GAUSS)%COORD(2) = GAUSS_COORD3(J)
           ELEMENT_LIB%INTE_COORD(N_GAUSS)%WEIGHT   = GAUSS_WEIGHT3(I) * GAUSS_WEIGHT3(J)
        ENDDO
     ENDDO
  ```

### 4. 形状函数修改
- **SHAPE_2D.F90**：
  - 修改了`NODE_SIGN`数组的定义，以正确的顺序和Fortran语法排列9个节点的坐标符号，确保Q4、Q8和Q9单元都能正确索引。角点（1-4）、Q8边中间节点（5-8）和Q9中心节点（9）现在按此顺序排列：
  ```fortran
  REAL(REAL_KIND), PARAMETER::    NODE_SIGN(2,9) = RESHAPE( &
                                (/ &
                                 -1.0, -1.0,  &  ! Node 1: (-1,-1)
                                  1.0, -1.0,  &  ! Node 2: ( 1,-1)
                                  1.0,  1.0,  &  ! Node 3: ( 1, 1)
                                 -1.0,  1.0,  &  ! Node 4: (-1, 1)
                                  0.0, -1.0,  &  ! Node 5: ( 0,-1) (Q8 midside)
                                  1.0,  0.0,  &  ! Node 6: ( 1, 0) (Q8 midside)
                                  0.0,  1.0,  &  ! Node 7: ( 0, 1) (Q8 midside)
                                 -1.0,  0.0,  &  ! Node 8: (-1, 0) (Q8 midside)
                                  0.0,  0.0   &  ! Node 9: ( 0, 0) (Q9 center)
                                 /), &
                                (/2, 9/))
  ```
  - 更新了8节点单元（Q8）部分的注释，以阐明`NODE_SIGN`数组的使用方式：
  ```fortran
      ! This section is for 8-node Serendipity element (Q8).
      ! NODE_SIGN array indices 1-8 now correctly correspond to Q8 node numbering (corners then midsides).
      ! Index 9 is for the Q9 center node.
  ```
  - 添加了9节点元素的形状函数计算，使用标准的Lagrangian多项式方法
  ```fortran
  ! FOR 9-NODE LAGRANGIAN ELEMENT (Q9)
  IF(NUM_NODE .EQ. 9) THEN
     ! Define 1D Quadratic Lagrangian Polynomials and their derivatives
     ! For nodes at zeta = -1, 0, 1
     L_neg1_ksi = 0.5 * ksi * (ksi - 1.0)
     L_zero_ksi = (1.0 - ksi**2)
     L_pos1_ksi = 0.5 * ksi * (ksi + 1.0)

     L_neg1_eta = 0.5 * eta * (eta - 1.0)
     L_zero_eta = (1.0 - eta**2)
     L_pos1_eta = 0.5 * eta * (eta + 1.0)

     dL_neg1_ksi = 0.5 * (2.0 * ksi - 1.0)
     dL_zero_ksi = -2.0 * ksi
     dL_pos1_ksi = 0.5 * (2.0 * ksi + 1.0)

     dL_neg1_eta = 0.5 * (2.0 * eta - 1.0)
     dL_zero_eta = -2.0 * eta
     dL_pos1_eta = 0.5 * (2.0 * eta + 1.0)

     ! Corner Nodes (1-4) - Standard Lagrangian
     ! Node 1: (-1,-1)
     SHAPES(1,I) = L_neg1_ksi * L_neg1_eta
     ! Node 2: (1,-1)
     SHAPES(2,I) = L_pos1_ksi * L_neg1_eta
     ! Node 3: (1,1)
     SHAPES(3,I) = L_pos1_ksi * L_pos1_eta
     ! Node 4: (-1,1)
     SHAPES(4,I) = L_neg1_ksi * L_pos1_eta

     ! Mid-side Nodes (5-8) - Standard Lagrangian
     ! Node 5: (0,-1)
     SHAPES(5,I) = L_zero_ksi * L_neg1_eta
     ! Node 6: (1,0)
     SHAPES(6,I) = L_pos1_ksi * L_zero_eta
     ! Node 7: (0,1)
     SHAPES(7,I) = L_zero_ksi * L_pos1_eta
     ! Node 8: (-1,0)
     SHAPES(8,I) = L_neg1_ksi * L_zero_eta

     ! Center node (9) - Standard Lagrangian
     SHAPES(9,I) = L_zero_ksi * L_zero_eta ! Equivalent to (1.0 - ksi**2) * (1.0 - eta**2)
  ENDIF
  ```

  - 添加了9节点元素的形状函数导数计算，同样使用标准的Lagrangian多项式方法
  ```fortran
  ! FOR 9-NODE LAGRANGIAN ELEMENT (Q9)
  IF(NUM_NODE .EQ. 9) THEN
     ! Define 1D Quadratic Lagrangian Polynomials and their derivatives again for safety / clarity
     ! For nodes at zeta = -1, 0, 1
     L_neg1_ksi = 0.5 * ksi * (ksi - 1.0)
     L_zero_ksi = (1.0 - ksi**2)
     L_pos1_ksi = 0.5 * ksi * (ksi + 1.0)

     L_neg1_eta = 0.5 * eta * (eta - 1.0)
     L_zero_eta = (1.0 - eta**2)
     L_pos1_eta = 0.5 * eta * (eta + 1.0)

     dL_neg1_ksi = 0.5 * (2.0 * ksi - 1.0)
     dL_zero_ksi = -2.0 * ksi
     dL_pos1_ksi = 0.5 * (2.0 * ksi + 1.0)

     dL_neg1_eta = 0.5 * (2.0 * eta - 1.0)
     dL_zero_eta = -2.0 * eta
     dL_pos1_eta = 0.5 * (2.0 * eta + 1.0)

     ! Corner Nodes (1-4) - Standard Lagrangian Derivatives
     ! Node 1: (-1,-1)
     D_SHAPE(1,1,I) = dL_neg1_ksi * L_neg1_eta
     D_SHAPE(2,1,I) = L_neg1_ksi * dL_neg1_eta
     ! Node 2: (1,-1)
     D_SHAPE(1,2,I) = dL_pos1_ksi * L_neg1_eta
     D_SHAPE(2,2,I) = L_pos1_ksi * dL_neg1_eta
     ! Node 3: (1,1)
     D_SHAPE(1,3,I) = dL_pos1_ksi * L_pos1_eta
     D_SHAPE(2,3,I) = L_pos1_ksi * dL_pos1_eta
     ! Node 4: (-1,1)
     D_SHAPE(1,4,I) = dL_neg1_ksi * L_pos1_eta
     D_SHAPE(2,4,I) = L_neg1_ksi * dL_pos1_eta

     ! Mid-side Nodes (5-8) - Standard Lagrangian Derivatives
     ! Node 5: (0,-1)
     D_SHAPE(1,5,I) = dL_zero_ksi * L_neg1_eta
     D_SHAPE(2,5,I) = L_zero_ksi * dL_neg1_eta
     ! Node 6: (1,0)
     D_SHAPE(1,6,I) = dL_pos1_ksi * L_zero_eta
     D_SHAPE(2,6,I) = L_pos1_ksi * dL_zero_eta
     ! Node 7: (0,1)
     D_SHAPE(1,7,I) = dL_zero_ksi * L_pos1_eta
     D_SHAPE(2,7,I) = L_zero_ksi * dL_pos1_eta
     ! Node 8: (-1,0)
     D_SHAPE(1,8,I) = dL_neg1_ksi * L_zero_eta
     D_SHAPE(2,8,I) = L_neg1_ksi * dL_zero_eta

     ! Center node (9) - Standard Lagrangian Derivatives
     D_SHAPE(1,9,I) = dL_zero_ksi * L_zero_eta ! Equivalent to -2.0 * ksi * (1.0 - eta**2)
     D_SHAPE(2,9,I) = L_zero_ksi * dL_zero_eta ! Equivalent to -2.0 * eta * (1.0 - ksi**2)
  ENDIF
  ```

### 5. 元素读取修改
- **READ_ELEMENTS.F90**：
  - 添加了对RECTANGLE9元素类型的识别
  ```fortran
  CASE ('RECTANGLE9')
     NUMBER = 5
  ```

### 6. 形状函数初始化修改
- **INIT_SHAPE.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6', 'RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        !-- ALLOCATE MEMORY FOR SHAPE FUNCTION --
  ```

### 7. 力学初始化修改
- **INIT_MECHANICAL.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(P_ELE_TYPE%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6', 'RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9') ! 2D ELEMENTS
        NUM_STRAIN = 3
  ```

### 8. 刚度矩阵计算修改
- **STIFFNESS.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6', 'RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        CALL STIFFNESS_PLANE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE, OPTION)
  ```

### 9. 线力计算修改
- **READ_LINE_FORCE.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELE_TYPE))
     CASE ('TRIANGLE3', 'RECTANGLE4')
        NUM_NODE = 2

     CASE ('TRIANGLE6', 'RECTANGLE8', 'RECTANGLE9')
        NUM_NODE = 3
  ```

### 10. 应变计算修改
- **GET_STRAIN.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6', 'RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        CALL STRAIN_GAUSS_PLANE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE)
  ```

  - 在节点应变计算中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6')
        CALL STRAIN_NODAL_TRIANGLE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE, ADJ_ELE)

     CASE ('RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        CALL STRAIN_NODAL_RECTANGLE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE, ADJ_ELE)
  ```

### 11. 应力计算修改
- **GET_STRESS.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6', 'RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        CALL STRESS_GAUSS_PLANE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_INTP, OPTION(2))
  ```

  - 在节点应力计算中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'TRIANGLE6')
        CALL STRESS_NODAL_TRIANGLE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE, ADJ_ELE, OPTION(2))

     CASE ('RECTANGLE4', 'RECTANGLE8', 'RECTANGLE9')
        CALL STRESS_NODAL_RECTANGLE(I_ELEMENT, ELEMENT_LIB(ELEMENTS(I_ELEMENT)%ELEMENT)%NUM_NODE, ADJ_ELE, OPTION(2))
  ```

### 12. 力计算修改
- **GET_FORCE.F90**：
  - 在元素类型选择中添加了RECTANGLE9
  ```fortran
  SELECT CASE(TRIM(ELEMENT_LIB(ELE_TYPE)%ELEMENT_TYPE))
     CASE ('TRIANGLE3', 'RECTANGLE4')
        CALL LINE_FORCE_PLANE(LINE_FORCE(I), 2)

     CASE ('TRIANGLE6', 'RECTANGLE8', 'RECTANGLE9')
        CALL LINE_FORCE_PLANE(LINE_FORCE(I), 3)
  ```

## 二、前处理程序修改

### 1. 元素类型处理修改
- **pre_process.m**：
  - 添加了对9节点元素的支持
  ```matlab
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
  ```

### 2. 调试信息修改
- **pre_process.m**：
  - 添加了对9节点元素的调试信息
  ```matlab
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
  ```

### 3. 力的应用修改
- **pre_process.m**：
  - 修改了力的应用部分以支持9节点元素
  ```matlab
  if length(t(1,:)) == 4 || length(t(1,:)) == 3
      fprintf(bnd,'             %6d    SIDE    %6d    %6d    F_N     %4.1f    %4.1f\r\n',...
      Fset(i,1), Fset(i,2), Fset(i,3), Pa, Pa);
  elseif  length(t(1,:)) == 6 || length(t(1,:)) == 8 || length(t(1,:)) == 9
      fprintf(bnd,'             %6d    SIDE    %6d    %6d    %6d    F_N     %4.1f    %4.1f    %4.1f\r\n',...
      Fset(i,1), Fset(i,2), Fset(i,3), Fset(i,4) ,Pa , Pa, Pa);
  end
  ```

## 三、后处理程序修改

### 1. 元素连接矩阵处理修改
- **post_process.m**：
  - 添加了对9节点元素的支持
  ```matlab
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
  ```

### 2. 添加网格信息显示功能
- **post_process.m**：
  - 添加了网格信息显示功能
  ```matlab
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
  ```

### 3. 添加程序说明
- **post_process.m**：
  - 添加了对9节点元素支持的说明
  ```matlab
  % post_processing program for FEMT
  %
  % input: PROGRAM.CTR, PROGRAM.OUT
  %
  % Last Edited by  GAO HEXUAN at 2016-2-25
  % Modified to support 9-node elements
  ```

## 总结

通过以上修改，我们成功地为FEMT有限元分析系统添加了对9节点网格的支持。这些修改涵盖了系统的三个主要部分：计算程序、前处理程序和后处理程序。

9节点四边形元素（RECTANGLE9）是一种高阶元素，相比于4节点和8节点元素，它能够更准确地模拟复杂的应力场和变形，特别是在应力集中区域。通过添加对9节点元素的支持，FEMT系统的分析能力得到了显著提升。