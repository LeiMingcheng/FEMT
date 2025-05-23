# FEMT - Finite Element Method Toolkit

## 项目简介

FEMT (Finite Element Method Toolkit) 是一个计算力学基础有限元分析系统，包含前处理、计算和后处理三个主要部分。该项目原本仅支持Intel Fortran编译器，本版本进行了修改，使其同时支持：

1. 九节点四边形单元 (RECTANGLE9)
2. GNU Fortran (gfortran) 编译器

## 主要特性

- 支持多种单元类型：
  - 三节点三角形单元 (TRIANGLE3)
  - 六节点三角形单元 (TRIANGLE6)
  - 四节点四边形单元 (RECTANGLE4)
  - 八节点四边形单元 (RECTANGLE8)
  - 九节点四边形单元 (RECTANGLE9) - 新增支持

- 完整的有限元分析流程：
  - 前处理：网格生成与边界条件设置
  - 计算程序：有限元分析求解
  - 后处理：结果可视化与分析

- 编译器兼容性：
  - 支持Intel Fortran编译器
  - 支持GNU Fortran (gfortran) 编译器

## 项目结构

- `src/`: 源代码目录
  - 计算程序 (Fortran)
  - 前处理程序 (MATLAB)
  - 后处理程序 (MATLAB)

## 九节点单元实现

为支持九节点四边形单元 (RECTANGLE9)，对程序进行了以下主要修改：

1. 计算程序：
   - 添加了RECTANGLE9元素类型定义
   - 实现了9节点Lagrangian形状函数及其导数
   - 修改了积分点初始化、刚度矩阵计算等核心功能

2. 前处理程序：
   - 添加了9节点单元的网格生成支持
   - 修改了边界条件应用方式

3. 后处理程序：
   - 添加了9节点单元的结果可视化支持
   - 优化了网格信息显示

## 使用说明

1. 前处理：使用MATLAB运行前处理程序生成网格和边界条件
2. 计算：使用Fortran编译器编译并运行计算程序
3. 后处理：使用MATLAB运行后处理程序可视化结果

## 编译说明

### 使用gfortran编译

```bash
gfortran -o femt *.f90
```

### 使用Intel Fortran编译

```bash
ifort -o femt *.f90
```

## 许可证

[待添加]
