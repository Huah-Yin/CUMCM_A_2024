clc;
clear;
close all;
%%%%%%%%%%%%问题一

%%附件4

% H_f = 3; %浮子圆柱部分高度 (m)
% h_f = 0.8; %浮子圆锥部分高度 (m)

% rz = 0.5; %振子半径 (m)
% h_z = 0.5; %振子高度 (m)

% L0 = 0.5; %弹簧原长 (m)
% k1 = 8890.7; %静水恢复力矩系数 (N·m)
%%%parmas 参数.
mf = 4866; %浮子质量:kg
rf = 1; %浮子底半径：m
mz = 2433; %振子质量 (kg)
rou = 1025; %海水的密度 (kg/m3)
g = 9.8;
k4 = 656.3616; %垂荡兴波阻尼系数 (N·s/m)
m0 = 1335.535; %垂荡附加质量
k2 = 80000; %弹簧劲度系数
A = 6250; %垂荡激励力振幅 (N)
omega = 1.4005; %入射波浪频率 (s-1)
%%%求解
syms x1(t) x2(t) X1 X2 s t
detK = 10000;
% F1=-rou*g*pi*(rf^2)*x1;
equ1 = (m0 + mf) * diff(x1, t, 2) + rou * g * pi * (rf ^ 2) * x1 + k2 * (x1 - x2) +detK * (diff(x1, t, 1) - diff(x2, t, 1)) + k4 * diff(x2, t, 1) + A * cos(omega * t) == 0;
equ2 = mz * diff(x2, t, 2) + k2 * (x2 - x1) + detK * (diff(x2, t, 1) - diff(x1, t, 1)) == 0;

%%%laplace
X1_s = laplace(equ1, t, s);
lap1 = {laplace(x1(t), t, s), diff(x1(t), t,1), x1(0)};
val0 = {X1, 0, 0};
X1_s = subs(X1_s, lap1, val0);
disp(X1_s);

X2_s = laplace(equ2, t, s);
lap2 = {laplace(x2(t), t, s), diff(x2(t), t), x2(0)};
val1 = {X2, 0, 0};
X2_s = subs(X2_s, lap2, val1);
disp(X2_s);

% 求解频域方程
[X1_solved, X2_solved] = solve([X1_s, X2_s], [X1, X2]);
% 进行拉普拉斯逆变换
x1_solved = ilaplace(X1_solved, s, t);
x2_solved = ilaplace(X2_solved, s, t);
% 显示结果
disp('解 x1(t):');
disp(x1_solved);
disp('解 x2(t):');
disp(x2_solved);
