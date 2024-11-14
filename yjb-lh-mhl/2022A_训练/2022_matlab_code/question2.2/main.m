clc;
clear ;
close all;
%遗传算法

%%问题定义
problem.costFunction = @(x) test(x);
problem.nVar = 10; %问题的变量个数
problem.Min = [0 0 0 0 0 0 0 0 0 0]; %变量的下界
problem.Max = [100000 100000 100000 100000 100000 100000 100000 100000 100000 100000]; %变量的上界
%%参数
params.n_pop = 100; %种群数量
params.MaxGen = 100; %种群最大迭代数目
params.Beta = 1; %e的负指数函数的参数
params.PCrossover = 0.81; %交叉概率
params.mu = 0.01; %变异率
params.MySigma = 0.01;
params.MyGama = 0.01;
%%开始计时
tic;
%%Run

out = RunGA(problem, params);

%%结果可视化
figure;
disp(out.bestsolu);
plot(out.bestcost);
semilogy(out.bestcost);
xlabel('迭代次数');
ylabel('最佳效益');
grid on;

% 结束计时并获取时间
elapsedTime = toc;
% 显示运行时间
disp(['遗传算法运行时间: ' num2str(elapsedTime) ' s']);
