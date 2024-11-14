%% 问题2-二维遗传算法
%% 初始化
clear;
close all;      % 清图
clc;
NP=100;         % 种群数量
L=30;           % 单个维度二进制数段长度
Dim=2;          % 维度
Pc=0.5;         % 交叉率
Pm=0.2;         % 变异率
G=200;          % 最大遗传次数
for miz=[0 1]
Xs=[100000,miz];              % 定义域上限
Xx=[0,0];            % 定义域下限
f=randi([0,1],NP,L*Dim);    % 随机获取初始种群

%% 遗传算法
tic
% 预声明变量
trace1(G)=0;
trace2(G)=0;
trace3(G)=0;
maxIndex(G,Dim)=0;

for k=1:G
    fprintf('正在计算问题2(%d)第%d代...',miz+1,k);
    %% 二进制解码为定义域内的十进制
    for i=1:NP
        U=f(i,:);   % 获取第i条染色体
        x(i,:)=b2d(U,Xx,Xs,L);
        Fit(i)=Q2_func(x(i,:));          % 计算每个样本适应度
    end
    maxFit=max(Fit);    % 最大适应度
    minFit=min(Fit);    % 最小适应度
    rr=find(maxFit==Fit);   % 找出最大适应度所在个体索引，返回数或数组
    fBest=f(rr(1,1),:);     % 历代最优个体
    xBest=x(rr(1,1),:);       % 最优个体对应的十进制染色体
    Fit1=(Fit-minFit)/(maxFit-minFit);   % 适应度归一化
    
    %% 基于轮盘赌选择法的复制操作
    sum_Fit=sum(Fit1);   % 计算种群中所有个体适应值的和
    fitvalue=Fit1/sum_Fit;   % 计算每个种群的选择概率
    fitvalue=(cumsum(fitvalue));    % 计算每个种群的累计概率
    ms=sort(rand(NP,1));    % 随机生成NP个(0,1)的值，并排序
    fiti=1;     % 旧种群当前指针
    newi=1;     % 新种群当前指针
    while newi<=NP      % 随机复制，并使适应度大的遗传下去
        if ms(newi) < fitvalue(fiti)
            nf(newi,:)=f(fiti,:);   % 复制
            newi=newi+1;
        else
            fiti=fiti+1;
            if fiti>NP
                break;
            end
        end
    end
    
    %% 基于概率的交叉操作
    for i=1:NP-1
        p=rand;     % 随机生成一个处于[0,1]的概率p
        if p<Pc     % 满足交叉条件
            q=randi([0,1],1,L*Dim);     % 随机生成要交叉的基因位置
            for j=1:L*Dim
                if q(j)==1  % 以下代码实现变量交换
                    temp=nf(i+1,j);
                    nf(i+1,j)=nf(i,j);
                    nf(i,j)=temp;
                end
            end
        end
    end
    
    %% 基于概率的变异操作
    i=1;
    while i<=round(NP*Pm)   % 控制变异染色体总数
        h=randi([1,NP],1,1);    % 随机选择一个需要变异的染色体索引
        for j=1:round(L*Dim*Pm) % 控制变异基因数
            g=randi([1,L*Dim],1,1); % 随机需要变异的基因索引
            nf(h,g)=~nf(h,g);   % 取反
        end
        i=i+1;
    end
    
    %% 下一代预备
    f=nf;   % 种群位置搬家
    f(1,:)=fBest;   % 保留上代最优个体
    trace1(k)=maxFit;       % 历代最优适应度
    trace2(k)=mean(Fit);    % 平均适应度
    trace3(k)=minFit;       % 最小适应度
    maxIndex(k,:)=xBest(:);    % 最大适应度索引
    fprintf('最大适应度:%f\n',maxFit);
end
toc

%% 绘制适应度过程图
figure
subplot(2,1,1);
hold on
plot(1:k,trace1(1:k),'^r-');
plot(1:k,trace2(1:k),'ob-');
plot(1:k,trace3(1:k),'vg-');
title(strcat('第',num2str(miz+1),'小问历代适应度'));
xlabel('进化代数');
ylabel('适应度');
legend({'最大值','平均值','最小值'},'Location','SouthEast');
hold off
subplot(2,1,2);
hold on
plot(1:G,trace1(1:G),'^r-');
plot([0,G],[trace1(G),trace1(G)],'LineWidth',1.5);
yl=ylim;
for i=1:G
    if trace1(i)==trace1(G)
        plot(i,trace1(G),'b*','LineWidth',1.5,'MarkerSize',6);
        text(i,trace1(G)-(yl(2)-yl(1))/30,strcat(num2str(i),'代'));
        break
    end
end
title('历代最大适应度');
xlabel('进化代数');
ylabel('适应度最大值');
legend({'历代最大值','全局最大值','初次达到最大值点'},'Location','SouthEast');
hold off

%% 输出结果
fprintf('==========问题2第%d小问==========\n',miz+1)
fprintf('最大值：');
disp(trace1(G));
fprintf('最大值点：');
disp(maxIndex(G,:));
fprintf('精度：x±')
disp((Xs-Xx)/(2^L-1));
end