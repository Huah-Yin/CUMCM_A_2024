%% 问题二多维函数遗传算法与非线性优化
%% 初始化
clear;
close all;      % 清图
clc;
NP=50;         % 种群数量
L=10;           % 单个维度二进制数段长度
Dim=6;          % 维度
Pc=0.5;         % 交叉率
Pm=0.2;         % 变异率
Mc=20;          % 交叉最多尝试次数
Mm=50;          % 变异最多尝试次数
G=100;         % 最大遗传次数
Xs=[350,15,20,6,8,8];              % 定义域上限
Xx=[0,5,10,2,2,2];            % 定义域下限
f=randi([0,1],NP,L*Dim);    % 随机获取初始种群
boxMax=20;      % 箱型图绘制代数
boxPage=20;     % 箱型图单图显示数据数
%% 遗传算法
tic
% 预声明变量
trace1(G)=0;
trace2(G)=0;
trace3(G)=0;
maxIndex(G,Dim)=0;
fitHis(NP,boxMax)=0;
if boxMax>G     % 纠错
    boxMax=g;
end
% 修改不符合约束条件的初始种群，直至全部符合条件为止
for i=1:NP
    while ~Q2_judge(b2d(f(i,:),Xx,Xs,L))
    f(i,:)=randi([0,1],1,L*Dim);
    end
end
for k=1:G
    fprintf('正在计算问题2第%d代...\n',k);
    %% 二进制解码为定义域内的十进制
    for i=1:NP
        tic
        U=f(i,:);   % 获取第i条染色体
        x(i,:)=b2d(U,Xx,Xs,L);
        Fit(i)=Q2_func_dim(x(i,:));          % 计算每个样本适应度
        fprintf('代数%d，种群序号%d，适应度%f\n',k,i,Fit(i))
        toc
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
            t=1;
            while t<Mc
                % 开始尝试交叉
                tnf=zeros(1,L*Dim);
                q=randi([0,1],1,L*Dim);     % 随机生成要交叉的基因位置
                for j=1:L*Dim
                    tnf(j)=nf(i+q(j),j);
                end
                if Q2_judge(b2d(tnf,Xx,Xs,L))  % 判断当前交叉结果是否符合约束条件
                    nf(i,:)=tnf(:);
                    break
                end
                t=t+1;
            end
        end
    end
    %% 基于概率的变异操作
    i=1;
    while i<=round(NP*Pm)   % 控制变异染色体总数
        h=randi([1,NP],1,1);    % 随机选择一个需要变异的染色体索引
        t=1;
        while t<Mm
            tnf=nf(h,:);
            for j=1:round(L*Dim*Pm) % 控制变异基因数
                g=randi([1,L*Dim],1,1); % 随机需要变异的基因索引
                tnf(g)=~tnf(g);   % 取反
            end
            if Q2_judge(b2d(tnf,Xx,Xs,L))  % 判断当前变异结果是否符合约束条件
                nf(i,:)=tnf(:);
                break
            end
            t=t+1;
        end
        i=i+1;
    end
    %% 下一代预备
    if k<=boxMax
        fitHis(:,k)=Fit';
    end
    f=nf;   % 种群位置搬家
    f(1,:)=fBest;   % 保留上代最优个体
    trace1(k)=maxFit;       % 历代最优适应度
    trace2(k)=mean(Fit);    % 平均适应度
    trace3(k)=minFit;       % 最小适应度
    maxIndex(k,:)=xBest(:);    % 最大适应度索引
end
toc
%% 绘制适应度过程图
figure
subplot(2,1,1);
hold on
plot(1:k,trace1(1:k),'^r-');
plot(1:k,trace2(1:k),'ob-');
plot(1:k,trace3(1:k),'vg-');
title('历代适应度');
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
%% 绘制条件范围内适应度箱型图
figure
page=ceil(boxMax/boxPage);
for i=1:page
    subplot(page,1,i);
    hold on
    boxplot(fitHis(:,(i-1)*boxPage+1:min([boxMax,i*boxPage])),(i-1)*boxPage+1:min([boxMax,i*boxPage]));
    title(strcat('第',num2str((i-1)*boxPage+1),'~',num2str(min([boxMax,i*boxPage])),'代适应度箱型图'));
    xlabel('进化代数');
    ylabel('适应度');
    plot([0,boxMax+1],[trace1(G),trace1(G)],'LineWidth',1.5);
    legend({strcat('最大值:',num2str(trace1(G)))},'Location','SouthEast');
    set(gca,'XTickLabelRotation',70);   % 倾斜x坐标轴数字
    hold off
end
%% 输出结果
fprintf('最大值：');
disp(trace1(G));
fprintf('最大值点：');
disp(maxIndex(G,:));
fprintf('精度：x±')
disp((Xs-Xx)/(2^L-1));