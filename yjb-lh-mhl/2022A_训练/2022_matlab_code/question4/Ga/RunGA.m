function out = RunGA(problem, params)
    %problem
    CostFunction = problem.costFunction; %需要优化的problem
    nVar = problem.nVar;
    VarMin = problem.Min;
    VarMax = problem.Max;
    VarSize = [1, nVar];

    %params
    MaxGener = params.MaxGen; %最大迭代次数
    nPop = params.n_pop; %种群数量
    MyBeta = params.Beta; %参数传递
    MySigma = params.MySigma;
    MyGama = params.MyGama;
    P_crossover = params.PCrossover; %交叉概率
    %crossover=1，保证2个父代产生两个子代
    N_crossover = round(P_crossover * nPop / 2) * 2; %矩阵乘积并四舍五入
    mu = params.mu; %变异率

    %%初始化
    empty_individual.position = [];
    empty_individual.Cost = [];
    %%Best Solution

    % 初始值矩阵
    initial_position = zeros(100, 2);

    % 设置范围和步长
    r_d = 10000:10000:100000; % 第一列的值
    c_d = 10000:10000:100000;
    % 生成所有组合
    [R, C] = meshgrid(r_d, c_d);

    % 将组合转换为列向量
    initial_position = [R(:), C(:)]; % 100x2的矩阵，包含所有组合

    best_solu.Cost = inf; %无穷大

    %填写空数组
    pop = repmat(empty_individual, nPop, 1);

    for i = 1:nPop

        if i <= nPop / 3
            idx_i = mod(i, size(initial_position, 1)) + 1;
            pop(i).position = initial_position(idx_i, :);
        else
            %%随机生成过程
            pop(i).position = unifrnd(VarMin, VarMax, VarSize);
        end

        %%评价方案
        pop(i).Cost = CostFunction(pop(i).position);

        %%寻找最佳方案
        if pop(i).Cost < best_solu.Cost
            best_solu = pop(i); %最佳基因个体
        end

    end

    %每次迭代结束时保留最佳成本的记录
    BestCost_matrix = zeros(MaxGener, 1); %创建空matrix

    %%%初始化 相同代数的计数器
    status_num = 0;

    %主循环
    for j = 1:MaxGener

        %选择概率,采取e的负次幂指数函数
        c = [pop.Cost];
        average = mean(c); %求此矩阵的均值
        %average后面要做分母得进行判0操作.
        if average ~= 0
            c = c / average; %标准化
        end

        probs = exp(-MyBeta * c); %负指数函数

        %初始化种群
        popc = repmat(empty_individual, N_crossover / 2, 2);

        %crossover
        for k = 1:N_crossover / 2
            %select 父母代分别各选择一次,轮盘赌法则

            p1 = pop(WheelSelection(probs));
            p2 = pop(WheelSelection(probs));

            %执行crossover
            [popc(k, 1).position, popc(k, 2).position] = UnformCrossover(p1.position, p2.position, MyGama);

        end

        %evaluate 评估
        %将popc转化为单列矩阵
        popc = popc(:);

        %mutation ->变异,N_crossover为产生子代的数量，必须保证每个子代都有变异的机会
        for l = 1:N_crossover
            %变异函数
            popc(l).position = mutation(popc(l).position, mu, MySigma);
            %检查上下界,如果比最小值小，替换为下界，
            pop(l).position = max(pop(l).position, VarMin);
            %相反，如果比最大值大,则替换为上界
            pop(l).position = min(pop(l).position, VarMax);
            %成本函数,evaluate
            popc(l).Cost = CostFunction(popc(l).position);

            %%最佳方案
            if popc(l).Cost < best_solu.Cost
                best_solu = popc(l);
            end

        end

        %排序
        pop = Sort([pop; popc]); %自己定义排序
        %删除额外的成员个体.
        pop = pop(1:nPop); %垂直连接矩阵
        %pop = pop(l:nPop);
        %%更新迭代的最佳成本
        BestCost_matrix(j) = best_solu.Cost;

        if j > 1
            % 检查是否连续5代最优结果相同
            if abs((BestCost_matrix(j) - BestCost_matrix(j - 1))) <= 0.01
                status_num = status_num + 1;
                fprintf('status_num = %d  mu=%0.2f MyGama=%0.2f\n', status_num, mu, MyGama);
            else
                status_num = 0; % 重置计数器
                mu = params.mu; %重置变异系数
                MyGama = params.MyGama; %重置交叉系数
                fprintf('status_num = %d  mu=%0.2f   MyGama=%0.2f\n', status_num, mu, MyGama);
            end

        end

        if (status_num >= 5)
            staus = round(status_num / 5) + status_num - 5;
            MyGama = get_crossover(staus, MyGama);
            mu = get_mutation(staus, mu);

        end

        %打印提示信息
        disp(['迭代次数：' num2str(j) ':最佳效益 = ' num2str(-BestCost_matrix(j))]);

    end

    %%结果
    out.pop = pop;
    out.bestsolu = best_solu; %将最佳个体换取到out.bestsolu
    out.bestcost = BestCost_matrix;

end
