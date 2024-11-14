function y = mutation(x, mu, MySigma) %mu为变异概率
    my_flag = (rand(size(x)) < mu); %my_flag为矩阵
    y = x;
    r = randn(size(x));
    %normal distribution 正态分布 N(0.sigma)
    y(my_flag) = x(my_flag) + MySigma * r(my_flag); %逻辑索引，只选择my_flag=1即为真的元素，将其进行取反操作
    %具有高斯步长的突变
end
