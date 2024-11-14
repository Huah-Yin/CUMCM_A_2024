function [y1, y2] = UnformCrossover(x1, x2, MyGama)
    %生成一个与变量 x1 相同尺寸的向量,在0到1之间均匀分布
    %alpha = rand(size(x1));
    alpha = unifrnd(-MyGama, 1 + MyGama, size(x1));
    %%点积
    y1 = alpha .* x1 + (1 - alpha) .* x2;
    y2 = alpha .* x2 + (1 - alpha) .* x1;
     % 在交叉后的子代中应用变异
    y1=mutation(y1,0.1,0.1);
    y2=mutation(y2,0.1,0.1);
end
