function i = WheelSelection(p) %p是每一个种群被选择的概率

    number = rand * sum(p); %防止p的求和大于1，从而使得生成的rand函数结果无法被随机选择
    c = cumsum(p); %将圆盘抽象为一个线性的直线，累计求其概率
    i = find(number <= c, 1, 'first'); %寻找从第一个开始符合要求的父本

end
