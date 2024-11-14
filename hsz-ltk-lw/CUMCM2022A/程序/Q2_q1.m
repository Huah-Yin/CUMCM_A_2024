% 问题2第1小问
%% 初始化
clear
clc

%% 主程序
omg=1.7152;
T=100*2*pi/omg;
global C;
e=0.001;
left_barrier=0;
right_barrier=100000;
dist=right_barrier-left_barrier;
while dist>e
    %% 左试验点
    left_try=right_barrier-0.618*dist;
    C(1)=left_try;
    C(2)=0;
    [~,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
    left_ans=sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T;
    %% 右试验点
    right_try=left_barrier+0.618*dist;
    C(1)=right_try;
    C(2)=0;
    [~,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
    right_ans=sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T;
    %% 边界更新
    if left_ans>right_ans
        right_barrier=right_try;
    elseif left_ans<right_ans
        left_barrier=left_try;
    else
        right_barrier=right_try;
        left_barrier=left_try;
    end
    dist=right_barrier-left_barrier;
    fprintf("当前长度=%.6f\n",dist)
end
if left_ans > right_ans
    max_num=left_ans;
    max_idx=left_try;
else
    max_num=right_ans;
    max_idx=right_try;
end
fprintf("============结果============\n")
disp("最大值"+max_num)
disp("最大值点"+max_idx)

% 计算结果：
% 最大值229.4678
% 最大值点37265.4646
