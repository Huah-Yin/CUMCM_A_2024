% 问题2第2小问
%% 初始化
clear
clc

%% main
barrier1=[60000,100000];
barrier2=[0,0.5];
area=get_area(barrier1,barrier2);
e=0.00001;
B=[(barrier1(2)-barrier1(1))/2,0,0];    % 优化因素2得到的最优点x,y,z坐标
flag=0; % 首次不划因素2
while area>e
    %% 优化因素2
    A=[B(1),0,0];
    [A(3),A(2)]=get_best2(B(1),barrier2);
    %% 因素2边界更新
    if flag
        if A(2)>B(2)
            barrier2(1)=B(2);
        else
            barrier2(2)=B(2);
        end
        area=get_area(barrier1,barrier2);
    else
        flag=1;
    end
    fprintf("当前面积=%.6f,区域[%.4f,%.4f]×[%.4f,%.4f],最大值=%.6f\n",area,barrier1(1),barrier1(2),barrier2(1),barrier2(2),A(3))
    %% 优化因素1
    B=[0,A(2),0];
    [B(3),B(1)]=get_best1(barrier1,A(2));
    %% 因素1边界更新
    if A(1)>B(1)
        barrier1(2)=A(1);
    else
        barrier1(1)=A(1);
    end
    area=get_area(barrier1,barrier2);
    fprintf("当前面积=%.6f,区域[%.4f,%.4f]×[%.4f,%.4f],最大值=%.6f\n",area,barrier1(1),barrier1(2),barrier2(1),barrier2(2),B(3))
end
[a,k]=meshgrid(barrier1,barrier2);
res=zeros(1,4);
for i=1:4
    res(i)=aim([a(i),k(i)]);
end
[max_num,idx]=max(res);
max_idx=[a(idx),k(idx)];
disp("最大值"+max_num)
disp(["alpha=","k="]+max_idx)

%% 目标函数
function result=aim(x)
    omg=1.7152;
    T=100*2*pi/omg;
    global C
    C(1)=x(1);
    C(2)=x(2);
    [~,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
    result=(sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T);
end

function area=get_area(barrier1,barrier2)
    area=(barrier1(2)-barrier1(1))*(barrier2(2)-barrier2(1));
end

function [max_num,max_idx]=get_best1(barrier1,val2)
    % 优化因素1
    omg=1.7152;
    T=100*2*pi/omg;
    global C;
    left_barrier=barrier1(1);
    right_barrier=barrier1(2);
    dist=right_barrier-left_barrier;
    e=0.001;
    while dist>e
        %% 左试验点
        left_try=right_barrier-0.618*dist;
        C(1)=left_try;
        C(2)=val2;
        [~,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
        left_ans=sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T;
        %% 右试验点
        right_try=left_barrier+0.618*dist;
        C(1)=right_try;
        C(2)=val2;
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
%         fprintf("当前长度=%.6f\n",dist)
    end
    if left_ans > right_ans
        max_num=left_ans;
        max_idx=left_try;
    else
        max_num=right_ans;
        max_idx=right_try;
    end
end

function [max_num,max_idx]=get_best2(val1,barrier2)
    % 优化因素2
    omg=1.7152;
    T=100*2*pi/omg;
    global C;
    left_barrier=barrier2(1);
    right_barrier=barrier2(2);
    dist=right_barrier-left_barrier;
    e=0.000001;
    while dist>e
        %% 左试验点
        left_try=right_barrier-0.618*dist;
        C(1)=val1;
        C(2)=left_try;
        [~,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
        left_ans=sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T;
        %% 右试验点
        right_try=left_barrier+0.618*dist;
        C(1)=val1;
        C(2)=right_try;
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
%         fprintf("当前长度=%.6f\n",dist)
    end
    if left_ans > right_ans
        max_num=left_ans;
        max_idx=left_try;
    else
        max_num=right_ans;
        max_idx=right_try;
    end
end