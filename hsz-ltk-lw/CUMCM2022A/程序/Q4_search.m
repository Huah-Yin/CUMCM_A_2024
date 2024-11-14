% 问题4
%% 初始化
clear
clc

%% 计算恒定的质心与转动惯量，存至全局变量之中
m1=4866;
R0=1;
h1=3;
h2=0.8;
global I1 rz1;
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(pi*R0*h1^2+2*pi*R0/h2*(1+R0/h2)*(-(h2^3)/6)+h1*pi*R0^2)/m1;
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;

%% main
barrier1=[30000,70000];
barrier2=[30000,100000];
area=get_area(barrier1,barrier2);
e=0.0001;
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
[beta,gamma]=meshgrid(barrier1,barrier2);
res=zeros(1,4);
for i=1:4
    res(i)=aim([beta(i),gamma(i)]);
end
[max_num,idx]=max(res);
max_idx=[beta(idx),gamma(idx)];
disp("最大值"+max_num)
disp(["beta=","gamma="]+max_idx)

%% 目标函数
function result=aim(x)
    omg=1.9806;
    T=100*2*pi/omg;
    global C
    C(1)=x(1);
    C(2)=x(2);
    [~,n]=ode45('Q4_func1',0:0.01:100+T,[0,0,0,0,0,0,0,0]);
    result=sum((C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^2 ...
        +C(2)*(abs(n(10000:end,8)-n(10000:end,7))).^2).*0.01)/T;
end

function area=get_area(barrier1,barrier2)
    area=(barrier1(2)-barrier1(1))*(barrier2(2)-barrier2(1));
end

function [max_num,max_idx]=get_best1(barrier1,val2)
    % 优化因素1
    global C;
    left_barrier=barrier1(1);
    right_barrier=barrier1(2);
    dist=right_barrier-left_barrier;
    e=0.001;
    while dist>e
        %% 左试验点
        left_try=right_barrier-0.618*dist;
        left_ans=aim([left_try,val2]);
        %% 右试验点
        right_try=left_barrier+0.618*dist;
        right_ans=aim([right_try,val2]);
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
    global C;
    left_barrier=barrier2(1);
    right_barrier=barrier2(2);
    dist=right_barrier-left_barrier;
    e=0.001;
    while dist>e
        %% 左试验点
        left_try=right_barrier-0.618*dist;
        left_ans=aim([val1,left_try]);
        %% 右试验点
        right_try=left_barrier+0.618*dist;
        right_ans=aim([val1,right_try]);
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