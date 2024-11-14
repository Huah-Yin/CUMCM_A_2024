% 问题2
%% 初始化
clear
clc

%%主程序
alpha_ls=37131:0.001:37133;
k_ls=0;
P=zeros(length(alpha_ls),length(k_ls));
omg=1.7152;
T=100*2*pi/omg;
global C;
alpha_idx=0;
for alpha=alpha_ls
    tic
    alpha_idx=alpha_idx+1;
    k_idx=0;
    for k=k_ls
        k_idx=k_idx+1;
        C(1)=alpha;
        C(2)=k;
        [t,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
        su=sum(C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^(C(2)+2).*0.01)/T;
        P(alpha_idx,k_idx)=su;
    end
    time=toc;
    fprintf("完成计算alpha=%f,进度%f,用时%fs,预计还需%.2fmin\n",alpha,alpha_idx/length(alpha_ls),time,(length(alpha_ls)-alpha_idx)*time/60)
end

%% 绘图
plot(alpha_ls,P',"o-")
% [x,y]=meshgrid(alpha_ls,k_ls);
% mesh(x,y,P')
% colorbar
xlabel("$$\alpha$$",Interpreter="latex")
ylabel("$$\overline{P_1}$$",Interpreter="latex")
