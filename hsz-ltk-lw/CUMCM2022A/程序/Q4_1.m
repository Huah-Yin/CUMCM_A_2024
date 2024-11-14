%% 初始化
clear
clc
close all

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

%% 主函数
beta_ls=0;%0:10000:100000;
gamma_ls=50000;%0:10000:100000;
P=zeros(length(beta_ls),length(gamma_ls));
omg=1.9806;
T=100*2*pi/omg;
beta_idx=0;
for beta=beta_ls
    tic
    beta_idx=beta_idx+1;
    gamma_idx=0;
    for gamma=gamma_ls
        gamma_idx=gamma_idx+1;
        global C
        C(1)=beta;
        C(2)=gamma;
        [~,n]=ode45('Q4_func1',0:0.01:100+T,[0,0,0,0,0,0,0,0]);
        p1=sum((C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^2).*0.01)/T;
        p2=sum((C(2)*(abs(n(10000:end,8)-n(10000:end,7))).^2).*0.01)/T;
        su=p1+p2;
        P(beta_idx,gamma_idx)=su;
    end
    time=toc;
    fprintf("完成计算beta=%f,进度%f,用时%fs,预计还需%.2fmin\n",beta,beta_idx/length(beta_ls),time,(length(beta_ls)-beta_idx)*time/60)
end

%% 绘图
[x,y]=meshgrid(beta_ls,gamma_ls);
mesh(x,y,P')
colorbar
xlabel("$$\beta$$",Interpreter="latex")
ylabel("$$\gamma$$",Interpreter="latex")
