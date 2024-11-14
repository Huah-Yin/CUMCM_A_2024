%% 初始化
clear
close all
clc

%% GA参数赋值
NP=50; % 种群数量
G=100;  % 最大遗传次数
L=20;   % 单维二进制数段长度
Dim=5;  % 维度
Max=[185,205,245,265,100];  % 定义域上限
Min=[165,185,225,245,60];    % 定义域下限
max_sta=inf;     % 最大稳定代数，稳定超过此代数将结束算法
max_err=30;     % 最大不符合约束条件次数

%% GA开始
f=zeros(NP,L*Dim);    % 初始种群
for i=1:NP
    f(i,:)=d2b([182,203,237,254,81.5-(i-1)/10],Min,Max,L);
end
max_fit_rec=zeros(1,G);

for k=1:G
    tic
    %% 开始迭代
    for i=1:NP
        U=f(i,:);
        x(i,:)=b2d(U,Min,Max,L);
        Fit(i)=aim(x(i,:));
    end
    [max_fit,max_fit_idx]=max(Fit);
    f_best=f(max_fit_idx,:);
    Fit1=(Fit-min(Fit))./(max_fit-min(Fit));

    %% 计算稳定代数
    for idx=k-1:-1:1
        if max_fit_rec(idx)~=max_fit
            break
        end
    end
    sta_num=k-idx;  % 稳定代数
    
    %% 判断是否稳定超max_sta代
    if sta_num>=max_sta
        break
    end

    %% 基于轮盘赌选择法的复制操作
    sum_Fit=sum(Fit1);
    fitvalue=Fit1/sum_Fit;
    fitvalue=(cumsum(fitvalue));
    ms=sort(rand(NP,1));
    fitidx=1;
    newidx=1;
    while newidx<=NP-1
        if ms(newidx) < fitvalue(fitidx)
            nf(newidx,:)=f(fitidx,:);   % 复制
            newidx=newidx+1;
        else
            fitidx=fitidx+1;
            if fitidx>NP
                break;
            end
        end
    end
    
    %% 基于概率的交叉操作
    for i=1:NP-2
        p=rand;
        if p<get_pc(sta_num)
            t=1;
            for trytime=1:max_err
                tnf=zeros(1,L*Dim);
                q=randi([0,1],1,L*Dim); % 多点交叉
                for j=1:L*Dim
                    tnf(j)=nf(i+q(j),j);
                end
                if judge(b2d(tnf,Min,Max,L))
                    nf(i,:)=tnf(:);
                    break
                end
            end
        end
    end
    
    %% 基于概率的变异操作
    for i=1:round((NP-1)*get_pm(sta_num))
        h=randi([1,NP-1],1,1);
        t=1;
        for trytime=1:max_err
            tnf=nf(h,:);
            for j=1:round(L*Dim*get_pm(sta_num))
                g=randi([1,L*Dim],1,1);
                tnf(g)=~tnf(g); % 变异
            end
            if judge(b2d(tnf,Min,Max,L))
                nf(i,:)=tnf(:);
                break
            end
        end
    end

    %% 下一代预备
    f=nf;
    f(NP,:)=f_best; % 保留精英
    max_fit_rec(k)=max_fit;
    usetime=toc;
    fprintf("完成计算第%d代，最大适应度%.2f，进度%.6f，用时%.2fs，预计还需%.3fmin\n",k,max_fit,k/G,usetime,usetime*(G-k)/60)
end

%% 收尾
for i=1:NP
    Fit(i)=aim(b2d(f(i,:),Min,Max,L));  % 计算每个样本适应度
end
[max_fit,max_fit_idx]=max(Fit);

%% 显示结果
fprintf("最大值:")
disp(max_fit)
fprintf("最大值点:")
disp(x(max_fit_idx,:))
fprintf('精度：x±')
disp((Max-Min)/(2^L-1));

function fit=aim(x)
    % 目标函数
    %% 设置炉温
    U_info=[ones(1,5)*x(1),x(2),x(3),ones(1,2)*x(4),ones(1,2)*25];
    % 炉前,1-4,5,6,7,8,9,10-11,炉后
    D_info=    [0.00014,ones(1,4)*0.00020,  0.00020,0.00020,0.00027,0.00020,0.00020,0.00015,0.00015];
    sigma_info=[0.012,  ones(1,4)*0.021,    0.024,  0.030,  0.034,  0.025,0.025,0.010,0.010];
    len=30.5;
    gap=5;
    global x0 y0
    x0=[0,25];
    y0=[25];
    D_ls=[D_info(1)];
    sigma_ls=[sigma_info(1)];
    for i=1:11
        x0(end+1)=x0(end)+len;
        x0(end+1)=x0(end)+gap;
        y0(end+1)=U_info(i);
        y0(end+1)=U_info(i);
        D_ls(end+1)=D_info(i+1);
        D_ls(end+1)=D_info(i+1);
        sigma_ls(end+1)=sigma_info(i+1);
        sigma_ls(end+1)=sigma_info(i+1);
    end
    x0(end)=435.5;
    y0(end+1)=25;
    D_ls(end+1)=D_ls(end);
    sigma_ls(end+1)=sigma_info(end);
    %% 计算
    v=x(5);
    d=0.15;     % 厚度(mm)
    j_step=0.005;   % 厚度步长(mm)
    k_step=0.005;     % 时间步长(s)
    vv=v/60;    % 传送带前进速度(cm/s)[不要修改参数]
    data=zeros(floor(435/(v/60)/k_step)+1,floor(d/j_step)+1);     % (k,j)->(时间s,空间cm);首行/首列为零时刻/位置
    [k_max,j_max]=size(data);
    data(1,:)=25;   % 初始条件
    key=ones(k_max,1)*25; % 中心温度记录
    for k=2:k_max
        t_now=k*k_step; % 当前时间(s)
        DD=k_step*get_D(t_now*vv,D_ls)/(j_step*j_step);
        sigma=get_sigma(t_now*vv,sigma_ls);
        for j=2:j_max-1
            data(k,j)=DD*(data(k-1,j+1)-2*data(k-1,j)+data(k-1,j-1))+data(k-1,j);
        end
        data(k,1)=-k_step*sigma*(data(k-1,1)-get_U(t_now*vv))+data(k-1,1);
        data(k,j_max)=-k_step*sigma*(data(k-1,j_max)-get_U(t_now*vv))+data(k-1,j_max);
        key(k)=data(k,d/2/j_step);
        if key(k)-key(k-1)<0    % 一下降就停止迭代
            break
        end
    end
    [~,k_idx_umax]=max(key);
    fit=0;
    for k=find(key>217)'
        if k>k_idx_umax
            break
        end
        fit=fit+((key(k)+key(k+1))/2-217)*k_step;
    end
    fit=-fit;
end

function res=judge(x)
    % 约束条件
    % 输入速度v(cm/s);输出1-符合约束,0-不符合约束
    %% 设置炉温
    U_info=[ones(1,5)*x(1),x(2),x(3),ones(1,2)*x(4),ones(1,2)*25];
    % 炉前,1-4,5,6,7,8,9,10-11,炉后
    D_info=    [0.00014,ones(1,4)*0.00020,  0.00020,0.00020,0.00027,0.00020,0.00020,0.00015,0.00015];
    sigma_info=[0.012,  ones(1,4)*0.021,    0.024,  0.030,  0.034,  0.025,0.025,0.010,0.010];
    len=30.5;
    gap=5;
    global x0 y0
    x0=[0,25];
    y0=[25];
    D_ls=[D_info(1)];
    sigma_ls=[sigma_info(1)];
    for i=1:11
        x0(end+1)=x0(end)+len;
        x0(end+1)=x0(end)+gap;
        y0(end+1)=U_info(i);
        y0(end+1)=U_info(i);
        D_ls(end+1)=D_info(i+1);
        D_ls(end+1)=D_info(i+1);
        sigma_ls(end+1)=sigma_info(i+1);
        sigma_ls(end+1)=sigma_info(i+1);
    end
    x0(end)=435.5;
    y0(end+1)=25;
    D_ls(end+1)=D_ls(end);
    sigma_ls(end+1)=sigma_info(end);
    %% 约束
    v=x(5);
    d=0.15;     % 厚度(mm)
    j_step=0.005;   % 厚度步长(mm)
    k_step=0.005;     % 时间步长(s)
    vv=v/60;    % 传送带前进速度(cm/s)[不要修改参数]
    data=zeros(floor(435/(v/60)/k_step)+1,floor(d/j_step)+1);     % (k,j)->(时间s,空间cm);首行/首列为零时刻/位置
    [k_max,j_max]=size(data);
    data(1,:)=25;   % 初始条件
    key=ones(k_max,1)*25; % 中心温度记录
    count1=0;   % 约束2计时器
    count2=0;   % 约束3计时器
    for k=2:k_max
        t_now=k*k_step; % 当前时间(s)
        DD=k_step*get_D(t_now*vv,D_ls)/(j_step*j_step);
        sigma=get_sigma(t_now*vv,sigma_ls);
        for j=2:j_max-1
            data(k,j)=DD*(data(k-1,j+1)-2*data(k-1,j)+data(k-1,j-1))+data(k-1,j);
        end
        data(k,1)=-k_step*sigma*(data(k-1,1)-get_U(t_now*vv))+data(k-1,1);
        data(k,j_max)=-k_step*sigma*(data(k-1,j_max)-get_U(t_now*vv))+data(k-1,j_max);
        key(k)=data(k,d/2/j_step);
        if key(k)-key(k-1)>3*k_step || key(k)-key(k-1)<-3*k_step    % 约束(1)
            res=0;
%             fprintf("v=%f,约束(1)跳出\n",v)
            return
        end
        if key(k)-key(k-1)>0 && key(k)>150 && key(k)<190
            count1=count1+1;
        elseif key(k)>217
            count2=count2+1;
        end
        if count1*k_step>120 || count2*k_step>90 || key(k)>250   % 约束(2)(3)(4)上界
%             fprintf("v=%f,约束(234)上界跳出\n",v)
            res=0;
            return
        end
        if key(k)-key(k-1)<0 && key(k)<217  % 下降，减少递推次数
            break
        end
    end
    if count1*k_step<60 || count2*k_step<40 || max(key)<240 % 约束(2)(3)(4)下界
%         fprintf("v=%f,约束(234)下界跳出\n",v)
        res=0;
        return
    end
%     fprintf("v=%f,符合条件\n",v)
    res=1;
end

function pc=get_pc(sta_num)
    % 交叉率，sta_num-稳定代数
    pc=min(0.2+sta_num*0.05,0.9);   % 防止pc超过0.9
end

function pm=get_pm(sta_num)
    % 变异率，sta_num-稳定代数
    pm=min(0.2+sta_num*0.05,0.9);   % 防止pm超过0.9
end

function d=b2d(str,Min,Max,L)
    %格雷解码并将二进制转为定义域内的十进制
    %b2d(str,Xx,Xs,L),str-二进制行向量，Min-下限，Max-上限，L单维度二进制串长度
    Dim=length(Min);
    n=zeros(Dim,L);
    for k=1:Dim
        n(k,1)=str((k-1)*L+1);  % 格雷编码解码
        for j=2:L
            n(k,j)=xor(n(k,j-1),str((k-1)*L+j));    % xor—异或运算
        end
    end
    for di=1:Dim
        m(di)=0;
        for j=1:L   % 二进制转十进制
            m(di)=m(di)+n(di,L-j+1)*(2^(j-1));
        end
    end
    d=Min+m.*(Max-Min)./(2^L-1);    % 十进制解码至定义域
end

function str = d2b (arr,Min,Max,L)
    %格雷编码
    %d2b(arr,Xx,Xs,L),arr-十进制数据数组（行向量），Min-下限，Max-上限，L单维度二进制串长度；输出str为二进制行向量
    Dim=length(Min);
    str=zeros(1,Dim*L);
    m=round((arr-Min).*(2^L-1)./(Max-Min)); % 分散至二进制定义域
    n=double(boolean(dec2bin(m')-'0'));   % 十进制转二进制
    [~,col]=size(n);
    while col<L     % 首位补0至长为L
        n=[zeros(length(arr),1),n];
        [~,col]=size(n);
    end   
    for k=1:Dim % 格雷编码
        str((k-1)*L+1)=n(k,1);  % 首位不变
        str((k-1)*L+2:k*L)=xor(n(k,1:end-1),n(k,2:end));    % 异或运算
    end
end

function U=get_U(x)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0 y0
    U=interp1(x0,y0,x);
end

function D=get_D(x,D_ls)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0
    D=interp1(x0,D_ls,x);
end

function sigma=get_sigma(x,sigma_ls)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0
    sigma=interp1(x0,sigma_ls,x);
end