function E_year = Q2_func_dim(arg)
% arg=[348.9736,8.9150,8.9443,2.0978,3.6481,3.4428];
%问题二目标函数   
%% 计算吸收塔坐标
xi=arg(1)*[cos(-77.3*pi/180),sin(-77.3*pi/180)];
%% 布阵，生成定日镜坐标
aff_num=-1;  % 当前行定日镜数目
lin_ind=1;  % 行号
distan=100; % 当前离吸收塔距离
coo_ls=[];  % 坐标列表
num_coo=0;  % 当前坐标数
while aff_num ~= 0  % 循环行
    %% 实现do-while循环
    aff_num=0;
    %% 生成相对于吸收塔的极坐标
    if lin_ind<=20
        lin_num=floor(distan*2*pi/arg(2));
    else
        lin_num=floor(distan*2*pi/arg(3));
    end
    the=2*pi/lin_num;
    now_the=0;
    for the_ind=1:lin_num
        %% 转直角坐标
        now_the=now_the+the;
        coo_now=distan.*[cos(now_the),sin(now_the)];
        %% 改变参考系
        coo_now=xi-coo_now;
        %% 判断是否在场地内
        if dot(coo_now,coo_now)<=122500
            num_coo=num_coo+1;
            coo_ls(num_coo,1)=coo_now(1);
            coo_ls(num_coo,2)=coo_now(2);
            aff_num=aff_num+1;
        end
    end
    if lin_ind<=20
        distan=distan+15;
    else
        distan=distan+0.5*lin_ind+5;
    end
    lin_ind=lin_ind+1;
end
coo_tran=coo_ls;
coo_tran(:,1)=coo_tran(:,1)-xi(1);
coo_tran(:,2)=coo_tran(:,2)-xi(2);
E_year=Q2_func(coo_tran,arg(5),arg(6),arg(4));
%% 出图
% close all
% figure
% hold on
% scatter(coo_ls(:,1), coo_ls(:,2))
% rectangle('Position', [0-350,0-350,2*350,2*350], 'Curvature', [1 1],'EdgeColor', 'r');
% plot(xi(1),xi(2),'b*')
% axis equal
% grid on
fprintf('计算完成，镜面数目%d',length(coo_ls))
end