function E_year = Q3_func_dim(arg)
%问题三目标函数   
arg=xlsread("Q3_result(1).xlsx");
ALL_NUM=3500;
%% 计算吸收塔坐标
xi=arg(5)*[cos(-77.3*pi/180),sin(-77.3*pi/180)];
%% 布阵，生成定日镜坐标
aff_num=-1;  % 当前行定日镜数目
lin_ind=1;  % 行号
distan=100; % 当前离吸收塔距离
coo_ls=[];  % 坐标列表
num_coo=0;  % 当前坐标数
while aff_num ~= 0 && num_coo<ALL_NUM  % 循环行
    %% 实现do-while循环
    aff_num=0;
    %% 生成相对于吸收塔的极坐标
    the=0;  % 下一角
    now_the=0;  % 当前转角
    while now_the+the<2*pi-the
        if num_coo>=ALL_NUM
            break
        end
        %% 转直角坐标
        now_the=now_the+the;
        coo_now=distan.*[cos(now_the),sin(now_the)];
        %% 改变参考系
        coo_now=xi-coo_now;
        %% 判断是否在场地内
        the_app=0;
        if dot(coo_now,coo_now)<=122500
            num_coo=num_coo+1;
            coo_ls(num_coo,1)=coo_now(1);
            coo_ls(num_coo,2)=coo_now(2);
            %% 计算安装高度
            d=abs(-coo_now(1)+0.48*coo_now(2)+388.75)/sqrt(1+0.48^2);
            coo_ls(num_coo,3)=arg(2)+arg(1)*d;
            aff_num=aff_num+1;
            if lin_ind<=20
                the_app=arg(3);
            else
                the_app=arg(4);
            end
        end
        the=asin((the_app+(arg(5+num_coo)+arg(6+num_coo))/2)/(2*distan));
    end
    if lin_ind<=20
        distan=distan+15;
    else
        distan=distan+0.5*lin_ind+5;
    end
    lin_ind=lin_ind+1;
end
len=length(coo_ls);
coo_tran=coo_ls;
coo_tran(:,1)=coo_tran(:,1)-xi(1);
coo_tran(:,2)=coo_tran(:,2)-xi(2);
coo_tran(:,5)=coo_tran(:,3);
coo_tran(:,3)=arg(6:5+len);
coo_tran(:,4)=arg(ALL_NUM+6:ALL_NUM+5+len);
E_year=Q3_func(coo_tran);
%% 出图
close all
figure
hold on
scatter3(coo_ls(:,1), coo_ls(:,2),coo_ls(:,3))
rectangle('Position', [0-350,0-350,2*350,2*350], 'Curvature', [1 1],'EdgeColor', 'r');
plot(xi(1),xi(2),'b*')
axis equal
grid on
fprintf('计算完成，镜面数目%d，',len)
end