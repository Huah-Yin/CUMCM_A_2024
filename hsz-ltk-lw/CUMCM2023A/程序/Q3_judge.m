function result = Q3_judge( arg )
%判断x是否符合约束条件，若符合，则返回1，若因宽高度不符，则自动更正，返回2与更正后十进制数组
ALL_NUM=3500;
% arg=rand(1,5+2*ALL_NUM);
% arg(1:5)=[0.003,2,10,15,100];
flag=1;
if 700*arg(1)+arg(2)>6
    result={0};
    fprintf('坡度错误')
    return
end
%% 计算吸收塔坐标
xi=arg(5)*[cos(0.4503),-sin(0.4503)];
%% 布阵，生成定日镜坐标
aff_num=-1;  % 当前行定日镜数目
lin_ind=1;  % 行号
distan=100; % 当前离吸收塔距离
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
            %% 计算安装高度
            d=abs(-coo_now(1)+0.48*coo_now(2)+388.75)/sqrt(1+0.48^2);
            if arg(2)+arg(1)*d <= arg(5+ALL_NUM+num_coo)/2
                result={0};
                return
            end
            %% 宽高度试验
            if arg(5+ALL_NUM+num_coo)>arg(5+num_coo)
                temp=arg(5+num_coo);
                arg(5+num_coo)=arg(5+ALL_NUM+num_coo);
                arg(5+ALL_NUM+num_coo)=temp;
                flag=2;
            end
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
if flag==2
    result={2,arg};
else
    result={flag};
end
end