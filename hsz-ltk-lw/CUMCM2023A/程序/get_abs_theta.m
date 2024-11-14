%% 初始化
clear
clc
close all

%% 读数据
data=xlsread("2023A-定日镜问题1结论.xlsx");
len=length(data);

%% 坐标转化
[theta,rho]=cart2pol(data(:,1),data(:,2));

%% 圈数确定
r_ls=[107.8, 121.3, 134.8, 148.3, 161.8, 175.3, 188.7, 202.2, 215.7, 229.2, 242.7, 256.2, ...
      269.7, 283.1, 296.6, 310.1, 323.6, 337.1];
loop_max_theta=zeros(1,length(r_ls));
loop_max_value=zeros(1,length(r_ls));
for idx=1:len
    loop=find(r_ls-rho(idx)<=2);
    %% 最大值更新
    if data(idx,3)>loop_max_value(loop)
        loop_max_value(loop)=data(idx,3);
        loop_max_theta(loop)=theta(idx);
    end
end

%% 统计
disp(loop_max_theta)
fprintf("均值=%f rad(%f°)",mean(loop_max_theta),mean(loop_max_theta)*180/pi)