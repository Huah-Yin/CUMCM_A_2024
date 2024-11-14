%% 求太阳方位角平均值
%% 初始化
clear
clc
close all
%% 参数赋值
fai=39.4;
D_ind=1;
for D=[-59,-28,0,31,61,92,122,153,184,214,245,275]
    for ST=[9,10.5,12,13.5,15]
    %% 计算太阳高度角与方位角
    de=asin(sin(2*pi*D/350)*sin(2*pi*23.45/360));
    omg=pi*(ST-12)/12;
    as=asin(cos(de)*cos(fai*2*pi/360)*cos(omg)+sin(de)*sin(fai*2*pi/360));          % 太阳高度角
    gs(D_ind)=acos((sin(de)-sin(as)*sin(fai*2*pi/360))/cos(as)*cos(fai*2*pi/360))*180/pi;  % 太阳方位角
    D_ind=D_ind+1;
    end
end
gs_bar=mean(gs)