% 问题三求解程序
%% 初始化
clear
clc
close all

%% 参数赋值
v=100;  % 龙头速度(cm/s)
return_radius=4.5*100;  % 掉头区间半径(cm)
%% 二分查找
b_left=30;
b_right=55;
figure
hold on
while b_right-b_left>1e-4
%% 进行尝试
    b=(b_left+b_right)/2   % 螺距(cm)
    s=get_s(b,16*2*pi);   % 自然坐标系龙头位置初值(cm)
    %% 时间递推
    t=1;    % t-时间(s)
    bk=0;   % 跳出flag
    while true
        t=t+1;
        xdata=zeros(223+1,1);  % x
        ydata=zeros(223+1,1);  % y
        vdata=ones(223+1,1);   % 速度(m/s)
        theta=zeros(1,223+1);
        r=zeros(1,233+1);
        %% 求解龙头位置
        if t~=0
            s=s-0.01*v;  % 龙头位置更新
        end
        theta(1)=inv_s(b,s);    % 龙头极角(rad)
        r(1)=get_curve(b,theta(1));   % 龙头极径(cm)
        [xdata(1),ydata(1)]=get_xy(r(1),theta(1));
        %% 递推龙身位置
        for idx=1:223
            if idx==1
                R=341-27.5*2;
            else
                R=220-27.5*2;
            end
            theta(idx+1)=get_body_theta(b,theta(idx),R);
            [xdata(idx+1),ydata(idx+1)]=get_xy(get_curve(b,theta(idx+1)),theta(idx+1));
            %% 求速度
            x0=xdata(idx);
            y0=ydata(idx);
            x1=xdata(idx+1);
            y1=ydata(idx+1);
            k=(y1-y0)/(x1-x0);  % 割线斜率
            k0=(sin(theta(idx))+theta(idx)*cos(theta(idx)))/(cos(theta(idx))-theta(idx)*sin(theta(idx)));  % 切线0斜率
            k1=(sin(theta(idx+1))+theta(idx+1)*cos(theta(idx+1)))/(cos(theta(idx+1))-theta(idx+1)*sin(theta(idx+1)));  % 切线1斜率
            the0=abs(abs(atan(k0))-abs(atan(k)));
            the1=abs(abs(atan(k))-abs(atan(k1)));
            vdata(idx+1)=vdata(idx)*cos(the0)/cos(the1);
        end
        %% 判断是否相撞
        for point_wood_idx=1:4
            if point_wood_idx==1
                chang=314/2;
            else
                chang=220/2;
            end
            kuan=30/2;
            point=[chang,chang,-chang,-chang;
                   kuan,-kuan,kuan,-kuan];
            zhuan=-atan((ydata(point_wood_idx+1)-ydata(point_wood_idx))/(xdata(point_wood_idx+1)-xdata(point_wood_idx)));
            T=[cos(zhuan),sin(zhuan);
               -sin(zhuan),cos(zhuan)];
            point_aft_T=T*point;
            move=[(xdata(point_wood_idx+1)+xdata(point_wood_idx))/2;
                (ydata(point_wood_idx+1)+ydata(point_wood_idx))/2];
            point_aft_move=point_aft_T+move; % 2x4 array
            for line_wood_idx=point_wood_idx+2:100
                x1=xdata(line_wood_idx);
                x2=xdata(line_wood_idx+1);
                y1=ydata(line_wood_idx);
                y2=ydata(line_wood_idx+1);
                % 平行线系数ABCp
                Ap=(y1-y2);
                Bp=-(x1-x2);
                Cp=(x1-x2)*y2-(y1-y2)*x2;
                dis_p=abs(Ap*point_aft_move(1,:)+Bp*point_aft_move(2,:)+Cp)/(Ap^2+Bp^2).^0.5;    % 到平行线距离
                if any(dis_p<15)
                    line_mid=[(xdata(line_wood_idx+1)+xdata(line_wood_idx))/2;
                              (ydata(line_wood_idx+1)+ydata(line_wood_idx))/2];
                    k=(ydata(line_wood_idx+1)-ydata(line_wood_idx))/(xdata(line_wood_idx+1)-xdata(line_wood_idx));
                    % 中垂线系数ABCc
                    Ac=1;
                    Bc=k;
                    Cc=-line_mid(1)-k*line_mid(2);
                    dis_c=abs(Ac*point_aft_move(1,:)+Bc*point_aft_move(2,:)+Cc)/(Ac^2+Bc^2).^0.5;    % 到中垂线距离
                    if any(dis_c<110)
                        can=0;
                        bk=1;
                        break
                    end
                end
            end
            if bk==1
                break
            end
        end
        if bk==1
            break
        end
        %% 判断是否进入调头区域
        if (xdata(1)^2+ydata(1)^2)^0.5<=return_radius
            can=1;
            break
        end
    end
    %% 二分更新
    if can==1
        plot(b,1,"b+",'MarkerSize',12,'LineWidth',2);
        % 可掉头,b太大
        b_right=b;
    else
        plot(b,0,"rx",'MarkerSize',12,'LineWidth',2);
        % 不可掉头,b太小
        b_left=b;
    end
end

%% 画图
% [x_std,y_std]=get_xy(get_curve(b,0:0.1:44*pi),0:0.1:44*pi);
% figure
% hold on
% plot(x_std,y_std,"-",'Color',[0.01,0.01,0.01])
% plot(xdata(:),ydata(:),"-",'Color',[0.1,0.3,0.9])
% plot(xdata(:),ydata(:),"*",'Color',[0.9,0.1,0.1])

%% 螺线方程
function r=get_curve(b,theta)
    % in:b-螺距(cm);theta-极角(rad)
    % out:r-极径(cm)
    r=b.*theta./(2.*pi);
end

%% 螺线方程转平面直角
function [x,y]=get_xy(r,theta)
    x=r.*cos(theta);
    y=r.*sin(theta);
end

%% 螺线弧长方程
function s=get_s(b,theta)
    % in:b-螺距(cm);theta-极角(rad)
    % out:s-弧长(cm)
    s=b./(4*pi).*(theta.*(1+theta.^2).^0.5+log(theta+(1+theta.^2).^0.5));
end

%% 反-螺线弧长方程,获取龙头极角
function theta=inv_s(b,s)
    % in:b-螺距(cm);s-弧长(cm)
    % out:theta-极角(rad)
    %% 二分求解
    left=0;
    right=17*2*pi;
    err=right-left;
    while err>1e-8
        mid=(right+left)/2;
        mid_s=get_s(b,mid);
        if mid_s>s
            right=mid;
        else
            left=mid;
        end
        err=right-left;
    end
    theta=(right+left)/2;
end

%% 递推龙身
function t=get_body_theta(b,theta_bef,R)
    [x0,y0]=get_xy(get_curve(b,theta_bef),theta_bef);
    p0=[x0,y0];
    ktemp=(sin(theta_bef)+theta_bef*cos(theta_bef))/(cos(theta_bef)-theta_bef*sin(theta_bef));  % 切线斜率
    k=-1/ktemp;
    gamma=atan(ktemp);
    rho=get_rho(b,theta_bef);
    point_E1=[x0+rho*sin(gamma),y0-rho*cos(gamma)];
    point_E2=[x0-rho*sin(gamma),y0+rho*cos(gamma)];
    if norm(point_E1)<norm(point_E2)
        point_E=point_E1;
    else
        point_E=point_E2;
    end
    t=get_cross_point(theta_bef,point_E,rho,p0,R,R,x0,y0);
end

%% 求曲率半径
function rho=get_rho(b,theta)
    rho=b/(2*pi)*(theta^2+1)^(3/2)/(theta^2+2);
end

%% 求两圆交点
function theta=get_cross_point(theta_bef,p1,r1,p2,r2,R,x0,y0)
    % p1,p2-输入点横向量;r1,r2-圆半径
    d=norm(p1-p2);
    a=(r1^2-r2^2+d^2)/(2*d);
    h=sqrt(r1^2-a^2);
    p=p1+a*(p2-p1)/d;
    x1=p(1)+h*(p2(2)-p1(2))/d;
    y1=p(2)-h*(p2(1)-p1(1))/d;
    x2=p(1)-h*(p2(2)-p1(2))/d;
    y2=p(2)+h*(p2(1)-p1(1))/d;
    if y1*y2<0
        if x0>0
            if y1>0
                y=y1;
                x=x1;
            else
                x=x2;
                y=y2;
            end
        else
            if y1<0
                y=y1;
                x=x1;
            else
                x=x2;
                y=y2;
            end
        end
    elseif x1*x2<0
        if y0>0
            if x1<0
                x=x1;
                y=y1;
            else
                x=x2;
                y=y2;
            end
        else
            if x1>0
                x=x1;
                y=y1;
            else
                x=x2;
                y=y2;
            end
        end
    else
        if y0<0
            x=max(x1,x2);
        elseif y0>0
            x=min(x1,x2);
        end
        if x0>0
            y=max(y1,y2);
        else
            y=min(y1,y2);
        end
    end
    R1=norm(p2);
    R2=norm([x,y]);
    dtheta=acos((R1^2+R2^2-R^2)/(2*R1*R2));
    theta=(dtheta)+theta_bef;
end