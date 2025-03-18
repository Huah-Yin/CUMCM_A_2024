% 问题五求解程序
%% 初始化
clear
clc
close all

%% 参数赋值
% v=100;  % 龙头速度(cm/s)
b=170;   % 螺距(cm)
s=get_s(b,8*2*pi);   % 自然坐标系龙头位置初值(cm)
head_s=-1;   % 龙头在掉头曲线自然坐标系下坐标(cm)
dt=0.1; % 采样时间(s)
%% 二分查找
v_left=15;
v_right=20;
while v_right-v_left>1e-4
%% 进行尝试
%     hold on
    return_radius=450;  % 掉头空间半径(cm)
    v=(v_left+v_right)/2
%     v=18.74
    boundary_theta=2*pi*return_radius/b;    % 进入掉头空间的临界极角(rad)
    rotate=pi/2-mod(boundary_theta,2*pi);
    t_in=-1;    % 进入掉头区域时间索引
    %% 时间递推
    t=-1;    % t-时间(s)
    bk=0;   % 跳出flag
    rec_x=zeros(223+1,3000);    % 每秒结果保存(m)
    rec_y=zeros(223+1,3000);    % 每秒结果保存(m)
    rec_v=zeros(223+1,30000);    % 每秒结果保存(m/s)
    state=zeros(233+1,1);   % 节点状态:0-盘入,1-掉头圆1,2-掉头圆2,3-盘出
    while true
        t=t+1;
        if t_in~=-1 && t>t_in+100/dt
            break
        end
        xdata=zeros(223+1,1);  % x(cm)
        ydata=zeros(223+1,1);  % y(cm)
        vdata=ones(223+1,1)*v/100;   % 速度(m/s)
        theta=zeros(1,223+1);   % 极角
        r=zeros(1,233+1);
        %% 求解龙头位置
        if t~=0 && (head_s==-1 || head_s==-2)
            s=s-dt*v;  % 龙头位置更新
        end
        theta(1)=inv_s(b,s);    % 龙头极角(rad)
        if theta(1)>boundary_theta && head_s == -1
            % 掉头前
            r(1)=get_curve(b,theta(1));   % 龙头极径(cm)
            [xdata(1),ydata(1)]=get_xy(r(1),theta(1));
        elseif (abs(theta(1))<boundary_theta || head_s > 0) && head_s<(calculate_return_length(return_radius)-dt*v) && head_s ~= -2
            % 掉头时
            state(1)=1;
            if head_s == -1
                % 初次进入
                head_s=(s+dt*v)-get_s(b,boundary_theta);  % 进入路程(cm)
            else
                head_s=head_s+dt*v;
            end
            if head_s>pi*return_radius*2/3
                state(1)=2;
            end
            if t_in==-1
                t_in=t;
            end
            [xdata(1),ydata(1)]=from_s_get_return_xy(head_s,return_radius,rotate);
        else
            % 掉头后
            if t_in==-1
                t_in=t;
            end
            state(1)=3;
            if s>0 % 初次离开
                s=-get_s(b,boundary_theta)-(head_s-calculate_return_length(return_radius)+dt*v);
                head_s=-2;
            end
            theta(1)=inv_s(b,abs(s));
            r(1)=abs(get_curve(b,theta(1)));   % 龙头极径(cm)
            [xdata(1),ydata(1)]=get_xy(r(1),theta(1));
            xdata(1)=-xdata(1);
            ydata(1)=-ydata(1);
        end
%         plot(xdata(1),ydata(1),"*r")
        %% 递推龙身位置
        for idx=1:223
            if idx==1
                L=341-27.5*2;
            else
                L=220-27.5*2;
            end
            T=[cos(rotate),sin(rotate);
               -sin(rotate),cos(rotate)];
            res=T*[0,               0;
                   return_radius/3, -return_radius*2/3];
            xo1=res(1,1);
            yo1=res(2,1);
            xo2=res(1,2);
            yo2=res(2,2);
            ro1=return_radius*2/3;
            ro2=return_radius/3;
            [xin,yin]=get_xy(return_radius,boundary_theta);
            if state(idx)==0
                % 前一节点在盘入
                theta(idx+1)=get_body_theta(b,theta(idx),L);
                [xdata(idx+1),ydata(idx+1)]=get_xy(get_curve(b,theta(idx+1)),theta(idx+1));
            elseif state(idx)==1
                % 前一节点在掉头圆1
                temp1=acos((xdata(idx)^2+ydata(idx)^2+return_radius^2-L^2)/(2*return_radius*(xdata(idx)^2+ydata(idx)^2)^0.5));
                temp2=acos((return_radius^2+xdata(idx)^2+ydata(idx)^2-(xdata(idx)-xin)^2-(ydata(idx)-yin)^2)/(2*return_radius*(xdata(idx)^2+ydata(idx)^2)^0.5));
                if temp1>temp2
                    % 此节点未进入掉头圆1
                    state(idx+1)=0;
                    theta(idx+1)=boundary_theta+temp1-temp2;
                    [xdata(idx+1),ydata(idx+1)]=get_xy(get_curve(b,theta(idx+1)),theta(idx+1));
                else
                    % 此节点已进入掉头圆1
                    state(idx+1)=1;
                    temp3=2*asin(L/(2*ro1));
                    Trans=[cos(temp3),-sin(temp3);
                           sin(temp3),cos(temp3)];  % 逆时针旋转
                    tem=Trans*[xdata(idx)-xo1;ydata(idx)-yo1];
                    xdata(idx+1)=tem(1)+xo1;
                    ydata(idx+1)=tem(2)+yo1;
                end
            elseif state(idx)==2
                % 前一节点在掉头圆2
                T=[cos(rotate),-sin(rotate);
                   sin(rotate),cos(rotate)];
                res=T*[xdata(idx);ydata(idx)];
                y_bef=res(2);
                if y_bef>=-L^2/(2*ro2)-return_radius/3
                    % 此节点未进入掉头圆2
                    state(idx+1)=1;
                    temp=acos(((xdata(idx)-xo1)^2+(ydata(idx)-yo1)^2+return_radius^2-L^2)/(2*return_radius*((xdata(idx)-xo1)^2+(ydata(idx)-yo1)^2)^0.5));
                    Trans=[cos(temp),-sin(temp);
                           sin(temp),cos(temp)];  % 逆时针旋转
                    res=Trans*[xdata(idx)-xo1;ydata(idx)-yo1];% +[xo1;yo1];
                    det_r=res;%-[xo1;yo1];
                    det_r=det_r/norm(det_r)*ro1;
                    xdata(idx+1)=det_r(1)+xo1;
                    ydata(idx+1)=det_r(2)+yo1;
                else
                    % 此节点已进入掉头圆2
                    state(idx+1)=2;
                    temp3=2*asin(L/(2*ro2));
                    Trans=[cos(temp3),sin(temp3);
                           -sin(temp3),cos(temp3)];  % 顺时针旋转
                    tem=Trans*[xdata(idx)-xo2;ydata(idx)-yo2];
                    xdata(idx+1)=tem(1)+xo2;
                    ydata(idx+1)=tem(2)+yo2;
                end
            else
                % 前一节点在盘出
                temp1=acos(((xdata(idx)-xo2)^2+(ydata(idx)-yo2)^2+ro2^2-L^2)/(2*ro2*((xdata(idx)-xo2)^2+(ydata(idx)-yo2)^2)^0.5));
                temp2=acos(((xdata(idx)-xo2)^2+(ydata(idx)-yo2)^2+ro2^2-(xdata(idx)+xin)^2-(ydata(idx)+yin)^2)/(2*ro2*((xdata(idx)-xo2)^2+(ydata(idx)-yo2)^2)^0.5));
                if temp1>temp2
                    % 此节点未进入盘出曲线
                    state(idx+1)=2;
                    xout=-xin;
                    yout=-yin;
                    temp3=temp1-temp2;
                    Trans=[cos(temp3),sin(temp3);
                           -sin(temp3),cos(temp3)];  % 顺时针旋转
                    tem=Trans*[xout-xo2;yout-yo2];
                    xdata(idx+1)=tem(1)+xo2;
                    ydata(idx+1)=tem(2)+yo2;
                else
                    % 进入盘出曲线
                    state(idx+1)=3;
                    theta(idx+1)=get_body_theta_out(b,theta(idx),L);
                    [xtemp,ytemp]=get_xy(get_curve(b,theta(idx+1)),theta(idx+1));
                    xdata(idx+1)=-xtemp;
                    ydata(idx+1)=-ytemp;
                end

            end
            %% 求速度
            x0=xdata(idx);
            y0=ydata(idx);
            x1=xdata(idx+1);
            y1=ydata(idx+1);
            k=(y1-y0)/(x1-x0);  % 割线斜率
            if state(idx)==0
                k0=(sin(theta(idx))+theta(idx)*cos(theta(idx)))/(cos(theta(idx))-theta(idx)*sin(theta(idx)));  % 切线0斜率
            elseif state(idx)==3
                k0=-(sin(theta(idx))+theta(idx)*cos(theta(idx)))/(cos(theta(idx))-theta(idx)*sin(theta(idx)));  % 切线0斜率
            elseif state(idx)==1
                k0=-1/((yo1-y0)/(xo1-x0));
            else
                k0=-1/((yo2-y0)/(xo2-x0));
            end
            if state(idx+1)==0
                k1=(sin(theta(idx+1))+theta(idx+1)*cos(theta(idx+1)))/(cos(theta(idx+1))-theta(idx+1)*sin(theta(idx+1)));  % 切线1斜率
            elseif state(idx+1)==3
                k1=-(sin(theta(idx+1))+theta(idx+1)*cos(theta(idx+1)))/(cos(theta(idx+1))-theta(idx+1)*sin(theta(idx+1)));  % 切线1斜率
            elseif state(idx+1)==1
                k0=-1/((yo1-y1)/(xo1-x1));
            else
                k0=-1/((yo2-y1)/(xo2-x1));
            end
            the0=abs(abs(atan(k0))-abs(atan(k)));
            the1=abs(abs(atan(k))-abs(atan(k1)));
            vdata(idx+1)=vdata(idx)*cos(the0)/cos(the1);
            %% 绘图
            if 0
                figure
                [x_std,y_std]=get_xy(get_curve(b,boundary_theta:0.1:44*pi),boundary_theta:0.1:44*pi);
                [x_circle,y_circle]=get_xy(return_radius,0:0.1:2*pi);
                S_x=[];
                S_y=[];
                for sss=0:0.1:calculate_return_length(return_radius)
                    [S_x(end+1),S_y(end+1)]=from_s_get_return_xy(sss,return_radius,rotate);
                end
                hold on
                plot(x_std,y_std,"-",'Color',"b")
                plot(-x_std,-y_std,"-",'Color',"g")
                plot(S_x,S_y,"-",'Color',"k")
                plot(x_circle,y_circle,"--",'Color',[0.1,0.7,0.1],'LineWidth',1)
                plot(xdata(:),ydata(:),"-",'Color',[0.1,0.3,0.9])
                plot(xdata(:),ydata(:),"*",'Color',[0.9,0.1,0.1])
                axis square
            end
        end
        
%         if abs(mod(t*dt,1))<1e-10
%             rec_x(:,round(t*dt)+1)=xdata(:,1)/100;
%             rec_y(:,round(t*dt)+1)=ydata(:,1)/100;
            rec_v(:,t+1)=vdata(:,1);
%         end
        if any(vdata>2)
            bk=1;
            break;
        end
    end
    %% 二分更新
    if bk==1
        % 可掉头,r太大
        v_right=v;
    else
        % 不可掉头,r太小
        v_left=v;
    end
end

%% 画图
figure
[x_std,y_std]=get_xy(get_curve(b,boundary_theta:0.1:44*pi),boundary_theta:0.1:44*pi);
[x_circle,y_circle]=get_xy(return_radius,0:0.1:2*pi);
S_x=[];
S_y=[];
for sss=0:0.1:calculate_return_length(return_radius)
    [S_x(end+1),S_y(end+1)]=from_s_get_return_xy(sss,return_radius,rotate);
end
hold on
plot(x_std,y_std,"-",'Color',"b")
plot(-x_std,-y_std,"-",'Color',"g")
plot(S_x,S_y,"-",'Color',"k")
plot(x_circle,y_circle,"--",'Color',[0.1,0.7,0.1],'LineWidth',1)
plot(xdata(:),ydata(:),"-",'Color',[0.1,0.3,0.9])
plot(xdata(:),ydata(:),"*",'Color',[0.9,0.1,0.1])
axis square

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

%% 盘入递推龙身
function t=get_body_theta(b,theta_bef,R)
    % in:b-螺距(cm);theta_bef-前节点极角(rad),R-龙身节点间长度(cm);rerurn_radius-掉头空间半径(cm)
    % out:t-下一节点极角(rad)
    [x0,y0]=get_xy(get_curve(b,theta_bef),theta_bef);
    p0=[x0,y0];
    ktemp=(sin(theta_bef)+theta_bef*cos(theta_bef))/(cos(theta_bef)-theta_bef*sin(theta_bef));  % 切线斜率
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

%% 盘出递推龙身
function t=get_body_theta_out(b,theta_bef,R)
    % in:b-螺距(cm);theta_bef-前节点极角(rad),R-龙身节点间长度(cm);rerurn_radius-掉头空间半径(cm)
    % out:t-下一节点极角(rad)
    [x0,y0]=get_xy(get_curve(b,theta_bef),theta_bef);
    p0=[x0,y0];
    ktemp=(sin(theta_bef)+theta_bef*cos(theta_bef))/(cos(theta_bef)-theta_bef*sin(theta_bef));  % 切线斜率
    gamma=atan(ktemp);
    rho=get_rho(b,theta_bef);
    point_E1=[x0+rho*sin(gamma),y0-rho*cos(gamma)];
    point_E2=[x0-rho*sin(gamma),y0+rho*cos(gamma)];
    if norm(point_E1)<norm(point_E2)
        point_E=point_E1;
    else
        point_E=point_E2;
    end
    t=get_cross_point_out(theta_bef,point_E,rho,p0,R,R,x0,y0);
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
        else
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

%% 求盘出两圆交点
function theta=get_cross_point_out(theta_bef,p1,r1,p2,r2,R,x0,y0)
    % p1,p2-输入点横向量;r1,r2-圆半径
    d=norm(p1-p2);
    a=(r1^2-r2^2+d^2)/(2*d);
    h=sqrt(r1^2-a^2);
    p=p1+a*(p2-p1)/d;
    x1=p(1)+h*(p2(2)-p1(2))/d;
    y1=p(2)-h*(p2(1)-p1(1))/d;
    x2=p(1)-h*(p2(2)-p1(2))/d;
    y2=p(2)+h*(p2(1)-p1(1))/d;
    if norm([x1,y1])<norm([x2,y2])
        x=x1;
        y=y1;
    else
        x=x2;
        y=y2;
    end
    R1=norm(p2);
    R2=norm([x,y]);
    dtheta=acos((R1^2+R2^2-R^2)/(2*R1*R2));
    theta=-(dtheta)+theta_bef;
end

%% 求掉头曲线s对应的两相对theta
function [theta1,theta2]=get_return_theta(s,return_radius)
    r1=return_radius*2/3;
    r2=return_radius/3;
    theta1=-1;
    theta2=-1;
    s1=pi*r1;
    s2=pi*r2;
    if s<s1 % 在前圆内
        theta1=pi*s/s1;
    else % 在后圆内
        theta2=pi*(s-s1)/s2;
    end
end

%% 掉头曲线s自然坐标转直角坐标
function [x,y]=from_s_get_return_xy(s,return_radius,rotate)
    % in:s-自然坐标(cm),return_radius-掉头空间半径(cm),rotate-旋转角度(rad)
    % out:x,y-坐标(cm)
    [theta1,theta2]=get_return_theta(s,return_radius);
    r1=return_radius*2/3;
    r2=return_radius/3;
    if theta1 >= 0   % 在前圆内
        x1=r1*sin(theta1);
        y1=r1*cos(theta1)+return_radius/3;
    else  % 在后圆内
        x1=-r2*sin(theta2);
        y1=r2*cos(theta2)-return_radius*2/3;
    end
    %% 旋转
    T=[cos(rotate),sin(rotate);
        -sin(rotate),cos(rotate)];
    res=T*[x1;y1];
    x=res(1);
    y=res(2);
end

function s=calculate_return_length(r)
    % in:r-掉头空间半径(cm)
    % out:s-掉头曲线长度(cm)
    s=pi*r;
end