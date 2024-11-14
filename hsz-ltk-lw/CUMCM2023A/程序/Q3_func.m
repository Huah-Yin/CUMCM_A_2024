function ret = Q3_func(data)
%问题三求解函数
% data为4列矩阵，列1、2：定日镜坐标，列3：镜面宽度，列4：镜面高度，列5：安装高度
%% 参数赋值
step=2;   % 网格生成步长
H=3;
h2=80;
h3=8;
R0=3.5;
Rm=100;
d=149597870700; % 日地距离
R1=6.963e8; % 太阳半径
fai=39.4;
q=0.05;   % 遮挡角度允许误差（角度）
%% 结果初始化
yita_sb_day=zeros(1,12);
yita_cos_day=zeros(1,12);
yita_at_day=zeros(1,12);
yita_trunc_day=zeros(1,12);
yita_ref_day=ones(1,12).*0.92;
%% 遍历时间
D_ind=0;
for D=[-59,-28,0,31,61,92,122,153,184,214,245,275]
    D_ind=D_ind+1;
    ST_ind=0;
    %% 结果初始化
    yita_sb_time=zeros(1,5);
    yita_cos_time=zeros(1,5);
    yita_at_time=zeros(1,5);
    yita_trunc_time=zeros(1,5);
    E_time=zeros(1,5);
    for ST=[9,10.5,12,13.5,15]
        ST_ind=ST_ind+1;
        E=0;    % 输出热功率初始化
        %% 计算太阳高度角与方位角
        de=asin(sin(2*pi*D/350)*sin(2*pi*23.45/360));
        omg=pi*(ST-12)/12;
        as=asin(cos(de)*cos(fai*pi/180)*cos(omg)+sin(de)*sin(fai*2*pi/360))*180/pi;          % 太阳高度角
        gs=acos((sin(de)-sin(as*pi/180)*sin(fai*2*pi/360))/cos(as*pi/180)*cos(fai*2*pi/360))*180/pi;  % 太阳方位角
        %% 计算DNI
        DNI=1.366*((0.4237-0.00821*(6-H)^2)+(0.5055+0.00595*(6.5-H)^2)*exp(-(0.2711+0.01858*(2.5-H)^2)/sin(as*pi/180)));
        %% 遍历定日镜
        len=length(data);
        in=cell(1,len);
        yita_sb=zeros(len,1);
        yita_cos=zeros(len,1);
        yita_at=zeros(len,1);
        yita_trunc=ones(len,1).*(-1);
        for num=1:len
            x=data(num,1);
            y=data(num,2);
            r=sqrt(x^2+y^2);
            %% 计算定日镜高度角
            theta1=atan((h2-data(num,5))/r)*180/pi;
            am=theta1+(as-theta1)/2;
            %% 计算定日镜方位角
            if y>0 && ( x>=0 || -atan(x/y)*180/pi+gs<=180 ) % 情况1
                beta=180*atan(y/x)/pi;
                theta2=beta+gs-sign(beta)*90;
                l=sqrt(r^2+d^2-2*r*d*cos(theta2*pi/180));
                theta3=180*acos((r^2+l^2-d^2)/(2*r*l))/pi;
                if x>=0
                    gm=180-(theta3/2-(90-beta));
                else
                    gm=-(theta3/2-(90-beta));
                end
            elseif y<0 && ( x<=0 || -atan(x/y)*180/pi+gs<=180 ) % 情况2
                beta=180*atan(x/y)/pi;
                theta4=beta+180-gs;
                l=sqrt(r^2+d^2-2*r*d*cos(theta4*pi/180));
                gm=beta+1/2*180*acos((r^2+l^2-d^2)/(2*r*l))/pi;
            elseif  y>=0 && ( x<=0 && -atan(x/y)*180/pi+gs>=180 ) % 情况3
                beta=180*abs(atan(x/y))/pi;
                theta5=(180-beta)+(180-gs);
                l=sqrt(r^2+d^2-2*r*d*cos(theta5*pi/180));
                gm=180-(beta-1/2*180/pi*acos((r^2+l^2-d^2)/(2*r*l)));
            else % 情况4
                beta=180*atan(x/y)/pi;
                theta6=180-gs+beta;
                l=sqrt(r^2+d^2-2*r*d*cos(theta6*pi/180));
                gm=180-1/2*180/pi*acos((r^2+l^2-d^2)/(2*r*l))+beta;
            end
            bla_num=0;
            unbla_num=0;
            in_num=0;
            unin_num=0;
            %% 反射光锥轴线向量计算
            s0=-[cos(as*pi/180)*sin(gs*pi/180);cos(as*pi/180)*cos(gs*pi/180);sin(as*pi/180)];
            n=[cos(am*pi/180)*sin(gm*pi/180);cos(am*pi/180)*cos(gm*pi/180);sin(am*pi/180)];
            s_dot=s0+2*abs(dot(s0,n))*n;
            s_dot=s_dot./norm(s_dot);
            %% 划分网格
            ind_max=floor(data(num,3)/step+1)*floor(data(num,4)/step+1);
            cosn=zeros(1,ind_max);
            ind=1;
            for i=-data(num,3)/2:step:data(num,3)/2
                for k=-data(num,4)/2:step:data(num,4)/2
                    x0=i;
                    y0=0;
                    z0=k;     
                    %% 变换
                    A=[cos(gm*pi/180),sin(gm*pi/180),0;-sin(gm*pi/180),cos(gm*pi/180),0;0,0,1];
                    B=[1,0,0;0,cos(am*pi/180),-sin(am*pi/180);0,sin(am*pi/180),cos(am*pi/180)];
                    t=A*B*[x0;0;z0];
                    x1=t(1)+x;
                    y1=t(2)+y;
                    z1=t(3)+data(num,5);
                    %% 遮挡判断
                    s_f=[cos(as*pi/180)*sin(gs*pi/180);cos(as*pi/180)*cos(gs*pi/180);sin(as*pi/180)];
                    if x1<0 && y1>0   % 仅判断第二象限是否被遮挡
                        jijiao=(atan(y1/x1)+pi)*180/pi; % 极角
                        theta_s=(180-gs)+90;
                        det_j=atan(R0/Rm)*180/pi+q;
                        if jijiao>theta_s-det_j && jijiao<theta_s+det_j
                            det=4*(x1*s_f(1)+y1*s_f(2))^2-4*(s_f(1)^2+s_f(2)^2)*(x1^2+y1^2-R0^2);
                            if det<0
                                black=1;  % 不被遮挡
                                unbla_num=unbla_num+1;
                            else
                                t_sol=(-2*(x1*s_f(1)+y1*s_f(2))-sqrt(det))/(2*(s_f(1)^2+s_f(2)^2));
                                z_sol=z1+s_f(3)*t_sol;
                                if z_sol>=0 && z_sol<=h2+h3/2
                                    black=0;  % 被遮挡
                                    bla_num=bla_num+1;
                                else
                                    black=1;  % 不被遮挡
                                    unbla_num=unbla_num+1;
                                end
                            end
                        else
                        black=1;
                        unbla_num=unbla_num+1;
                        end
                    else
                        black=1;
                        unbla_num=unbla_num+1;
                    end
                    %% 截断效率计算
                    if black==1
                        det=4*(x1*s_dot(1)+y1*s_dot(2))^2-4*(s_dot(1)^2+s_dot(2)^2)*(x1^2+y1^2-R0^2);
                        if det<0
                            unin_num=unin_num+1;    % 不能射入
                        else
                            t_sol=(-2*(x1*s_dot(1)+y1*s_dot(2))-sqrt(det))/(2*(s_dot(1)^2+s_dot(2)^2));
                            z_sol=z1+s_dot(3)*t_sol;
                            
                            if z_sol>=h2-h3/2 && z_sol<=h2+h3/2
                                in_num=in_num+1;    % 能射入
                                in{num}(ind)=0;
                            else
                                unin_num=unin_num+1;    % 不能射入
                                in{num}(ind)=1;
                            end
                        end
                    end
                    %% 计算余弦效率
                    n_f=[cos(am*pi/180)*sin(gm*pi/180);cos(am*pi/180)*cos(gm*pi/180);sin(am*pi/180)];
                    cosn(ind)=real(dot(s_f,n_f));
                    ind=ind+1;
                end
            end
            %% 统计阴影遮挡效率
            yita_sb(num)=1-bla_num/(bla_num+unbla_num);
            %% 计算单板平均余弦效率
            yita_cos(num)=mean(cosn);
            %% 计算大气透射率
            dHR=sqrt(x^2+y^2+(h2-data(num,5))^2);
            yita_at(num)=0.99321-0.0001176*dHR+1.97e-8*dHR^2;
            %% 计算截断效率
            if yita_sb(num)~=0
                yita_trunc(num)=1-in_num/(in_num+unin_num);
            else
                yita_trunc(num)=0;
            end
            yita=yita_sb(num)*yita_cos(num)*yita_at(num)*yita_trunc(num)*yita_ref_day(D_ind);
            E=E+data(num,3)*data(num,4)*yita;
            
        end
        E=E*DNI;
        yita_sb_time(ST_ind)=mean(yita_sb);
        yita_cos_time(ST_ind)=mean(yita_cos);
        yita_at_time(ST_ind)=mean(yita_at);
        yita_trunc_time(ST_ind)=mean(yita_trunc);
        E_time(ST_ind)=E/(sum(data(:,3).*data(:,4)));
        E_all_time(ST_ind)=E;
    end
    fprintf('\t完成计算%d月...\n',D_ind);
    yita_sb_day(D_ind)=mean(yita_sb_time);
    yita_cos_day(D_ind)=mean(yita_cos_time);
    yita_at_day(D_ind)=mean(yita_at_time);
    yita_trunc_day(D_ind)=mean(yita_trunc_time);
    E_day=mean(E_time);
    E_all_day(D_ind)=mean(E_all_time);
    %% 结果整理
    res{D_ind,1}=yita_sb_day(D_ind)*yita_cos_day(D_ind)*yita_at_day(D_ind)*yita_trunc_day(D_ind)*yita_ref_day(D_ind);
    res{D_ind,2}=yita_cos_day(D_ind);
    res{D_ind,3}=yita_sb_day(D_ind);
    res{D_ind,4}=yita_trunc_day(D_ind);
    res{D_ind,5}=E_day;
    E_year=E_day;
end
E_all_year=mean(E_all_day(D_ind))+6e3;
fprintf('年平均输出热功率%d，单位面积功率%d\n',E_all_year,E_year)
ret=E_all_year^E_year;