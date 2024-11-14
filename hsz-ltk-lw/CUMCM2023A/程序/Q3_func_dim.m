function E_year = Q3_func_dim(arg)
%������Ŀ�꺯��   
arg=xlsread("Q3_result(1).xlsx");
ALL_NUM=3500;
%% ��������������
xi=arg(5)*[cos(-77.3*pi/180),sin(-77.3*pi/180)];
%% �������ɶ��վ�����
aff_num=-1;  % ��ǰ�ж��վ���Ŀ
lin_ind=1;  % �к�
distan=100; % ��ǰ������������
coo_ls=[];  % �����б�
num_coo=0;  % ��ǰ������
while aff_num ~= 0 && num_coo<ALL_NUM  % ѭ����
    %% ʵ��do-whileѭ��
    aff_num=0;
    %% ����������������ļ�����
    the=0;  % ��һ��
    now_the=0;  % ��ǰת��
    while now_the+the<2*pi-the
        if num_coo>=ALL_NUM
            break
        end
        %% תֱ������
        now_the=now_the+the;
        coo_now=distan.*[cos(now_the),sin(now_the)];
        %% �ı�ο�ϵ
        coo_now=xi-coo_now;
        %% �ж��Ƿ��ڳ�����
        the_app=0;
        if dot(coo_now,coo_now)<=122500
            num_coo=num_coo+1;
            coo_ls(num_coo,1)=coo_now(1);
            coo_ls(num_coo,2)=coo_now(2);
            %% ���㰲װ�߶�
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
%% ��ͼ
close all
figure
hold on
scatter3(coo_ls(:,1), coo_ls(:,2),coo_ls(:,3))
rectangle('Position', [0-350,0-350,2*350,2*350], 'Curvature', [1 1],'EdgeColor', 'r');
plot(xi(1),xi(2),'b*')
axis equal
grid on
fprintf('������ɣ�������Ŀ%d��',len)
end