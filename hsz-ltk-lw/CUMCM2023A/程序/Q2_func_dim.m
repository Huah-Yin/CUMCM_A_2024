function E_year = Q2_func_dim(arg)
% arg=[348.9736,8.9150,8.9443,2.0978,3.6481,3.4428];
%�����Ŀ�꺯��   
%% ��������������
xi=arg(1)*[cos(-77.3*pi/180),sin(-77.3*pi/180)];
%% �������ɶ��վ�����
aff_num=-1;  % ��ǰ�ж��վ���Ŀ
lin_ind=1;  % �к�
distan=100; % ��ǰ������������
coo_ls=[];  % �����б�
num_coo=0;  % ��ǰ������
while aff_num ~= 0  % ѭ����
    %% ʵ��do-whileѭ��
    aff_num=0;
    %% ����������������ļ�����
    if lin_ind<=20
        lin_num=floor(distan*2*pi/arg(2));
    else
        lin_num=floor(distan*2*pi/arg(3));
    end
    the=2*pi/lin_num;
    now_the=0;
    for the_ind=1:lin_num
        %% תֱ������
        now_the=now_the+the;
        coo_now=distan.*[cos(now_the),sin(now_the)];
        %% �ı�ο�ϵ
        coo_now=xi-coo_now;
        %% �ж��Ƿ��ڳ�����
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
%% ��ͼ
% close all
% figure
% hold on
% scatter(coo_ls(:,1), coo_ls(:,2))
% rectangle('Position', [0-350,0-350,2*350,2*350], 'Curvature', [1 1],'EdgeColor', 'r');
% plot(xi(1),xi(2),'b*')
% axis equal
% grid on
fprintf('������ɣ�������Ŀ%d',length(coo_ls))
end