function result = Q3_judge( arg )
%�ж�x�Ƿ����Լ�������������ϣ��򷵻�1�������߶Ȳ��������Զ�����������2�������ʮ��������
ALL_NUM=3500;
% arg=rand(1,5+2*ALL_NUM);
% arg(1:5)=[0.003,2,10,15,100];
flag=1;
if 700*arg(1)+arg(2)>6
    result={0};
    fprintf('�¶ȴ���')
    return
end
%% ��������������
xi=arg(5)*[cos(0.4503),-sin(0.4503)];
%% �������ɶ��վ�����
aff_num=-1;  % ��ǰ�ж��վ���Ŀ
lin_ind=1;  % �к�
distan=100; % ��ǰ������������
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
            %% ���㰲װ�߶�
            d=abs(-coo_now(1)+0.48*coo_now(2)+388.75)/sqrt(1+0.48^2);
            if arg(2)+arg(1)*d <= arg(5+ALL_NUM+num_coo)/2
                result={0};
                return
            end
            %% ��߶�����
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