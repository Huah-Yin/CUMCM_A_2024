%% ��ʼ��
clear
close all
clc

%% ������ֵ
m1=4866;
mf=1028.876;
R0=1;
h1=3;
h2=0.8;
m2=2433;
R1=0.5;
h3=0.5;
rou=1025;
g=9.8;
k=80000;
d0=0.5;
f=3640;
omg=1.7152;
lmd=683.4558;
mu=654.3383;
beta=10000;
gamma=1000;
syms d;

%% ��������
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(2*pi*R0*h1+pi*R0*h2)/m1;
rz2=d+h3/2;
rz=(m1*rz1+m2*rz2)/(m1+m2);

%% ����ת������
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;
I2=m2*R0^2/4+m2*h3^2/12+m2*(d+h3/2)^2;
t0=0;
tf=40*2*pi/omg;

%% ���΢�ַ����顢��ͼ����������
global wc;
index=1;
for wc=[-0.01,0,0.01]
    [t,n{index}]=ode45('test_func',t0:0.2:tf,[0,0,0,0,0,0,0,0]);
    index=index+1;
end
for i=1:8
    cha(i,1)=mean(n{2}(:,i)-n{1}(:,i));
    cha(i,2)=mean(n{3}(:,i)-n{1}(:,i));
end
titl={'����λ��','����λ��','�����ٶ�','�����ٶ�','���ӽ�λ��','���ӽ�λ��','���ӽ��ٶ�','���ӽ��ٶ�'};
ylab={'λ��(m)','λ��(m)','�ٶ�(m/s)','�ٶ�(m/s)','��λ��(rad)','��λ��(rad)','���ٶ�(rad/s)','���ٶ�(rad/s)'};
figure
for i=1:8
    subplot(2,4,i)
    hold on
    plot(t(:),n{1}(:,i),'r-');
    plot(t(:),n{2}(:,i),'g-');
    plot(t(:),n{3}(:,i),'b-');
    legend('�Ŷ�-1%','���Ŷ�','�Ŷ�+1%')
    xlabel('ʱ��(s)');
    ylabel(ylab{i});
    title(titl{i})
    axis square
    hold off
end
figure
titl2={'����','ҡ��'};
ylab2={{'ƽ��λ�Ʋ�(m)','ƽ���ٶȲ�(m/s)'},{'ƽ����λ�Ʋ�(rad)','ƽ�����ٶȲ�(rad/s)'}};
xlab2={{'����λ��','����λ��','�����ٶ�','�����ٶ�'},{'���ӽ�λ��','���ӽ�λ��','���ӽ��ٶ�','���ӽ��ٶ�'}};
for i=1:2
    subplot(1,2,i)
    hold on
    [AX,H1,H2]=plotyy(1:2:8,cha(4*(i-1)+1:4*(i-1)+4,1),2:2:8,cha(4*(i-1)+1:4*(i-1)+4,2),'bar','bar');
    set(H1,'BarWidth',0.3,'FaceColor',[0.9 0.9 0.9])
    set(H2,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
    set(AX(1),'xtick',1.5:2:8.5,'xlim',[0 9]);
    set(AX(1),'XTickLabel',xlab2{i});
    set(AX(2),'xtick',1.5:2:8.5,'xlim',[0 9]);
    set(AX(2),'XTickLabel','')
    set(get(AX(1),'ylabel'),'string',ylab2{i}{1});
    set(get(AX(2),'ylabel'),'string',ylab2{i}{2});
    title(titl2{i})
    legend('�Ŷ�-1%','�Ŷ�+1%','Location','SouthWest')
    hold off
end