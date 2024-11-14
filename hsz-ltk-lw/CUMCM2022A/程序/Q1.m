%% ��ʼ��
clear
close all
clc

%% ������ֵ
m1=4866;
mf=1335.535;
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
f=6250;
omg=1.4005;
lmd=656.3616;
t=0;
t0=0;
tf=40*2*pi/omg;

%% ���΢�ַ�����
beta=10000;
P12=-m1-mf;
P11=beta-lmd;
P21=-beta;
P10=rou*g*pi*R0*R0+k;
P20=-k;
P=m1*g-rou*g*pi*R0*R0*h1-rou*g*pi*R0*R0*h2/3+f*cos(omg*t)-k*(h1-h3-d0);
Q22=-m2;
Q11=-beta;
Q21=beta;
Q10=-k;
Q20=k;
Q=m2*g+k*(h1-h3-d0);

%% (1)
[t,m]=ode45('Q1_func1',t0:0.2:tf,[0,0,0,0]);
figure
subplot(1,2,1)
hold on
plot(t(:),m(:,1),'r-');
plot(t(:),m(:,2),'b-');
legend('����λ��','����λ��')
title('��1С��λ��-ʱ��ͼ��')
xlim([0,180])
xlabel('ʱ��(s)');
ylabel('λ��(m)');
hold off
subplot(1,2,2)
hold on
plot(t(:),m(:,3),'r-');
plot(t(:),m(:,4),'b-');
legend('�����ٶ�','�����ٶ�')
title('��1С���ٶ�-ʱ��ͼ��')
xlim([0,180])
xlabel('ʱ��(s)');
ylabel('�ٶ�(m/s)');
hold off
A=t;
A(:,2:5)=[m(:,1),m(:,3),m(:,2),m(:,4)];
xlswrite('����1(1)����.xlsx',A,1);

%% (2)
[t,n]=ode45('Q1_func2',t0:0.2:tf,[0,0,0,0]);
figure
subplot(1,2,1)
hold on
plot(t(:),n(:,1),'r-');
plot(t(:),n(:,2),'b-');
legend('����λ��','����λ��')
title('��2С��λ��-ʱ��ͼ��')
xlim([0,180])
xlabel('ʱ��(s)');
ylabel('λ��(m)');
hold off
subplot(1,2,2)
hold on
plot(t(:),n(:,3),'r-');
plot(t(:),n(:,4),'b-');
legend('�����ٶ�','�����ٶ�')
title('��2С���ٶ�-ʱ��ͼ��')
xlim([0,180])
xlabel('ʱ��(s)');
ylabel('�ٶ�(m/s)');
hold off
B=t;
B(:,2:5)=[n(:,1),n(:,3),n(:,2),n(:,4)];;
xlswrite('����1(2)����.xlsx',B,1);
ind=2;
C={'ʱ��','(1)����λ��','(1)�����ٶ�','(1)����λ��','(1)�����ٶ�','(2)����λ��','(2)�����ٶ�','(2)����λ��','(2)�����ٶ�'};
for i=1:length(A(:,1))
    if ismember(A(i,1),[10 20 40 60 100])
        for j=1:5
            C{ind,j}=A(i,j);
        end
        for j=6:9
            C{ind,j}=B(i,j-4);
        end
        ind=ind+1;
   end
end

%% ����Ҷ����
L=256;
% ����λ��
% A1=fft(m(end-area+1:end,1));
% f=1/0.2/area*(-area/2:area/2-1);
% plot(f,abs(fftshift(A1)))
Fs=t(end)-t(end-L+1);
f = Fs/L*(0:(L/2));

Y = fft(m(end-L+1:end,1));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")