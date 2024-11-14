%% 初始化
clear
close all
clc

%% 参数赋值
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

%% 求解微分方程组
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
legend('浮子位移','振子位移')
title('第1小问位移-时间图像')
xlim([0,180])
xlabel('时间(s)');
ylabel('位移(m)');
hold off
subplot(1,2,2)
hold on
plot(t(:),m(:,3),'r-');
plot(t(:),m(:,4),'b-');
legend('浮子速度','振子速度')
title('第1小问速度-时间图像')
xlim([0,180])
xlabel('时间(s)');
ylabel('速度(m/s)');
hold off
A=t;
A(:,2:5)=[m(:,1),m(:,3),m(:,2),m(:,4)];
xlswrite('问题1(1)结论.xlsx',A,1);

%% (2)
[t,n]=ode45('Q1_func2',t0:0.2:tf,[0,0,0,0]);
figure
subplot(1,2,1)
hold on
plot(t(:),n(:,1),'r-');
plot(t(:),n(:,2),'b-');
legend('浮子位移','振子位移')
title('第2小问位移-时间图像')
xlim([0,180])
xlabel('时间(s)');
ylabel('位移(m)');
hold off
subplot(1,2,2)
hold on
plot(t(:),n(:,3),'r-');
plot(t(:),n(:,4),'b-');
legend('浮子速度','振子速度')
title('第2小问速度-时间图像')
xlim([0,180])
xlabel('时间(s)');
ylabel('速度(m/s)');
hold off
B=t;
B(:,2:5)=[n(:,1),n(:,3),n(:,2),n(:,4)];;
xlswrite('问题1(2)结论.xlsx',B,1);
ind=2;
C={'时间','(1)浮子位移','(1)浮子速度','(1)振子位移','(1)振子速度','(2)浮子位移','(2)浮子速度','(2)振子位移','(2)振子速度'};
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

%% 傅里叶分析
L=256;
% 浮子位移
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