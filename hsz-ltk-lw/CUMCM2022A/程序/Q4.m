%% ����4-��ά�Ŵ��㷨
%% ��ʼ��
clear;
close all;      % ��ͼ
clc;
NP=100;         % ��Ⱥ����
L=30;           % ����ά�ȶ��������γ���
Dim=2;          % ά��
Pc=0.5;         % ������
Pm=0.2;         % ������
G=150;          % ����Ŵ�����
Xs=[100000,100000];              % ����������
Xx=[0,0];            % ����������
f=randi([0,1],NP,L*Dim);    % �����ȡ��ʼ��Ⱥ

%% ����㶨��������ת������������ȫ�ֱ���֮��
m1=4866;
R0=1;
h1=3;
h2=0.8;
global I1 rz1;
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(2*pi*R0*h1+pi*R0*h2)/m1;
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;

%% �Ŵ��㷨
tic
% Ԥ��������
trace1(G)=0;
trace2(G)=0;
trace3(G)=0;
maxIndex(G,Dim)=0;

for k=1:G
    fprintf('���ڼ�������4��%d��...',k);
    %% �����ƽ���Ϊ�������ڵ�ʮ����
    for i=1:NP
        U=f(i,:);   % ��ȡ��i��Ⱦɫ��
        x(i,:)=b2d(U,Xx,Xs,L);
        Fit(i)=Q4_func(x(i,:));          % ����ÿ��������Ӧ��
    end
    maxFit=max(Fit);    % �����Ӧ��
    minFit=min(Fit);    % ��С��Ӧ��
    rr=find(maxFit==Fit);   % �ҳ������Ӧ�����ڸ���������������������
    fBest=f(rr(1,1),:);     % �������Ÿ���
    xBest=x(rr(1,1),:);       % ���Ÿ����Ӧ��ʮ����Ⱦɫ��
    Fit1=(Fit-minFit)/(maxFit-minFit);   % ��Ӧ�ȹ�һ��
    
    %% �������̶�ѡ�񷨵ĸ��Ʋ���
    sum_Fit=sum(Fit1);   % ������Ⱥ�����и�����Ӧֵ�ĺ�
    fitvalue=Fit1/sum_Fit;   % ����ÿ����Ⱥ��ѡ�����
    fitvalue=(cumsum(fitvalue));    % ����ÿ����Ⱥ���ۼƸ���
    ms=sort(rand(NP,1));    % �������NP��(0,1)��ֵ��������
    fiti=1;     % ����Ⱥ��ǰָ��
    newi=1;     % ����Ⱥ��ǰָ��
    while newi<=NP      % ������ƣ���ʹ��Ӧ�ȴ���Ŵ���ȥ
        if ms(newi) < fitvalue(fiti)
            nf(newi,:)=f(fiti,:);   % ����
            newi=newi+1;
        else
            fiti=fiti+1;
            if fiti>NP
                break;
            end
        end
    end
    
    %% ���ڸ��ʵĽ������
    for i=1:NP-1
        p=rand;     % �������һ������[0,1]�ĸ���p
        if p<Pc     % ���㽻������
            q=randi([0,1],1,L*Dim);     % �������Ҫ����Ļ���λ��
            for j=1:L*Dim
                if q(j)==1  % ���´���ʵ�ֱ�������
                    temp=nf(i+1,j);
                    nf(i+1,j)=nf(i,j);
                    nf(i,j)=temp;
                end
            end
        end
    end
    
    %% ���ڸ��ʵı������
    i=1;
    while i<=round(NP*Pm)   % ���Ʊ���Ⱦɫ������
        h=randi([1,NP],1,1);    % ���ѡ��һ����Ҫ�����Ⱦɫ������
        for j=1:round(L*Dim*Pm) % ���Ʊ��������
            g=randi([1,L*Dim],1,1); % �����Ҫ����Ļ�������
            nf(h,g)=~nf(h,g);   % ȡ��
        end
        i=i+1;
    end
    
    %% ��һ��Ԥ��
    f=nf;   % ��Ⱥλ�ð��
    f(1,:)=fBest;   % �����ϴ����Ÿ���
    trace1(k)=maxFit;       % ����������Ӧ��
    trace2(k)=mean(Fit);    % ƽ����Ӧ��
    trace3(k)=minFit;       % ��С��Ӧ��
    maxIndex(k,:)=xBest(:);    % �����Ӧ������
    fprintf('%d�������Ӧ��:%f\n',k,maxFit);
end
toc

%% ������Ӧ�ȹ���ͼ
figure
subplot(2,1,1);
hold on
plot(1:G,trace1(1:G),'^r-');
plot(1:G,trace2(1:G),'ob-');
plot(1:G,trace3(1:G),'vg-');
title(strcat('����4������Ӧ��'));
xlabel('��������');
ylabel('��Ӧ��');
legend({'���ֵ','ƽ��ֵ','��Сֵ'},'Location','SouthEast');
hold off
subplot(2,1,2);
hold on
plot(1:G,trace1(1:G),'^r-');
plot([0,G],[trace1(G),trace1(G)],'LineWidth',1.5);
yl=ylim;
for i=1:G
    if trace1(i)==trace1(G)
        plot(i,trace1(G),'b*','LineWidth',1.5,'MarkerSize',6);
        text(i,trace1(G)-(yl(2)-yl(1))/30,strcat(num2str(i),'��'));
        break
    end
end
title('���������Ӧ��');
xlabel('��������');
ylabel('��Ӧ�����ֵ');
legend({'�������ֵ','ȫ�����ֵ','���δﵽ���ֵ��'},'Location','SouthEast');
hold off

%% ������
fprintf('==========����4==========\n')
fprintf('���ֵ��');
disp(trace1(G));
fprintf('���ֵ�㣺');
disp(maxIndex(G,:));
fprintf('���ȣ�x��')
disp((Xs-Xx)/(2^L-1));