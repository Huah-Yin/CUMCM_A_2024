%% �������ά�����Ŵ��㷨��������Ż�
%% ��ʼ��
clear;
close all;      % ��ͼ
clc;
NP=50;         % ��Ⱥ����
L=10;           % ����ά�ȶ��������γ���
Dim=6;          % ά��
Pc=0.5;         % ������
Pm=0.2;         % ������
Mc=20;          % ������ೢ�Դ���
Mm=50;          % ������ೢ�Դ���
G=100;         % ����Ŵ�����
Xs=[350,15,20,6,8,8];              % ����������
Xx=[0,5,10,2,2,2];            % ����������
f=randi([0,1],NP,L*Dim);    % �����ȡ��ʼ��Ⱥ
boxMax=20;      % ����ͼ���ƴ���
boxPage=20;     % ����ͼ��ͼ��ʾ������
%% �Ŵ��㷨
tic
% Ԥ��������
trace1(G)=0;
trace2(G)=0;
trace3(G)=0;
maxIndex(G,Dim)=0;
fitHis(NP,boxMax)=0;
if boxMax>G     % ����
    boxMax=g;
end
% �޸Ĳ�����Լ�������ĳ�ʼ��Ⱥ��ֱ��ȫ����������Ϊֹ
for i=1:NP
    while ~Q2_judge(b2d(f(i,:),Xx,Xs,L))
    f(i,:)=randi([0,1],1,L*Dim);
    end
end
for k=1:G
    fprintf('���ڼ�������2��%d��...\n',k);
    %% �����ƽ���Ϊ�������ڵ�ʮ����
    for i=1:NP
        tic
        U=f(i,:);   % ��ȡ��i��Ⱦɫ��
        x(i,:)=b2d(U,Xx,Xs,L);
        Fit(i)=Q2_func_dim(x(i,:));          % ����ÿ��������Ӧ��
        fprintf('����%d����Ⱥ���%d����Ӧ��%f\n',k,i,Fit(i))
        toc
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
            t=1;
            while t<Mc
                % ��ʼ���Խ���
                tnf=zeros(1,L*Dim);
                q=randi([0,1],1,L*Dim);     % �������Ҫ����Ļ���λ��
                for j=1:L*Dim
                    tnf(j)=nf(i+q(j),j);
                end
                if Q2_judge(b2d(tnf,Xx,Xs,L))  % �жϵ�ǰ�������Ƿ����Լ������
                    nf(i,:)=tnf(:);
                    break
                end
                t=t+1;
            end
        end
    end
    %% ���ڸ��ʵı������
    i=1;
    while i<=round(NP*Pm)   % ���Ʊ���Ⱦɫ������
        h=randi([1,NP],1,1);    % ���ѡ��һ����Ҫ�����Ⱦɫ������
        t=1;
        while t<Mm
            tnf=nf(h,:);
            for j=1:round(L*Dim*Pm) % ���Ʊ��������
                g=randi([1,L*Dim],1,1); % �����Ҫ����Ļ�������
                tnf(g)=~tnf(g);   % ȡ��
            end
            if Q2_judge(b2d(tnf,Xx,Xs,L))  % �жϵ�ǰ�������Ƿ����Լ������
                nf(i,:)=tnf(:);
                break
            end
            t=t+1;
        end
        i=i+1;
    end
    %% ��һ��Ԥ��
    if k<=boxMax
        fitHis(:,k)=Fit';
    end
    f=nf;   % ��Ⱥλ�ð��
    f(1,:)=fBest;   % �����ϴ����Ÿ���
    trace1(k)=maxFit;       % ����������Ӧ��
    trace2(k)=mean(Fit);    % ƽ����Ӧ��
    trace3(k)=minFit;       % ��С��Ӧ��
    maxIndex(k,:)=xBest(:);    % �����Ӧ������
end
toc
%% ������Ӧ�ȹ���ͼ
figure
subplot(2,1,1);
hold on
plot(1:k,trace1(1:k),'^r-');
plot(1:k,trace2(1:k),'ob-');
plot(1:k,trace3(1:k),'vg-');
title('������Ӧ��');
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
%% ����������Χ����Ӧ������ͼ
figure
page=ceil(boxMax/boxPage);
for i=1:page
    subplot(page,1,i);
    hold on
    boxplot(fitHis(:,(i-1)*boxPage+1:min([boxMax,i*boxPage])),(i-1)*boxPage+1:min([boxMax,i*boxPage]));
    title(strcat('��',num2str((i-1)*boxPage+1),'~',num2str(min([boxMax,i*boxPage])),'����Ӧ������ͼ'));
    xlabel('��������');
    ylabel('��Ӧ��');
    plot([0,boxMax+1],[trace1(G),trace1(G)],'LineWidth',1.5);
    legend({strcat('���ֵ:',num2str(trace1(G)))},'Location','SouthEast');
    set(gca,'XTickLabelRotation',70);   % ��бx����������
    hold off
end
%% ������
fprintf('���ֵ��');
disp(trace1(G));
fprintf('���ֵ�㣺');
disp(maxIndex(G,:));
fprintf('���ȣ�x��')
disp((Xs-Xx)/(2^L-1));