function d=b2d(str,Xx,Xs,L)
%���׽��벢��������תΪ�������ڵ�ʮ����
%b2d(str,Xx,Xs,L),str-��������������Xx-���ޣ�Xs-���ޣ�L��ά�ȶ����ƴ�����
    Dim=length(Xx);
    m(Dim)=0;
    n=zeros(1,L*Dim);
    for k=1:Dim
        n((k-1)*L+1)=str((k-1)*L+1);  % ���ױ������
        for j=2:L
            n((k-1)*L+j)=xor(n((k-1)*L+j-1),str((k-1)*L+j));      % xor���������
        end
    end
    for di=1:Dim
        for j=1:L   % ������תʮ����
            m(di)=m(di)+n((di-1)*L+j)*2^(j-1);
        end
        d(di)=Xx(di)+m(di)*(Xs(di)-Xx(di))/(2^L-1);  % ʮ���ƽ�����������
        m(di)=0;
    end
end