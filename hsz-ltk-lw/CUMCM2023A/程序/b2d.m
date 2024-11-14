function d=b2d(str,Xx,Xs,L)
%格雷解码并将二进制转为定义域内的十进制
%b2d(str,Xx,Xs,L),str-二进制行向量，Xx-下限，Xs-上限，L单维度二进制串长度
    Dim=length(Xx);
    m(Dim)=0;
    n=zeros(1,L*Dim);
    for k=1:Dim
        n((k-1)*L+1)=str((k-1)*L+1);  % 格雷编码解码
        for j=2:L
            n((k-1)*L+j)=xor(n((k-1)*L+j-1),str((k-1)*L+j));      % xor—异或运算
        end
    end
    for di=1:Dim
        for j=1:L   % 二进制转十进制
            m(di)=m(di)+n((di-1)*L+j)*2^(j-1);
        end
        d(di)=Xx(di)+m(di)*(Xs(di)-Xx(di))/(2^L-1);  % 十进制解码至定义域
        m(di)=0;
    end
end