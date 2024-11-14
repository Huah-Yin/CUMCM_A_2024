function str = d2b (arr,Xx,Xs,L)
%格雷编码
%d2b(arr,Xx,Xs,L),arr-十进制数据数组，Xx-下限，Xs-上限，L单维度二进制串长度；输出str为二进制行向量
Dim=length(Xx);
str=zeros(1,Dim*L);
m=zeros(1,Dim);
n=zeros(1,Dim);
for k=1:Dim
    m(k)=round((arr(k)-Xx(k))*(2^L-1)/(Xs(k)-Xx(k)));   % 分散至二进制定义域
    for i=1:L   % 十进制转二进制
        n((k-1)*L+i)=mod(m(k),2);
        m(k)=floor(m(k)/2);
    end
end
for k=1:Dim
    str((k-1)*L+1)=n((k-1)*L+1);
    for i=2:L   % 格雷编码
        str((k-1)*L+i)=xor(n((k-1)*L+i-1),n((k-1)*L+i));
    end
end
end