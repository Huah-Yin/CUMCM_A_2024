function flag = Q2_judge( x )
%判断x是否符合约束条件，若符合，则返回1
flag=0;
if x(2)>5+x(5) && x(3)>x(2) && x(4)>x(6)/2 && x(5)>x(6)
    flag=1;
    return
end
end