clc,clear,close all;
t = 100;
dt = 0.1;
dx = 1.7;
n = round(16*0.55/dx);
flag =1;
for v = 1:0.0001:2
  [~, ~, ~, v0] = new_prepare(t, dx, dt, n,v);
 for j = 1:224
     for i = 1:1001
         if v0(i,j)>=2
          flag = 0;
         end
     end
 end
 if flag == 0
         disp(v);
         break;
end
end