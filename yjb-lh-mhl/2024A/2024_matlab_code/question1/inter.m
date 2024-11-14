syms x a v
f = (a ./ v) .* sqrt(x ^ 2 + 1); % 定义符号表达式
I = int(f, x); % 计算不定积分

% 使用 simplify 简化表达式
simplified_I = simplify(I)

% % 使用 rewrite 将表达式转换为对数形式
log_form = rewrite(simplified_I, 'log')

% 使用 matlabFunction 将表达式转换为 MATLAB 函数，并保存到文件
matlabFunction(log_form, 'File', 'generated_function.m');
