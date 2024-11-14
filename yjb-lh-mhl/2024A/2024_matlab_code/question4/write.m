clc; clear; close all;

% 读取 Excel 文件
filename = 'F:\EXCEL\result6.xlsx'; % 替换为你的文件名
data = readtable(filename);

% 将表格数据转换为矩阵以便操作
data_matrix = table2array(data);
[r, c] = size(data);
% % 使用循环将左上角的值赋给右下角
for i = 178:c % % % %列的循环

    for j = 46:1:r

        if data_matrix(j, i) == 0 && data_matrix(j - 1, i - 2) ~= 0
            data_matrix(j, i) = data_matrix(j - 1, i - 2);
        end

    end

end

% for j = 18:33

%     for i = 158:c

%         if data_matrix(j, i) == 0
%             data_matrix(j, i) = data_matrix(j, 236 - i);
%         else continue;

%         end

%     end

% end

% 将矩阵转换回表格格式
updated_data = array2table(data_matrix, 'VariableNames', data.Properties.VariableNames);

% 保存到 Excel 文件
% 确保写入时路径是正确的
writetable(updated_data, 'F:/EXCEL/result6.xlsx');
disp("数据已成功写入到 Excel 文件");
