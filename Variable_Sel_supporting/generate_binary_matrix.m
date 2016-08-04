function binary_matrix = generate_binary_matrix(X_num, row,predefine)
%+++ Binary matrix sampling (BMS) in the variable space  
%+++ Input:  X_num: the number of vairbales in the X matrix
%+++         row: the number of the row of binary matrix
%+++         predefine:the mean number of '1' of every row
%+++ Output: binary_matrix: the matrix contains 0 and 1.
%+++ Yonghuan Yun, Dec.15, 2013. yunyonghuan@foxmail.com

control = 0;
while control == 0
    alpha=round(predefine*row/X_num)/row;
    num_one=round(alpha*row);
    Initial_binary=[ones(num_one,X_num);zeros(row-num_one,X_num)];
    binary_matrix=zeros(row,X_num);
    for i=1:X_num
        RandNumber=randperm(row);
        binary_matrix(:,i)=Initial_binary(RandNumber,i);
    end

    for m = 1:row
        if sum(binary_matrix(m,:)) <= 1   % Guarantee that every row has at least two '1'
            control = 0;
            break; 
        else
            control = 1;
        end
    end
end
