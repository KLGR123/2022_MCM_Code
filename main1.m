clear;
clc;
% ---------------
% 模型一的参数接口
type = 4; %
SI = 1;
%k = 25;  %
Q = 0.33; %
size = 10000;
w = 17/20;
year = 200;
grow_type = 0;
cut_type = 0; %
shape = 0; %
tree_pro = 0.4; %

% ---------------
t = [1:year];
% [B, P, C, S, tree, tree_num] = Forest(type, SI, k, Q, size, w, year, grow_type);
% [B, P, C, S, tree, tree_num, tree_type, C_max] = Forest_update(type, SI, k, Q, size, w, year, grow_type, cut_type, shape, tree_pro);

Clist = zeros(9, year);
%k = [1:80];
k = [5, 10, 15, 20, 25, 30, 35, 50, 80];
for ii = 1:9
    [B, P, C, S, tree, tree_num, tree_type, C_max] = Forest_update(type, SI, k(ii), Q, size, w, year, grow_type, cut_type, shape, tree_pro); 
    Clist(ii, 1:year) = C(1:year); 
end
plot(Clist(1))
