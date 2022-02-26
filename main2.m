clear;
clc;
% ---------------
% 本闭环模拟函数只考虑简单的 Management Plan
% 即只有 k 和 Q 两个标量组成
k = 27; % 初始化
Q = [0.33 0];

% ---------------
% 和管理有关的参数 组成更复杂的 Management Plan
type = 1;
cut_type = 0;
shape = 0;
tree_pro = 1;
weight = [0.5 0.25 0.25]; % 由需求 needs 决定

% sensitivity = delta value / delta weight;

% ---------------
% 和管理无关的参数
SI = 1;
size = 100;
w = 17/20;
year = 500;
grow_type = 0;
rain = 30; % 有待确定
humid = 50; 
wind = 30;

% ---------------
% 深度搜索
forcast_depth = 5;

tic
[B, P, C, S, tree, tree_num, tree_type, v1, v2, v3, value, klist, Qlist, q, P_mean, B_mean] = Forest_decision(type, SI, k, Q, size, w, year, grow_type, cut_type, shape, tree_pro, rain, humid, wind, forcast_depth, weight);
toc

t = toc;

B = B * 10000;
P = P * 10000;
C = C * 10000;
v1 = v1*10;
v2 = v2*10;
v3 = v3*100;
value = value*10 + B(1:year);
 
% ------------------------
% 可视化
subplot(3,1,1);
year = 500;
x = 1 : year;
plot(x, P(1:year), x, C(1:year), x, B(1:year), x, P_mean(1:year)*10000, x, B_mean(1:year)*10000); % 三指标时序曲线

subplot(3,1,2); 
plot(x, klist(1:year), x, Qlist(1:year)*100);

subplot(3,1,3);
plot(x, v1(1:year), x, v2(1:year), x, v3(1:year), x, B(1:year), x, value(1:year));

[mean(v1(2:year)) mean(v2(2:year)) mean(v3(2:year)) mean(value(2:year))]

