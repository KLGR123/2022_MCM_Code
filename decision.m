function [k_new, Q_new] = decision(tree, tree_num, rain, humid, wind, tree_type, w, weight, forcast_depth, size, c, shape, cut_type, p, P, T, TP, B, year)
    % 决策器
    % ------------------------
    % 深度是指基于一个 plan 往后看 depth 个深度

    value_max = 0;
    k_new = 0;
    Q_new = 0;
    for k = [15 27 50]
        for Q = [0.33 0.5]
            [tree_new, tree_num_new, tree_type_new, Bdiff, Pdiff] = forest_step(k, Q, tree, tree_num, tree_type, w, size, c, shape, cut_type, p, P, T, TP, B, forcast_depth, year);
            [~,~,~,value] = value_score(Bdiff, tree_new, tree_num_new, rain, humid, wind, Pdiff, tree_type_new, 1-w, weight);
            if value > value_max % 目标函数value
                k_new = k;
                Q_new = Q;
                value_max = value;
            end
        end
    end
end

function [tree_new, tree_num_new, tree_type_new, Bdiff, Pdiff] = forest_step(k, Q, tree, tree_num, tree_type, w, size, c, shape, cut_type, p, P, T, TP, B, forcast_depth, year)
    % 只迭代一定长度的森林函数
    
    Q = [Q 0]; % 只考虑纯种林，有待改进
    q = 1; 
    fan = 1; 
    for ii = year : (k * forcast_depth + year)
        if mod(ii, k) == 0 % 判断到达采伐年
            % ------------------------
            % 首先进行采伐活动 更新各个数组
            mark = 0; % 已采伐数量  
            count1 = 0;
            count2 = 0;
            for s = 1 : size
                count1 = count1 + (tree(s) > c(tree_type(s)+1, 2)); % 大龄树木的数量
                count2 = count2 + (tree(s) <= c(tree_type(s)+1, 1)); % 幼龄树木的数量
            end
            
            if shape == 0
                r = randperm(size); % 随机选取位置采伐（可修改，更优空间方案）
                resize = size;
            else 
                resize = size / 4;
                switch fan % 对应空间上四个不同区域
                    case 1
                        r = randperm(size / 4);
                        fan = 2;
                    case 2
                        r = randperm(size / 4) + size / 4;
                        fan = 3;
                    case 3
                        r = randperm(size / 4) + size / 2;
                        fan = 4;
                    case 4
                        r = randperm(size / 4) + 3 * size / 4;
                        fan = 1;
                end
            end
    
            jj = 1;
            while jj <= resize
                if mark == floor(count1 * Q(1) + count2 * Q(2)) % 到达采伐强度要求
                    break; % 停止采伐
                end
                if  tree(r(jj)) > c(tree_type(jj)+1, 2) || (tree(r(jj)) < c(tree_type(jj)+1, 1) && cut_type == 1)
                    p(1, q) = p(1, q) + T(tree_type(jj)+1, tree(r(jj)));
                    p(2, q) = 1;
                    tree(r(jj)) = 0; % 采伐后位置为新生树
                    mark = mark + 1; % 已采伐增加
                end  
                jj = jj + 1;
            end
            q = q + 1; % 采伐次数增加          
        end
        
        tree = tree + 1; % 年龄所有增加
    
        for ss = 1 : size
            tree_num(3, ii) = tree_num(3, ii) + (tree(ss) > c(tree_type(ss)+1, 2)); % 大龄树木的数量
            tree_num(1, ii) = tree_num(1, ii) + (tree(ss) <= c(tree_type(ss)+1, 1)); % 幼龄树木的数量
        end
        tree_num(2, ii) = size - tree_num(1, ii) - tree_num(3, ii);
    
        % 更新森林、生产品生物量指标
        pp = p(2, :) ~= 0;
        p(2, :) = p(2, :) + pp;
        for s = 1: size
            if p(2, s) ~= 0
                P(ii) = P(ii) + p(1, s) * (1-w) * TP(1, p(2, s)) + p(1, s) * w * TP(2, p(2, s));              
            end
            B(ii) = B(ii) + T(tree_type(s)+1, tree(s));
        end 
    end
    tree_new = tree;
    tree_num_new = tree_num; 
    tree_type_new = tree_type;
    Bdiff = B(k * forcast_depth + year) - B(year);
    Pdiff = P(k * forcast_depth + year) - P(year);
end