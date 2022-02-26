function [B, P, C, S, tree, tree_num] = Forest(type, SI, k, Q, size, w, year, grow_type)
    % 本函数是第一问的森林模型函数的基础版本 
    % ------------------------
    % 输入变量
    % 树种 type
    % 立地条件 SI 归一化后
    % 采伐频率 k
    % 采伐强度 Q
    % 微分生长率选项，默认为分段 grow_type=0
    % 森林面积 size 默认为矩形
    % 长周期生产品占比系数 w 一般为 3/20
    % year 迭代年数
    % ------------------------
    % 输出变量
    % B Biomass 森林总生物量时间序列
    % P Product 生产品生物量时间序列
    % C Carbon Sequestration 系统碳固存总量时间序列
    % S Soil 土壤碳量时间序列
    % tree 森林树木矩阵，数值为年龄
    % tree_num 树龄统计时间序列

    % ------------------------
    % 根据树种、立地条件不同判断A0和alpha
    % 数据需要根据文献修改
    if type == 1 % A类树
        A0 = 80; % 0岁的初始生物量 目前以 kg 考虑
        c = [25, 50]; % c1 c2 树木成熟期和老龄期时间指标
        if grow_type == 0
            alpha = [0.033 0.017] * SI; % 生长率离散（可修改为函数）
        else
            alpha = [0 0] * SI; % 可修改，生长率连续
        end
    end

    % ------------------------
    % 初始化森林矩阵，初始年龄
    tree = rand(sqrt(size)); % 二维森林矩阵
    first_year = 10; 
    tree = round(first_year * tree) + 1; % 随机分布在 1 到 first_year + 1 之间

    % ------------------------
    % 初始化森林生物量指标 Biomass
    tree_num = zeros(3, year); % 统计每年三种树龄对应树的数目
    B = zeros(1, year); % 森林总生物量

    T = zeros(1, year * 2); % 一棵树生物量和树龄关系的树生长表，时间长度取 year 的两倍
    T(1) = A0; % 初始化生长表
    for t = 1 : year * 2
        if t <= c(1)
            T(t+1) = T(t)*(1 + alpha(1));
        elseif t <= c(2)
            T(t+1) = T(t)*(1 + alpha(2));
        else
            T(t+1) = T(t); % 维持不变
        end
    end
    for s = 1: size
        B(1) = B(1) + T(tree(s)); % 初始化第 0 年森林总生物量
    end
    % ------------------------
    % 初始化生产品指标，采伐指标 Produce Update
    beta = [0.02 0.3]; % 长短周期生产品损耗率
    P = zeros(1, year); % 统计每年产品总生物量
    
    % cut_num = ceil(year/k); % 采伐次数 = 总年数 / 采伐频率
    p = zeros(2, size); % 采伐变量表，记录每次采伐的生物量（第一行）和间隔tau（第二行）
    TP = zeros(2, year * 2); % 生产品使用衰减表，横坐标=2是长短周期两种生产品
    TP(:, 1) = 1;
    for t = 1 : year * 2 % 初始化衰减表TP
        TP(1, t+1) = TP(1, t)*(1 - beta(1));
        TP(2, t+1) = TP(2, t)*(1 - beta(2));
    end
    TP(TP < 0.01) = 0; % 过小则视为完全排放
    
    % ------------------------
    % 初始化土壤参数
    bound = year * 2;
    Sa = 2600 * sqrt(size) * SI; % deep > 10cm
    Sb = 6200 * sqrt(size) * SI; % deep < 10cm
    s_t = [5 30]; % 两个转折时间点
    
    % ------------------------
    % 计算土壤含碳量
    % 分为 10cm 以上和 10cm 以下考虑
    gamma1 = zeros(1,bound); % 浅于10变化率
    gamma2 = zeros(1,bound); % 深于10变化率
    gamma1(1:s_t(2)) = 0.0086;
    gamma2(1:s_t(1)) = -0.0346;
    for ii = s_t(2) + 1 : bound
        gamma1(ii) = gamma1(ii-1)-0.0001;
        if gamma1(ii) <= 0
            break;
        end
    end
    for ii = s_t(1) + 1 : bound
        gamma2(ii) = gamma2(ii - 1) + 0.001;
        if gamma2(ii)>=0
            break;
        end
    end
    soil = zeros(2, bound); % 10cm为界
    soil(1, 1) = Sa;
    soil(2, 1) = Sb;
    for ii = 1:(bound - 1)
        soil(1, ii + 1)= soil(1, ii) * (1 + gamma2(ii));
        soil(2, ii + 1)= soil(2, ii) * (1 + gamma1(ii));
    end
    S = soil(1, :) + soil(2, :);
    S = S(1:year);
    
    % ------------------------
    % 蒙特卡罗方法，迭代 year 次
    q = 1;
    for ii = 1 : year
        if mod(ii, k) == 0 % 每 k 年采伐一次
            count = sum(sum(tree > c(2))); % 大龄树木的数量
            mark = 0; % 已采伐数量  
            r = randperm(size); % 随机选取位置采伐（可修改，更优空间方案）

            jj = 1;
            while jj <= size
                if mark == floor(count * Q) % 到达采伐强度要求
                    break; % 停止采伐
                end
                if  tree(r(jj)) > 50 % 符合采伐年龄
                    p(1, q) = p(1, q) + T(tree(r(jj)));
                    p(2, q) = 1;
                    tree(r(jj)) = 0; % 采伐后位置为新生树
                    mark = mark + 1; % 已采伐增加
                end
                jj = jj + 1;
            end
            q = q + 1;
        end
        tree = tree + 1; % 年龄所有增加

        tree_num(3, ii) = sum(sum(tree > c(2)));
        tree_num(2, ii) = sum(sum(tree <= c(2) & tree > c(1)));
        tree_num(1, ii) = sum(sum(tree <= c(1)));

        % 更新森林、生产品生物量指标
        pp = p(2, :) ~= 0;
        p(2, :) = p(2, :) + pp;
        for s = 1: size
            if p(2, s) ~= 0
                P(ii) = P(ii) + p(1, s) * (1-w) * TP(1, p(2, s)) + p(1, s) * w * TP(2, p(2, s));              
            end
            B(ii) = B(ii) + T(tree(s));
        end  
    end

    % ------------------------
    % 最终模型（暂未考虑 Soil 土壤）
    C = B + P + S;

    % ------------------------
    % 可视化
    x = 1 : year;
    plot(x, P, x, C, x, B, x, S); % 三指标时序曲线 
end



