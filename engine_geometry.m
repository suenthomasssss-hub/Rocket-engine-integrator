% 20吨级分级燃烧引擎基本火箭公式计算（含效率公式）
% 目标推力：20吨 = 196,200 N
% 推进剂：液氧/煤油，初始假设 Isp = 320 s

% 输入参数
g0 = 9.80665; % 重力加速度 (m/s²)
F_target = 20 * 1000 * g0; % 目标推力 (N)
Isp = 320; % 比冲 (s)
pa = 101325; % 大气压力 (Pa)
p_c = 100*101325;
p_e = 101325; % 出口压力 (Pa，假设完全膨胀)
epsilon = 10; % 喷管扩张比
rho_c = 100; % 燃烧室密度 (kg/m³，经验值)
v_c = 50; % 燃烧室流速 (m/s，经验值)
delta_t = 160; % 燃烧时间 (s，典型任务)
Q = 1.15e7; % 推进剂化学能 (J/kg，液氧/煤油)
v_flight = 0; % 飞行速度 (m/s，地面测试，任务中可调整)

% 1. 基本火箭公式
ve = Isp * g0; % 出口速度 (m/s)
m_dot = F_target / ve; % 质量流率 (kg/s，假设 pe = pa)
Ae = 0.01; % 初始出口面积 (m²)
F = m_dot * ve + (p_e - pa) * Ae; % 实际推力
Isp_calculated = ve / g0; % 验证比冲
m_prop = m_dot * delta_t; % 推进剂总质量 (kg)
I_total = F * delta_t; % 总冲 (N·s)

% 2. 初步几何尺寸
Ac = m_dot / (rho_c * v_c); % 燃烧室截面积 (m²)
Dc = sqrt(4 * Ac / pi); % 燃烧室直径 (m)
At = Ac / 5; % 喉部面积 (m²，经验值)
Dt = sqrt(4 * At / pi); % 喉部直径 (m)
Ae = epsilon * At; % 出口面积 (m²)
De = sqrt(4 * Ae / pi); % 出口直径 (m)

% 3. 新增公式：推进剂化学能到动能的转换效率
P_int = 0.5 * m_dot * ve^2; % 喷气动能功率 (W)
P_chem = m_dot * Q; % 化学能功率 (W)
eta_energy = P_int / P_chem; % 转换效率

% 4. 新增公式：引擎功率
P_engine = F * v_flight; % 引擎功率 (W)

% 5. 新增公式：推进效率
eta_p = (2 * (v_flight / ve)) / (1 + (v_flight / ve)); % 推进效率

% 输出结果
fprintf('=== 20吨级引擎基本计算结果（含效率） ===\n');
fprintf('推力: %.2f N\n', F);
fprintf('比冲: %.2f s\n', Isp_calculated);
fprintf('总冲: %.2e N·s\n', I_total);
fprintf('质量流率: %.2f kg/s\n', m_dot);
fprintf('推进剂总质量: %.2f kg\n', m_prop);
fprintf('燃烧室直径: %.4f m\n', Dc);
fprintf('喉部直径: %.4f m\n', Dt);
fprintf('出口直径: %.4f m\n', De);
fprintf('化学能到动能转换效率: %.2f%%\n', eta_energy * 100);
fprintf('引擎功率: %.2e W (飞行速度 %.1f m/s)\n', P_engine, v_flight);
fprintf('推进效率: %.2f%%\n', eta_p * 100);

% 可视化初步几何
x = [0, 0.5, 0.7]; % 燃烧室、喉部、出口的轴向位置 (m)
D = [Dc, Dt, De]; % 对应直径
figure;
subplot(2, 1, 1);
plot(x, D, 'b-o', 'LineWidth', 2); hold on;
plot(x, -D, 'b-o', 'LineWidth', 2); % 对称喷管
xlabel('轴向位置 (m)'); ylabel('直径 (m)');
title('初步推力室几何');
grid on;

% 可视化推进效率随飞行速度的变化
v_flight_range = linspace(0, ve, 100);
eta_p_range = (2 * (v_flight_range / ve)) ./ (1 + (v_flight_range / ve));
subplot(2, 1, 2);
plot(v_flight_range, eta_p_range * 100, 'r-', 'LineWidth', 2);
xlabel('飞行速度 (m/s)'); ylabel('推进效率 (%)');
title('推进效率随飞行速度变化');
grid on;

% 保存结果，为后续模块提供输入
save('rocket_params.mat', 'F', 'Isp', 'm_dot', 'Dc', 'Dt', 'De', 'eta_energy', 'eta_p','p_c','p_e');