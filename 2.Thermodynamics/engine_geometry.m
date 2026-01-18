% 20吨级分级燃烧引擎基本火箭公式计算（含效率公式）
% 目标推力：20吨 = 196,200 N
% 推进剂：液氧/煤油，初始假设 Isp = 320 s
% 整合热力学模块结果

% 输入参数
g0 = 9.80665; % 重力加速度 (m/s²)
F_target = 20 * 1000 * g0; % 目标推力 (N)
Isp = 320; % 比冲 (s)
pa = 101325; % 大气压力 (Pa)

% 热力学参数（假设从 thermodynamic 模块加载）
load('rocket_params_module2.mat', 'p_c', 'p_e', 'gamma', 'T_c','ve'); % 加载热力学结果
if exist('p_c_thermo', 'var') && exist('p_e_thermo', 'var')
    p_c = p_c_thermo; % 更新燃烧室压力
    p_e = p_e_thermo; % 更新出口压力
end

gamma = 1.4; % 比热比（若未加载则使用默认值）
T_c = 3300; % 总温 (K)，基于你提到的燃烧室温度

% 1. 基本火箭公式
m_dot = F_target / ve; % 质量流率 (kg/s，假设 pe ≈ pa 初始，后调整)
Ae = m_dot * ve / (p_c * (2 / (gamma + 1)) ^ ((gamma + 1) / (2 * (gamma - 1)))); % 更精确出口面积
F = m_dot * ve + (p_e - pa) * Ae; % 实际推力
Isp_calculated = ve / g0; % 验证比冲
delta_t = 160; % 燃烧时间 (s)
m_prop = m_dot * delta_t; % 推进剂总质量 (kg)
I_total = F * delta_t; % 总冲 (N·s)

% 3. 新增公式：推进剂化学能到动能的转换效率
P_int = 0.5 * m_dot * ve^2; % 喷气动能功率 (W)
Q = 1.15e7; % 推进剂化学能 (J/kg)
P_chem = m_dot * Q; % 化学能功率 (W)
eta_energy = P_int / P_chem; % 转换效率

% 4. 新增公式：引擎功率
v_flight = 0; % 飞行速度 (m/s，地面测试)
P_engine = F * v_flight; % 引擎功率 (W)

% 5. 新增公式：推进效率
eta_p = (2 * (v_flight / ve)) / (1 + (v_flight / ve)); % 推进效率

% 输出结果
fprintf('=== 20吨级引擎基本计算结果（含热力学调整） ===\n');
fprintf('推力: %.2f N\n', F);
fprintf('比冲: %.2f s\n', Isp_calculated);
fprintf('总冲: %.2e N·s\n', I_total);
fprintf('质量流率: %.2f kg/s\n', m_dot);
fprintf('推进剂总质量: %.2f kg\n', m_prop);
fprintf('燃烧室直径: %.4f m\n', Dc);
fprintf('喉部直径: %.4f m\n', Dt);
fprintf('化学能到动能转换效率: %.2f%%\n', eta_energy * 100);
fprintf('引擎功率: %.2e W (飞行速度 %.1f m/s)\n', P_engine, v_flight);
fprintf('推进效率: %.2f%%\n', eta_p * 100);

% 保存结果，为后续模块提供输入
save('rocket_params.mat', 'F', 'Isp', 'm_dot', 'Dc', 'Dt', 'eta_energy', 'eta_p', 'p_c', 'p_e');