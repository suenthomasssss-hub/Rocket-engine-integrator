% 主脚本：模块2 - 热力学计算
% 加载模块1结果文件，调用 thermodynamics_functions.m 计算

% 确保 thermodynamics_functions.m 在路径中（根据需要取消注释）
% addpath('path_to_your_functions_directory'); % 指向包含 thermodynamics_functions.m 的目录

% 加载模块1结果文件
if isfile('rocket_params.mat')
    load('rocket_params.mat', 'p_c', 'p_e');
else
    p_c = 100 * 101325; % 默认燃烧室压力 (Pa，约100 atm)
    p_e = 101325; % 默认出口压力 (Pa，约1 atm)
    fprintf('未找到 rocket_params.mat，使用默认值：p_c = %.2e Pa, p_e = %.2e Pa\n', p_c, p_e);
end

% 比热比和其他参数（假设模块1未保存，手动指定）
gamma = 1.4; % 比热比（理想气体假设，LOX/LCH4燃烧产物典型值）
T_c = 3300; % 燃烧室温度 (K)
R_M = 287; % 气体常数 R/M (J/(kg·K))

% 调用 thermodynamics_functions.m 计算
% 1. 计算压强比、温度比、密度比
[p_ratio, T_ratio, rho_ratio] = thermodynamics_functions('isentropic_ratios', p_c, p_e, gamma);
% 2. 计算出口速度、出口马赫数、出口温度、出口声速
[ve, Me, Te, a_e, p_t_ratio, T_t_ratio, rho_t_ratio] = thermodynamics_functions('velocity_mach', p_ratio, T_c, gamma, R_M);

% 输出结果
fprintf('=== 模块2：等熵流动计算结果 ===\n');
fprintf('压强比 (p_c / p_e): %.4f\n', p_ratio);
fprintf('温度比 (T_c / T_e): %.4f\n', T_ratio); % 修正为 T_e / T_c
fprintf('密度比 (rho_c / rho_e): %.4f\n', rho_ratio); % 修正为 rho_e / rho_c
fprintf('出口速度 (m/s): %.4f\n', ve);
fprintf('出口声速 (m/s): %.4f\n', a_e);
fprintf('出口温度 (K): %.4f\n', Te);

fprintf('===喉部参数===\n')
fprintf('喉部压强比: %.6f\n',p_t_ratio)
fprintf('喉部温度比: %.6f\n',T_t_ratio)
fprintf('喉部密度比: %.6f\n',rho_t_ratio)
% 马赫数范围：0 到 Me
M = linspace(0, Me, Me*100);

% 计算压强比、温度比、密度比随马赫数的变化
[p_ratios, T_ratios, rho_ratios] = thermodynamics_functions('ratios_vs_mach', M, gamma);

% temperature and ratio change
p = p_c * p_ratios;
T = T_c * T_ratios;
a = sqrt(gamma*R_M*T);

% 绘制曲线
subplot(3,1,1)
plot(M,p,'b-', 'LineWidth', 2, 'DisplayName','pressure');
xlabel('Mach Number (M)');
ylabel('pressure (Pa)');
hold on;
grid on;
subplot(3,1,2)
plot(M,T,'r-', 'LineWidth', 2, 'DisplayName','temperature');
xlabel('Mach Number (M)');
ylabel('temperature (K)');
hold on;
grid on;
subplot(3,1,3)
plot(M,a,'g-.', 'LineWidth', 2, 'DisplayName','sonic_speed')
xlabel('Mach Number (M)');
ylabel('sonic_speed (ms-1)')
grid on;
figure;
plot(M, p_ratios, 'b-', 'LineWidth', 2, 'DisplayName', 'Pressure Ratio (p/p_c)');
hold on
plot(M, T_ratios, 'r-', 'LineWidth', 2, 'DisplayName', 'Temperature Ratio (T/T_c)');
hold on
plot(M, rho_ratios, 'g-.', 'LineWidth', 2, 'DisplayName', 'Density Ratio (\rho/\rho_c)');
xlabel('Mach Number (M)');
ylabel('Ratios');
title('Isentropic Flow: Pressure, Temperature, and Density Ratios vs Mach Number');
legend('Location', 'best');
grid on;

save_dir = pwd; % 当前工作目录
save_file = fullfile(save_dir, ['rocket_params_module2', '.mat']);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('创建保存目录：%s\n', save_dir);
end
save(save_file, 'params', 'a','p','T','ve','a_e','Te','p_t_ratio','T_t_ratio','rho_t_ratio');
fprintf('模块2 结果已保存到：%s\n', save_file);