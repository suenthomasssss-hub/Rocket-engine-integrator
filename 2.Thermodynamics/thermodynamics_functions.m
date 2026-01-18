% 主函数：整合热力学计算功能
% 输入：
%   mode - 计算模式 ('isentropic_ratios', 'velocity_mach', 'ratios_vs_mach')
%   varargin - 根据模式传递的参数（支持结构体 params 或单独参数）
% 输出：
%   varargout - 根据模式返回不同的结果，包括喉部参数

function varargout = thermodynamics_functions(mode, varargin)
    % 根据 mode 提取参数
    switch mode
        case 'isentropic_ratios'
            % 提取参数
            if isstruct(varargin{1})
                params = varargin{1};
                p_c = params.p_c;
                p_e = params.p_e;
                gamma = params.gamma;
            else
                p_c = varargin{1};
                p_e = varargin{2};
                gamma = varargin{3};
            end
            [p_ratio, T_ratio, rho_ratio, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_isentropic_ratios(p_c, p_e, gamma);
            varargout = {p_ratio, T_ratio, rho_ratio, p_t_ratio, T_t_ratio, rho_t_ratio};
            
        case 'velocity_mach'
            % 提取参数
            if isstruct(varargin{1})
                params = varargin{1};
                p_ratio = varargin{2};
                T_c = params.T_c;
                gamma = params.gamma;
                R_M = params.R_M;
            else
                p_ratio = varargin{1};
                T_c = varargin{2};
                gamma = varargin{3};
                R_M = varargin{4};
            end
            [ve, Me, Te, a_e, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_velocity_mach(p_ratio, T_c, gamma, R_M);
            varargout = {ve, Me, Te, a_e, p_t_ratio, T_t_ratio, rho_t_ratio};
            
        case 'ratios_vs_mach'
            % 提取参数
            M = varargin{1};
            if isstruct(varargin{2})
                params = varargin{2};
                gamma = params.gamma;
            else
                gamma = varargin{2};
            end
            [p_ratios, T_ratios, rho_ratios, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_ratios_vs_mach(M, gamma);
            varargout = {p_ratios, T_ratios, rho_ratios, p_t_ratio, T_t_ratio, rho_t_ratio};
            
        otherwise
            error('未知模式：%s', mode);
    end
end

% 子函数：计算压强比、温度比、密度比及喉部参数
function [p_ratio, T_ratio, rho_ratio, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_isentropic_ratios(p_c, p_e, gamma)
    % 计算压强比
    p_ratio = p_c / p_e;
    
    % 计算温度比
    T_ratio = (p_ratio)^((gamma - 1) / gamma);
    
    % 计算密度比
    rho_ratio = (p_ratio)^(1 / gamma);
    
    % 计算喉部参数
    p_t_ratio = (2 / (gamma + 1))^(gamma / (gamma - 1)); % p_t / p_c
    T_t_ratio = 2 / (gamma + 1); % T_t / T_c
    rho_t_ratio = (2 / (gamma + 1))^(1 / (gamma - 1)); % rho_t / rho_c
end

% 子函数：计算出口速度、出口马赫数、出口温度、出口声速及喉部参数
function [ve, Me, Te, a_e, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_velocity_mach(p_ratio, T_c, gamma, R_M)
    % 计算温度比
    T_ratio = (p_ratio)^((gamma - 1) / gamma);
    
    % 计算出口温度
    Te = T_c / T_ratio;
    
    % 计算出口速度
    ve = sqrt((2 * gamma / (gamma - 1)) * R_M * T_c * (1 - (1/p_ratio)^((gamma - 1) / gamma)));
    
    % 计算出口马赫数
    Me = sqrt((2 / (gamma - 1)) * ((p_ratio)^((gamma - 1) / gamma) - 1));
    
    % 计算出口声速
    a_e = sqrt(gamma * R_M * Te); % 出口声速
    
    % 验证马赫数：Me = ve / a_e
    Me_check = ve / a_e;
    fprintf('验证：通过马赫数定义计算的 Me = %.4f\n', Me_check);
    
    % 计算喉部参数
    p_t_ratio = (2 / (gamma + 1))^(gamma / (gamma - 1)); % p_t / p_c
    T_t_ratio = 2 / (gamma + 1); % T_t / T_c
    rho_t_ratio = (2 / (gamma + 1))^(1 / (gamma - 1)); % rho_t / rho_c
end

% 子函数：计算压强比、温度比、密度比随马赫数的变化及喉部参数
function [p_ratios, T_ratios, rho_ratios, p_t_ratio, T_t_ratio, rho_t_ratio] = calculate_ratios_vs_mach(M, gamma)
    % 初始化数组
    p_ratios = zeros(size(M));
    T_ratios = zeros(size(M));
    rho_ratios = zeros(size(M));
    
    % 计算每个马赫数对应的比值
    for i = 1:length(M)
        factor = 1 + ((gamma - 1) / 2) * M(i)^2;
        p_ratios(i) = factor^(-gamma / (gamma - 1));
        T_ratios(i) = 1 / factor;
        rho_ratios(i) = factor^(-1 / (gamma - 1));              

    end
    
    % 计算喉部参数
    p_t_ratio = (2 / (gamma + 1))^(gamma / (gamma - 1)); % p_t / p_c
    T_t_ratio = 2 / (gamma + 1); % T_t / T_c
    rho_t_ratio = (2 / (gamma + 1))^(1 / (gamma - 1)); % rho_t / rho_c
end