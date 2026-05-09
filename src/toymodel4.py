import numpy as np
import matplotlib.pyplot as plt

def generate_masses(n, m_min=500, m_max=1e5):
    """逆变换采样生成 dN/dM ~ M^-2 的分子云质量分布"""
    u = np.random.uniform(0, 1, n)
    return 1.0 / ((1.0 - u) / m_min + u / m_max)

def weighted_binned_statistic(x, z, w, bins):
    """计算质量加权的区间平均值 Z"""
    sum_wz, _ = np.histogram(x, bins=bins, weights=z * w)
    sum_w, _ = np.histogram(x, bins=bins, weights=w)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_z = np.where(sum_w > 0, sum_wz / sum_w, np.nan)
    return 0.5 * (bins[:-1] + bins[1:]), mean_z

def simulate_hunter_systematic_forward_model(n_clouds=30000):
    R_sun = 8.15
    V_0 = 236.0
    np.random.seed(42)

    # =======================================================
    # 1. 生成真实空间分布 (Outer disk: R = 9~16 kpc, 聚焦二三象限)
    # =======================================================
    l_true_deg = np.random.uniform(90, 160, n_clouds)
    l_true = np.radians(l_true_deg)
    d_true = np.random.uniform(2.0, 14.0, n_clouds)
    R_true = np.sqrt(R_sun**2 + d_true**2 - 2 * R_sun * d_true * np.cos(l_true))
    
    # 限制真实 R 在 9.0 到 16.0 之间
    valid_R = (R_true >= 9.0) & (R_true <= 16.0)
    l_true, d_true, R_true, l_true_deg = l_true[valid_R], d_true[valid_R], R_true[valid_R], l_true_deg[valid_R]
    n_valid = len(R_true)

    # 生成真实的连续质量分布
    masses = generate_masses(n_valid)

    # =======================================================
    # 2. 生成真实的物理形态 (波浪盘 vs 平坦盘)
    # =======================================================
    # 设定真实波浪：振幅 200 pc, 波长 4.5 kpc
    wave_amplitude = 0.200 
    wave_length = 4.5
    Z_true_wave = wave_amplitude * np.sin(2 * np.pi * (R_true - 9.0) / wave_length)
    Z_true_flat = np.zeros(n_valid)
    
    # 引入本征的气体盘厚度 (Flaring: 随 R 变厚，且质量越大的云越沉淀在中心)
    intrinsic_disp = 0.04 + 0.02 * (R_true - 9.0) 
    thickness_noise = np.random.normal(0, intrinsic_disp * (500/masses)**0.1)
    
    Z_true_wave += thickness_noise
    Z_true_flat += thickness_noise

    # 计算真实的银纬 b
    b_true_wave = np.arcsin(Z_true_wave / np.sqrt(d_true**2 + Z_true_wave**2))
    b_true_flat = np.arcsin(Z_true_flat / np.sqrt(d_true**2 + Z_true_flat**2))

    # =======================================================
    # 3. 核心物理引擎：注入 Hunter 2024 的系统性交替非圆运动
    # =======================================================
    # 真实的纯圆周运动视向速度
    V_true_wave = V_0 * (R_sun / R_true - 1.0) * np.sin(l_true) * np.cos(b_true_wave)
    V_true_flat = V_0 * (R_sun / R_true - 1.0) * np.sin(l_true) * np.cos(b_true_flat)

    # [大杀器] 模拟 Hunter 2024 中绿色和粉紫色交替的非圆运动条带！
    # 假设流体运动在径向和方位角上呈现周期性交替，振幅高达 12 km/s
    V_systematic = 12.0 * np.sin(2 * np.pi * (R_true - 9.0) / 3.0) + \
                    5.0 * np.cos(2 * np.pi * (l_true_deg - 90) / 40.0)
    
    # 加上气体本身自带的随机湍流弥散 (5 km/s)
    V_random = np.random.normal(0, 5.0, n_valid)
    
    # 最终望远镜“观测”到的速度 = 真实速度 + 系统交替偏差 + 随机噪声
    V_obs_wave = V_true_wave + V_systematic + V_random
    V_obs_flat = V_true_flat + V_systematic + V_random

    # =======================================================
    # 4. 逆推运动学距离 (几何联动：V -> R, d, Z)
    # =======================================================
    def solve_kinematic(V_obs, l_rad, b_rad):
        # 视向速度反解公式
        term = V_obs / (V_0 * np.sin(l_rad) * np.cos(b_rad)) + 1.0
        R_obs = np.zeros_like(V_obs) * np.nan
        valid = term > 0
        R_obs[valid] = R_sun / term[valid]
        
        # 外盘几何解
        valid = valid & (R_obs > R_sun * np.abs(np.sin(l_rad)))
        d_obs = np.zeros_like(R_obs) * np.nan
        d_obs[valid] = R_sun * np.cos(l_rad[valid]) + np.sqrt(R_obs[valid]**2 - (R_sun * np.sin(l_rad[valid]))**2)
        
        # 重新计算投影高度 Z
        Z_obs = d_obs * np.tan(b_rad)
        return R_obs, Z_obs

    R_obs_wave, Z_obs_wave = solve_kinematic(V_obs_wave, l_true, b_true_wave)
    R_obs_flat, Z_obs_flat = solve_kinematic(V_obs_flat, l_true, b_true_flat)

    # =======================================================
    # 5. 绘图：准备封杀审稿人的 2x2 面板
    # =======================================================
    fig, axs = plt.subplots(2, 2, figsize=(15, 11))
    R_bins = np.arange(9.5, 15.5, 0.4)

    def plot_panel(ax, R_data, Z_data, M_data, title, color_map):
        # 过滤无效值
        mask = ~np.isnan(R_data) & ~np.isnan(Z_data) & (R_data > 9) & (R_data < 16)
        R_c, Z_c, M_c = R_data[mask], Z_data[mask], M_data[mask]
        
        # 绘制底部散点 (大量低质量云被距离误差拉扯，形成 Noisy 的视觉效果)
        ax.scatter(R_c, Z_c * 1000, s=3, c='gray', alpha=0.15, label='Individual Clouds')
        
        # 绘制质量加权平均线 (大质量云稳住阵脚)
        bin_centers, mean_Z = weighted_binned_statistic(R_c, Z_c, M_c, R_bins)
        ax.plot(bin_centers, mean_Z * 1000, color='red', lw=3.5, marker='o', 
                markersize=8, markeredgecolor='white', label='Mass-weighted Mean')
        
        ax.axhline(0, color='black', lw=1.5, ls='--')
        ax.set_title(title, fontsize=15, pad=10)
        ax.set_ylim(-600, 600)
        ax.set_xlim(9, 16)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=11, loc='upper right')

    # 第一行：真实的波浪盘能否抗住系统性误差？
    plot_panel(axs[0, 0], R_true, Z_true_wave, masses, 
               'Model A: True Corrugated Disk\n($A=200$ pc, $\lambda=4.5$ kpc)', 'Blues')
    plot_panel(axs[0, 1], R_obs_wave, Z_obs_wave, masses, 
               f'Model A Observed:\nWith Hunter 2024 Systematic Streaming ($12$ km/s)', 'Blues')
    
    # 第二行：平坦盘能否被系统性误差捏造出波浪？
    plot_panel(axs[1, 0], R_true, Z_true_flat, masses, 
               'Model B: True Flat Disk\n(Null Hypothesis, $Z=0$)', 'Greys')
    plot_panel(axs[1, 1], R_obs_flat, Z_obs_flat, masses, 
               f'Model B Observed:\nWith Hunter 2024 Systematic Streaming ($12$ km/s)', 'Greys')

    for ax in axs[:, 0]: ax.set_ylabel("Vertical Height $Z$ (pc)", fontsize=14)
    for ax in axs[1, :]: ax.set_xlabel("Galactocentric Radius $R$ (kpc)", fontsize=14)
    
    plt.tight_layout()
    plt.savefig('Rigorous_Hunter2024_Forward_Model.pdf', dpi=300)
    plt.show()

# 运行终极模拟
simulate_hunter_systematic_forward_model(n_clouds=30000)
