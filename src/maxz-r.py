import numpy as np
import matplotlib.pyplot as plt

# Define theta grid for scanning (same as HI and CO)
THETA_DEG = np.linspace(-30, 170, 200)

def warp_cheng_max_over_theta(R, R1, R2, R_w, alpha):
    """
    Calculate maximum Cheng2020 warp amplitude by scanning over theta
    """
    h = np.zeros_like(R)
    
    mask1 = R <= R1
    mask2 = (R > R1) & (R <= R2)
    mask3 = R > R2
    
    h[mask1] = 0
    h[mask2] = (R_w / (R2 - R1)) * (R[mask2] - R1)**alpha
    
    h_R2 = (R_w / (R2 - R1)) * (R2 - R1)**alpha
    dh_dR_at_R2 = (R_w * alpha / (R2 - R1)) * (R2 - R1)**(alpha - 1)
    h[mask3] = h_R2 + dh_dR_at_R2 * (R[mask3] - R2)
    
    return h

def warp_poggio_max_over_theta(R, R_w, h_w0, alpha_w, phi_w=0):
    """
    Calculate maximum Poggio2025 warp amplitude by scanning over theta
    """
    R = np.asarray(R)
    h_max = np.zeros_like(R, dtype=float)
    mask = R > R_w
    if np.any(mask):
        amp = h_w0 * (R[mask] - R_w)**alpha_w
        s = np.sin(np.deg2rad(THETA_DEG)[None, :] - phi_w)
        h_vals = amp[:, None] * s
        h_max[mask] = np.max(np.abs(h_vals), axis=1)
    return h_max

def chen2019_max_over_theta(R, R_w=7.72, a=0.060, b=1.33, phi_w=17.5):
    """
    Calculate maximum Chen2019 warp amplitude by scanning over theta
    """
    R = np.asarray(R)
    z_max = np.zeros_like(R, dtype=float)
    mask = R > R_w
    if np.any(mask):
        amp = a * (R[mask] - R_w)**b
        s = np.sin(np.deg2rad(THETA_DEG)[None, :] - np.radians(phi_w))
        z_vals = amp[:, None] * s
        z_max[mask] = np.max(np.abs(z_vals), axis=1)
    return z_max

def khanna2025_max_over_theta(R, R_w, h_w0, a_w, phi_w=178):
    """
    Calculate maximum Khanna2025 warp amplitude by scanning over theta
    """
    R = np.asarray(R)
    h_max = np.zeros_like(R, dtype=float)
    mask = R > R_w
    if np.any(mask):
        amp = h_w0 * (R[mask] - R_w)**a_w
        s = np.sin(np.deg2rad(THETA_DEG)[None, :] - np.radians(phi_w))
        h_vals = amp[:, None] * s
        h_max[mask] = np.max(np.abs(h_vals), axis=1)
    return h_max

def hi_warp_model(rg):
    """
    Calculate the maximum warp amplitude for a given Galactocentric radius (R_GC).
    """
    if rg < 10:
        return 0

    w1 = 9 + 197 * (rg - 10) - 3.1 * (rg - 10)**2
    w0 = -66 + 150 * (rg - 15) - 0.47 * (rg - 15)**2
    w2 = -70 + 171 * (rg - 15) - 5.3 * (rg - 15)**2

    theta = np.linspace(-30, 170, 200)

    if w2 >= 150 and rg > 15:
        w_val = w0 + w1 * np.sin(np.deg2rad(theta)) + w2 * np.sin(2 * np.deg2rad(theta))
    elif w2 < 150 and rg > 15:
        w_val = w0 + w1 * np.sin(np.deg2rad(theta))
    else:
        w_val = w1 * np.sin(np.deg2rad(theta))

    return np.max(np.abs(w_val)) / 1000

# CO m=1 parameters
a1    = 0.09363
Rw1   = 8.568
bw1   = 1.050
phi_w1 = -0.7050

def co_warp_max_over_theta(R, a1, Rw1, bw1, phi_w1, theta_deg=THETA_DEG):
    """
    Max CO m=1 warp amplitude at each R by scanning theta in [-30°, 170°].
    """
    R = np.asarray(R)
    zw_max = np.zeros_like(R, dtype=float)
    mask = R >= Rw1
    if np.any(mask):
        amp = a1 * (R[mask] - Rw1)**bw1
        s = np.sin(np.deg2rad(theta_deg)[None, :] - phi_w1)
        zw = amp[:, None] * s
        zw_max[mask] = np.max(np.abs(zw), axis=1)
    return zw_max

# CO m=1,2复合模型参数
a1_m12 = 0.1192
Rw1_m12 = 9.015
phi_w1_m12 = -3.538  # 弧度
a2_m12 = -0.01241
Rw2_m12 = 12.71
bw2_m12 = 2.031
phi_w2_m12 = -20.45  # 弧度

def co_m12_warp_max_over_theta(R, a1, Rw1, phi_w1, a2, Rw2, bw2, phi_w2, theta_deg=THETA_DEG):
    """
    Calculate maximum CO m=1,2复合翘曲模型振幅 by scanning over theta
    Using: Zw(R ≥ Rw1,≥ Rw2) = a1(R−Rw1)sin(φ −φw1)+a2(R−Rw2)^(bw2)sin(2(φ −φw2))
    """
    R = np.asarray(R)
    zw_max = np.zeros_like(R, dtype=float)
    
    # 将角度转换为弧度
    theta_rad = np.deg2rad(theta_deg)
    
    for i, r in enumerate(R):
        if r < Rw1:
            zw_max[i] = 0
        else:
            # 计算m=1部分
            if r >= Rw1:
                m1_amp = a1 * (r - Rw1)
                m1_vals = m1_amp * np.sin(theta_rad - phi_w1)
            else:
                m1_vals = np.zeros_like(theta_rad)
            
            # 计算m=2部分
            if r >= Rw2:
                m2_amp = a2 * (r - Rw2)**bw2
                m2_vals = m2_amp * np.sin(2 * (theta_rad - phi_w2))
            else:
                m2_vals = np.zeros_like(theta_rad)
            
            # 组合m=1和m=2
            zw_vals = m1_vals + m2_vals
            zw_max[i] = np.max(np.abs(zw_vals))
    
    return zw_max

# Model parameters
R1_cheng = 8.87
R2_cheng = 17.78
R_w_cheng = 1.20
alpha_cheng = 1.53

R_w_young = 5.5
h_w0_young = 0.012
alpha_w_young = 1.9
phi_w_young = 9.9

R_w_ceph = 7.7
h_w0_ceph = 0.057
alpha_w_ceph = 1.3
phi_w_ceph = 14.0

R_w_chen_power = 7.72
a_chen_power = 0.060
b_chen_power = 1.33
phi_w_chen = 17.5

R_w_rc3 = 8.79
h_w0_rc3 = 0.083
a_w_rc3 = 2.0
phi_w_rc3 = 178

# Create Galactocentric radius array - 调整为8-20 kpc
R = np.linspace(7, 20, 200)

# Calculate warp amplitudes
h_R_rc3 = khanna2025_max_over_theta(R, R_w_rc3, h_w0_rc3, a_w_rc3, phi_w_rc3)
h_R_cheng = warp_cheng_max_over_theta(R, R1_cheng, R2_cheng, R_w_cheng, alpha_cheng)
h_R_hi = np.array([hi_warp_model(r) for r in R])
h_R_young_giants = warp_poggio_max_over_theta(R, R_w_young, h_w0_young, alpha_w_young, np.radians(phi_w_young))
h_R_cepheids = warp_poggio_max_over_theta(R, R_w_ceph, h_w0_ceph, alpha_w_ceph, np.radians(phi_w_ceph))
h_R_chen_power = chen2019_max_over_theta(R, R_w_chen_power, a_chen_power, b_chen_power, phi_w_chen)
h_R_co_max = co_warp_max_over_theta(R, a1, Rw1, bw1, phi_w1)
# CO m=1,2复合模型
h_R_co_m12 = co_m12_warp_max_over_theta(R, a1_m12, Rw1_m12, phi_w1_m12, a2_m12, Rw2_m12, bw2_m12, phi_w2_m12)

# 创建优化的图形 - 缩小图形尺寸
fig, ax1 = plt.subplots(1, 1, figsize=(6.5, 4.5))

# =============================================================================
# 视觉优化方案：更新图例标签
# =============================================================================

# 1. 老恒星群体 - 使用暖色调和实线
ax1.plot(R, h_R_rc3, color='#8B4513', linewidth=2.5, label='Red clump (Khanna+2025)', 
         alpha=0.8, linestyle='-')
ax1.plot(R, h_R_cheng, color='#D2691E', linewidth=2.5, label='Gaia DR2/APOGEE (Cheng+2020)', 
         alpha=0.8, linestyle='-')

# 2. 年轻恒星群体 
ax1.plot(R, h_R_chen_power, color='#1E90FF', linewidth=2.5, label='Cepheids (Chen+2019)', 
         alpha=0.8, linestyle='-')  # 改为实线

ax1.plot(R, h_R_young_giants, color='#4682B4', linewidth=2.5, label='Young Giants (Poggio+2025)', 
         alpha=0.8, linestyle='--')

ax1.plot(R, h_R_cepheids, color='b', linewidth=2.5, label='Cepheids (Poggio+2025)', 
         alpha=0.8, linestyle='--')

# 3. 原子气体 - 使用绿色系和点划线
ax1.plot(R, h_R_hi, color='#2E8B57', linewidth=2.5, label='HI (Levine+2006)', 
         alpha=0.8, linestyle='-.')

# 4. 分子气体 - 特别突出显示（我们的结果）
ax1.plot(R, h_R_co_max, color='r', linewidth=3.5, label='CO m=1 - This work', 
         alpha=0.9, linestyle='-')
#ax1.plot(R, h_R_co_m12, color='#FF4500', linewidth=4.5, label='CO m=1,2 - This work', 
#         alpha=1.0, linestyle='-', marker='')

# 1. 添加单一浅灰色背景填充
ax1.fill_between(R, 0, 2.5, alpha=0.05, color='#f5f5f5', edgecolor='none')

# 设置标签
ax1.set_xlabel('R (kpc)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Maximum Z (kpc)', fontsize=11, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9, loc='upper left', framealpha=0.95, ncol=2)
ax1.set_ylim(-0.1, 2.5)
ax1.set_xlim(9.5, 16.6)  # 调整x轴范围

# 调整群体标注框位置
ax1.text(14, 1.3, 'OLD \nPOPULATIONS', fontsize=10, color='#8B4513', 
         ha='center', va='center', weight='bold',
         bbox=dict(boxstyle="round,pad=0.3", facecolor='#FAEBD7', alpha=0.8))

# 年轻恒星群体标注 
ax1.text(14, 0.2, 'YOUNG \nPOPULATIONS', fontsize=10, color='#1E90FF', 
         ha='center', va='center', weight='bold',
         bbox=dict(boxstyle="round,pad=0.3", facecolor='#F0F8FF', alpha=0.8))

plt.tight_layout()
ax1.text(-0.08, 1., "a", fontsize=11,transform=ax1.transAxes, va="top", ha="right",
         fontweight="bold")
plt.savefig("multi_warp_model.png", dpi=300, bbox_inches="tight")
# 打印分组统计信息，特别突出CO结果
print("=" * 100)
print("GALACTIC WARP AMPLITUDES BY POPULATION TYPE")
print("=" * 100)
print("\n--- OLD STELLAR POPULATIONS ---")
print(f"Red clump (Khanna+2025):      Max = {np.max(h_R_rc3):.3f} kpc at R = {R[np.argmax(h_R_rc3)]:.1f} kpc")
print(f"Gaia DR2/APOGEE (Cheng+2020): Max = {np.max(h_R_cheng):.3f} kpc at R = {R[np.argmax(h_R_cheng)]:.1f} kpc")

print("\n--- YOUNG STELLAR POPULATIONS ---")
print(f"Young Giants (Poggio+2025):   Max = {np.max(h_R_young_giants):.3f} kpc at R = {R[np.argmax(h_R_young_giants)]:.1f} kpc")
print(f"Cepheids (Poggio+2025):       Max = {np.max(h_R_cepheids):.3f} kpc at R = {R[np.argmax(h_R_cepheids)]:.1f} kpc")
print(f"Cepheids (Chen+2019):         Max = {np.max(h_R_chen_power):.3f} kpc at R = {R[np.argmax(h_R_chen_power)]:.1f} kpc")

print("\n--- GASEOUS COMPONENTS ---")
print(f"HI (Levine+2006):             Max = {np.max(h_R_hi):.3f} kpc at R = {R[np.argmax(h_R_hi)]:.1f} kpc")
print(f"CO m=1 - This work:           Max = {np.max(h_R_co_max):.3f} kpc at R = {R[np.argmax(h_R_co_max)]:.1f} kpc")
print(f"CO m=1,2 - This work:         Max = {np.max(h_R_co_m12):.3f} kpc at R = {R[np.argmax(h_R_co_m12)]:.1f} kpc ★")

# 特别比较CO m=1,2与其他模型
print(f"\n--- CO m=1,2 MODEL COMPARISON ---")
co_m12_max_R = R[np.argmax(h_R_co_m12)]
co_m12_max_val = np.max(h_R_co_m12)
print(f"CO m=1,2 warp peaks at R = {co_m12_max_R:.1f} kpc with amplitude {co_m12_max_val:.3f} kpc")

# 打印CO m=1,2模型参数
print(f"\n--- CO m=1,2 MODEL PARAMETERS ---")
print(f"a1 = {a1_m12}, Rw1 = {Rw1_m12} kpc, φw1 = {phi_w1_m12} rad")
print(f"a2 = {a2_m12}, Rw2 = {Rw2_m12} kpc, bw2 = {bw2_m12}, φw2 = {phi_w2_m12} rad")

plt.show()
