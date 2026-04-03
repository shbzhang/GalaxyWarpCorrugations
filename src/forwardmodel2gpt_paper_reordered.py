
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from forwardmodel2gpt import load_real_data, weighted_binned_statistic, a5_vlsr_from_distance
from Distance import Distance, A5


def direct_los_streaming_field(x, y, pitch_angle_deg=12.0, v_sys_amp=6.0, m_arms=2.0, R_ref=9.0):
    R = np.sqrt(x**2 + y**2)
    phi = np.arctan2(x, y)
    pitch = np.radians(pitch_angle_deg)
    phase = m_arms * (phi - np.log(np.clip(R, 1e-6, None) / R_ref) / np.tan(pitch))
    return v_sys_amp * np.sin(phase), phase


def plot_direct_streaming_faceon(
    x_true, y_true, v_systematic,
    pitch_angle_deg=12.0, v_sys_amp=6.0,
    save_prefix='Empirical_DirectLOS_Paper'
):
    xg = np.linspace(-16.0, 16.0, 500)
    yg = np.linspace(-16.0, 16.0, 500)
    XG, YG = np.meshgrid(xg, yg)
    VG, _ = direct_los_streaming_field(XG, YG, pitch_angle_deg=pitch_angle_deg, v_sys_amp=v_sys_amp)

    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    im = ax.imshow(
        VG, origin='lower', extent=[xg.min(), xg.max(), yg.min(), yg.max()],
        cmap='coolwarm', vmin=-max(12, 1.2*v_sys_amp), vmax=max(12, 1.2*v_sys_amp)
    )
    ax.scatter(x_true, y_true, c=v_systematic, s=8, cmap='coolwarm',
               vmin=-max(12, 1.2*v_sys_amp), vmax=max(12, 1.2*v_sys_amp),
               edgecolors='none', alpha=0.55)
    ax.scatter([0.0], [A5['Ro']], s=130, c='gold', edgecolors='black', marker='*', zorder=5, label='Sun')
    ax.scatter([0.0], [0.0], s=60, c='black', marker='x', zorder=5, label='GC')
    ax.set_xlabel('x (kpc)', fontsize=14)
    ax.set_ylabel('y (kpc)', fontsize=14)
    ax.set_title('Face-on distribution of injected LOS streaming', fontsize=15)
    ax.set_xlim(-16, 16)
    ax.set_ylim(-16, 16)
    ax.legend(loc='upper right', framealpha=0.9)
    cbar = plt.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label(r'Injected $\Delta v_{\rm sys,LOS}$ (km s$^{-1}$)', fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{save_prefix}_StreamingFaceOn.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{save_prefix}_StreamingFaceOn.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def _solve_kinematic_with_A5(l_true_deg, l_true_rad, V_obs, b_rad):
    R_sun = A5['Ro']
    b_deg = np.degrees(b_rad)
    skyc = SkyCoord(l=l_true_deg * u.deg, b=b_deg * u.deg, frame='galactic')
    farnear = np.zeros(len(V_obs), dtype=int)
    dist_solver = Distance(skyc=skyc, v_lsr=V_obs, farnear=farnear, unit=A5)
    _, d_obs, _, _ = dist_solver.Kdist()
    d_obs = np.asarray(d_obs, dtype=float)

    Z_obs = d_obs * np.tan(b_rad)
    d_proj_obs = d_obs * np.cos(b_rad)
    x_obs = d_proj_obs * np.sin(l_true_rad)
    y_obs = R_sun - d_proj_obs * np.cos(l_true_rad)
    R_obs = np.sqrt(R_sun**2 + d_proj_obs**2 - 2.0 * R_sun * d_proj_obs * np.cos(l_true_rad))
    phi_obs_deg = np.degrees(np.arctan2(x_obs, y_obs))
    return R_obs, phi_obs_deg, Z_obs, d_obs


def _build_mock(df_real, pitch_angle_deg=12.0, v_sys_amp=6.0, v_rand_amp=3.0,
                intrinsic_sigma_z_kpc=0.06, radial_wavelength_kpc=4.5, random_seed=42):
    rng = np.random.default_rng(random_seed)

    l_true_deg = df_real['l_cen'].values
    l_true_rad = np.radians(l_true_deg)
    d_true = df_real['distance'].values
    masses = df_real['mass'].values
    x_true = df_real['xx'].values
    y_true = df_real['yy'].values
    R_true = df_real['gal_ridus'].values
    phi_true_deg = np.degrees(np.arctan2(x_true, y_true))

    intrinsic_noise = rng.normal(0.0, intrinsic_sigma_z_kpc, len(df_real))
    Z_true_radial = 0.200 * np.sin(2.0 * np.pi * (R_true - 9.0) / radial_wavelength_kpc) + intrinsic_noise
    Z_true_flat = intrinsic_noise

    b_true_rad = np.arcsin(Z_true_radial / np.sqrt(d_true**2 + Z_true_radial**2))
    b_true_flat = np.arcsin(Z_true_flat / np.sqrt(d_true**2 + Z_true_flat**2))

    V_systematic, phase = direct_los_streaming_field(x_true, y_true, pitch_angle_deg=pitch_angle_deg, v_sys_amp=v_sys_amp)
    V_random_rad = rng.normal(0.0, v_rand_amp, len(df_real))
    V_random_flat = rng.normal(0.0, v_rand_amp, len(df_real))

    def make_obs(b_true, v_rand):
        V_model = a5_vlsr_from_distance(l_true_deg, b_true, d_true, unit=A5)
        V_obs = V_model + V_systematic + v_rand
        R_obs, phi_obs_deg, Z_obs, d_obs = _solve_kinematic_with_A5(l_true_deg, l_true_rad, V_obs, b_true)
        return dict(R_obs=R_obs, phi_obs_deg=phi_obs_deg, Z_obs=Z_obs, d_obs=d_obs, V_obs=V_obs)

    return {
        'l_true_deg': l_true_deg, 'd_true': d_true, 'masses': masses,
        'x_true': x_true, 'y_true': y_true, 'R_true': R_true, 'phi_true_deg': phi_true_deg,
        'Z_true_radial': Z_true_radial, 'Z_true_flat': Z_true_flat,
        'radial': make_obs(b_true_rad, V_random_rad),
        'null': make_obs(b_true_flat, V_random_flat),
        'V_systematic': V_systematic,
    }


def _plot_slice(ax, X_obs, Z_obs, X_true, Z_true, mass, mask_obs, mask_true, bins, panel_label, title_text):
    X_c, Z_c, M_c = X_obs[mask_obs], Z_obs[mask_obs], mass[mask_obs]
    sizes = np.clip(np.sqrt(M_c) / 3.0, 5, 120)
    ax.scatter(X_c, Z_c * 1000.0, s=sizes, c='gray', alpha=0.24, edgecolors='none', label='Mock data')
    bin_c, mean_Z, std_Z = weighted_binned_statistic(X_c, Z_c, M_c, bins)
    valid = ~np.isnan(mean_Z) & ~np.isnan(std_Z)
    ax.fill_between(bin_c[valid], (mean_Z[valid]-std_Z[valid])*1000.0, (mean_Z[valid]+std_Z[valid])*1000.0,
                    color='red', alpha=0.14, label=r'1$\sigma$ mass-weighted dispersion')
    ax.plot(bin_c, mean_Z*1000.0, color='red', lw=3.0, marker='o', markersize=7.5, markeredgecolor='white', label='Recovered mean')
    bin_c_true, mean_Z_true, _ = weighted_binned_statistic(X_true[mask_true], Z_true[mask_true], mass[mask_true], bins)
    ax.plot(bin_c_true, mean_Z_true*1000.0, color='blue', lw=2.5, ls='--', label='Injected true signal')
    ax.axhline(0.0, color='black', lw=1.4, ls='--')
    ax.text(-0.10, 1.01, panel_label, transform=ax.transAxes, fontsize=20, fontweight='bold', va='bottom', ha='right')
    ax.text(0.50, 0.95, title_text, transform=ax.transAxes, fontsize=13.0, ha='center', va='top',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.86, edgecolor='gray'))
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(-700, 700)
    ax.tick_params(top=True, right=True, labelsize=12.5, direction='in', length=5)
    ax.set_xlabel('Galactocentric Radius $R$ (kpc)', fontsize=14)


def make_paper_figure(
    df_real,
    baseline_kwargs=None,
    stress_kwargs=None,
    radial_sector_phi=(29.0, 38.0),
    radial_bins=np.arange(9.0, 16.0, 0.6),
    save_prefix='Empirical_DirectLOS_Paper',
):
    if baseline_kwargs is None:
        baseline_kwargs = dict(pitch_angle_deg=12.0, v_sys_amp=6.0, v_rand_amp=3.0, intrinsic_sigma_z_kpc=0.06, radial_wavelength_kpc=4.5, random_seed=42)
    if stress_kwargs is None:
        stress_kwargs = dict(pitch_angle_deg=12.0, v_sys_amp=12.0, v_rand_amp=5.0, intrinsic_sigma_z_kpc=0.06, radial_wavelength_kpc=4.5, random_seed=42)

    base = _build_mock(df_real, **baseline_kwargs)
    stress = _build_mock(df_real, **stress_kwargs)

    phi_min, phi_max = radial_sector_phi
    plt.rcParams.update({'font.family': 'sans-serif', 'axes.linewidth': 1.4, 'xtick.major.width': 1.3, 'ytick.major.width': 1.3, 'xtick.direction': 'in', 'ytick.direction': 'in'})
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.2))

    panels = [
        ('a', axs[0,0], stress, 'radial', 'Injected radial corrugation\nStress test: $v_{\\rm sys}=12$, $\\sigma=5$ km s$^{-1}$'),
        ('b', axs[0,1], base, 'radial', 'Injected radial corrugation\nBaseline: $v_{\\rm sys}=6$, $\\sigma=3$ km s$^{-1}$'),
        ('c', axs[1,0], stress, 'null', 'Null hypothesis: flat disk\nStress test: $v_{\\rm sys}=12$, $\\sigma=5$ km s$^{-1}$'),
        ('d', axs[1,1], base, 'null', 'Null hypothesis: flat disk\nBaseline: $v_{\\rm sys}=6$, $\\sigma=3$ km s$^{-1}$'),
    ]
    for panel_label, ax, dat, key, title_text in panels:
        obs = dat[key]
        z_true = dat['Z_true_radial'] if key == 'radial' else dat['Z_true_flat']
        mask_obs = (obs['phi_obs_deg'] >= phi_min) & (obs['phi_obs_deg'] <= phi_max) & ~np.isnan(obs['R_obs'])
        mask_true = (dat['phi_true_deg'] >= phi_min) & (dat['phi_true_deg'] <= phi_max)
        _plot_slice(ax, obs['R_obs'], obs['Z_obs'], dat['R_true'], z_true, dat['masses'], mask_obs, mask_true,
                    radial_bins, panel_label, title_text + f'\nSector: $\\phi \\in [{phi_min:.0f}^\\circ, {phi_max:.0f}^\\circ]$')
        if ax in [axs[0,0], axs[1,0]]:
            ax.set_ylabel('$\mathbf{\Delta}$Z (pc)', fontsize=14)
        if panel_label == 'a':
            ax.legend(fontsize=10.2, loc='lower left', framealpha=0.88, edgecolor='black')

    plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.08, hspace=0.22, wspace=0.16)
    plt.savefig(f'{save_prefix}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{save_prefix}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    plot_direct_streaming_faceon(base['x_true'], base['y_true'], base['V_systematic'],
                                 pitch_angle_deg=baseline_kwargs['pitch_angle_deg'], v_sys_amp=baseline_kwargs['v_sys_amp'],
                                 save_prefix=f'{save_prefix}_Baseline')
    plot_direct_streaming_faceon(stress['x_true'], stress['y_true'], stress['V_systematic'],
                                 pitch_angle_deg=stress_kwargs['pitch_angle_deg'], v_sys_amp=stress_kwargs['v_sys_amp'],
                                 save_prefix=f'{save_prefix}_Stress')

    diagnostics = {
        'baseline': {'params': baseline_kwargs,
                     'median_abs_distance_error_radial_kpc': float(np.nanmedian(np.abs(base['radial']['d_obs'] - base['d_true']))),
                     'median_abs_distance_error_null_kpc': float(np.nanmedian(np.abs(base['null']['d_obs'] - base['d_true'])))},
        'stress': {'params': stress_kwargs,
                   'median_abs_distance_error_radial_kpc': float(np.nanmedian(np.abs(stress['radial']['d_obs'] - stress['d_true']))),
                   'median_abs_distance_error_null_kpc': float(np.nanmedian(np.abs(stress['null']['d_obs'] - stress['d_true'])))},
        'note': 'Azimuthal injection removed from figures by design; mention only briefly in text.'
    }
    return diagnostics


if __name__ == '__main__':
    df_catalog = load_real_data()
    diagnostics = make_paper_figure(df_catalog)
    print(diagnostics)
