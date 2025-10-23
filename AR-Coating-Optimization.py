#!/usr/bin/env python3
"""
Optimize AR Coating for CaF₂ Filter
====================================
Use optimization algorithm to find best AR coating thicknesses
to minimize ripple in 250-600nm while maximizing transmission.

Goal: 
- Minimize ripple in passband (250-600nm)
- Maximize average transmission
- Maintain good performance at critical wavelengths
"""

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Import optical data
CAF2_DATA = {
    'wavelength_nm': np.array([
        150, 160, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250,
        260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600
    ]),
    'n': np.array([
        1.50, 1.48, 1.47, 1.46, 1.46, 1.45, 1.45, 1.45, 1.45, 1.44, 
        1.44, 1.44, 1.43, 1.43, 1.43, 1.43, 1.43, 1.43, 1.43, 1.43, 
        1.43, 1.43, 1.43, 1.43, 1.43, 1.43
    ]),
    'k': np.array([
        0.002, 5e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 
        2e-7, 1e-7, 5e-8, 2e-8, 1e-8, 5e-9, 2e-9, 1e-9, 5e-10, 2e-10, 
        1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10
    ])
}

MGF2_DATA = {
    'wavelength_nm': np.array([
        160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 280,
        300, 320, 340, 360, 380, 400, 450, 500, 550, 600
    ]),
    'n': np.array([
        1.430, 1.422, 1.416, 1.411, 1.407, 1.404, 1.401, 1.399, 1.397, 
        1.395, 1.393, 1.390, 1.388, 1.386, 1.385, 1.383, 1.382, 1.381, 
        1.379, 1.378, 1.377, 1.376
    ]),
    'k': np.array([1e-9] * 22)
}

# Create interpolators
caf2_n_interp = interp1d(CAF2_DATA['wavelength_nm'], CAF2_DATA['n'], 
                         kind='cubic', fill_value='extrapolate')
caf2_k_interp = interp1d(CAF2_DATA['wavelength_nm'], CAF2_DATA['k'],
                         kind='cubic', fill_value='extrapolate')
mgf2_n_interp = interp1d(MGF2_DATA['wavelength_nm'], MGF2_DATA['n'],
                         kind='cubic', fill_value='extrapolate')

def get_n(material, wavelength_nm):
    """Get refractive index"""
    if material == 'air':
        return 1.0 + 0j
    elif material == 'caf2':
        n = float(caf2_n_interp(wavelength_nm))
        k = float(caf2_k_interp(wavelength_nm))
        return n - 1j*k
    elif material == 'mgf2':
        n = float(mgf2_n_interp(wavelength_nm))
        return n - 1j*1e-9
    return 1.0 + 0j

def transfer_matrix(n, d_nm, wl_nm):
    """2x2 transfer matrix for single layer"""
    delta = 2 * np.pi * n * d_nm / wl_nm
    cos_d = np.cos(delta)
    sin_d = np.sin(delta)
    
    m11 = cos_d
    m12 = 1j * sin_d / n
    m21 = 1j * n * sin_d
    m22 = cos_d
    
    return np.array([[m11, m12], [m21, m22]])

def calculate_transmission(d1_nm, d2_nm, wavelength_nm):
    """
    Calculate transmission with given AR coating thicknesses
    
    Structure: Air | MgF2(d1) | MgF2(d2) | CaF2(substrate)
    """
    # Get refractive indices
    n0 = get_n('air', wavelength_nm)
    n1 = get_n('mgf2', wavelength_nm)
    n2 = get_n('mgf2', wavelength_nm)
    ns = get_n('caf2', wavelength_nm)
    
    # Total transfer matrix
    M = transfer_matrix(n1, d1_nm, wavelength_nm) @ \
        transfer_matrix(n2, d2_nm, wavelength_nm)
    
    m11, m12 = M[0,0], M[0,1]
    m21, m22 = M[1,0], M[1,1]
    
    # Reflection coefficient
    r = (n0*m11 + n0*ns*m12 - m21 - ns*m22) / \
        (n0*m11 + n0*ns*m12 + m21 + ns*m22)
    
    # AR coating transmission (no substrate absorption here)
    T_coating = 1 - np.abs(r)**2
    
    # Add substrate absorption
    k_sub = -np.imag(ns)
    if k_sub > 1e-10:
        T_substrate = np.exp(-4*np.pi*k_sub*1.5e-3/(wavelength_nm*1e-9))
    else:
        T_substrate = 1.0
    
    return T_coating * T_substrate

def fitness_function(params):
    """
    Fitness function to minimize
    
    params: [d1_nm, d2_nm] - thicknesses of two MgF2 layers
    
    Objectives:
    1. Minimize ripple in 250-600nm
    2. Maximize average transmission in 200-600nm
    3. Good performance at key wavelengths (200, 220, 250nm)
    """
    d1, d2 = params
    
    # Wavelengths to evaluate
    wavelengths = np.linspace(200, 600, 81)  # 5nm steps
    
    T = np.array([calculate_transmission(d1, d2, wl) for wl in wavelengths])
    
    # Region 1: 250-600nm (should be flat and high)
    idx_250_600 = wavelengths >= 250
    T_passband = T[idx_250_600]
    
    ripple = (T_passband.max() - T_passband.min()) / T_passband.mean()
    avg_transmission = T_passband.mean()
    
    # Region 2: 200-250nm (critical UV region)
    idx_200_250 = (wavelengths >= 200) & (wavelengths < 250)
    T_uv = T[idx_200_250].mean()
    
    # Key wavelengths
    T_200 = calculate_transmission(d1, d2, 200)
    T_220 = calculate_transmission(d1, d2, 220)
    T_400 = calculate_transmission(d1, d2, 400)
    
    # Fitness (lower is better)
    # Penalize: ripple, low transmission, low UV
    fitness = (
        10 * ripple +                    # Minimize ripple (weight=10)
        2 * (1 - avg_transmission) +     # Maximize passband transmission
        3 * (1 - T_uv) +                 # Maximize UV transmission (weight=3)
        5 * (1 - T_200) +                # Critical: 200nm (weight=5)
        2 * (1 - T_220) +                # Important: 220nm
        1 * (1 - T_400)                  # Reference: 400nm
    )
    
    return fitness

def main():
    print("=" * 80)
    print("AR COATING OPTIMIZATION FOR CaF₂ FILTER")
    print("=" * 80)
    print()
    print("Current design:")
    print("  Layer 1 (outer): 72.4 nm MgF₂")
    print("  Layer 2 (inner): 61.6 nm MgF₂")
    print("  Substrate: 1.5 mm CaF₂")
    print()
    print("Goals:")
    print("  1. Minimize ripple in 250-600nm passband")
    print("  2. Maximize transmission in 200-250nm critical region")
    print("  3. Maintain high overall transmission")
    print()
    
    # Run optimization
    print("🚀 Running optimization (differential evolution)...")
    print("   Search space: 20-150 nm per layer")
    print("   Strategy: Adaptive, multi-core")
    print()
    
    bounds = [(20, 150), (20, 150)]  # Thickness bounds for d1, d2
    
    result = differential_evolution(
        fitness_function,
        bounds,
        strategy='best1bin',
        maxiter=200,
        popsize=15,
        tol=0.0001,
        mutation=(0.5, 1.5),
        recombination=0.7,
        seed=42,
        workers=-1,  # Use all CPU cores
        updating='deferred',
        disp=True
    )
    
    d1_opt, d2_opt = result.x
    
    print()
    print("=" * 80)
    print("OPTIMIZATION RESULTS")
    print("=" * 80)
    print()
    print("Optimized AR coating:")
    print(f"  Layer 1 (outer): {d1_opt:.2f} nm MgF₂")
    print(f"  Layer 2 (inner): {d2_opt:.2f} nm MgF₂")
    print()
    
    # Compare performance
    wavelengths = np.linspace(160, 600, 441)
    
    # Original design
    T_original = np.array([calculate_transmission(72.4, 61.6, wl) for wl in wavelengths])
    
    # Optimized design
    T_optimized = np.array([calculate_transmission(d1_opt, d2_opt, wl) for wl in wavelengths])
    
    # Calculate metrics
    idx_250_600 = wavelengths >= 250
    
    print("Performance comparison:")
    print("-" * 80)
    print(f"{'Metric':<30} {'Original':<15} {'Optimized':<15} {'Change'}")
    print("-" * 80)
    
    for wl in [175, 200, 220, 250, 300, 400, 500]:
        idx = np.argmin(np.abs(wavelengths - wl))
        T_orig = T_original[idx] * 100
        T_opt = T_optimized[idx] * 100
        delta = T_opt - T_orig
        print(f"T({wl}nm) [%]               {T_orig:>6.2f}%        {T_opt:>6.2f}%        {delta:+6.2f}%")
    
    print()
    
    # Ripple
    ripple_orig = (T_original[idx_250_600].max() - T_original[idx_250_600].min()) / \
                  T_original[idx_250_600].mean() * 100
    ripple_opt = (T_optimized[idx_250_600].max() - T_optimized[idx_250_600].min()) / \
                 T_optimized[idx_250_600].mean() * 100
    
    print(f"Ripple (250-600nm) [%]    {ripple_orig:>6.2f}%        {ripple_opt:>6.2f}%        {ripple_opt-ripple_orig:+6.2f}%")
    
    # Average transmission
    avg_orig = T_original[idx_250_600].mean() * 100
    avg_opt = T_optimized[idx_250_600].mean() * 100
    
    print(f"Avg T (250-600nm) [%]     {avg_orig:>6.2f}%        {avg_opt:>6.2f}%        {avg_opt-avg_orig:+6.2f}%")
    
    print("-" * 80)
    
    # Plot comparison
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Full spectrum
    ax1.plot(wavelengths, T_original*100, 'b-', linewidth=2.5, 
            label='Original (72.4, 61.6 nm)', alpha=0.7)
    ax1.plot(wavelengths, T_optimized*100, 'r-', linewidth=2.5,
            label=f'Optimized ({d1_opt:.1f}, {d2_opt:.1f} nm)')
    
    ax1.axvline(175, color='gray', linestyle='--', alpha=0.5)
    ax1.axvline(250, color='gray', linestyle='--', alpha=0.5)
    ax1.axhline(95, color='gray', linestyle=':', alpha=0.5)
    
    ax1.set_xlabel('Wavelength (nm)', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Transmittance (%)', fontsize=13, fontweight='bold')
    ax1.set_title('AR Coating Optimization: Before vs After',
                 fontsize=15, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(160, 600)
    ax1.set_ylim(0, 102)
    
    # Zoom on passband
    idx_zoom = wavelengths >= 250
    ax2.plot(wavelengths[idx_zoom], T_original[idx_zoom]*100, 'b-', 
            linewidth=2.5, label='Original', alpha=0.7)
    ax2.plot(wavelengths[idx_zoom], T_optimized[idx_zoom]*100, 'r-',
            linewidth=2.5, label='Optimized')
    
    ax2.axhline(95, color='gray', linestyle=':', alpha=0.5, label='95% target')
    
    ax2.set_xlabel('Wavelength (nm)', fontsize=13, fontweight='bold')
    ax2.set_ylabel('Transmittance (%)', fontsize=13, fontweight='bold')
    ax2.set_title(f'Passband Detail: Ripple {ripple_orig:.1f}% → {ripple_opt:.1f}%',
                 fontsize=15, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(250, 600)
    ax2.set_ylim(80, 102)
    
    plt.tight_layout()
    plt.savefig('ar_optimization_result.png', dpi=300, bbox_inches='tight')
    plt.savefig('ar_optimization_result.pdf', bbox_inches='tight')
    print()
    print("✅ Plots saved: ar_optimization_result.png / .pdf")
    
    print()
    print("=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print()
    if ripple_opt < ripple_orig * 0.8:
        print("✅ Significant improvement in ripple reduction!")
    else:
        print("⚠️ Limited improvement - current design is already near-optimal")
    print()
    print("Note: 175nm blocking and transition steepness are determined by")
    print("      CaF₂ substrate properties and CANNOT be improved by AR coating.")
    print("=" * 80)

if __name__ == '__main__':
    main()

