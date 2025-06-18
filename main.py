import numpy as np
import math

def exp_calc(m_flask, m_cap, m_aqueous_phase, m_liquid_phase, feed_pnsh_frac, equilibrium_pnsh_frac):

    mass_aqueous_phase = m_aqueous_phase - m_flask
    feed_frac_aqueous_phase = mass_aqueous_phase / mass_aqueous_phase * 100

    mass_liquid_phase = m_liquid_phase - mass_aqueous_phase - m_flask - m_cap

    mass_pnsh = (feed_pnsh_frac / 100) * mass_liquid_phase
    feed_frac_pnsh = mass_pnsh / mass_liquid_phase * 100

    mass_oct = (1 - (feed_pnsh_frac / 100)) * mass_liquid_phase
    feed_frac_oct = mass_oct / mass_liquid_phase * 100

    pnsh_frac_liquid_phase = equilibrium_pnsh_frac * (76.16 / 32.065)

    equilibrium_liquid_phase_pnsh_mass = mass_liquid_phase * (pnsh_frac_liquid_phase / 100)
    equilibrium_liquid_phase_mass = equilibrium_liquid_phase_pnsh_mass + mass_oct
    oct_frac_liquid_phase = mass_oct / equilibrium_liquid_phase_mass * 100

    equilibrium_aqueous_phase_pnsh_mass = mass_pnsh - equilibrium_liquid_phase_pnsh_mass
    equilibrium_aqueous_phase_mass = equilibrium_aqueous_phase_pnsh_mass + mass_aqueous_phase
    pnsh_frac_aqueous_phase = equilibrium_aqueous_phase_pnsh_mass / mass_aqueous_phase * 100
    water_frac_aqueous_phase = mass_aqueous_phase / equilibrium_aqueous_phase_mass * 100

    print('-----------#Исходнй состав-----------')
    print(f'-----------Водная фаза - {mass_aqueous_phase:.5g}, грамм:')
    print(f'Вода: {mass_aqueous_phase:.5g} грамм,  {feed_frac_aqueous_phase:.5g} % мас.')
    print(f'-----------Органическая фаза фаза - {mass_liquid_phase:.5} грамм:')
    print(f'Октан {mass_oct:.5g}, грамм,  {feed_frac_oct:.5g} % мас.')
    print(f'Пропантиол {mass_pnsh:.5g}, грамм,  {feed_frac_pnsh:.5g} % мас.')

    print('-----------#Равновесный состав-----------')
    print(f'-----------Водная фаза - {equilibrium_aqueous_phase_mass:.5g}, грамм:')
    print(f'Вода {mass_aqueous_phase:.5g} грамм, {water_frac_aqueous_phase:.5g} % мас.')
    print(f'Пропантиол {equilibrium_aqueous_phase_pnsh_mass:.5g} грамм, {pnsh_frac_aqueous_phase:.5g} % мас.')

    print(f'-----------Органическая фаза фаза - {equilibrium_liquid_phase_mass:.5g} грамм:')
    print(f'Октан {mass_oct:.5g} грамм, {oct_frac_liquid_phase:.5g} % мас.')
    print(f'Пропантиол {equilibrium_liquid_phase_pnsh_mass:.5g} грамм, {pnsh_frac_liquid_phase:.5g} % мас.')

    return pnsh_frac_aqueous_phase, pnsh_frac_liquid_phase, equilibrium_liquid_phase_pnsh_mass, equilibrium_aqueous_phase_pnsh_mass

def m_calc(Na, x, K_a, K_w, m0):
    return x / ((x / m0) + (2 * K_a * (x / m0)) / (np.sqrt(Na**2 + 4 * (K_w + K_a * (x / m0))) - Na))

def k_w(T):
    pKw = 14.94355 - 0.0429299 * T + 0.00021447938* T**2 - 5.6625156*10**(-7) * T**3
    return 10**-pKw

def pH_calc(m0, m_d, pKa):
    return pKa + np.log10(m0 / m_d - 1)

def m_d_calc(m0, pH, pKa):
    return m0 / (1 + 10**(pH - pKa))


def solve_equilibrium(
    C_NaOH0,          # моль/л, нач. NaOH в воде
    C_o0_HA,          # моль/л, нач. C3H7SH в органике
    m0,               # коэффициент распределения нейтральной формы
    Ka,               # константа диссоциации C3H7SH в воде (моль/л)
    Kw,               # ионное произведение воды
    Vo,               # л, объём органической фазы
    Vw,               # л, объём водной фазы
    xtol=1e-20,       # точность по H+
    max_iter=10000    # максимум итераций
):

    R = Vw / Vo

    def f(x):
        x = np.array(x, dtype=float)
        if x <= 0.0:
            return 1.0
        Cw_HA = C_o0_HA / (m0 + R * (1.0 + Ka / x))
        return x**2 + C_NaOH0 * x - Ka * Cw_HA - Kw

    a, b = 1e-20, 1e-1
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise RuntimeError("Bisection: root not bracketed on [1e-20, 1e-1]")

    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < xtol or (b - a) / 2 < xtol:
            x = c
            break
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    else:
        raise RuntimeError("Bisection did not converge")

    x = np.maximum(x, 1e-30)
    Cw_HA = C_o0_HA / (m0 + R * (1.0 + Ka / x))
    Cw_A  = Ka * Cw_HA / x
    Co_HA = m0 * Cw_HA
    D = m0 / (1.0 + Ka / x)

    return {
        "pH": -np.log10(x),
        "m_D": D,
        "H_plus": x,
        "OH_minus": Kw / x,
        "Cw_HA": Cw_HA,
        "Cw_A": Cw_A,
        "Co_HA": Co_HA
    }

def equilibrium_from_Co_eq(
    C_NaOH0,
    Co_eq_HA,
    P,
    Ka,
    Kw
):

    S = Ka * Co_eq_HA / P + Kw
    disc = C_NaOH0**2 + 4.0 * S
    x = (-C_NaOH0 + math.sqrt(disc)) / 2.0
    if x <= 0.0:
        raise ValueError("Calculated [H+] is non‑positive; check inputs.")

    pH = -math.log10(x)
    Cw_HA = Co_eq_HA / P
    A_minus = Ka * Cw_HA / x
    D = Co_eq_HA / (Cw_HA + A_minus)

    return {
        "H_plus / mol L⁻¹": x,
        "pH": pH,
        "m_D": D,
        "Cw_HA / mol L⁻¹": Cw_HA,
        "Cw_A⁻ / mol L⁻¹": A_minus,
        "OH_minus / mol L⁻¹": Kw / x
    }