import numpy as np


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


T = 40 + 273
def rho_calc(TC, PC, ZC, MW):
    Tr = T / TC
    exponent = 1 + ((1 - Tr) ** (2 / 7))
    Vm = (8.314 * TC * (ZC ** exponent)) / PC
    rho = MW / Vm / 1e3
    return rho


print(rho_calc(647.0960000, 2.20640000E+7, 0.2290000000, 18.01528000), '992,298297429942')
print(rho_calc(568.7000000, 2.49000000E+6, 0.2560000000, 114.2309200), '690,323049539957')
print(rho_calc(536.6000000, 4.63000000E+6, 0.2640000000, 76.16252000), '819,328703204514')
