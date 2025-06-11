print(f'Концентрация NaOH 10.18 % масс')
print(f'40')
exp10_50 = main.exp_calc(16.4635, 1.4885, 27.7406, 36.3746, feed_pnsh_frac, 0.0227)

print('-----------##Равновесный состав по hysys-----------')
equilibrium_aqueous_phase_mass_hysys_10_50, frac_aqueous_phase_hysys_10_50 = 0.1480, 1.392985
print(f'{equilibrium_aqueous_phase_mass_hysys_10_50:.5g} грамм, {frac_aqueous_phase_hysys_10_50:.5g} % масс пропантиола в воде')

equilibrium_liquid_phase_mass_hysys_10_50, frac_liquid_phase_hysys_10_50 = 0.001831, 0.026209
print(f'{equilibrium_liquid_phase_mass_hysys_10_50:.5g} грамм, {frac_liquid_phase_hysys_10_50:.5g} % масс пропантиола в октане')

NaOH_frac_aqueous_phase_10_50 = 10.18 / 100
rho_NaOH_10_50 = 1107
C_NaOH_on_aqueous_phase_10_50 = NaOH_frac_aqueous_phase_10_50 * rho_NaOH_10_50 / MW_NaOH
print(f'Концентрация NaOH в воде: {C_NaOH_on_aqueous_phase_10_50:.5g} моль/л')

pnsh_frac_aqueous_phase_10_50 = exp10_50[0] / 100
rho_aqueous_phase_10_50 = 1097
C_pnsh_on_aqueous_phase_10_50 = pnsh_frac_aqueous_phase_10_50 * rho_aqueous_phase_10_50 / MW_pnsh
print(f'Концентрация пропантиола в воде: {C_pnsh_on_aqueous_phase_10_50:.5g} моль/л')

pnsh_frac_liquid_phase_10_50 = exp10_50[1] / 100
rho_liquid_phase_10_50 = 665.8
C_pnsh_on_liquid_phase_10_50 = pnsh_frac_liquid_phase_10_50 * rho_liquid_phase_10_50 / MW_pnsh
print(f'Концентрация пропантиола в октане: {C_pnsh_on_liquid_phase_10_50:.5g} моль/л')

m_10_50 = C_pnsh_on_liquid_phase_10_50 / C_pnsh_on_aqueous_phase_10_50
print(f'Коэффициент распределения: {m_10_50:.2g}')

pnsh_frac_aqueous_phase_10_50_hysys = frac_aqueous_phase_hysys_10_50 / 100
C_pnsh_on_aqueous_phase_10_50_hysys = pnsh_frac_aqueous_phase_10_50_hysys * rho_aqueous_phase_10_50 / MW_pnsh
print(f'Концентрация пропантиола в воде: {C_pnsh_on_aqueous_phase_10_50_hysys:.5g} моль/л (hysys)')

pnsh_frac_liquid_phase_10_50_hysys = frac_liquid_phase_hysys_10_50 / 100
C_pnsh_on_liquid_phase_10_50_hysys = pnsh_frac_liquid_phase_10_50_hysys * rho_liquid_phase_10_50 / MW_pnsh
print(f'Концентрация пропантиола в октане: {C_pnsh_on_liquid_phase_10_50_hysys:.5g} моль/л (hysys)')

m_10_50_hysys = C_pnsh_on_liquid_phase_10_50_hysys / C_pnsh_on_aqueous_phase_10_50_hysys
print(f'Коэффициент распределения: {m_10_50_hysys:.2g} (hysys)')