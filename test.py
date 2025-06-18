import numpy as np
import matplotlib.pyplot as plt
import main

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

if __name__ == "__main__":

    Exp_m_d_30 = np.array([5.40157032e+01, 7.94231577e-01, 2.73707924e-01, 1.25100363e-01,
 8.68999735e-02, 3.73046710e-02, 2.08875696e-02, 2.14223231e-02,
 1.98380262e-02, 1.84384422e-02])
    Exp_NaOH_frac = np.array([0.0, 0.25201827, 0.55192001, 0.77330606, 1.0868788,  1.59761583,
 2.34094472, 2.82770001, 3.65645509, 3.81696673])
    Exp_m_d_40 = np.array([5.71423715e+01, 9.72374218e-01, 3.79969890e-01, 1.66252148e-01,
 1.03471738e-01, 5.06351594e-02, 3.05205460e-02, 2.34641209e-02,
 2.27509544e-02, 1.96770724e-02])
    Exp_m_d_50 = np.array([8.86573121e+01, 1.18066936e+00, 5.48059538e-01, 2.49549777e-01,
 1.38884926e-01, 8.61727527e-02, 3.65393596e-02, 2.94443417e-02,
 3.26783734e-02, 2.15940498e-02])

    base_params = dict(
        C_o0_HA = 0.1958,
        m0      = 89,
        Ka      = 10**-10.4,
        Kw      = main.k_w(50),
        Vo      = 0.010,   # 10 мл органика
        Vw      = 0.010    # 10 мл вода
    )

    C_NaOH0_array = np.linspace(0.0, 4, 57)

    pH_list = []
    D_list = []

    for C_NaOH0 in C_NaOH0_array:
        params = {**base_params, "C_NaOH0": C_NaOH0}
        result = solve_equilibrium(**params)
        print(f"C_NaOH0 = {C_NaOH0:.2f} моль/л:")
        for k, v in result.items():
            print(f"  {k:12s}: {v:.5g}")
        print()

        # собираем данные
        pH_list.append(result["pH"])
        D_list.append(result["m_D"])

    #plt.figure()
    #plt.plot(C_NaOH0_array, pH_list, color='Black', linestyle='-')
    #plt.ylim(12.5, 14.5)
    #plt.xlabel('C_NaOH0, моль/л')
    #plt.ylabel('pH')
    #plt.title('Зависимость pH от начальной концентрации NaOH')
    #plt.grid(True)
    #plt.show()

    plt.figure()
    plt.plot(C_NaOH0_array, D_list, color='Black', linestyle='-')
    #plt.scatter(Exp_NaOH_frac, Exp_m_d_30,label='30', color='Red')
    #plt.scatter(Exp_NaOH_frac, Exp_m_d_40,label='40', color='Green')
    plt.scatter(Exp_NaOH_frac, Exp_m_d_50,label='50', color='Blue')
    plt.ylim(0, 1.25)
    plt.xlabel('C_NaOH0, моль/л')
    plt.ylabel('m_D (коэффициент распределения)')
    plt.title('Зависимость коэффициента распределения D от C_NaOH0')
    plt.grid(True)
    plt.show()

