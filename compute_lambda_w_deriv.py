import sympy as sp

if __name__ == '__main__':
    s, s_rw, s_ro, mu_w, mu_o = sp.symbols('s, s_rw, s_ro, mu_w, mu_o')

    k_rw = ((s - s_rw)**2 / (1 - s_rw)**2)**2
    k_ro = (1 - s / (1 - s_ro))**2

    lam = (k_ro / mu_o) + (k_rw / mu_w)
    lam_w = k_rw / (mu_w * lam)

    k_rw_deriv = sp.Derivative(k_rw, s).doit()
    k_ro_deriv = sp.Derivative(k_ro, s).doit()
    lam_deriv = k_ro_deriv / mu_o + k_rw_deriv / mu_w

    lam_w_deriv = (mu_w * lam * k_rw_deriv - k_rw * mu_w * lam_deriv) / (mu_w * lam)**2

    print lam.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:1})
    print sp.Derivative(lam_w, s).doit().subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:1})
    print lam_w_deriv.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:1})
