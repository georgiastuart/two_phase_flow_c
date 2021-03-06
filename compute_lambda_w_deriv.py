import sympy as sp

if __name__ == '__main__':
    s, s_rw, s_ro, mu_w, mu_o, eta, k = sp.symbols('s, s_rw, s_ro, mu_w, mu_o, eta, k')

    k_rw = ((s - s_rw)**2 / (1 - s_rw)**2)
    k_ro = (1 - s / (1 - s_ro))**2

    lam = (k_ro / mu_o) + (k_rw / mu_w)
    lam_w = k_rw / (mu_w * lam)
    lam_o = k_ro / (mu_o * lam)

    k_rw_deriv = sp.Derivative(k_rw, s).doit()
    k_ro_deriv = sp.Derivative(k_ro, s).doit()
    lam_deriv = k_ro_deriv / mu_o + k_rw_deriv / mu_w

    lam_w_deriv = (mu_w * lam * k_rw_deriv - k_rw * mu_w * lam_deriv) / (mu_w * lam)**2

    zeta = s_ro**2 * (1 - s_ro - s_rw)**(-2)
    p_c = eta * ((s - s_rw)**(-2) - zeta * (1 - s)**(-2))
    diff = k * lam * lam_o * lam_w * sp.Derivative(p_c, s).doit()
    diff1 = diff.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.84, eta:3000, k:10**(-11)})
    diff2 = diff.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21, eta:3000, k:10**(-11)})
    diff_eff = 2 * diff1 * diff2 / (diff1 + diff2)
    beta = -400 / diff_eff
    xi = 2 * diff1 / 400
    A = xi / (beta * xi - 1)
    phi_h_dt = 0.2 * 400 / 86400
    sum_A = 7.306955e-12
    sum_A_R = -1.027624e-07
    g = 5.171493e-06

    print lam.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0})
    print sp.Derivative(lam_w, s).doit().subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21})
    print lam_w_deriv.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21})
    print k_rw_deriv.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21})
    print k_ro_deriv.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21})
    print sp.Derivative(p_c, s).doit().subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.21, eta:3000})
    print diff.subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, s:0.84, eta:3000, k:10**(-11)})
    print 'beta: %e' % (beta)
    print 'xi: %e' % (xi)
    print 'A: %e' % (A)
    print 'phi_h_dt: %e' % (phi_h_dt)
    print 'Sat: %s' % repr((g*400 + sum_A_R + phi_h_dt * 0.21) / (phi_h_dt + sum_A))
    # sp.plotting.plot(sp.Derivative(p_c, s).doit().subs({s_rw:0.2, s_ro: 0.15, mu_w:0.5, mu_o:10, eta:3000}), (s, 0.21, 0.84))
