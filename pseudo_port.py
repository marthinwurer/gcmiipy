
def aflux(u, v, p):
    """
    THIS SUBROUTINE CALCULATES THE HORIZONTAL AIR MASS FLUXES
    AND VERTICAL AIR MASS FLUXES AS DETERMINED BY U, V AND P.

    Args:
        u:
        v:
        p:

    Returns:

    """
    # indexes into the arrays
    i = 0
    j = 0
    l = 0 # layer

    # COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
    spa[i, j, l] = u[i, j, l] + u[i, j+1, l]
    avrx(spa)

    pu[i, j, l] = .25 * dyp[j]*spa[i, j, l]*(p[i, j]+p[i+1, j])

    # COMPUTE PV, THE SOUTH-NORTH MASS FLUX
    pv[i, j, l] = 0.25 * dxv[j] * (v[i, j, l] + v[i-1, j, l]) * (p(i, j)+p(i, j-1))

    # CONTINUITY EQUATION
    # COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
    conv[i, j, l] = (pu[i-1, j, l] - pu[i, j, l] + pv[i, j, l] - pv[i, j+1, l]) * dsig(l)

    # END OF HORIZONTAL ADVECTION LAYER LOOP

    # COMPUTE PIT, THE PRESSURE TENDENCY
    pit[i, j] += conv[i, j, l]

    # COMPUTE SD, SIGMA DOT
    sd[i, j, l] = sd[i, j, l+1] + conv[i, j, l+1] - dsig[l+1] * pit[i, j]


    return pit, conv, sd

def advec_m(p, pa, dt):
    """
    THIS SUBROUTINE CALCULATES UPDATED COLUMN PRESSURES AS
    DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES.

    Args:
        p:
        pa: the new surface pressure
        dt:

    Returns:

    """
    # COMPUTE PA, THE NEW SURFACE PRESSURE
    pa[i, j] = p[i, j] + (dt * pit[i, j] / dxyp[j])

def advec_v(pa, ut, vt, pb, u, v, p, dt):
    """
    THIS SUBROUTINE ADVECTS MOMENTUM (INCLUDING THE CORIOLIS FORCE)
    AS DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES.

    Args:
        pa:
        ut:
        vt:
        pb:
        u:
        v:
        p:
        dt:

    Returns:

    """
    sha = rgas/kapa
    kapap1 = kapa + 1

    # SCALE QT.  UT AND VT MAY THEN BE PHYSICALLY INTERPRETED AS
    # MOMENTUM COMPONENTS, TT AS HEAT CONTENT, AND QT AS WATER CONTENT
    fd[i, j] = pa[i, j] * dxyp[j]

    dut[i, j, l] = 0
    ut[i, j, l] = ut[i, j, l] * (.25 * fd[i, j] + fd[i+1, j] + fd[i, j-1] + fd[i+1, j-1])

    # HORIZONTAL ADVECTION OF MOMENTUM

    # CONTRIBUTION FROM THE WEST-EAST MASS FLUX
    # todo

    # CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX
    # todo

    # CONTRIBUTION FROM THE SOUTHWEST-NORTHEAST MASS FLUX
    # todo

    # CONTRIBUTION FROM THE SOUTHEAST-NORTHWEST MASS FLUX
    # todo

    # VERTICAL ADVECTION OF MOMENTUM
    # todo

    # ADD ADVECTION INCREMENTS TO UT AND VT, CALL DIAGNOSTICS
    # todo

    # CORIOLIS FORCE
    # TODO

    # ADD CORIOLIS FORCE INCREMENTS TO UT AND VT, CALL DIAGNOSTICS
    # todo

    # UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
    fd[i, j] = pb[i, j] * dxyp[j]

    ut[i, j, l] = ut[i, j, l] * 4 / (fd[i, j] + fd[i+1, j] + fd[i, j-1] + fd[i+1, j-1])

    pass

def pgf(ut, vt, pb, u, v, t, p, dt):
    """
    THIS SUBROUTINE ADDS TO MOMENTUM THE TENDENCIES DETERMINED BY
    THE PRESSURE GRADIENT FORCE

    Args:
        ut:
        vt:
        pb:
        u:
        v:
        t:
        p:
        dt:

    Returns:

    """
    sha = rgas/kapa
    kapap1 = kapa + 1

    # VERTICAL DIFFERENCING

    # this is that one thing with the air density and stuff, iirc
    spa[i, j, l] = sig[l] * p[i, j] * rgas * t(i, j, l) * pkdn / pdn

    theta = thbar(t[i, j, l+1], t[i, j, l])

    # CALCULATE THE DIFFERENCE IN PHI BETWEEN ODD LEVELS LP1 AND L
    phi[i, j, l+1] = sha * theta * (pkdn - pkup)
    phi[i, j, 1] = fdata[i, j, 1] + sum1 - sum2

    # TODO there's a lot of stuff in here

    # PRESSURE GRADIENT FORCE

    # NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
    flux[i, j, l] = dt / 4 * ((p[i, j] + p[i, j-1]) * (phi[i, j, l] - phi[i, j-1, l]) +
                              (spa[i, j, l] + spa[i, j-1, l]) * (p[i, j] - p[i, j-1])) * dxv(j)
    dvt[i, j, l] -= flux
    dvt[i-1, j, l] -= flux

    # SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
    pu[i, j, l] = (p[i+1, j] + p[i, j]) * (phi[i+1, j, l] - phi[i, j, l]) + \
                  (spa[i+1, j, l] - spa[i, j, l]) * (p[i+1, j] - p[i, j])
    avrx(pu)

    dut[i, j, l] = dut[i, j, l] - dt/4 * dyv[j] * (pu[i, j, l] + pu[i, j-1, l])

    # UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
    # TODO







def advec_t(pa, tt, pb, t, dt):
    """
    THIS SUBROUTINE ADVECTS POTENTIAL TEMPERATURE

    Args:
        pa:
        tt:
        pb:
        t:
        dt:

    Returns:

    """
    sha = rgas/kapa
    kapap1 = kapa + 1

    # SCALE QT.  UT AND VT MAY THEN BE PHYSICALLY INTERPRETED AS
    # MOMENTUM COMPONENTS, TT AS HEAT CONTENT, AND QT AS WATER CONTENT
    fd[i, j] * dxyp[j]
    tt[i, j, l] *= fd[i, j]

    # HORIZONTAL ADVECTION OF HEAT

    # WEST-EAST ADVECTION OF HEAT
    fluxq[i] = dt/2 * pu(i, j, l) * (t[i, j, l] + t[i+1, j, l])
    tt(i, j, l) += (fluxq[i - 1] - fluxq[i])


