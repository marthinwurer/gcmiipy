"""
Discretizing the full primitive equations on a 1d staggered grid

grid is:
   i h ip
j  P U P

U is effectively at iph
"""
from coordinates_1d import *
from constants import *


# def advect_upwind(q, u, dx):




def advect_v_u(u, dx):
    # get staggered velocity
    u_iph = iph(u)
    u_imh = imh(u)

    du_ip = (ip(u) - u)
    du_im = (u - im(u))

    finite = (u_iph * du_ip + u_imh * du_im) / dx
    return finite




def advect_forward_euler(rho, u, dx, dt):
    rho_next = rho - advect_rho(rho, u, dx) * dt
    u_next = u
    return rho_next, u_next




# need to try formula at end of https://climate.ucdavis.edu/AOSS605-NumericalMethodsLectures.pdf


# need to copy the dang code of GCMII

def aflux(u, p, dy):
    """
    Calculate air mass fluxes
    Args:
        u: x component velocity for each layer
        p: surface pressure
        dy: change in y for each latitude

    Returns:
        pu: flux at x edge
        conv: horizontal mass convergence
        pit: pressure tendency
    """
    # SPA(I,J,L)=U(I,J,L)+U(I,J+1,L) # B grid code to get u in c grid position
    # PU(I,J,L)=.25*DYP(J)*SPA(I,J,L)*(P(I,J)+P(IP1,J)) # What is the dy for? turn it into FV scheme?
    pu = iph(p) * u * dy

    # pu is at iph, j

    # CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))*DSIG(L)
    conv = im(pu) - pu

    # EQUIVALENCE (CONV,PIT) # pit starts as first layer of conv[0]
    # PIT(I,J)=PIT(I,J)+CONV(I,J,L) # pit is the sum of conv
    pit = conv  # this is fine for a single layer

    return pu, conv, pit


def advecm(p, pit, dt, area):
    """
    Calculate new surface pressures - Advect mass
    Args:
        p: surface pressure
        pit: pressure tendency
        dt: change in time
        area: area of grid squares

    Returns:
        pa: New surface pressure
    """
    # PA(I,J)=P(I,J)+(DT1*PIT(I,J)/DXYP(J)+PTRUNC)
    pa = p + (dt * pit / area)
    return pa


def scaling(pa, q, dx):
    """
    Scale the quantity by the pressure in the grid cell
    Args:
        pa: current surface pressure
        q: quantity to scale
        dx:

    Returns:
        qq: the scaled quantity
    """
    qq = pa * q * dx * dx
    return qq


def unscaling(pb, qq, dx):
    """
    Undo the scaling of the quantity using the new pressure for the grid cell
    Args:
        pb: new pressure
        qq: scaled quantity
        dx:

    Returns:
        q: the unscaled quantity
    """
    q = qq / (pb * dx * dx)
    return q


def advecv(ut, pu, p, pa, u, dt, dx):
    # TODO: add u_prev so matsuno can be done
    """
    i  h  ip   h
       u  uph  up
       pu

    """
    # FD(I,J)=PA(I,J)*DXYP(J)
    # UT(I,J,L)=UT(I,J,L)*FDU # fd at u
    ut_s = scaling(p, ut, dx)

    # pu is at iph, j
    """
    Wait, their grid is backwards, with positive j towards north
    so U is at  jmh, not jph
        i  iph ip  i2h
    jm     4       2
    jmh    u   uph,f
    j      pu       1
    jph    up
    jp
    j2h
    """
    # FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L)) # b grid to get pu at the right spot to advect b grid u
    # FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))  # i is their x coordinate

    # (U(IP1,J,L)+U(I,J,L))
    uph = iph(u)

    # getting pu at the right spot to compute fluxes of u
    puph = iph(pu)

    fluxu = dt * puph * uph  # I could move the dt into the next calc
    # DUT(I,J+1,L)=DUT(I,J+1,L)+FLUXU
    # DUT(I,J,L)=DUT(I,J,L)-FLUXU

    dut = (im(fluxu) - fluxu)  # dt that feels like it should be here is in fluxu calc
    # UT(I,J,L)=(UT(I,J,L)+DUT(I,J,L))*(1/PB(I,J)*DXYP(J))
    ut_next = ut_s + dut
    u_next = unscaling(pa, ut_next, dx)

    return u_next


def thbar(t1, t2):
    """
    FUNCTION THBAR (X,Y)
c  ** TH-mean used for vertical differencing (Arakawa)
c  ** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
c  **              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
c  **      g(x) is replaced by a rational function
c  **           (a+bx+cxx+dxxx+cx**4)/(e+fx+gxx)
c  **      approx.error <1.E-6 for x between .9 and 1.7
c  **
    """
    x = t1 / t2
    g = np.log(x) / (x-1)
    return t1 * g


def pgf(u, p, pa, t, dt, dx):
    # SHA=RGAS/KAPA
    sha = Rd / kappa
    # KAPAP1=KAPA+1.
    kappa_p1 = kappa + 1

    # SP=P(I,J)
    sp = p
    # PDN=SIG(1)*SP+PTOP # ptop is pressure at top of atmosphere
    # PUP=SIG(LP1)*SP+PTOP
    # sig is the sigma of the layer - sigma coordinates of middle of layer. will be 1 for this.
    # pdn is the pressure of bottom of the layer?
    # ptop is pressure at 0 sigma, top of model. For this test it will be zero
    # PKUP=EXPBYK(PUP)
    pdn = sp
    pkdn = pdn.m ** kappa  # strip units to deal with kappa - see notes below
    pup = p_mesopause
    pkup = pup.m ** kappa  # same

    # This is wonky because they use a different formulation of potential temperature.
    """
    earlier you have this: 
    REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
    T(I,J,L)=T(I,J,L)/EXPBYK(SIG(L)*P(I,J)+PTOP)
    
    They use P_0 = 1, which lets them ignore that calculation
    """
    # SPA(I,J,L)=SIG(L)*SP*RGAS*T(I,J,L)*PKDN/PDN # is this computing that specific gravity?
    # this is the scaling term for the pressure gradient force
    # usually, this is pressure over density, but we have to compute that density at our pressure level.

    spa = 1 * sp * Rd * t * pkdn / pdn

    # THETA=THBAR(T(I,J,LP1),T(I,J,L))
    theta = thbar(t, t_mesopause)
    # SUM1=SUM1+SPA(I,J,L)*DSIG(L)
    # PHI(I,J,LP1)=SHA*THETA*(PKDN-PKUP)
    # is this layer thickness? because we then add on the height map
    # yeah, I think this is where we compute geopotential
    # it looks like this is just computing the difference between the layers.
    # for our uses with a single layer, it's good enough.
    phi = sha * theta * (pkdn - pkup)
    # SUM1=SUM1+SPA(I,J,L)*DSIG(L)
    # SUM2=SUM2+SIGE(LP1)*PHI(I,J,LP1)

    # PHI(I,J,1)=FDATA(I,J,1)+SUM1-SUM2
    # first layer phi is geopotential plus density sum minus layer thickness sum



    # PU(I,J,L)=(P(IP1,J)+P(I,J))*(PHI(IP1,J,L)-PHI(I,J,L))
    #            +(SPA(IP1,J,L)+SPA(I,J,L))*(P(IP1,J)-P(I,J))
    # for some ungodly reason they're reusing pu
    # probably because PGF is done after all the other advection code in the main loop
    # this is a mixture of the geopotential and the pgf terms

    # don't divide by dx here, look below for why
    dp = (ip(p) - p)
    dphi = (ip(phi) - phi)
    geo = iph(p) * dphi
    pg = iph(spa) * dp



    # DUT(I,J,L)=DUT(I,J,L)-DT4*DYV(J)*(PU(I,J,L)+PU(I,J-1,L)) # putting B pu at u post
    dut = (geo + pg) * dt * dx

    # oh, they don't divide by dx because it would be canceled out by multiplying by the area.
    # that's why they multiply by dy for the difference term

    # FD(I,J)=PB(I,J)*DXYP(J)
    # RFDU=4./(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
    # UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)*RFDU
    # at the end they just scale the acceleration by

    paph = iph(pa)

    # u_next = u + dut / paph #unscaling(paph, dut, dx)
    u_next = u + unscaling(paph, dut, dx)

    return spa, theta, phi, geo, pg, u_next


def advect(pu, pa, tt, pb, t, dt, dx):
    # FD(I,J)=PA(I,J)*DXYP(J)
    # TT(I,J,L)=TT(I,J,L)*FD(I,J)
    tt_s = scaling(pa, tt, dx)

    # FLUXQ(I)=DT2*PU(I,J,L)*(T(I,J,L)+T(IP1,J,L))
    fluxq = pu * iph(t) * dt

    # TT(I,J,L)=TT(I,J,L)+(FLUXQ(IM1)-FLUXQ(I))
    tt_s_next = tt_s + im(fluxq) - fluxq

    # UNDO SCALING PERFORMED AT BEGINNING OF DYNAM
    tt_next = unscaling(pb, tt_s_next, dx)
    return tt_next


def advecq(pu, pa, qt, pb, q, dt, dx):
    # SCALE QT.
    qt_s = scaling(pa, qt, dx)

    # FLUXQ(I)=DT2*PU(I,J,L)*(Q(I,J,L)+Q(IP1,J,L))
    fluxq = pu * iph(q) * dt

    # IF(FLUXQ(I).GT..5*QT(I,J,L)) FLUXQ(I)=.5*QT(I,J,L)
    # IF(FLUXQ(I).LT.-.5*QT(IP1,J,L)) FLUXQ(I)=-.5*QT(IP1,J,L)
    # limit flux
    half = qt_s / 2

    fluxq_limited = np.maximum(np.minimum(fluxq, half), -ip(half)) * fluxq.u

    # QT(I,J,L)=QT(I,J,L)+(FLUXQ(IM1)-FLUXQ(I))
    qt_s_next = qt_s + im(fluxq_limited) - fluxq_limited

    qt_next = unscaling(pb, qt_s_next, dx)

    return qt_next


def dynam_matsuno(u, p, t, q, dt, dx):
    pu, conv, pit = aflux(u, p, dx)

    pa = advecm(p, pit, dt, dx * dx)

    u_next = advecv(u, pu, p, pa, u, dt, dx)
    t_star = advect(pu, p, t, pa, t, dt, dx)
    q_star = advecq(pu, p, q, pa, q, dt, dx)
    spa, theta, phi, geo, pg, u_star = pgf(u_next, p, pa, t, dt, dx)
    p_star = pa

    # second pass
    pu, conv, pit = aflux(u_star, p_star, dx)
    pa = advecm(p, pit, dt, dx * dx)

    u_next = advecv(u, pu, p, pa, u_star, dt, dx)
    t_next = advect(pu, p, t, pa, t_star, dt, dx)
    q_next = advecq(pu, p, q, pa, q_star, dt, dx)
    spa, theta, phi, geo, pg, u_next = pgf(u_next, p_star, pa, t, dt, dx)
    p_next = pa

    return u_next, p_next, t_next, q_next

# ok, flux form should work great. input - output.
# flux form of drho/dt
def advect_rho(rho, u, dx):
    rho_ph = iph(rho)

    urho = u * rho_ph

    finite = (im(urho) - urho) / dx
    return finite


def advect_u_scaled(ut, u, p, pa, dt, dx):
    # TODO fix this
    ut_scaled = iph(p) * ut

    u_at_h = imh(u)

    adv_val = u_at_h * u_at_h * p

    adv_diff = (adv_val - ip(adv_val)) / dx

    geo_val = p * p * G / 2

    geo_diff = (geo_val - ip(geo_val)) / dx

    output = ((ut * p) - (adv_diff + geo_diff) * dt) / pa

    return output

    uph = iph(u)
    uph_s = iph(ut_scaled)
    """
    These are cells centered on U
    u is scaled, then advected.
    """
    suu = uph * uph_s

    dut = (suu - im(suu)) / dx * dt

    grad = 1/2 * G * p * p

    dgrad = (grad - ip(grad)) / dx * dt

    ut_s_next = ut_scaled - dut + dgrad
    u_next = ut_s_next / iph(pa)
    return u_next




def advect_matsumo(rho, u, dt, dx):
    rho_star = rho - advect_rho(rho, u, dx) * dt
    u_star = u
    rho_next = rho - advect_rho(rho_star, u_star, dx) * dt
    u_next = u

    return rho_next, u_next


def shallow_water_matsuno(h, u, dt, dx):
    h_star = h - advect_rho(h, u, dx) * dt
    u_star = advect_u_scaled(u, u, h, h_star, dt, dx)
    u_star[-1] = 0 * u.u
    h_next = h - advect_rho(h_star, u_star, dx) * dt
    u_next = advect_u_scaled(u, u_star, h_star, h_next, dt, dx)
    u_next[-1] = 0 * u.u

    return h_next, u_next


def advect_maccormack(rho, u, dt, dx):
    rho_star = rho - advect_rho(rho, u, dx) * dt
    u_star = u
    rho_next = ((rho + rho_star) - advect_rho(rho_star, u_star, dx) * dt) / 2
    u_next = u

    return rho_next, u_next


def shallow_water_maccormack(h, u, dt, dx):
    # https://en.wikipedia.org/wiki/MacCormack_method
    pass
















