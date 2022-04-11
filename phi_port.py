import numpy as np
import constants


def PGF(T, P, geom):
    """
C**** THIS SUBROUTINE ADDS TO MOMENTUM THE TENDENCIES DETERMINED BY     2503.
C**** THE PRESSURE GRADIENT FORCE                                       2504.
    """

    # localizing the missing arrays and data structures
    PTOP = geom.ptop
    RGAS = constants.Rd
    KAPA = constants.kappa
    # T and P must be transposed first
    IM, JM, LM = T.shape

    SIG = geom.sig.flatten()
    DSIG = geom.dsig.flatten()
    SIGE = geom.sige.flatten()

    FDATA = np.transpose(geom.heightmap)

    SHA = RGAS/KAPA

    # get the units right for this
    PDN = SIG[0] * P[0, 0] + PTOP
    PKDN = EXPBYK(PDN)
    SPA = np.zeros_like(T) * (SIG[0] * P[0, 0] * RGAS * T[0, 0, 0] * PKDN / PDN).u
    PHI = np.zeros_like(T) * (SHA * THBAR(T[0, 0, 0], T[0, 0, 0]) * (PKDN - PKDN)).u
    PHI_final = np.zeros_like(T)

    """
C****                                                                   2521.
C**** VERTICAL DIFFERENCING                                             2522.
C****                                                                   2523.
      IMAX=1                                                            2524.
      DO 3071 J=1,JM                                                    2525.
      IF(J.EQ.JM) IMAX=1                                                2526.
      DO 3070 I=1,IMAX                                                  2527.
      SUM1=0.                                                           2528.
      SUM2=0.                                                           2529.
      SP=P(I,J)                                                         2530.
      PDN=SIG(1)*SP+PTOP                                                2531.
      PKDN=EXPBYK(PDN)                                                  2532.
      DO 3040 L=1,LM-1                                                  2533.
      LP1=L+1                                                           2534.
    """

    IMAX = 1
    for J in range(JM):
        if J == JM - 1:
            IMAX = 1
        for I in range(IMAX):
            SUM1 = 0.
            SUM2 = 0.
            SP = P[I, J]
            # switching index 1 to 0 for numpy array
            PDN = SIG[0] * SP + PTOP
            PKDN = EXPBYK(PDN)
            for L in range(LM - 1):
                LP1 = L + 1
                """
C**** CALCULATE SPA                                                     2535.
      SPA(I,J,L)=SIG(L)*SP*RGAS*T(I,J,L)*PKDN/PDN                       2536.
      SUM1=SUM1+SPA(I,J,L)*DSIG(L)                                      2537.
      PUP=SIG(LP1)*SP+PTOP                                              2538.
      PKUP=EXPBYK(PUP)                                                  2539.
      THETA=THBAR(T(I,J,LP1),T(I,J,L))                                  2540.
                """
                # getting magnitude to deal with units stuff
                SPA[I, J, L] = (SIG[L] * SP * RGAS * T[I, J, L] * PKDN / PDN)
                SUM1 = SUM1 + SPA[I, J, L] * DSIG[L]
                PUP = SIG[LP1] * SP + PTOP
                PKUP = EXPBYK(PUP)
                # this goes NAN with equal t. Setting it to the average for now
                # THETA = THBAR(T[I, J, LP1], T[I, J, L])
                THETA = (T[I, J, LP1] + T[I, J, L]) / 2
                """
C**** CALCULATE THE DIFFERENCE IN PHI BETWEEN ODD LEVELS LP1 AND L      2541.
      PHI(I,J,LP1)=SHA*THETA*(PKDN-PKUP)                                2542.
      SUM2=SUM2+SIGE(LP1)*PHI(I,J,LP1)                                  2543.
      PDN=PUP                                                           2544.
 3040 PKDN=PKUP                                                         2545.
                """
                PHI[I, J, LP1] = SHA * THETA * (PKDN - PKUP)
                SUM2 = SUM2 + SIGE[LP1] * PHI[I, J, LP1]
                PDN = PUP
                PKDN = PKUP
            """
      SPA(I,J,LM)=SIG(LM)*SP*RGAS*T(I,J,LM)*PKDN/PDN                    2546.
      SUM1=SUM1+SPA(I,J,LM)*DSIG(LM)                                    2547.
 3050 PHI(I,J,1)=FDATA(I,J,1)+SUM1-SUM2                                 2548.
            """
            # adding -1 for 0 indexing
            SPA[I, J, LM-1] = SIG[LM-1] * SP * RGAS * T[I, J, LM-1] * PKDN / PDN

            # ok, so SPA is the sigma*pi/rho part of the PGF term.
            # They don't use a real reference temp for their virtual temperature

            SUM1 = SUM1 + SPA[I, J, LM-1] * DSIG[LM-1]
            # got rid of the index to fit with geom height map
            # also switched to zero indexing
            PHI[I, J, 0] = FDATA[I, J] + SUM1 - SUM2
            """
          DO 3070 L=2,LM                                                    2549.
     3070 PHI(I,J,L)=PHI(I,J,L)+PHI(I,J,L-1)                                2550.
            """
            # switched to zero indexing
            for L in range(1, LM):
                PHI[I, J, L] = PHI[I, J, L] + PHI[I, J, L - 1]

    return PHI


def EXPBYK(X):
    return X ** constants.kappa


def THBAR(X, Y):
    """
c  **
c  ** TH-mean used for vertical differencing (Arakawa)
c  ** THBAR(T1,T2) = (ln(T1) - ln(T2))/(1/T2 - 1/T1)
c  **              = T1*g(x) with x=T1/T2 , g(x)=ln(x)/(x-1)
c  **      g(x) is replaced by a rational function
c  **           (a+bx+cxx+dxxx+cx**4)/(e+fx+gxx)
c  **      approx.error <1.E-6 for x between .9 and 1.7
c  **
    """

    # we live in the 21st century, we can just do the damn math
    def g(x):
        return np.log(x) / (x - 1)

    return X * g(X/Y)
