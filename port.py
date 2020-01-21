
from constants import units
from constants import units as u


JM = 16
IM = 16
LM = 16


def stencil_generator(im):
    for i in range(im):
        yield i, i - 1, (i + 1) % im


def aflux(P, U, V, PU, PV, DYP, DXV, CONV, DSIG, PIT):
    """
    THIS SUBROUTINE CALCULATES THE HORIZONTAL AIR MASS FLUXES
    AND VERTICAL AIR MASS FLUXES AS DETERMINED BY U, V AND P.
    """
    # starts at line 1869 in Mjal
    # L=LM
    for L in range(LM, 0, -1):
        """
        COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
        ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
        """
        """
 2150 DO 2154 J=2,JM-1                                                  1026.
      DO 2154 I=1,IM                                                    1027.
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)                                    1028.
      CALL AVRX (SPA(1,1,L))                                            1029.
      I=IM                                                              1030.
      DO 2166 IP1=1,IM                                                  1031.
      DO 2165 J=2,JM-1                                                  1032.
 2165 PU(I,J,L)=.25*DYP(J)*SPA(I,J,L)*(P(I,J)+P(IP1,J))                 1033.
 2166 I=IP1
        """
        # I'm pretty sure this is just calculating PU from U and P
        # with smoohthing from the call to avrx.
        for J in range(JM):
            JP1 = (J + 1) % JM  # I'm making this toroidal
            for I in range(IM):
                IP1 = (I + 1) % IM  # I'm making this toroidal
                PU[I,J,L] = .25 * DYP[J] * (U[I,J,L] + U[I,JP1,L]) * \
                            (P[I,J] + P[IP1,J])

        """
C**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX                             1035.
      IM1=IM                                                            1036.
      DO 2172 I=1,IM                                                    1037.
      DO 2170 J=2,JM                                                    1038.
 2170 PV(I,J,L)=.25*DXV(J)*(V(I,J,L)+V(IM1,J,L))*(P(I,J)+P(I,J-1))
        """
        for I in range(IM):
            for J in range(JM):
                PV[I,J,L] = .25 * DXV[J] * (V[I,J,L] + V[I-1,J,L]) * \
                            (P[I,J] + P[I,J-1])

        # Then there's some BS for the poles which I'm skipping

        """
C**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE                     1075.
      DO 1510 J=2,JM-1                                                  1076.
      IM1=IM                                                            1077.
      DO 1510 I=1,IM                                                    1078.
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))*DSIG(L) 1079.
 1510 IM1=I 
        """
        for J in range(JM):
            JP1 = (J + 1) % JM  # I'm making this toroidal
            for I in range(IM):
                IP1 = (I + 1) % IM  # I'm making this toroidal
                CONV[I,J,L] = (PU[I-1,J,L]-PU[I,J,L]+PV[I,J,L]-PV[I,JP1,L])*DSIG[L]

    """
      L=L-1                                                             1083.
      IF(L.GE.1) GO TO 2150                                             1084.
C****                                                                   1085.
C**** END OF HORIZONTAL ADVECTION LAYER LOOP                            1086.
    """
    """
C**** COMPUTE PIT, THE PRESSURE TENDENCY                                1088.
C     PIT(I,J)=CONV(I,J,1)                                              1089.
      DO 2420 LX=2,LM                                                   1090.
      L=2+LM-LX 
      DO 2420 J=2,JM-1                                                  1094.
      DO 2420 I=1,IM                                                    1095.
 2420 PIT(I,J)=PIT(I,J)+CONV(I,J,L)
    """
    # no clue why they do the L loop that way; they're just summing up CONV
    # Ok, so they're reusing the first layer of CONV for PIT to save on computation.
    # I can ignore that because I'm making it its own array.
    for L in range(LM):
        for J in range(JM):
            JP1 = (J + 1) % JM  # I'm making this toroidal
            for I in range(IM):
                IP1 = (I + 1) % IM  # I'm making this toroidal
                PIT[I,J]=PIT[I,J]+CONV[I,J,L]

    """
C**** COMPUTE SD, SIGMA DOT                                             1097.
      SD(1, 1,LM-1)=CONV(1, 1,LM)-DSIG(LM)*PIT(1, 1)                    1098.
      SD(1,JM,LM-1)=CONV(1,JM,LM)-DSIG(LM)*PIT(1,JM)                    1099.
      DO 2430 J=2,JM-1                                                  1100.
      DO 2430 I=1,IM                                                    1101.
 2430 SD(I,J,LM-1)=CONV(I,J,LM)-DSIG(LM)*PIT(I,J)                       1102.
      DO 2440 LX=2,LM-1                                                 1103.
      L=LM-LX                                                           1104.
      SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)-DSIG(L+1)*PIT(1, 1)        1105.
      SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)-DSIG(L+1)*PIT(1,JM)        1106.
      DO 2440 J=2,JM-1                                                  1107.
      DO 2440 I=1,IM                                                    1108.
 2440 SD(I,J,L)=SD(I,J,L+1)+CONV(I,J,L+1)-DSIG(L+1)*PIT(I,J)
    """
    # TODO - This is for sigma dot



















