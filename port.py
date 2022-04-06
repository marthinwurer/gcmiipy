import numpy as np

from constants import units
from constants import units as u


from port_BAjal import *
# time for some good old global arrays

"""
      COMMON/WORK1/PIT(IM,JM),SD(IM,JM,LM-1),PU(IM,JM,LM)               1007.
     *      ,PV(IM,JM,LM)                                               1008.
"""
PIT = np.zeros((IM, JM))
SD = np.zeros((IM, JM, LM - 1))
PU = np.zeros((IM, JM, LM))
PV = np.zeros((IM, JM, LM))

"""
      COMMON/WORK3/PHI(IM,JM,LM),SPA(IM,JM,LM)                          1009.
"""
PHI = np.zeros((IM, JM, LM))
SPA = np.zeros((IM, JM, LM))

"""
      COMMON/WORK4/FD(IM,0:JM+1),FLUXQ(IM),DUMMYS(IM),DUMMYN(IM)        1010.
"""
FD = np.zeros((IM, JM+2)) # some weird shit here with indexing
FLUXQ = np.zeros((IM,))
DUMMYS = np.zeros((IM,))
SUMMYN = np.zeros((IM,))

"""
      COMMON/WORK5/DUT(IM,JM,LM),DVT(IM,JM,LM)                          1011.
"""
DUT = np.zeros((IM, JM, LM))
DVT = np.zeros((IM, JM, LM))


"""
      COMMON/WORK6/WORKX(IM,JM,LM),UT(IM,JM,LM),VT(IM,JM,LM),     
     *  TT(IM,JM,LM),PT(IM,JM),QT(IM,JM,LM)                      
      COMMON/WORK2/UX(IM,JM,LM),VX(IM,JM,LM),TX(IM,JM,LM),PX(IM,JM)
      DIMENSION PA(IM,JM),PB(IM,JM),PC(IM,JM)
"""
WORKX = np.zeros((IM, JM, LM))
UT = np.zeros((IM, JM, LM))
VT = np.zeros((IM, JM, LM))
TT = np.zeros((IM, JM, LM))
PT = np.zeros((IM, JM))
QT = np.zeros((IM, JM, LM))

UX = np.zeros((IM, JM, LM))
VX = np.zeros((IM, JM, LM))
TX = np.zeros((IM, JM, LM))
PX = np.zeros((IM, JM))

PA = np.zeros((IM, JM))
PB = np.zeros((IM, JM))
PC = np.zeros((IM, JM))


def stencil_generator(im):
    for i in range(im):
        yield i, i - 1, (i + 1) % im

def MOD(A, P):
    return A - (int(A / P) * P)


def main():
    """
    I think this starts at line 283 in mjal
    """

    # main loop starts at 348
    pass


def GEOM():
    """
C**** CALCULATE SPHERICAL GEOMETRY  (for 4x5 or 7.8x10)                  402.
C**** This is as in Model II'
    """


def GEOM_8x10():
    """
C**** CALCULATE SPHERICAL GEOMETRY (True 8x10)
    """


def INPUT():
    """
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE   503.
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS  504.
    """


def DYNAM():
    """
    INTEGRATE THE DYNAMICS TERMS
    starts line 1763
    """
    DTFS = DT * 2./3.
    DTLF = 2.*DT
    NDYNO = MOD(NDYN, 2)

    for L in range(LM*3+1):  # ???? Why the *3?
        for J in range(JM):
            for I in range(IM):
                UX[I, J, L] = U[I, J, L]
                UT[I, J, L] = U[I, J, L]

    for L in range(LM):
        for J in range(JM):
            for I in range(IM):
                QT[I, J, L] = Q[I, J, L]

    for J in range(JM):
        for I in range(IM):
            PA[I, J] = P[I, J]
            PB[I, J] = P[I, J]
            PC[I, J] = P[I, J]

    """
C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)                         97.   
      NS=0                                                                98.   
      MRCH=0                                                              99.   
    """
    NS = 0
    MRCH = 0

    """
C     CALL DYNAM (UX,VX,TX,PX,Q,U,V,T,P,Q,DTFS)                          100.   
      CALL AFLUX (U,V,P)
      CALL ADVECM (P,PB,DTFS)
      CALL ADVECV (P,UX,VX,PB,U,V,P,DTFS)
      CALL ADVECT (P,TX,PB,T,DTFS)
      CALL PGF (UX,VX,PB,U,V,T,P,DTFS)
    """
    AFLUX(U, V, P)
    ADVECM(P, PB, DTFS)
    ADVECV(P, UX, VX, PB, U, V, P, DTFS)
    ADVECT(P, TX, PB, T, DTFS)
    PGF(UX, VX, PB, U, V, T, P, DTFS)

    """
      IF(NDYNO.EQ.1) GO TO 320                                           101.   
C**** INITIAL BACKWARD STEP IS ODD, QT = QT + DT*F(QX)                   102.   
      MRCH=-1                                                            103.   
C     CALL DYNAM (UT,VT,TT,PT,QT,UX,VX,TX,PX,Q,DT)                       104.   
      CALL AFLUX (UX,VX,PB)
      CALL ADVECM (P,PA,DT)
      CALL ADVECV (P,UT,VT,PA,UX,VX,PB,DT)
      CALL ADVECT (P,TT,PA,TX,DT)
      CALL ADVECQ (P,QT,PA,Q,DT)
      CALL PGF (UT,VT,PA,UX,VX,TX,PB,DT)
    """


def AFLUX(U: np.ndarray, V, P):
    """
    THIS SUBROUTINE CALCULATES THE HORIZONTAL AIR MASS FLUXES
    AND VERTICAL AIR MASS FLUXES AS DETERMINED BY U, V AND P.

    starts at line 1869 in Mjal
    """
    # EQUIVALENCE(CONV, PIT)
    CONV = np.zeros_like(U)
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

    # there's some weird stuff with CONV being a union of two other arrays.
    PIT[:, :] = CONV[:, :, 0]
    SD[:, :, :] = CONV[:, :, 1:]

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
    # this means that their SD gets overwritten with the the CONV stuff
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
    # I think this is just for the top layer
    for J in range(JM):
        for I in range(IM):
            SD[I, J, LM - 1] = CONV[I, J, LM] - DSIG[LM] * PIT[I, J]

    for LX in range(LM - 1):
        L = LM - (LX + 1)
        for J in range(JM):
            for I in range(IM):
                # not quite sure how this works with the aliasing
                SD[I, J, L] = SD[I, J, L + 1] + CONV[I, J, L + 1] - DSIG[L + 1] * PIT[I, J]

    return


def ADVECM(P, PA, DT1):
    """
C**** THIS SUBROUTINE CALCULATES UPDATED COLUMN PRESSURES AS            1503.
C**** DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES.                1504.
    """

    """
C**** COMPUTE PA, THE NEW SURFACE PRESSURE                              1518.
      PA(1,1)=P(1,1)+(DT1*PIT(1,1)/DXYP(1)+PTRUNC)                      1519.
         IF(PA(1,1).GT.1150.) WRITE(6,991) L,L,MRCH,P(1,1),PA(1,1),     1520.
     *     (FDATA(1,1,K),K=1,22),(T(1,1,L),Q(1,1,L),L=1,LM)             1521.
      PA(1,JM)=P(1,JM)+(DT1*PIT(1,JM)/DXYP(JM)+PTRUNC)                  1522.
         IF(PA(1,JM).GT.1150.) WRITE(6,991) L,JM,MRCH,P(1,JM),PA(1,JM), 1523.
     *     (FDATA(1,JM,K),K=1,22),(T(1,JM,L),Q(1,JM,L),L=1,LM)          1524.
      DO 2424 I=2,IM                                                    1525.
      PA(I,1)=PA(1,1)                                                   1526.
 2424 PA(I,JM)=PA(1,JM)                                                 1527.
      DO 2426 J=2,JM-1                                                  1528.
      DO 2426 I=1,IM                                                    1529.
      PA(I,J)=P(I,J)+(DT1*PIT(I,J)/DXYP(J)+PTRUNC)                      1530.
         IF(PA(I,J).GT.1150.) WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),    1531.
     *     (FDATA(I,J,K),K=1,22),(U(I-1,J,L),U(I,J,L),U(I-1,J+1,L),     1532.
     *     U(I,J+1,L),V(I-1,J,L),V(I,J,L),V(I-1,J+1,L),V(I,J+1,L),      1533.
     *     T(I,J,L),Q(I,J,L),L=1,LM)                                    1534.
 2426 CONTINUE                                                          1535.
    """
    for J in range(JM):
        for I in range(IM):
            PA[I, J] = P[I, J] + (DT1 * PIT[I, J]/DXYP(J))
    return


def ADVECV(PA, UT, VT, PB, U, V, P, DT1):
    """
C**** THIS SUBROUTINE ADVECTS MOMENTUM (INCLUDING THE CORIOLIS FORCE)   2003.
C**** AS DETERMINED BY DT1 AND THE CURRENT AIR MASS FLUXES.             2004.
    """
    SHA = RGAS/KAPA
    KAPAP1 = KAPA + 1
    JMM2 = JM - 2
    DT2=DT1/2.
    DT4=DT1/4.
    DT8=DT1/8.
    DT12=DT1/12.
    DT24=DT1/24.

    """
C****                                                                   2026.
C**** SCALE QT.  UT AND VT MAY THEN BE PHYSICALLY INTERPRETED AS        2027.
C**** MOMENTUM COMPONENTS, TT AS HEAT CONTENT, AND QT AS WATER CONTENT  2028.
C****                                                                   2029.
      DO 101 J=2,JM-1                                                   2030.
      DO 101 I=1,IM                                                     2031.
  101 FD(I,J)=PA(I,J)*DXYP(J)                                           2032.
      FDSP=PA(1,1)*DXYP(1)                                              2033.
      FDNP=PA(1,JM)*DXYP(JM)                                            2034.
      FDSP=FDSP+FDSP                                                    2035.
      FDNP=FDNP+FDNP                                                    2036.
      DO 120 I=1,IM                                                     2037.
      FD(I,1)=FDSP                                                      2038.
  120 FD(I,JM)=FDNP                                                     2039.
      I=IM                                                              2040.
      DO 140 IP1=1,IM                                                   2041.
      DO 130 J=2,JM                                                     2042.
      FDU=.25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))                 2043.
      DO 130 L=1,LM+LM                                                  2044.
      DUT(I,J,L)=0.                                                     2045.
  130 UT(I,J,L)=UT(I,J,L)*FDU                                           2046.
  140 I=IP1                                                             2047.
    """
    for J in range(JM):
        for I in range(IM):
            FD[I, J] = PA[I, J] * DXYP(J)

    for J in range(JM):
        for I in range(IM):
            FDU = .25 * (FD[I, J] + FD[I+1, J] + FD[I, J-1] + FD[I+1, J-1])
            for L in range(LM):
                DUT[I, J, L] = 0
                UT[I, J, L] *= FDU

    """
C**** HORIZONTAL ADVECTION OF MOMENTUM                                  2051.
C****                                                                   2052.
 1400 I=IM                                                              2053.
      DO 2320 IP1=1,IM                                                  2054.
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX                         2055.
      DO 2310 J=2,JM                                                    2056.
      FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))       2057.
      FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))                                  2058.
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU                                   2059.
      DUT(I,J,L)=DUT(I,J,L)-FLUXU                                       2060.
      FLUXV=FLUX*(V(IP1,J,L)+V(I,J,L))                                  2061.
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV                                   2062.
 2310 DVT(I,J,L)=DVT(I,J,L)-FLUXV                                       2063.
    """
    for L in range(LM):
        for I in range(IM):
            IP1 = (I + 1) % IM  # I'm making this toroidal
            for J in range(JM):
                JP1 = (J + 1) % JM  # I'm making this toroidal
                FLUX = DT12 * (PU[IP1, J, L] + PU[IP1, J-1, L] + PU[I, J, L] + PU[I, J-1, L])
                FLUXU = FLUX * (U[IP1, J, L] + U[I, J, L])
                DUT[IP1, J, L] = DUT[IP1, J, L] + FLUXU
                DUT[I, J, L] = DUT[I, J, L] - FLUXU
                FLUXV = FLUX * (V[IP1, J, L] + V[I, J, L])
                DVT[IP1, J, L] = DUT[IP1, J, L] + FLUXV
                DVT[I, J, L] = DUT[I, J, L] - FLUXV
                # south-north
                # SW-NE
                # SE-NW
                # TODO the rest of these

    """
C**** VERTICAL ADVECTION OF MOMENTUM                                    2093.
C****                                                                   2094.
      DO 2470 L=1,LM-1                                                  2095.
      LP1=L+1                                                           2096.
      DO 2470 J=2,JM                                                    2097.
      I=IM                                                              2098.
      DO 2470 IP1=1,IM                                                  2099.
      SDU=DT2*((SD(I,J-1,L)+SD(IP1,J-1,L))*RAVPN(J-1)+                  2100.
     *  (SD(I,J,L)+SD(IP1,J,L))*RAVPS(J))                               2101.
      SDUDN=SDU/DSIG(L)                                                 2102.
      SDUUP=SDU/DSIG(LP1)                                               2103.
      DUT(I,J,L)  =DUT(I,J,L)  +SDUDN*(U(I,J,L)+U(I,J,LP1))             2104.
      DUT(I,J,LP1)=DUT(I,J,LP1)-SDUUP*(U(I,J,L)+U(I,J,LP1))             2105.
      DVT(I,J,L)  =DVT(I,J,L)  +SDUDN*(V(I,J,L)+V(I,J,LP1))             2106.
      DVT(I,J,LP1)=DVT(I,J,LP1)-SDUUP*(V(I,J,L)+V(I,J,LP1))             2107.
 2470 I=IP1                                                             2108.
    """
    for L in range(LM-1):
        LP1 = L + 1
        for I in range(IM):
            IP1 = (I + 1) % IM  # I'm making this toroidal
            for J in range(JM):
                SDU = DT2 * ((SD[I, J-1, L] + SD[IP1, J-1, L])*RAVPN[J-1] +
                             (SD[I, J, L] + SD[IP1, J, L])*RAVPS[J])
                SDUDN = SDU / DSIG[L]
                SDUUP = SDU / DSIG[LP1]
                DUT[I, J, L] = DUT[I, J, L] + SDUDN * (U[I, J, L] + U[I, J, LP1])
                DUT[I, J, LP1] = DUT[I, J, LP1] - SDUUP * (U[I, J, L] + U[I, J, LP1])
                DVT[I, J, L] = DVT[I, J, L] + SDUDN * (V[I, J, L] + V[I, J, LP1])
                DVT[I, J, LP1] = DVT[I, J, LP1] - SDUUP * (V[I, J, L] + V[I, J, LP1])

    """
C**** ADD ADVECTION INCREMENTS TO UT AND VT, CALL DIAGNOSTICS           2109.
         IF(MODD5K.LT.MRCH) CALL DIAG5A (4,MRCH)                        2110.
         IF(MODD5K.LT.MRCH) CALL DIAG9D (1,DT1,U,V)                     2111.
      DO 2900 L=1,LM+LM                                                 2112.
      DO 2900 J=2,JM                                                    2113.
      DO 2900 I=1,IM                                                    2114.
      UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)                                    2115.
 2900 DUT(I,J,L)=0.                                                     2116.
    """
    # I don't know why they don't actually modify VT here
    for L in range(LM):
        for I in range(IM):
            for J in range(JM):
                UT[I, J, L] = UT[I, J, L] + DUT[I, J, L]
                DUT[I, J, L] = 0

    """
    Coriolis force stuff
    """
    # TODO





def PGF(UT,VT,PB,U,V,T,P,DT1):
    """
C**** THIS SUBROUTINE ADDS TO MOMENTUM THE TENDENCIES DETERMINED BY     2503.
C**** THE PRESSURE GRADIENT FORCE                                       2504.
    """

def ADVECT(PA,TT,PB,T,DT1):
    """
C**** THIS SUBROUTINE ADVECTS POTENTIAL TEMPERATURE                     3003.
    """

def ADVECQ(PA,QT,PB,Q,DT1):
    """
C**** THIS SUBROUTINE ADVECTS HUMIDITY AS DETERMINED BY DT1 AND THE     3503.
C**** CURRENT AIR MASS FLUXES
    """

def AVRX(PU):
    """
C**** THIS SUBROUTINE SMOOTHES THE ZONAL MASS FLUX AND GEOPOTENTIAL     1803.
C**** GRADIENTS NEAR THE POLES TO HELP AVOID COMPUTATIONAL INSTABILITY. 1804.
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.    1805.
    """

def SDRAG():
    """
C**** THIS SUBROUTINE PUTS A DRAG ON THE WINDS ON THE TOP LAYER OF      7003.
C**** THE ATMOSPHERE                                                    7004.
    """

def FILTER():
    """
C**** THIS SUBROUTINE PERFORMS AN 8-TH ORDER SHAPIRO FILTER ON          7503.
C**** SELECTED PROGNOSTIC QUANTITIES IN THE ZONAL DIRECTION             7504.
C****                                                                   7505.
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER                7506.
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER    7507.
C****        3  SMOOTH P AND T                                          7508.
    """

def SHAP1D(NORDER):
    """
C**** THIS SUBROUTINE SMOOTHES THE ARRAY X IN THE ZONAL DIRECTION       7803.
C**** USING AN N-TH ORDER SHAPIRO FILTER.  N MUST BE EVEN.              7804.
C**** (USES ONLY IM,JM,JM-1,IM, AND JM FROM COMMON BLOCK)               7804.1
    """

def DAILY():
    """
C**** THIS SUBROUTINE PERFORMS THOSE FUNCTIONS OF THE PROGRAM WHICH     8003.
C**** TAKE PLACE AT THE BEGINNING OF A NEW DAY.                         8004.
    """

def CHECKT(N):
    """
C**** THIS SUBROUTINE CHECKS WHETHER THE TEMPERATURES ARE REASONABLE    9003.
C**** FOR DEBUGGING PURPOSES. IT IS TURNED ON BY SETTING IDACC(11)      9004.
C**** TO BE POSITIVE.  REMEMBER TO SET IDACC(11) BACK TO ZERO AFTER     9005.
C**** THE ERRORS ARE CORRECTED.                                         9006.
    """

















