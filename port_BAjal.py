import numpy as np

# From BA94jalC9.COM
IM = 36
JM = 24
LM = 9
KTD = 6
KAIJ = 80
KAJK = 50


"""
C**** THERE ARE 100 INTEGER PARAMETERS IN COMMON (JC-ARRAY)
      COMMON /IPARMB/IM0,JM0,LM0,JMM1,LMM1,   LS1,LTM,LBLM,LMCM,LSSM,
     *  KOCEAN,KDISK,KEYCT,KACC0,KCOPY,  IRAND,IJRA,MFILTR,NDYN,NCNDS,
     *  NRAD,NSURF,NGRND,NFILTR,NDAA,   NDA5D,NDA5K,NDA5S,NDA4,NDASF,
     *  MLAST,MDYN,MCNDS,MRAD,MSURF,    MDIAG,MELSE,MODRD,MODD5K,MODD5S,
     *  IYEAR,IDAY,IDAY0,JYEAR,JYEAR0,  JDAY,JDATE,JDATE0,NSTEP,MRCH,
     *  IDUM(4),NDZERO(13),NDPRNT(13),  IJD6(2,4),IDACC(12)
"""
IM0 = IM
JM0 = JM
LM0 = LM
JMM1 = 0
LMM1 = 0
LS1 = 8
LTM = 0
LBLM = 2


NDYN = 0


"""
      COMMON /RPARMB/
     *  TAU,TAU0,TOFDAY,TOFDY0,DT,      TAUP,TAUI,TAUE,TAUT,TAUO,
     *  TWOPI,SDAY,LHE,LHM,LHS,         RADIUS,GRAV,RGAS,KAPA,OMEGA,
     *  CCMCX,ETA,S0X,CO2,SRCOR,        PTOP,PSF,PSL,PTRUNC,AREAG,
     *  XCDNST(2),XINT,DLAT,DLON,       SKIPSE,USESLP,USEP,USET,FIM,
     *  RSDIST,SIND,COSD,DOPK,SIG(36),SIGE(37),RDM2(44)
"""
TAU = 0
TAU0 = 0
TOFDAY = 0
TOFDY0 = 0
DT = 900.

SDAY = 86400.
LHE = 2500000.
LHM = 334000.
LHS = 2834000.
RADIUS = 6375000.
GRAV = 9.81
RGAS = 287.
KAPA = 0.286
OMEGA = 0

PTOP = 10.

SIG = np.zeros((36,))
SIGE = np.zeros((37,))
RDM2 = np.zeros((44,))

"""
      COMMON /GEOMCB/ RAPVS(JM),RAPVN(JM),RAVPS(JM),RAVPN(JM),F(JM),
     *  DXYP(JM),DXP(JM),DYP(JM),DXYS(JM),SINP(JM),LAT(JM),
     *  DXYV(JM),DXV(JM),DYV(JM),DXYN(JM),COSP(JM),COSV(JM),lat_dg(jm,2)
"""
RAPVS = np.zeros((JM,))
RAPVN = np.zeros((JM,))
RAVPS = np.zeros((JM,))
RAVPN = np.zeros((JM,))
F = np.zeros((JM,))
DXYP = np.zeros((JM,))
DXP = np.zeros((JM,))
DYP = np.zeros((JM,))
DXYS = np.zeros((JM,))
SINP = np.zeros((JM,))
LAT = np.zeros((JM,))
DXYV = np.zeros((JM,))
DXV = np.zeros((JM,))
DYV = np.zeros((JM,))
DXYN = np.zeros((JM,))
COSP = np.zeros((JM,))
COSV = np.zeros((JM,))
lat_dg = np.zeros((JM, 2))

"""
      COMMON /BNDYCB/ FDATA(IM,JM,3),ODATA(IM,JM,5),GDATA(IM,JM,14),
     *  BLDATA(IM,JM,8),VDATA(IM,JM,10),
     *  Z1O(IM,JM),Z12O(IM,JM)
"""
FDATA = np.zeros((IM, JM, 3))

"""
      COMMON /LAYACB/ DSIG(37),DSIGO(36)
"""
DSIG = np.zeros((37,))
DSIGO = np.zeros((36,))

"""
      DIMENSION U(IM,JM,LM),V(IM,JM,LM),T(IM,JM,LM),P(IM,JM),Q(IM,JM,LM)
"""
U = np.zeros((IM, JM, LM))
V = np.zeros((IM, JM, LM))
T = np.zeros((IM, JM, LM))
P = np.zeros((IM, JM))
Q = np.zeros((IM, JM, LM))





