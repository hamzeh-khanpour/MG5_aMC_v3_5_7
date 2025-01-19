ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_50 = MDL_COMPLEXI*MDL_GW*MDL_SW
      GC_56 = MDL_COMPLEXI*MDL_GW__EXP__2*MDL_SW__EXP__2
      GC_143 = (MDL_EE__EXP__2*MDL_FM0*MDL_COMPLEXI*MDL_V__EXP__2)
     $ /2.000000D+00+(MDL_CW__EXP__2*MDL_EE__EXP__2*MDL_FM2
     $ *MDL_COMPLEXI*MDL_V__EXP__2)/MDL_SW__EXP__2
      GC_144 = (MDL_EE__EXP__2*MDL_FM1*MDL_COMPLEXI*MDL_V__EXP__2)
     $ /8.000000D+00+(MDL_CW__EXP__2*MDL_EE__EXP__2*MDL_FM3
     $ *MDL_COMPLEXI*MDL_V__EXP__2)/(4.000000D+00*MDL_SW__EXP__2)
      END
