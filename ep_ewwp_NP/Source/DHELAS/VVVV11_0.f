C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1)*Metric(3,4) -
C      P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV11_0(V1, V2, V3, V4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 TMP1
      COMPLEX*16 TMP12
      COMPLEX*16 TMP4
      COMPLEX*16 TMP7
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 V2(*)
      COMPLEX*16 V3(*)
      COMPLEX*16 V4(*)
      COMPLEX*16 VERTEX
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      TMP1 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP12 = (V3(3)*V4(3)-V3(4)*V4(4)-V3(5)*V4(5)-V3(6)*V4(6))
      TMP4 = (P1(0)*V2(3)-P1(1)*V2(4)-P1(2)*V2(5)-P1(3)*V2(6))
      TMP7 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      TMP9 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      VERTEX = COUP*TMP12*(-CI*(TMP4*TMP7)+CI*(TMP1*TMP9))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1)*Metric(3,4) -
C      P(-1,1)*P(-1,2)*Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE VVVV11_3_0(V1, V2, V3, V4, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP1
      COMPLEX*16 COUP2
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 V1(*)
      COMPLEX*16 V2(*)
      COMPLEX*16 V3(*)
      COMPLEX*16 V4(*)
      COMPLEX*16 TMP
      COMPLEX*16 VERTEX
      CALL VVVV11_0(V1,V2,V3,V4,COUP1,VERTEX)
      CALL VVVV3_0(V1,V2,V3,V4,COUP2,TMP)
      VERTEX = VERTEX + TMP
      END


