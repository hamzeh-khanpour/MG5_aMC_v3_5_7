c/* ********************************************************* */
c/*  Equivalent Photon Approximation (EPA) Functions         * */
c/*  Using Improved Weizsäcker-Williams Formula             * */
c/*  Reference: V.M. Budnev et al., Phys.Rep. 15C (1975) 181 * */
c/* ********************************************************* */
c/*   Provided by Tomasz Pierzchala - UCL                   * */
c/* ********************************************************* */

      real*8 function epa_electron(x, q2max)
      implicit none
      integer i
      real*8 x, alpha, xin, q2min, PI, f
      real*8 q2max

      ! Constants
      data PI /3.14159265358979323846/   ! Value of π
      data xin /0.511d-3/                ! Electron mass in GeV
      alpha = 0.0072992701               ! Fine-structure constant

C     // x = omega/E = (E-E')/E
      if (x < 1.d0) then
         ! Compute minimum virtuality for the photon
         q2min = xin * xin * x * x / (1.d0 - x)

         ! Check if q2min < q2max, otherwise set flux to zero
         if (q2min < q2max) then
            f = alpha / (2.d0 * PI) * &
     &          (2.d0 * xin * xin * x * (-1.d0 / q2min + 1.d0 / q2max) + &
     &           (2.d0 - 2.d0 * x + x * x) / x * dlog(q2max / q2min))
         else
            f = 0.d0
         endif
      else
         f = 0.d0
      endif

C     Ensure the photon flux is non-negative
      if (f < 0.d0) f = 0.d0
      epa_electron = f

      return
      end

      real*8 function epa_proton(x, q2max)
      implicit none
      integer i
      real*8 x, alpha, xin, qz, qmi, PI, f
      real*8 q2max
      real*8 phi_f

      ! Constants
      data PI /3.14159265358979323846/   ! Value of π
      data xin /0.938d0/                 ! Proton mass in GeV
      alpha = 0.0072992701               ! Fine-structure constant
      qz = 0.71                          ! Proton structure parameter (GeV^2)

C     // x = omega/E = (E-E')/E
      if (x < 1.d0) then
         ! Compute minimum virtuality for the photon
         qmi = xin * xin * x * x / (1.d0 - x)

         ! Check if qmi < q2max, otherwise set flux to zero
         if (qmi < q2max) then
            f = alpha / PI * &
     &          (phi_f(x, q2max / qz) - phi_f(x, qmi / qz)) * (1.d0 - x) / x
         else
            f = 0.d0
         endif
      else
         f = 0.d0
      endif

C     Ensure the photon flux is non-negative
      if (f < 0.d0) f = 0.d0
      epa_proton = f

      return
      end

      real*8 function phi_f(x, qq)
      implicit none
      real*8 x, qq, y, qq1, f
      real*8 a, b, c

      ! Parameters for the proton structure function
      a = 7.16
      b = -3.96
      c = 0.028

      ! Derived variables
      qq1 = 1.d0 + qq
      y = x * x / (1.d0 - x)

      ! Compute the structure function φ_f(x, qq)
      f = (1.d0 + a * y) * &
     &    (-dlog(qq1 / qq) + 1.d0 / qq1 + 1.d0 / (2.d0 * qq1 * qq1) + &
     &     1.d0 / (3.d0 * qq1 * qq1 * qq1))
      f = f + (1.d0 - b) * y / (4.d0 * qq * qq1 * qq1 * qq1)
      f = f + c * (1.d0 + y / 4.d0) * &
     &    (dlog((qq1 - b) / qq1) + b / qq1 + &
     &     b * b / (2.d0 * qq1 * qq1) + &
     &     b * b * b / (3.d0 * qq1 * qq1 * qq1))

      phi_f = f

      return
      end
