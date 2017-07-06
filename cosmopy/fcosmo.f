c     ---------------------------------------------
c     To create the python module
c     f2py -c -m fcosmo fcosmo.f -DNUMARRAY --build-dir AA
c     For nympy (Aug 2009) -- use this one!!!
c     f2py -c -m fcosmo fcosmo.f --build-dir AA
c     
c     F Menanteau, JHU Nov 2004
c     Rutgers, Dec 2011
c     ---------------------------------------------
c
c     ******************************************************************

c     Routine to compute the comoving distance DC as defined in eq. 16
c     from Hogg (1999) astro-ph 9905116 between 0 and z for a flat
c     cosmology. It use the fuction inv_ez(z) = 1/E(z) IT DOESN'T
c     INCLUDE THE c/Ho
c     
      subroutine dist_comov(z,dc,n)
      integer n,i
      real*8 z(n)
      real*8 dc(n)
      real*8 inv_ez
      external inv_ez
      real*8 sum,z1,z2
C     F2PY atributes
Cf2py intent(in) z
Cf2py intent(out) dc
Cf2py intent(hide) n
      do i = 1,n
         sum = 0.0
         z1 = 0.0
         z2 = z(i)
         call qtrap(inv_ez,z1,z2,sum)
         dc(i) = sum
      enddo
      return
      end


c     ******************************************************************
c     Routine to compute the integral of the lookback time from eq 30
c     from Hogg (1998) astro-ph 9905116 between 0 and z. It uses the
c     fuction dt(z) = 1/(1+z)E(z) 
c     IT DOESN'T INCLUDE THE 1/Ho factor
c     
      subroutine lkbt(z,t,n)
      integer n,i
      real*8 z(n)
      real*8 t(n)
      real*8 dt
      external dt
      real*8 sum,z1,z2
C     F2PY atributes
Cf2py intent(in) z
Cf2py intent(out) t
Cf2py intent(hide) n
      do i = 1,n
         sum = 0.0
         z1 = 0.0
         z2 = z(i)
         call qtrap(dt,z1,z2,sum)
         t(i) = sum
      enddo
      return
      end

c     ******************************************************************
c     ROUTINE to compute the integral of the time for eq 7.79 from
c     Weimberg (1972) between z and oo. It integrates the fuction
c     dt(z)=1/(1+z)E(z) 
c     It uses the function age_z over a do-loop
c     IT DOESN'T INCLUDE THE 1/Ho factor
c     
      subroutine time(z,t,n)
      integer n,i
      real*8 z(n)
      real*8 t(n)
      real*8 age_z
      external age_z
C     F2PY atributes
Cf2py intent(in) z
Cf2py intent(out) t
Cf2py intent(hide) n
      do i = 1,n
         t(i) = age_z(z(i))
      enddo
      return
      end


c     ******************************************************************
c     FUNCTION to compute the integral of the time for eq 7.79 from
c     Weimberg (1972) between z and oo. It integrates the fuction
c     dt(z)=1/(1+z)E(z) 
c     IT DOESN'T INCLUDE THE 1/Ho factor
c     
      real*8 function age_z(z)
      real*8 z
      real*8 dt
      external dt
      real*8 sum,z1,z2
      real*8 A,B,C
      real*8 Om,OL,Ok
      common /pars/ OL,Om,Ok
C     F2PY atributes
Cf2py intent(in) z
      sum = 0.0
      z1 = z
      z2 = 1000.0
      if (Ok.eq.0.and.OL.ne.0) then
         A = sqrt(Om*(1.+ z)**3 + OL) + sqrt(OL)
         B = sqrt(Om*(1.+ z)**3)
         C = 2./3./sqrt(OL)
         age_z = C*log(A/B)
      else
         call qtrap(dt,z1,z2,sum)
         age_z = sum
      endif
      return
      end


c     *******************************************************************
c     To get the redhift solution for an array of ages. It call n times
c     the function zx(age) the
      subroutine get_z(age,z,n)
      integer n,i
      real*8 age(n)
      real*8 z(n)
      real*8 zx
      external zx
C     F2PY atributes
Cf2py intent(in) age
Cf2py intent(out) z
Cf2py intent(hide) n
      do i=1,n
         z(i) = zx(age(i))
      enddo
      return
      end


c     *************************************************************
c     Adapted from from IDL's red cosmology routines
c     originaly adapted from ignacio ferreras' c code
c     /* secant method used to compute the root of age-get_age */
c     /* see numerical recipes cv2 pg.354 for details */
      real*8 function zx(age)
      implicit none
      integer j,niter
      real*8 age
      real*8 age_z
      external age_z
      real*8 tu,t
      real*8 zo,x1,x2,rts
      real*8 dx,xl,f,fl,swap
      parameter (niter=50)
C     F2PY atributes
Cf2py intent(in) age
      
      zo = 0.0
      tu = age_z(zo)
      if (age.gt.tu) then 
         zx = -1.0
         return  
      endif
      if ((tu - age).lt.1.0e-4) then 
         zx = 0.0
         return
      endif

      zx = 0.1
      t = age_z(zx)
      do while (age.lt.t)
         t = age_z(zx)
         zx = zx + 0.1
      enddo
      x2 = zx
      x1 = zx - 0.1
      fl = age - age_z(x1)
      f  = age - age_z(x2)
      if (abs(fl).lt.abs(f)) then 
         rts  = x1
         xl   = x2
         swap = fl
         fl   = f
         f    = swap
      else 
         xl  = x1
         rts = x2
      endif

      do j=1,niter 
         dx  = (xl-rts)*f/(f-fl)
         xl  = rts
         fl  = f
         rts = rts + dx
         f = age - age_z(rts)
         if (abs(dx).lt.0.001.or.abs(f).lt.1.0e-7) then 
            zx = rts
            return
         endif
      enddo
      write(*,*) 'Convergence error finding z for age'
      zx = -1.0
      return 
      end 


c     *******************************************************
c     Simple function that computes 1/E(z)
c     
c     O_matter and O_lamda are common block 'pars' that must be
c     set up from python like this
c     module.pars.ol = 0.7
c     module.pars.om = 0.3
c     *******************************************************
      real*8 function inv_ez(z)
      real*8 z
      real*8 Om,OL,Ok
      common /pars/ OL,Om,Ok
Cf2py intent(in) z    
      inv_ez = 0.0
      inv_ez = 1.0/sqrt(Om * (1.0 + z)**3. + Ok*(1.+z)**2 + OL)
      return
      end


c     ****************************************************************
c     Simple function that computes 1/(1+z)E(z) the lookbatime expresion
c     at a redshift z. O_matter, O_lamda and O_kappa are common block
c     'pars' that must be set up from python like this
c     
c     module.pars.ol = 0.7
c     module.pars.om = 0.3
c     *******************************************************
      real*8 function dt(z)
      real*8 z
      real*8 Om,OL,Ok
      common /pars/ OL,Om,Ok
Cf2py intent(in) z    
      dt = (1.+ z)*sqrt(Om * (1.0 + z)**3. + Ok*(1.+z)**2 + OL)
      dt = 1.0/dt
      return
      end

      
c     ********************************************
c     Trapezoidal integration of a funcion, func
c     Taken from Numerical Recipes
c     ********************************************
      subroutine trapzd(func,a,b,n,s)
      integer n
      integer it,j
      real*8 a,b,s
      real*8 del,sum,tnm,x
      real*8 func
      external func
cf2py intent(in) func,a,b,n
cf2py intent(out) s
      s = 0.0
      if (n.eq.1) then
         s=0.5*(b-a)*(func(a)+func(b))
      else
         it=2**(n-2)
         tnm=float(it)
         del=(b-a)/tnm
         x=a+0.5*del
         sum=0.
         do 11 j=1,it
            sum=sum+func(x)
            x=x+del
 11      continue
         s=0.5*(s+(b-a)*sum/tnm)
c     Added to fix factor 2 missing (?)
         s = 2.0*s
      endif
      return
      end

c     *****************************************************************
c     A robust integration routine using trapezoidal method. Runs a
c     trapezoidal integration until it converges to the desired
c     resolution designated by eps. It's supposed to be more efficient
c     than just running TRAPZD alone. Taken from Numerical Recipes and
c     hacked to work with Python's f2py and numarray. The original QTRAP
c     routine calls the TRAPZD, but f2py doesn't like to call 'func'
c     twice inside qtrap and trapz, so we write out the routine trapzd
c     rather than calling it. It includes the factor (2) missing in
c     trapzd
c     
c     F. Menanteau JHU, Nov 2004
c     *****************************************************************
      subroutine qtrap(func,a,b,s)
      integer jmax
      integer it,j
      real*8 a,b,s,eps,olds
      real*8 del,sum,tnm,x
      real*8 func
      external func
      parameter (eps=1.0e-3, jmax=20)
cf2py intent(in) func,a,b
cf2py intent(out) s
cf2py parameter eps,jmax
      olds = -1.e30
      do j=1,jmax
c     --------------------------------------------
c     This is just the trapzd function plugged in
c     to make it work with f2py
        s = 0.0
        if (j.eq.1) then
           s=(b-a)*(func(a)+func(b))
        else
           it=2**(j-2)
           tnm=float(it)
           del=(b-a)/tnm
           x=a+0.5*del
           sum=0.0
           do i=1,it
              sum=sum+func(x)
              x=x+del
           enddo
           s=(s+(b-a)*sum/tnm)
        endif
c     --------------------------------------------
        if (abs(s-olds).lt.eps*abs(olds)) return
        if (s.eq.0..and.olds.eq.0..and.j.gt.6) return
        olds=s
      enddo
      write(*,*) '# warning:too many steps in QTRAP -- did not converge'
      return
      end

      
C c     *****************************************************************
C c     Routine to compute the luminosity distance at z for a flat
C c     cosmology and lambda not zero. The redshift (array), is the only
C c     input. It call the DM function, 1/E(z). O_matter (Om) and O_lamda
C c     (OL) and O_kappa (Ok) are a common block 'pars' that must be set
C c     up from python as defined in the inv_ez function.  IT DOESN'T
C c     INCLUDE THE c/Ho FACTOR
C c     
C       subroutine dist_lum(z,dl,n)
C       integer n,i
C       real*8 z(n)
C       real*8 dl(n)
C       real*8 inv_ez
C       external inv_ez
C       real*8 sum,z1,z2
C C     F2PY atributes
C Cf2py intent(in) z
C Cf2py intent(out) dl
C Cf2py intent(hide) n
C       do i=1,n
C          sum = 0.0
C          z1  = 0.0
C          z2  = z(i)
C c     To use the trapzd method
C c     call trapzd(inv_ez,z1,z2,12,sum)
C c     To use the qtrap method -- more efficient
C          call qtrap(inv_ez,z1,z2,sum)
C          dl(i) = (1. + z(i))*sum
C       enddo
C       return
C       end

C c     ******************************************************************
C c     Routine to compute the comoving volume between 0 and z for a flat
C c     cosmology and lambda not zero. 
C c     IT DOESN'T INCLUDE THE c/Ho FACTOR or 4*pi/3
C c     
C       subroutine vol_comov(z,v,n)
C c     ****************************
C       integer n,i
C       real*8 z(n)
C       real*8 v(n)
C       real*8 inv_ez
C       external inv_ez
C       real*8 sum,z1,z2
C C     F2PY atributes
C Cf2py intent(in) z
C Cf2py intent(out) v
C Cf2py intent(hide) n
C       do i = 1,n
C          sum = 0.0
C          z1 = 0.0
C          z2 = z(i)
C          call qtrap(inv_ez,z1,z2,sum)
C          v(i) = sum**3
C       enddo
C       return
C       end


C c     ************************************************************************
C c     Routine to compute the comoving volume element at redshift z for a
C c     flat cosmology and lambda not zero.
C c     IT DOESN'T INCLUDE THE (c/Ho)**3 FACTOR
C c     
C       subroutine dvol_comov(z,dv,n)
C c     ****************************
C       integer n,i
C       real*8 z(n)
C       real*8 dv(n)
C       real*8 dl(n)
C       real*8 inv_ez
C       external inv_ez
C C     F2PY atributes
C Cf2py intent(in) z
C Cf2py intent(out) dv
C Cf2py intent(hide) n
C       call dist_lum(z,dl,n)         
C       do i = 1,n
C          dv(i) = inv_ez(z(i))*dl(i)**2 / (1.+z(i))**2
C       enddo
C       return
C       end

C c     ******************************************************************
C c     Routine to compute the comoving volume between redshifts z1 and z2
C c     for a flat cosmology and lambda not zero.  
C c     IT DOESN'T INCLUDE THE c/Ho FACTOR or 4*pi/3
C c     
C       subroutine dvolume(z1,z2,v,n)
C c     ****************************
C       integer n,i
C       real*8 z1(n)
C       real*8 z2(n)
C       real*8 v(n)
C       real*8 inv_ez
C       external inv_ez
C       real*8 sum
C C     F2PY atributes
C Cf2py intent(in) z1,z2
C Cf2py intent(out) v
C Cf2py intent(hide) n
C       do i = 1,n
C          sum = 0.0
C          call qtrap(inv_ez,z1(i),z2(i),sum)
C          v(i) = sum**3
C       enddo
C       return
C       end

      
