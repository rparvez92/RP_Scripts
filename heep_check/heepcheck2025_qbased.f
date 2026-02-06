      program heepcheck2025
      implicit none
C     Heepcheck - updated to use hcana/simc reconstruction and modernize the fortran a little
C
C=======================================================================
C  Q-BASED UPDATE (to match hcana/recon_hcana q-frame)
C
C  recon_hcana builds a rotation with:
C    rot_to_q.SetZAxis( q_vec, eprime_vec ).Invert();
C  This implies in the q-frame:
C    +Z is along q
C    +X lies in the (q, e') plane and points toward e'
C    +Y is OUT-OF-PLANE (normal to the (q, e') plane), right-handed
C
C  Here we reproduce the same basis in Fortran:
C    zhat = unit(qvec)
C    xhat = unit( eprime - (eprime·zhat) zhat )   ! e' projected ⟂ to q
C    yhat = unit( zhat × xhat )                   ! out-of-plane
C    xhat = unit( yhat × zhat )                   ! re-orthonormalize
C
C  Missing momentum convention matched to recon_hcana:
C    p_miss = -bq  and  bq = Q - X  (target has 0 3-mom)
C    => p_miss = X - Q
C
C  Outputs:
C    pmz = p_miss · zhat   (parallel to +q)
C    pmy = p_miss · yhat   (out-of-plane, "Pmy" in recon_hcana)
C
C  Compile/run (example):
C    gfortran -O2 -Wall -o heepcheck2025_qbased heepcheck2025_qbased.f
C    ./heepcheck2025_qbased
C=======================================================================

      real*8 par(5),diff(4),derv(4,5)
      real*8 am/938.30/
      real*8 radfac/0.017453293/
      real*8 unity
      real*8 e0,pe0,the0,thq0,q0,w0,em0,pmy0,pmz0
      real*8 de,dth,dpe,dthp,dpp
      real*8 w,angle,side3,pmy,pmz,em,wrecon
      integer i,j

      write(6,600)
 600  format(' program heepcheck, version Apr 28, 2025 (q-frame)')
   1  continue
      write(6,601)
 601  format(' enter electron energy(MeV) and angle; E0=0/: stop')     
      read(5,*)e0, the0
      if(e0.lt.0.1)stop
      call heepkin(e0,the0,pe0,thq0,q0)
      w0=w(e0,pe0,the0)
      write(6,602)e0,the0,pe0,q0,thq0,w0

 602  format(/'Nominal (correct) values:'
     1/'    e_e     th_e''    e_e''     q     th_q     W'/6f8.1)
      write(6,603)
 603  format(/'variations in W, E_m and p_m for changes in parms',
     1/'     dE   dth_e''  dp_e''  dthp  dp_p  ',
     2'    dW    dE_m   dpm_z   dpm_y'/)

      call heeprecon(e0,pe0,the0,thq0,q0,w0,em0,pmy0,pmz0)

      unity=1.0
C   variation of +0.1% in beam energy
      de=0.001*e0
      call heeprecon(e0+de,pe0,the0,thq0,q0,wrecon,em,pmy,pmz)
      derv(1,1)=wrecon-w0
      derv(2,1)=em-em0
      derv(3,1)=pmz-pmz0
      derv(4,1)=pmy-pmy0
      write(6,611)unity,(derv(i,1),i=1,4)
 611  format(f7.1,28x,4f8.2)

C   variation of +1 mrad in electron angle
      dth=0.001      
      call heeprecon(e0,pe0,the0+dth/radfac,thq0,q0,wrecon,em,pmy,pmz)
      derv(1,2)=wrecon-w0
      derv(2,2)=em-em0
      derv(3,2)=pmz-pmz0
      derv(4,2)=pmy-pmy0
      write(6,612)unity,(derv(i,2),i=1,4)
 612  format(7x,f7.1,21x,4f8.2)

C   variation of +0.1% in scattered electron energy
      dpe=0.001*pe0
      call heeprecon(e0,pe0+dpe,the0,thq0,q0,wrecon,em,pmy,pmz)
      derv(1,3)=wrecon-w0
      derv(2,3)=em-em0
      derv(3,3)=pmz-pmz0
      derv(4,3)=pmy-pmy0
      write(6,613)unity,(derv(i,3),i=1,4)
 613  format(14x,f7.1,14x,4f8.2)

C   variation of +1 mrad in proton angle
      dthp=0.001
      call heeprecon(e0,pe0,the0,thq0+dthp/radfac,q0,wrecon,em,pmy,pmz)
      derv(1,4)=wrecon-w0
      derv(2,4)=em-em0
      derv(3,4)=pmz-pmz0
      derv(4,4)=pmy-pmy0
      write(6,614)unity,(derv(i,4),i=1,4)
 614  format(21x,f7.1,7x,4f8.2)

C   variation of +0.1% in proton momentum
      dpp=0.001*q0
      call heeprecon(e0,pe0,the0,thq0,q0+dpp,wrecon,em,pmy,pmz)
      derv(1,5)=wrecon-w0
      derv(2,5)=em-em0
      derv(3,5)=pmz-pmz0
      derv(4,5)=pmy-pmy0
      write(6,615)unity,(derv(i,5),i=1,4)
 615  format(28x,f7.1,4f8.2)

      write(6,610)
 610  format(/' now give your own de0, dthe, dpe, dthp, dpp;',
     1/' the units are 0.1% for momenta and 1 mrad for angles',
     2/' the input can be repeated; back to E0 with de0>=9/')     

  10  read(5,*)par
      if(par(1).gt.8.0) goto 1
      do i=1,4
         diff(i)=0.0
      enddo
      do i=1,4
         do j=1,5
            diff(i)=diff(i)+par(j)*derv(i,j)
         enddo
      enddo
      write(6,620)par, diff
 620  format(5f7.1,4f8.2)
      goto 10
      end

      real*8 function side3(x,y,t)
      implicit none
      real*8 x,y,t,arg
      real*8 radfac/0.017453293/
      arg=x**2+y**2-2.0*x*y*cos(t*radfac)
      if (arg.lt.0.0) arg = 0.0
      side3=sqrt(arg)
      end

      real*8 function angle(x,y,z)
      implicit none
      real*8 x,y,z,arg
      real*8 radfac/0.017453293/
      arg=(y**2+z**2-x**2)/(2.0*y*z)
      if (arg.lt.-1.0) arg=-1.0
      if (arg.gt. 1.0) arg= 1.0
      angle=acos(arg)/radfac
      end

      subroutine heepkin(e,the,pe,thq,q)
      implicit none
      real*8 e,the,pe,thq,q
      real*8 am/938.30/
      real*8 radfac/0.017453293/
      real*8 angle,side3

      pe=e/(1.0+2.0*e/am*sin(the*radfac/2.0)**2)
      q=side3(e,pe,the)
      thq=angle(pe,q,e)
      return
      end

      real*8 function w(e,e1,th)
      implicit none
      real*8 e,e1,th
      real*8 am/938.30/
      real*8 radfac/0.017453293/
      w=sqrt(2.0*am*(e-e1)+am**2-4.0*e*e1*sin(th*radfac/2.0)**2)
      return
      end

      subroutine heeprecon(e,pe,the,thq,q,w,em,pmy,pmz)
      implicit none
      real*8 am/938.30/
      real*8 radfac/0.017453293/
      real*8 e,pe,the,thq,q
      real*8 w,Em,pmy,pmz
      real*8 Ep,nu,q1,q2,q3

C----- q-frame / vector vars
      real*8 kvec(3),kpvec(3),qvec(3),pvec(3),pmvec(3)
      real*8 xhat(3),yhat(3),zhat(3),xtmp(3),tmp(3)
      real*8 sgn,proj
      real*8 dot3


C----- 4-vector style pieces (same as original)
      nu=e-pe
      q1=0.0D0
      q2=-pe*dsin(the*radfac)
      q3=e-pe*dcos(the*radfac)

      w=dsqrt((am+nu)**2-q1**2-q2**2-q3**2)

      Ep=dsqrt(am**2+q**2)
      Em = e+am-pe-Ep

C=======================================================================
C  recon_hcana-like q-frame basis:
C    zhat = qhat
C    xhat in (q,e') plane pointing toward e'
C    yhat = zhat x xhat  (out-of-plane, right-handed)
C=======================================================================

C-- incoming electron 3-momentum (beam along +z)
      kvec(1)=0.0D0
      kvec(2)=0.0D0
      kvec(3)=e

C-- outgoing electron 3-momentum k' (lab)
      kpvec(1)=0.0D0
      kpvec(2)=pe*dsin(the*radfac)
      kpvec(3)=pe*dcos(the*radfac)

C-- q-vector (lab)
      qvec(1)=q1
      qvec(2)=q2
      qvec(3)=q3
      call unit3(qvec,zhat)

C-- xhat = unit( e' - (e'·zhat) zhat )
      proj = dot3(kpvec,zhat)
      xtmp(1)=kpvec(1) - proj*zhat(1)
      xtmp(2)=kpvec(2) - proj*zhat(2)
      xtmp(3)=kpvec(3) - proj*zhat(3)
      call unit3(xtmp,xhat)

C-- yhat = unit( zhat × xhat )   (OUT-OF-PLANE)
      call cross3(zhat,xhat,tmp)
      call unit3(tmp,yhat)

C-- Re-orthonormalize xhat = unit(yhat × zhat)
      call cross3(yhat,zhat,tmp)
      call unit3(tmp,xhat)

C=======================================================================
C  Proton 3-momentum in LAB from magnitude q and polar angle thq.
C  (thq in this program is treated as polar angle wrt +z beam axis.)
C  Assign sign so nominal elastic gives pm≈0.
C=======================================================================
      if (qvec(2) .ge. 0.0D0) then
         sgn = 1.0D0
      else
         sgn = -1.0D0
      endif

      pvec(1)=0.0D0
      pvec(2)=sgn*q*dsin(thq*radfac)
      pvec(3)=q*dcos(thq*radfac)

C=======================================================================
C  Missing momentum convention matched to recon_hcana:
C    p_miss = X - Q
C=======================================================================
      pmvec(1)=pvec(1) - qvec(1)
      pmvec(2)=pvec(2) - qvec(2)
      pmvec(3)=pvec(3) - qvec(3)

C-- Components in q-frame (like recon):
      pmz = dot3(pmvec,zhat)
      pmy = dot3(pmvec,yhat)

      return
      end

C=======================================================================
C  Vector utilities (F77-friendly)
C=======================================================================
      subroutine cross3(a,b,c)
      implicit none
      real*8 a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

      subroutine unit3(a,ahat)
      implicit none
      real*8 a(3),ahat(3),n
      n=dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      if (n .lt. 1.0D-30) then
         write(6,*) 'ERROR: cannot normalize near-zero vector in unit3'
         stop
      endif
      ahat(1)=a(1)/n
      ahat(2)=a(2)/n
      ahat(3)=a(3)/n
      return
      end

      real*8 function dot3(a,b)
      implicit none
      real*8 a(3),b(3)
      dot3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      return
      end
