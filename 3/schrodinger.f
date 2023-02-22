      program schrodinger
      
      implicit none
      
      integer Np,i,j,k,t
      real*8 h,s,st,k0t,lambda,L,pi,sum,nciclos
      real*8 Vt(0:1000),espx(0:1000),espp(0:1000),espx2(0:1000)
      real*8 espp2(0:1000),deltax(0:1000),deltap(0:1000),E(0:1000)
      complex*16 phi(0:1000,0:1000),chi(0:1000,0:1000)
      complex*16 alpha(0:1000),beta(0:1000,0:1000),omega(0:1000)
      complex*16 icomp
      parameter(icomp=(0.d0,1.d0),pi=dacos(-1.d0)) 
      !tomamos la longitud de la red como la unidad
      
      open(12,file="prob.dat")
      open(13,file="re_im.dat")
      open(14,file="espx_espp.dat")
          
      write(*,*) "Introduzca el valor de discretizacion espacial (N)"
      read(*,*) Np
      
      write(*,*) "Introduzca el factor de proporcionalidad
     & de la altura del potencial (lambda*k0²)"
      read(*,*) lambda
      
      write(*,*) "Introduzca el numero de ciclos en L"
      read(*,*) nciclos
      
      chi=0.d0
      phi=0.d0 !inicializando incluimos las condiciones de contorno phi(0,s)=phi(N,s)=0
      alpha=0.d0
      beta=0.d0
      omega=0.d0
      Vt=0.d0
      espx=0.d0
      espp=0.d0
      espx2=0.d0
      espp2=0.d0
      deltax=0.d0
      deltap=0.d0
      E=0.d0
            
C     CALCULO A PARTIR DE LA N EL VALOR DEL PASO ESPACIAL
      h=1.d0
      L=h*dfloat(Np)
      
C     PRIMEROS CALCULOS
      k0t=(2.d0*pi*nciclos)/(dfloat(Np))
      st=1.d0/(4.d0*(k0t**2.d0))      

      
C     ESTABLECEMOS EL VALOR FIJO DE V(j)
      do i=0,Np
      
       if((dfloat(i).lt.(2.d0*dfloat(Np)/5.d0)).OR.
     & (dfloat(i).gt.(3.d0*dfloat(Np)/5.d0))) then
     
        Vt(i)=0.d0
     
       else
        
        Vt(i)=lambda*(k0t**2.d0)
       
       end if
       
c       write(*,*) Vt(i),i
      
      end do 
      
C     AHORA LA FUNCIÓN DE ONDA EN TODAS LAS POSICIONES EN EL INSTANTE INICIAL (menos j=0 y j=N)
      do j=1,Np-1


        !FUNCION DE ONDA DEL GUION
       phi(j,0)=exp(icomp*k0t*dfloat(j))*
     & dexp((-8.d0*(dfloat(4*j-Np))**2.d0)/(Np**2.d0))

        !FUNCION DE ONDA PRUEBA
C        phi(j,0)=dsin(k0t*j)

     
c      write(*,*) phi(j,0)
      
      end do
      
c     Subrutina para normalizar
      call norm(phi,0,Np)
      
      
C     CALCULAMOS LAS ALPHAS INDEPENDIENTES DEL TIEMPO (fijas)
      do i=Np-2,0,-1 !empezamos en N-2 ya que una condición es que alpha(N-1) es 0 y ya está inicializado a ese valor
      
       alpha(i)=-1.d0/(-2.d0+(2.d0*icomp)/st-Vt(i+1)+alpha(i+1))
c       write(*,*) alpha(i),i
       
      end do
      
C     AHORA FIJAMOS LAS OMEGAS
      do j=0,Np
       omega(j)=1.d0/(-2.d0+2.d0*icomp/st-Vt(j)+alpha(j))
C       write(*,*) omega(j)
      end do
      
      
      
c  COMENZAMOS EL BUCLE*********

      do t=0,999
            
       call evolucion(phi,omega,beta,alpha,chi,t,Np,st)
          
      end do 
      
      
C     IMPRIMIMOS LAS POSIBILIDADES
      do t=0,1000
       do j=0,Np
       
        write(12,*) j,abs(phi(j,t))**2      
      
       end do
       !dejamos dos renglones en blanco
       write(12,*)
       write(12,*) 
      end do
      
C     AHORA LAS PARTES REALES E IMAGINARIAS
      do t=0,1000
       do j=0,Np
       
        write(13,*) j,real(phi(j,t)),imag(phi(j,t))      
      
       end do
       !dejamos dos renglones en blanco
       write(13,*)
       write(13,*) 
      end do
      
C     CALCULO TODOS LOS VALORES ESPERADOS DE LA POSICION Y MOMENTO PARA TODO T
      write(14,*) "#t|espx|espp|espx2|espp2|deltax|deltap|E" 
      do t=0,1000
       do j=1,Np-1 !la funcion de onda en 0 y N es 0, asi que no contribuye
       
        espx(t)=espx(t)+conjg(phi(j,t))*j*phi(j,t)
        espp(t)=espp(t)-icomp*conjg(phi(j,t))*(phi(j+1,t)-phi(j,t)) !aprox. de derivada tomando h=1 ( fi(j+1,t)-fi(j,t) )
        espx2(t)=espx2(t)+conjg(phi(j,t))*j*j*phi(j,t)
        espp2(t)=espp2(t)-conjg(phi(j,t))*
     &  (phi(j+1,t)-2*phi(j,t)+phi(j-1,t))
         
       
       end do
       
       E(t)=espp2(t) !SOLO CUANDO V=0 Y CON EL REESCALAMIENTO USADO
       deltax(t)=dsqrt(espx2(t)-espx(t)**2.d0)
       deltap(t)=dsqrt(espp2(t)-espp(t)**2.d0)
       
       write(14,*) t,espx(t),espp(t),espx2(t),espp2(t),
     & deltax(t),deltap(t),E(t)
       
      end do       
      
C     VAMOS A COMPROBAR SI SE HA CONSERVADO LA NORMA
      sum=0.d0
      do j=0,Np
        sum = sum + abs(phi(j,1000))**2.d0
      end do
      
      write(*,*) "*************************"
      write(*,*) "Cte de normalizacion en t=1000:",sum
      write(*,*) "*************************" 
      
       
      close(12)
      close(13)
      close(14)
      
      stop
      end
      
      subroutine evolucion(phi,omega,beta,alpha,chi,t,Np,st)
      integer t,j,Np
      real*8 st
      complex*16 phi(0:1000,0:1000),chi(0:1000,0:1000)
      complex*16 alpha(0:1000),beta(0:1000,0:1000),omega(0:1000)
      complex*16 icomp
      
      icomp=(0.d0,1.d0)
      
      do j=Np-2,0,-1 !aqui empezamos otra vez por N-2 
       beta(j,t)=omega(j+1)*(((4.d0*icomp*phi(j+1,t))/st)-beta(j+1,t))
c        write(11,*) j,t,beta(j,t)
       end do
      
c      LAS CHIS
       do j=1,Np
        chi(j,t)=alpha(j-1)*chi(j-1,t)+beta(j-1,t)
c        write(11,*) j,t,chi(j,t)
       end do
       
C      YA PUEDO CALCULAR LAS FUNCIONES DE ONDA
       do j=1,Np-1
        phi(j,t+1)=chi(j,t)-phi(j,t)
c        write(11,*) j,t,phi(j,t)
       end do 
      
      return
      end
      
      
      
      
      subroutine norm(phi,t,Np)
      complex*16 phi(0:1000,0:1000)
      real*8 sum
      integer t,Np,j
      
      sum=0.d0
       do j=0,Np
         sum = sum + abs(phi(j,t))**2.d0
       end do
       
       do j=0,Np
         phi(j,0)=phi(j,t)/dsqrt(sum)
       end do
      
      return 
      end
