      program cohete
      
      implicit none
      
      integer i,j,ki,n,it
      real*8 t,h,theta,pi,Rt,dTL,v,omega,cte,G,Mt,Ml
      parameter(pi=acos(-1.d0))
      real*8 f
      real*8 y(4),yp(4),k(4,4)
                        ! n,i
      
C     y(1) sera r, y(2) phi, y(3) pr, y(4) pphi
C     yp(n) sera la velocidad correspondiente a la funcion con indice n
C     las f(n) corresponderan a la ecuacion que determinan la velocidad de la funcion con indice n
      
      Rt=6.378160E6
      Mt=5.9736E24
      Ml=0.07349E24
      dTL=3.844E8
      G=6.67E-11
      omega=2.6617E-6     
      y=0.d0
      yp=0.d0
      h=0.5
      t=0.d0
      it=0
      
      open(12,file="test.dat")
      open(11,file="datos.dat")
      write(11,*) "#t|x_c|y_c|x_l|y_l|x_t|y_t|H'-omega*pfi"
      
C     ******CONDICIONES INICIALES*******
      !se dan las de r,phi,pr y phi y como consecuencia tenemos las f iniciales
      v=7903.75/(dTL*1.d0)
      theta=0.d0
      y(2)=pi/2.d0
      y(1)=Rt/dTL
      y(3)=v*dcos(theta-y(2))
      y(4)=y(1)*v*dsin(theta-y(2))
      
        cte=y(3)**2.d0*dTL**2.d0/2.d0+y(4)**2.d0/(2.d0*y(1))-             !para escribir H'-omega*Pphi inicialmente
     &  G*Mt/(y(1)*2.d0)-
     &  G*Ml/(dTl*dsqrt( (y(1)*dcos(y(2))-dcos(omega*t))**2.d0+
     &  (y(1)*dsin(y(2))-dsin(omega*t))**2.d0))
      
c     Imprimo las condiciones iniciales
      write(11,*) t,y(1)*dcos(y(2)),y(1)*dsin(y(2)),1.d0,0.d0,0.d0,0.d0
     &,(cte-omega*y(4))*1E-20

      do ki=0,300000
      
       do j=1,4 !calculo todas las k1                 
        k(j,1)=h*f(j,y(1),y(2),y(3),y(4),t)                    
       end do
       
       do j=1,4 !todas las k2
        k(j,2)=h*f(j,y(1)+k(1,1)/2.d0,y(2)+k(2,1)/2.d0,
     &  y(3)+k(3,1)/2.d0,y(4)+k(4,1)/2.d0,t+h/2.d0)
       end do
       
       do j=1,4 !todas las k3
        k(j,3)=h*f(j,y(1)+k(1,2)/2.d0,y(2)+k(2,2)/2.d0,
     &  y(3)+k(3,2)/2.d0,y(4)+k(4,2)/2.d0,t+h/2.d0)
       end do
       
       do j=1,4 !todas las k4
        k(j,4)=h*f(j,y(1)+k(1,3),y(2)+k(2,3),y(3)+k(3,3)
     &  ,y(4)+k(4,3),t+h)
     
       !ya tengo todas las k ya puedo encontrar y(t+h)
       y(j)=y(j)+(k(j,1)+2.d0*k(j,2)+2.d0*k(j,3)+k(j,4))/6.d0
        
       end do 
         
       it=it+1 !cuento las iteraciones
       
       if(it.eq.300) then  !cada 300 iteraciones escribo sobre el fichero

        cte=y(3)**2.d0*dTL**2.d0/2.d0+y(4)**2.d0/(2.d0*y(1))-
     &  G*Mt/(y(1)*2.d0)-
     &  G*Ml/(dTl*dsqrt( (y(1)*dcos(y(2))-dcos(omega*t))**2.d0+
     &  (y(1)*dsin(y(2))-dsin(omega*t))**2.d0))

        write(11,*) t,y(1)*dcos(y(2)),y(1)*dsin(y(2)),
     & dcos(omega*t),dsin(omega*t),0.d0,0.d0,(cte-omega*y(4))*1E-20
     
        it=0
        
       end if
       
       t=t+h 
      
      end do 
      
      close(11)
      close(12)
      
      stop
      end
      
      
      
      real*8 function f(i,r,phi,pr,pphi,t)
      
      integer i
      real*8 r,phi,pr,pphi,t
      real*8 delta,G,mu,Mt,Ml,omega,Rt,Rl,dTL
      
      G=6.67E-11
      Mt=5.9736E24
      Ml=0.07349E24
      dTL=3.844E8
      omega=2.6617E-6
      delta=G*Mt/dTL**3.d0
      Rt=6.378160E6
      Rl=1.7374E6
      mu=Ml/Mt
      
      if(i.eq.1) then
       
       f=pr
       
      else if(i.eq.2) then
      
       f=pphi/r**2.d0
      
      else if(i.eq.3) then
      
       f=pphi**2.d0/r**3.d0-delta*
     & (1.d0/r**2.d0+(mu*(r-dcos(phi-omega*t)))/
     & (1.d0+r**2.d0-2.d0*r*dcos(phi-omega*t))**1.5)
      
      else
      
       f=dsin(phi-omega*t)*(-delta*mu*r/(1.d0+r**2.d0-2.d0*
     & r*dcos(phi-omega*t))**1.5)
       
      end if
      
      
      return
      end
