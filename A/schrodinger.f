      program schrodinger
      
      implicit none
      
      integer Np,i,j,k,t,mt,nd,textra
      integer encontrada,intento
      real*8 h,s,st,k0t,lambda,L,pi,sum,nciclos,dran_u
      real*8 Vt(0:2000)
      real*8 espp2(0:2000)
      complex*16 phi(0:2000,0:2000),chi(0:2000,0:2000)
      complex*16 alpha(0:2000),beta(0:2000,0:2000),omega(0:2000)
      complex*16 icomp
      parameter(icomp=(0.d0,1.d0),pi=dacos(-1.d0)) 
      
          
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
      
      nd=150 ! PASOS TEMPORALES QUE EVOLUCIONA EL SISTEMA **********
      mt=0

            
C     CALCULO A PARTIR DE LA N EL VALOR DEL PASO ESPACIAL
      h=1.d0
      L=h*dfloat(Np)
      
C     PRIMEROS CALCULOS
      k0t=(2.d0*pi*nciclos)/(dfloat(Np)) !energias
      st=1.d0/(4.d0*(k0t**2.d0)) !paso temporal     

      
C     ESTABLECEMOS EL VALOR FIJO DE V(j)
      do i=0,Np
      
       if((dfloat(i).lt.(2.d0*dfloat(Np)/5.d0)).OR.
     & (dfloat(i).gt.(3.d0*dfloat(Np)/5.d0))) then
     
        Vt(i)=0.d0
     
       else
        
        Vt(i)=lambda*(k0t**2.d0)
       
       end if
      
      end do 
            
      
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
      
      
      !***COMENZAMOS LOS EXPERIMENTOS ***

      do k=1,1000 !experimentos, calculo las funciones de onda en cada uno
      
       call dran_ini(857235*k) !inicializamos numeros aleatorios diferentes para cada experimento
      
       do j=1,Np-1

        !PULSO GAUSSIANO INICIAL DEL GUION DEL OBLIGATORIO
        phi(j,0)=exp(icomp*k0t*dfloat(j))*
     &  dexp((-8.d0*(dfloat(4*j-Np))**2.d0)/(Np**2.d0))
      
       end do
      
       call norm(phi,0,Np) !normalizo en t=0
       
       
       
       !******funcion de onda en tiempos posteriores******* ALGORITMO OBLIGATORIO
      
       intento=0
       encontrada=0
       textra=0
      
      
      do t=0,nd !tiempo de evolucion
            
       call evolucion(phi,omega,beta,alpha,chi,t,Np,st) !evoluciona la funcion de onda desde 0 a nd
          
      end do
      
      
       !****YA TENEMOS LA FUNCION DE ONDA EVOLUCIONADA, SEGUIMOS CON EL ALGORITMO DEL VOLUNTARIO
       do while(encontrada.eq.0)
       
       
       !aqui vamos a evolucionar la funcion de onda solo si no se ha encontrado anteriormente
       if(intento.gt.0) then
       
        textra=nd+2*intento !añadimos un tiempo extra a nd para que evolucione un poco el sistema y volvemos a comprobar
               
        do t=nd,nd+textra !tiempo de evolucion
            
         call evolucion(phi,omega,beta,alpha,chi,t,Np,st) !evoluciona la funcion de onda desde nd a nd+textra
          
        end do       
       end if
       
       
       !calculo la probabilidad de que la particula este en el detector derecho
       sum=0.d0
       do j=floor(0.8*Np),Np !0.8=4/5
        sum=sum+abs(phi(j,nd+textra))**2.d0 !en t=nd+textra el ultimo instante
       end do
       
       
       !comparo la probalidad con un numero aleatorio entre 0 y 1
       if(dran_u().le.sum) then
       
        mt=mt+1 !hemos detectado la partícula, empezamos de cero
        encontrada=1
        
       else
       
       
        !no se ha detectado la particula, hacemos 0 la funcion de onda en el detector de la derecha y renormalizamos
        do j=floor(0.8*Np),Np
         phi(j,nd+textra)=0.d0
        end do
        
        
        call norm(phi,nd+textra,Np) !normaliza en t=nd+textra
        
       
        !probabilidad de que la particula este en el detector izquierdo
        sum=0.d0
        do j=0,floor(0.2*Np) !0.2=1/5
         sum=sum+abs(phi(j,nd+textra))**2.d0 !en t=nd el ultimo instante
        end do
       
        if(dran_u().le.sum) then
        
         !se ha detectado la particula en el detector de la izquierda, no hacemos nada, y repetimos.
         encontrada=1
         
         else
         
         
         !no se ha detectado la particula, hacemos 0 la funcion de onda en el detector de la izquierda y renormalizamos
         do j=0,floor(0.2*Np)
          phi(j,nd+textra)=0.d0
         end do
        
         call norm(phi,nd+textra,Np) !normaliza en t=nd
         
         intento=intento+1
         
         
         !****A PARTIR DE AQUI VAMOS TENEMOS QUE DEJAR EVOLUCIONAR UN POCO EL SISTEMA Y VOLVER A COMPROBAR*** 
        
        end if !comprobacion izquierda
        
       end if !comprobacion derecha
      
      end do !do while
      
      write(*,*) mt,k
      
      end do !experimentos
      
      write(*,*) "***Se ha encontrado en la derecha (veces):",mt
      

      stop
      end
      
      
      
      
      
      subroutine evolucion(phi,omega,beta,alpha,chi,t,Np,st)
      integer t,j,Np
      real*8 st
      complex*16 phi(0:2000,0:2000),chi(0:2000,0:2000)
      complex*16 alpha(0:2000),beta(0:2000,0:2000),omega(0:2000)
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
      complex*16 phi(0:2000,0:2000)
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
      
      include 'dranxor2.f'
