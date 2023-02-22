      program ising
      
      implicit none
      integer i,j,k,l,d,N,x,y,PMC,sumprom1
      real*8 NN,f(200)
      parameter(N=64,NN=10000.d0)
      integer s(0:N+1,0:N+1)
      
      real*8 T,p,dran_u,E,exponencial,ji,sumMn,sumeN,sumcN,sumM
      real*8 sumprom2
      real*8 prom1,prom2,prom3,prom4
      integer*4 i_dran
      
      call dran_ini(857235)
      
      f=0.d0
      s=0
C     INICIALIZAMOS EL CONTADOR DE ITERACIONES
      PMC=0
           
      write(*,*) "Introduzca una temperatura entre 0 y 5"
      read(*,*) T
      
C      open(10,file="matrix.dat")
       open(11,file="correlacion.dat")
       open(120,file="todo.dat")
       open(130,file="nicaso.dat")

C     COMENZAMOS EN UNA CONFIGURACION ORDENADA
      do i=1,N
       do j=1,N
      
       s(i,j)=1

       end do
      end do 
      

c     IMPONEMOS LAS CONDICIONES DE CONTORNO PERIODICAS*******
      
      do j=1,N
       s(0,j)=s(N,j)
       s(N+1,j)=s(1,j)
       s(j,0)=s(j,N)
       s(j,N+1)=s(j,1)
      end do               
     
      
C     COMENZAMOS EL ALGORITMO DE METROPOLIS ********

      sumMn=0.d0
      sumeN=0.d0
      sumcN=0.d0
      sumM=0.d0

      do l=1,1000000 !pasos montecarlo
      do k=1,N*N !iteraciones por paso montecarlo
      
       x=i_dran(N)
       y=i_dran(N)
       E=2.d0*s(x,y)*(s(x+1,y)+s(x-1,y)+s(x,y+1)+s(x,y-1))
       
       exponencial=exp(-E/T)
       
       if(1.lt.exponencial) then
        p=1.d0
       else
        p=exponencial
       end if
      
       ji=dran_u()
       
       if(ji.lt.p) then
        s(x,y)=-1*s(x,y)
       end if
       !**********************************!
       
      
      end do
      
       PMC=PMC+1 !añadimos un paso montecarlo
       if(PMC.eq.100) then
       
        sumprom1=0
        sumprom2=0.d0
        
        do i=1,N
         do j=1,N
         
          sumprom1=sumprom1+s(i,j) !para la magnetizacion promedio y suscept.
          sumprom2=sumprom2-0.5*(s(i,j)*
     &    (s(i,j+1)+s(i,j-1)+s(i-1,j)+s(i+1,j))) !energia para En y Cn
           
         end do
        end do
        
        !funcion de correlacion
        do d=1,floor(dfloat(N)/2.d0)
         do i=1,N
          do j=1,N
          
           f(d)=f(d)+s(i,j)*s(i+d,j) !suma de la funcion a la distancia d de todos los PMC pero sin promediar ni dividir entre N*N
           
          end do
         end do
         
         
        end do
        
               
        sumMn=sumMn+dfloat(sumprom1) !la suma total del primer promedio
        sumeN=sumeN+sumprom2 !la suma total del segundo promedio
        sumcN=sumcN+sumprom2**2.d0 !la suma de las energias al cuadrado en 100 PMC de diferencia
        sumM=sumM+sumprom1**2.d0 !la suma de las magnetizaciones al cuadrado
                
        PMC=0
       
       end if
       
      end do
      
      prom1=(1.d0/N**2.d0)*(sumMn/NN) !magnetizacion promedio
      prom2=sumeN/(NN*N**2.d0) !energia media
      prom3=(sumcN/NN-(sumeN/NN)**2.d0)/(N**2.d0*T**2.d0) !cN
      prom4=(sumM/NN-(sumMn/NN)**2.d0)/(N**2.d0*T) !susceptibilidad magnetica
      
      !promediamos y calculamos la funcion final para cada d
      do d=1,floor(dfloat(N)/2.d0)
      
       write(130,*) d,f(d)/(10**4.d0*N*N)
       write(11,*) f(d)/(10**4.d0*N*N)
             
      end do
      
     
     
C     ***************
      write(*,*) "Magnetizacion promedio:",prom1
      write(*,*) "Energia media:",prom2
      write(*,*) "Calor especifico:",prom3
      write(*,*) "Susceptibilidad magnetica:",prom4

C     *************** 
      write(120,*) prom1
      write(120,*) prom2
      write(120,*) prom3
      write(120,*) prom4

     
C      close(10)
       close(11)
       close(120)
       close(130)
      
      stop
      end
      
      include 'dranxor2.f'
