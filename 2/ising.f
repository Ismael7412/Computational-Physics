      program ising
      
      implicit none
      integer i,j,k,aux,N,it,x,y
      parameter(N=100)
      integer s(0:N+1,0:N+1)
      real*8 T,p,dran_u,E,exponencial,ji
      integer*4 i_dran
      
      call dran_ini(857235)
      
C     INICIALIZAMOS EL CONTADOR DE ITERACIONES
      it=0 
           
      write(*,*) "Introduzca una temperatura entre 0 y 5"
      read(*,*) T
      
C     INICIALIZAMOS LA RED CON UN NUMERO ALEATORIO DE 1 Y -1 PARA ELLO OBTENEMOS NUMEROS ENTRE 1 Y 2 Y ASOCIAMOS AL 1 CON EL -1 Y EL 2 CON      EL    1

      
      open(10,file="matrix.dat")

      do i=1,N
       do j=1,N
      
        aux=i_dran(2)
        if(aux.eq.1) then      
         s(i,j)=-1      
        else
         s(i,j)=1
        end if
      
       end do
      end do 
      

c     IMPONEMOS LAS CONDICIONES DE CONTORNO PERIÓDICAS*******
      
      do j=1,N
       s(0,j)=s(N,j)
       s(N+1,j)=s(1,j)
       s(j,0)=s(j,N)
       s(j,N+1)=s(j,1)
      end do               
     
C     IMPRIMIMOS MATRIZ INICIAL*********
            
      do i=1,N       
       write(10,*) (s(i,j),j=1,N)             
      end do 
       write(10,*) !hay que dejar dos lineas en blanco
       write(10,*)
      
C     COMENZAMOS EL ALGORITMO DE METROPOLIS ********
       
      do k=1,2000000
      
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
        s(x,y)=-1.d0*s(x,y)            
       end if
       
       
       it=it+1
       
       if(it.eq.N**2) then !vemos si se ha superado un paso montecarlo
         
             
        do i=1,N       
                   
         write(10,*) (s(i,j),j=1,N) !imprimir la matriz 
              
        end do 
        
        write(10,*) 
        write(10,*)     
       
        it=0       
       end if
       
      
      end do
      
     
      close(10)
      
      stop
      end
      
      include 'dranxor2.f'
