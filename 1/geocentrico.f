      program sistema_solar
      implicit none

      real*8 tn(10),x(10),y(10)
      real*8 t,xt,yt
      integer i,j,k,N,TR
      parameter(N=140001) !importante el numero de interacciones que se dio en sistemasolar.f
      character(len=20) str


C     ABRIMOS UN ARCHIVO DE DATOS POR CADA CUERPO
      do i=1,10
       open(10+i,file=trim(str(i))//'.dat')
       open(20+i,file=trim(str(i))//'_geo.dat')

       write(20+i,*) "#t  |  x  |  y  |  vx  |  vy  |  ax  |  ay  |
     &alfax  |  alfay  |"

      end do

       tn=0.d0
       x=0.d0
       y=0.d0
       xt=0.d0
       yt=0.d0
       TR=0


       do i=1,10
         if(i.ne.4) then
         
	   do k=1,N
             
C          LEO DEL PLANETA i
           read(10+i,*) tn(i),x(i),y(i)
C          AHORA DE LA TIERRA
           read(14,*) tn(4),xt,yt
           write(20+i,*) tn(i),x(i)-xt,y(i)-yt

           if(TR.eq.0) then
            write(24,*) 0.d0,0.d0,0.d0
           end if
           
           
           end do
           
           TR=1
           rewind(14) !para que relea el archivo de la tierra una y otra vez

          end if
       end do












C     CERRAMOS TODOS LOS ARCHIVOS
      do i=1,10
       close(10+i)
       close(20+i)
      end do



      stop
      end

C     FUNCION NECESARIA PARA ABRIR LOS ARCHIVOS .dat CON EL CICLO DO
      character(len=20) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str
