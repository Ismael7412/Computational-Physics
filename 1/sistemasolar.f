      program sistema_solar
      implicit none
      
      real*8 x(10),y(10),vx(10),vy(10),ax(10),ay(10),m(10),alfax(10)
      real*8 alfay(10),L(10),P(10),x0(10),E(10)
      real*8 t,h,hp,G,Ms,c
      integer i,j,k
      integer TR(10)
      character(len=20) str
      
      
c     INICIALIZACION DE CTES:
      Ms=1.989d30
      G=6.674d-11
      c=149.6d9
      h=59000
      hp=h*sqrt(G*Ms/c**3.d0)
      t=0.d0
      
c     *****INICILIZAMOS LOS VALORES INICIALES**** 
C     SOL:1 MERCURIO:2 VENUS:3 TIERRA:4 MARTE:5 JUPITER:6 SATURNO:7 URANO:8 NEPTUNO:9 PLUTON:10

      
C     DEFINIMOS LOS VALORES DE LAS MASAS REESCALADAS
      m(1)=1.d0
      m(2)=1.6591d-07
      m(3)=2.44846d-06
      m(4)=3.0015d-06
      m(5)=3.2277d-07
      m(6)=9.5475113d-04
      m(7)=2.8557063d-04
      m(8)=4.364002d-05
      m(9)=5.128205d-05
      m(10)=6.28d-09
      
C     ASIGNAMOS A LA COORDENADA X EL RADIO (CIRCULAR) REESCALADO DEL PLANETA i PARA COMENZAR LA SIMULACION EN EL EJE X.
      x(1)=0.d0
      x(2)=0.3870320856     
      x(3)=0.7232620321
      x(4)=1.d0
      x(5)=1.523395722
      x(6)=5.204545455
      x(7)=9.582219251
      x(8)=19.20120321
      x(9)=30.04745989
      x(10)=39.23796791      
      
C     INICIALIZAMOS TODAS LAS COORDENADAS 'Y' Y LAS VELOCIDADES INICIALES EN 'X' A 0
C     LA VARIABLE TR SERVIRÁ COMO UNA ESPECIE DE INDICADOR DE SI YA SE HA PASADO DE NUEVO POR LA POSICION INICIAL
C     LA VARIABLE x0 GUARDARÁ LAS POSICIONES INICIALES DE LA COORD X EN t=0 (cuando y=0)
      do i=1,10
       y(i)=0.d0
       vx(i)=0.d0
       ax(i)=0.d0
       ay(i)=0.d0
       alfax(i)=0.d0
       alfay(i)=0.d0
       TR(i)=0
       x0(i)=x(i)
       E(i)=0.d0
  
      end do
     
      write(*,*) "MOMENTOS ANGULARES INICIALES: ****************"
C     INICIALIZAMOS LA VELOCIDAD INICIAL EN 'Y' DETERMINANDO ORBITAS CIRCULARES
C     ADEMÁS CALCULAREMOS EL MOMENTO ANGULAR REESCALADO INICIAL DE CADA PLANETA SUPONIENDO LAS ORBITAS CIRCULARES Y EL PERIODO
      vy(1)=0.d0    
      do i=2,10
      
       vy(i)=dsqrt(1.d0/x(i)) !reescalando y con las unidades usadas va**2=1/x(i)
       L(i)=m(i)*x(i)*vy(i)    !ri=xi î, vi=vyi ĵ entonces sin(90)=1
       write(*,*) i,":",L(i)
       
      end do 
      
      
C     CALCULAMOS LAS ENERGÍAS INICIALES
      write(*,*) "ENERGÍAS INICIALES: **************************"
      do i=2,10
      E(i)=0.5*m(i)*(vx(i)**2.d0+vy(i)**2.d0)   
       do j=1,10
       
        if(j.ne.i) then
        
         E(i)=E(i)-((m(i)*m(j))/((x(i)-x(j))**2.d0
     &+(y(i)-y(j))**2.d0)**0.5)   
        
        end if
       
       end do 
       write(*,*) i,":",E(i)    
      end do 
    
      
C     CALCULAMOS LAS ACELERACIONES INICIALES Y LAS ALFAS INICIALES

      do i=2,10 !para todos los planetas menos el sol
      
      
       do j=1,10 !debe estar incluida la interacción con el sol
       
        if(j.ne.i) then !excluimos al planeta i
         
         ax(i)=ax(i)-((m(j)*(x(i)-x(j)))/
     &(((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0)**1.5))
     
         ay(i)=ay(i)-((m(j)*(y(i)-y(j)))/
     &(((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0)**1.5))
     
        end if
       
       end do
       
       
c      Ya tenemos las aceleraciones iniciales, calculamos las alfas iniciales para cada coordenada 
       alfax(i)=vx(i)+(hp/2.d0)*ax(i)
       alfay(i)=vy(i)+(hp/2.d0)*ay(i)
      
      end do 
      
      
      
C     ABRIMOS UN ARCHIVO DE DATOS POR CADA CUERPO  
      do i=1,10
       open(10+i,file=trim(str(i))//'.dat')
       
C       write(10+i,*) "#t  |  x  |  y  |  vx  |  vy  |  ax  |  ay  |
C     &alfax  |  alfay  |"
       write(10+i,*) t,x(i),y(i),vx(i),vy(i),ax(i),ay(i),alfax(i),
     &alfay(i)
      end do
        
      
C     *************APLICAMOS EL ALGORITMO DE VERLET PARA t>0*******************
          
      do k=1,140000
      
          

       do i=2,10
                       
C       CALCULAMOS x(0+h) y(0+h) PARA EL PLANETA i
	x(i)=x(i)+hp*alfax(i)
	y(i)=y(i)+hp*alfay(i)
	
	
C  AQUI COMPROBAREMOS SI SE HA PASADO YA POR LA POSICION INICIAL Y SI SI ACTIVAMOS LA VARIABLE TR Y GUARDAMOS LA DIF. DE TIEMPOS (PERIODO)

      if(i.lt.5) then
      
        if((TR(i).eq.0).and.(dabs(x(i)-x0(i)).lt.0.0001)
     &  .and.(t.gt.10)) then      
         
         TR(i)=1   !Ya ha sido guardado el periodo del planeta i
         P(i)=t+hp  ! A esta altura del ciclo aun no se ha actualizado el tiempo, así que le añadimos h'
                    ! Es así porque como nuestro tiempo de referencia es t=0, P(i)=(t+h')-0            
        end if
      
      else
      
        if((TR(i).eq.0).and.(dabs(x(i)-x0(i)).lt.0.01)
     &  .and.(t.gt.100)) then      
         
         TR(i)=1   !Ya ha sido guardado el periodo del planeta i
         P(i)=t+hp  ! A esta altura del ciclo aun no se ha actualizado el tiempo, así que le añadimos h'
                    ! Es así porque como nuestro tiempo de referencia es t=0, P(i)=(t+h')-0            
        end if
      
      end if	
      
      
C       PARA QUE NO SUMA LOS VALORES DE LA INTERACION k-1 VOLVEMOS A PONER A 0 LAS ACELERACIONES
        ax(i)=0.d0
        ay(i)=0.d0
	
C	CON LOS DATOS CALCULADOS DE x e y CALCULAMOS ax y ay en t=0+h NUEVAMENTE ES NECESARIO BARRER LOS PLANETAS
	do j=1,10
                           
         if(j.ne.i) then
         
         ax(i)=ax(i)-((m(j)*(x(i)-x(j)))/
     &(((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0)**1.5))
     
         ay(i)=ay(i)-((m(j)*(y(i)-y(j)))/
     &(((x(i)-x(j))**2.d0+(y(i)-y(j))**2.d0)**1.5))
     
         end if
         
        end do
        
C       CALCULO LAS VELOCIDADES EN t=0+h
        vx(i)=alfax(i)+(hp/2.d0)*ax(i)
        vy(i)=alfay(i)+(hp/2.d0)*ay(i)
        
C       LAS ALFAS PARA LA SIG. INTERACION 
        alfax(i)=vx(i)+(hp/2.d0)*ax(i)
        alfay(i)=vy(i)+(hp/2.d0)*ay(i) 
        
        t=t+hp
        
        write(10+i,*) t,x(i),y(i)
        
       end do
      
C     EN CADA INTERACIÓN DEL ALGORITMO EL SOL QUEDA INVARIANTE Y LO REFLEJAMOS AÑADIENDO UNA LINEA DE CEROS 
      write(11,*) t,0.d0,0.d0
      
      
      
      
      end do
      
C     VAMOS A ESCRIBIR EN PANTALLA LOS MOMENTOS ANGULARES FINALES DESPUES DE 140.000 ITERACIONES

      write(*,*) "MOMENTOS ANGULARES FINALES: ******************"
      do i=2,10
       L(i)=m(i)*(x(i)*vy(i)-vx(i)*y(i))
       write(*,*) i,":",L(i)
       E(i)=0.d0 !volvemos a poner a 0 las energías para calcularlas nuevamente
      end do 
      
C     CALCULAMOS LAS ENERGÍAS FINALES
      write(*,*) "ENERGÍAS FINALES: ****************************"
      do i=2,10
      E(i)=0.5*m(i)*(vx(i)**2.d0+vy(i)**2.d0)
       do j=1,10
               
        if(j.ne.i) then
        
         E(i)=E(i)-((m(i)*m(j))/((x(i)-x(j))**2.d0
     &+(y(i)-y(j))**2.d0)**0.5)   
        end if      
  
       end do 
       write(*,*) i,":",E(i)    
      end do 
      
      write(*,*) "PERIODOS ([s] reescalados): ******************"
      do i=2,10
       write(*,*) i,":",P(i)
      end do 
      
      
C     CERRAMOS TODOS LOS ARCHIVOS 
      do i=1,10
       close(10+i)
      end do
      
      
      
      stop
      end
      
C     FUNCION NECESARIA PARA ABRIR LOS ARCHIVOS .dat CON EL CICLO DO 
      character(len=20) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str
