!lluvia: Metodologias para interpolar lluvia sobre las cuencas
!Copyright (C) <2014>  <Nicolas Velasquez Giron>

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
!Descripcion
!-----------------------------------------------------------------------
!lluvia es un modulo escrito en fortran y diseñado para ser operado
!tanto desde fortran como desde python, para usar desde fortran 
!compilar: "gfortran -c cuencas.f90", para usar desde python debe 
!ser compilado con: "f2py -c -m nombre_modulo cuencas.f90". para usar 
!esta funcionalidad se debe tener instalado el modulo "numpy" de python
!y un compilador de c y fortran instalado, si se trabaja en linux 
!basta con tener instalado gcc y gfortran junto con numpy.
!
!lluvia se ha disenado para interpolar datos puntuales de lluvia y obtener
!a partri de estos campos de lluvia en el tiempo sobre una cuenca determinada.
!actualmente cuenta con los metodos de interpolacion IDW y TIN. 
!El modulo cuenta con herramientas para la lectura de datos binarios de 
!lluvia obtenidos a partir de radar y el modelo wrf (formatos propios)
!El modulo debe ser usado en conjunto con la cuenca trazada a partir del
!modulo "cuencas".
!
!Gran parte de los resultados son arrojados como arreglos de fortran o
!numpy segun el caso.
!
!Gran parte de los resultados del modulo pueden ser asociados a los modulos:
!lluvia, modelacion. Estos pueden ser compilados en conjunto.
!-----------------------------------------------------------------------
!Esrito por: Nicolas Velasquez Giron
!email: nicolas.velasquezgiron@gmail.com

module lluvia

!-----------------------------------------------------------------------
!Variableses globales
!-----------------------------------------------------------------------
!Variables de uso generico
integer Num_reg !Numero de registros 
!Variables de lluvia puntuales 
integer Num_est !Numero de estaciones
real noData
real, allocatable :: Punt_coord(:,:) !Vector con las coordenadas X,Y de las estaciones
integer, allocatable :: Punt_ides(:) !Vector con los ides de las estaciones
real, allocatable :: Punt_lluvia(:,:) !Vector con los registros de las estaciones
integer, allocatable :: tin(:,:) !Matriz con los indices de los triangulos que forman las estaciones
integer, allocatable :: tin_perte(:,:) !Matriz con las pertenencias de cada celda a cada triangulo

contains

!-----------------------------------------------------------------------
!Herramientas varias
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!Herramientas de lectura de datos
!-----------------------------------------------------------------------
subroutine read_rain_point(rutaCoord,rutaLluvia)
    !Define variables
    character*255, intent(in) :: rutaCoord,rutaLluvia
    !f2py intent(in) :: rutaCoord,rutaLluvia
    integer i,j
    !Lee las coordenadas
    open(10,file=rutaCoord,status='old')
	read(10,*)
	read(10,*) Num_est			
	do i=1,2; read(10,*); enddo
	allocate(Punt_coord(2,Num_est),Punt_ides(Num_est))
	do i=1,Num_est
	    read(10,*) Punt_ides(i),Punt_coord(2,i),Punt_coord(1,i)				
	enddo
    close(10)
    !Lee la lluvia de cada estacion
    open(10,file=rutaLluvia,status='old')
	read(10,'(20X,I6)') Num_reg
	read(10,'(21X,I7)') Num_est
	read(10,'(9X,F4.0)') dt
	read(10,'(9X,f7.2)') noData
	read(10,'(14X,3(I2,1X))') dia,mes,ano
	read(10,'(13X,3(I2,1X))') horas,minutos,seg
	read(10,'(9X,I3)') ventana
	read(10,*)
	allocate(Punt_lluvia(Num_est,Num_reg))
	read(10,*) (Punt_ides(i),i=1,Num_est)
	do i=1,Num_reg
	    read(10,*) (Punt_lluvia(j,i),j=1,Num_est)
	enddo
    close(10)
end subroutine
subroutine read_rain_lina(ruta)
    !Variables de entrada
    character*255 ruta
    !f2py intent(in) :: ruta
    integer i,j
    open(10,file=ruta,status='old')
	read(10,*) Num_est,Num_reg
	allocate(Punt_coord(2,Num_est),Punt_ides(Num_est),Punt_lluvia(Num_est,Num_reg))    
	read(10,*) (Punt_coord(2,i),i=1,Num_est)
	read(10,*) (Punt_coord(1,i),i=1,Num_est)
	do i=1,Num_reg
	    read(10,*) (Punt_lluvia(j,i),j=1,Num_est)
	enddo
    close(10)
end subroutine
subroutine read_rain_radar(ruta_radar,N_reg,ndatos,mat)
    character*255, intent(in) :: ruta_radar
    integer, intent(in) :: ndatos,N_reg
    real, intent(out) :: mat(ndatos,N_reg)
    !f2py intent(in) :: ruta_radar
    !f2py intent(in) :: ndatos
    !f2py intent(out) :: mat
    open(10,file=ruta_radar,status='old',form='unformatted',access='direct',recl=4*N_reg)
	do i=1,ndatos
	    read(10,rec=i) mat(i,:)
	enddo
    close(10)
end subroutine
subroutine read_rain_radar2(ruta_radar,N_reg,ndatos,mat)
    character*255, intent(in) :: ruta_radar
    integer, intent(in) :: ndatos,N_reg
    real, intent(out) :: mat(N_reg,ndatos)
    !f2py intent(in) :: ruta_radar
    !f2py intent(in) :: ndatos
    !f2py intent(out) :: mat
    print *, 'hola'
    open(10,file=ruta_radar,status='old',form='unformatted',access='direct',recl=4*ndatos)
	do i=1,n_reg
	    read(10,rec=i) mat(i,:)
	enddo
    close(10)
end subroutine
subroutine read_radar_byline(ruta_radar,line,nceldas,vec,read_stat)
    character*255, intent(in) :: ruta_radar
    integer, intent(in) :: nceldas,line
    real, intent(out) :: vec(nceldas)
    integer, intent(out) :: read_stat
    !f2py intent(in) :: ruta_radar,nceldas,line
    !f2py intent(out) :: vec,read_stat
    open(10,file=ruta_radar,status='old',form='unformatted',access='direct',recl=4*nceldas)
	read(10,rec=line,iostat=read_stat) vec(:)
	if (read_stat.ne.0) then
	    print*, 'Error en lectura'
	endif
    close(10)
end subroutine

!-----------------------------------------------------------------------
!Herramientas de interpolacion
!-----------------------------------------------------------------------
subroutine interpolation_idw(Rain,x,y,N_reg,N_est,nceldas,pp,ruta,correccion) !Interpola a partir de la metodologia IDW para todos los intervalos
    !Variables de entrada
    integer, intent(in) :: nceldas,N_reg,N_est
    real, intent(in) :: pp,x(nceldas),y(nceldas)
    character*255, intent(in), optional :: ruta
    character*3, intent(in), optional :: correccion
    !Variables de salida
    real, intent(out) :: Rain(N_reg,nceldas)
    !f2py intent(in) :: nceldas,x,y,pp,N_reg,N_est
    !f2py intent(in), optional :: ruta
    !f2py intent(out) :: Rain
    integer i,j,cont
    real W(N_est),WR,med_r,diferencia,med
    !Lo hace para cada una de las celdas
    Rain=0.0
    do i=1,nceldas
	!calcula los pesos
	do j=1,N_est
	    W(j)=1.0/(sqrt(((Punt_coord(1,j)-X(i))**2+(Punt_coord(2,j)-Y(i))**2)))**pp
	end do
	!Calcula la lluvia en todos los intervalos
	do j=1,N_reg
	    if (count(Punt_lluvia(:,j).ge.0).gt.0) then
		!calcula la suma del peso por la lluvia
		Wr=sum(W*Punt_lluvia(:,j),mask=Punt_lluvia(:,j).ge.0.0)
		!calcula la lluvia
		Rain(j,i)=Wr/sum(W,mask=Punt_lluvia(:,j).ge.0.0)
	    endif
	enddo	    
    enddo
    !correccion de magnitud: si el usuario lo desea se corrige la magnitud de la interpolacion con lo obs
    if (correccion.eq.'yes') then
	!calcula la precipitacion media y comienza a iterar por Rain para obtener las diferencias y corregir
	do i=1,N_reg
	    med_r=sum(punt_lluvia(:,i),mask=punt_lluvia(:,i).ne.nodata)/count(punt_lluvia(:,i).ne.nodata)
	    med=sum(Rain(i,:))/nceldas
	    if (med.ne.0.0) then
		diferencia=med_r/med
		Rain(i,:)=Rain(i,:)*diferencia
	    endif
	enddo
    endif
    !Guarda resultados
    if (len_trim(ruta).gt.0) then
	open(10,file=ruta,status='replace',form='unformatted',access='direct',recl=4*nceldas)
	    do i=1,N_reg
		write(10,rec=i) Rain(i,:)
	    enddo
	close(10)
    endif
end subroutine
subroutine interpolation_idw_one(Rain,x,y,intervalo,nceldas,pp,ruta,correccion) !Interpola a partir de la metodologia IDW para el intervalo definido
    !Variables de entrada
    integer, intent(in) :: nceldas,intervalo
    real, intent(in) :: x(nceldas),y(nceldas)
    real, intent(in), optional :: pp
    character*255, intent(in), optional :: ruta
    character*3, intent(in), optional :: correccion
    !Variables de salida
    real, intent(out) :: Rain(nceldas)
    !f2py intent(in) :: nceldas,x,y,intervalo
    !f2py intent(in), optional :: ruta,pp,correccion
    !f2py intent(out) :: Rain
    integer i,j,cont
    real W(Num_est),WR,med_r,diferencia,med,pp_interno
    !Evalua si se ingreso el valor de pp
    if (present(pp)) then
	    if (pp.lt.1.0) then
			pp_interno=1.0
	    else
			pp_interno=pp
	    endif
	endif
    !Lo hace para cada una de las celdas
    Rain=0.0
    do i=1,nceldas
		!calcula los pesos
		do j=1,Num_est
		    W(j)=1.0/(sqrt(((Punt_coord(1,j)-X(i))**2+(Punt_coord(2,j)-Y(i))**2)))**pp_interno
		end do
		!Calcula la lluvia en el intervalo definido
		if (count(Punt_lluvia(:,intervalo).ge.0).gt.0) then
		    !calcula la suma del peso por la lluvia
		    Wr=sum(W*Punt_lluvia(:,intervalo),mask=Punt_lluvia(:,intervalo).ge.0.0)
		    !calcula la lluvia
		    Rain(i)=Wr/sum(W,mask=Punt_lluvia(:,intervalo).ge.0.0)
		endif
    enddo
    !correccion de magnitud: si el usuario lo desea se corrige la magnitud de la interpolacion con lo obs
    if (present(correccion)) then
	    if (correccion.eq.'yes') then
		!calcula la precipitacion media y comienza a iterar por Rain para obtener las diferencias y corregir
		med_r=sum(punt_lluvia(:,intervalo),mask=punt_lluvia(:,intervalo).ne.nodata)/count(punt_lluvia(:,intervalo).ne.nodata)
		med=sum(Rain)/nceldas
			if (med.ne.0.0) then
			    diferencia=med_r/med
			    Rain=Rain*diferencia
			endif
	    endif
	endif
    !Guarda resultados
	if (present(ruta)) then
		if (len_trim(ruta).gt.0) then
			open(10,file=ruta,status='replace',form='unformatted',access='direct',recl=4*nceldas)
				write(10,rec=intervalo) Rain
			close(10)
		endif
	endif
end subroutine
subroutine interpolation_tin(Rain,x,y,tin_perte,N_reg,nceldas,ruta,correccion) !Interpola a partir de la metodologia IDW
    !Variables de entrada
    integer, intent(in) :: nceldas,tin_perte(nceldas),N_reg
    real, intent(in) :: x(nceldas),y(nceldas)
    character*255, intent(in), optional :: ruta
    character*3, intent(in), optional :: correccion
    !Variables de salida
    real, intent(out) :: Rain(N_reg,nceldas)
    !f2py intent(in) :: nceldas,x,y,tin_perte,N_reg
    !f2py intent(in), optional :: ruta
    !f2py intent(out) :: Rain
    integer i,j,cont
    real ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,coef1,coef2,det1,det2,det3,det4
    real med,med_r
    !Lo hace para cada una de las celdas
    Rain=0.0
    do i=1,nceldas
	!Determina la pertenencia como un id
	id_p=tin_perte(i)
	if (id_p.ne.noData) then
	    !Obtiene las coordenadas de los puntos que conforman el triangulo para la celda analizada
	    ax=punt_coord(1,tin(1,id_p)) ; ay=punt_coord(2,tin(1,id_p))
	    bx=punt_coord(1,tin(2,id_p)) ; by=punt_coord(2,tin(2,id_p))
	    cx=punt_coord(1,tin(3,id_p)) ; cy=punt_coord(2,tin(3,id_p))
	    dx=x(i) ; dy=y(i)
	    !Resume calculos para no hacerlos tantas veces
	    det1=(cx-ax)*(dy-ay); det2=(by-ay)*(dx-ax); det3=(cy-ay)*(dx-ax); det4=(dy-ay)*(bx-ax)
	    !Calcula el coeficiente para la celda
	    coef1=(bx-ax)*(cy-ay)-(cx-ax)*(by-ay)
	    !Calcula la lluvia en todos los intervalos
	    do j=1,N_reg
		!Determina el eje z como la lluvia de las 3 estaciones
		az=punt_lluvia(tin(1,id_p),j); bz=punt_lluvia(tin(2,id_p),j); cz=punt_lluvia(tin(3,id_p),j)
		if (az.lt.0) az=0
		if (bz.lt.0) bz=0
		if (cz.lt.0) cz=0
		!Calcula el coeficiente del determinante 
		coef2=det1*(bz-az)+det2*(cz-az)-det3*(bz-az)-det4*(cz-az)
		!Calcula la lluvia
		Rain(j,i)=az-coef2/coef1
		if (az-coef2/coef1 .lt. 0.0) Rain(j,i)=0
	    enddo
	else
	    Rain(:,i)=0.0
	endif
    enddo
    !correccion de magnitud: si el usuario lo desea se corrige la magnitud de la interpolacion con lo obs
    if (correccion.eq.'yes') then
	!calcula la precipitacion media y comienza a iterar por Rain para obtener las diferencias y corregir
	do i=1,N_reg
	    med_r=sum(punt_lluvia(:,i),mask=punt_lluvia(:,i).ne.nodata)/count(punt_lluvia(:,i).ne.nodata)	    
	    med=sum(Rain(i,:))/nceldas
	    if (med.ne.0.0) then
		diferencia=med_r/med
		Rain(i,:)=Rain(i,:)*diferencia
	    endif
	enddo
    endif
    !Guarda resultados
    if (len_trim(ruta).gt.0) then
	open(10,file=ruta,status='replace',form='unformatted',access='direct',recl=4*nceldas)
	    do i=1,N_reg
		write(10,rec=i) Rain(i,:)
	    enddo
	close(10)
    endif
end subroutine
subroutine interpolation_tin_one(Rain,x,y,intervalo,nceldas,ruta,correccion) !Interpola a partir de la metodologia TIN para un intervalo
    !Variables de entrada
    integer, intent(in) :: nceldas,intervalo
    real, intent(in) :: x(nceldas),y(nceldas)
    character*255, intent(in), optional :: ruta
    character*3, intent(in), optional :: correccion
    !Variables de salida
    real, intent(out) :: Rain(nceldas)
    !f2py intent(in) :: nceldas,x,y,intervalo
    !f2py intent(in), optional :: ruta
    !f2py intent(out) :: Rain
    integer i,j,cont
    real ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,coef1,coef2,det1,det2,det3,det4
    real med,med_r
    !Lo hace para cada una de las celdas
    Rain=0.0
    do i=1,nceldas
		!Determina la pertenencia como un id
		id_p=tin_perte(1,i)
		if (id_p.ne.noData) then
		    !Obtiene las coordenadas de los puntos que conforman el triangulo para la celda analizada
		    ax=punt_coord(1,tin(1,id_p)) ; ay=punt_coord(2,tin(1,id_p))
		    bx=punt_coord(1,tin(2,id_p)) ; by=punt_coord(2,tin(2,id_p))
		    cx=punt_coord(1,tin(3,id_p)) ; cy=punt_coord(2,tin(3,id_p))
		    dx=x(i) ; dy=y(i)
		    !Resume calculos para no hacerlos tantas veces
		    det1=(cx-ax)*(dy-ay); det2=(by-ay)*(dx-ax); det3=(cy-ay)*(dx-ax); det4=(dy-ay)*(bx-ax)
		    !Calcula el coeficiente para la celda
		    coef1=(bx-ax)*(cy-ay)-(cx-ax)*(by-ay)
		    !Calcula la lluvia en el intervalo
		    !Determina el eje z como la lluvia de las 3 estaciones
		    az=punt_lluvia(tin(1,id_p),intervalo)
		    bz=punt_lluvia(tin(2,id_p),intervalo)
		    cz=punt_lluvia(tin(3,id_p),intervalo)
		    if (az.lt.0) az=0
		    if (bz.lt.0) bz=0
		    if (cz.lt.0) cz=0
		    !Calcula el coeficiente del determinante 
		    coef2=det1*(bz-az)+det2*(cz-az)-det3*(bz-az)-det4*(cz-az)
		    !Calcula la lluvia
		    Rain(i)=az-coef2/coef1
		    if (az-coef2/coef1 .lt. 0.0) Rain(i)=0
		else
		    Rain(i)=0.0
		endif
    enddo
    !correccion de magnitud: si el usuario lo desea se corrige la magnitud de la interpolacion con lo obs
    if (present(correccion)) then
	    if (correccion.eq.'yes') then
			!calcula la precipitacion media y comienza a iterar por Rain para obtener las diferencias y corregir
			med_r=sum(punt_lluvia(:,intervalo),mask=punt_lluvia(:,intervalo).ne.nodata)/count(punt_lluvia(:,intervalo).ne.nodata)	    
			med=sum(Rain)/nceldas
			if (med.ne.0.0) then
			    diferencia=med_r/med
			    Rain=Rain*diferencia
			endif
	    endif
	endif
    !Guarda resultados
	if (present(ruta)) then	
		if (len_trim(ruta).gt.0) then
			open(10,file=ruta,status='replace',form='unformatted',access='direct',recl=4*nceldas)
				write(10,rec=intervalo) Rain
			close(10)
		endif
	endif
end subroutine

subroutine field_correction(Rain_in,nceldas,N_reg,Rain_out,ruta) !Corrige un campo de precipitacion cualquiera
    !variables de entrada
    integer, intent(in) :: N_reg,nceldas    
    real, intent(in) :: Rain_in(N_reg,nceldas)
    character*255, intent(in), optional :: ruta
    !Variables de salida
    real, intent(out) :: Rain_out(N_reg,nceldas)
    !f2py intent(in) :: N_reg,nceldas,Rain_in,ruta
    !f2py intent(out) Rain_out
    !Variables locales
    integer i
    real med,med_r
    !calcula la precipitacion media y comienza a iterar por Rain para obtener las diferencias y corregir
    Rain_out=Rain_in
    do i=1,N_reg
	med_r=sum(punt_lluvia(:,i),mask=punt_lluvia(:,i).ne.nodata)/count(punt_lluvia(:,i).ne.nodata)	    
	med=sum(Rain_in(i,:))/nceldas
	if (med.ne.0.0) then
	    diferencia=med_r/med
	    Rain_out(i,:)=Rain_in(i,:)*diferencia
	endif
    enddo
    !Guarda resultados
    if (len_trim(ruta).gt.0) then
	open(10,file=ruta,status='replace',form='unformatted',access='direct',recl=4*nceldas)
	    do i=1,N_reg
		write(10,rec=i) Rain_out(i,:)
	    enddo
	close(10)
    endif
end subroutine


!-----------------------------------------------------------------------
!Herramientas de pre-interpolacion
!-----------------------------------------------------------------------
subroutine pre_idw_points_select(basin_f,nrows,xll,yll,dx,N_est,xy_stat,idw_stations,resultado,N_est_cel,nceldas) !Arroja una matriz con las estaciones que se van a emplear para cada una de las celdas
    !Variables de entrada
    integer, intent(in) :: nceldas,N_est,N_est_cel
    real, intent(in) :: xll,yll,dx
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: xy_stat(2,N_est)
    !Variables de salida
    integer, intent(out) :: idw_stations(N_est_cel,nceldas),resultado
    !f2py intent(in) :: nceldas,N_est,N_est_cel,xll,yll,dx,basin_f
    !f2py intent(out) :: idw_stations,resultado
    !Variables internas
    integer i,j,min_pos,numeros(N_est_cel)
    real max_dist
    real distancias(N_est)
    !Evalua si se ha seleccionado o no un numero de estaciones mayor a las que hay
    if (N_est_cel.lt.N_est) then
	!Pasa por todas las celdas de la cuenca para asignar estaciones
	do i=1,nceldas
	    !Calcula la posicion X, Y de la celda
	    X=xll+(basin_f(2,i)-0.5)*dx
	    Y=yll+dx*((nrows-basin_f(3,i))+0.5)
	    !Calcula la distancia de esa celda con respecto a todas las estaciones
	    distancias=0
	    tomados=0
	    distancias=sqrt((xy_stat(1,:)-X)**2+(xy_stat(2,:)-Y)**2)
	    !itera sobre las distancias hasta encontrar las primeras N_est_cel
	    max_dist=maxval(distancias)
	    do j=1,N_est_cel
		min_pos=minloc(distancias,dim=1)
		idw_stations(j,i)=min_pos
		distancias(min_pos)=max_dist
	    enddo
	enddo
	resultado=0 !Se ha n asignado las primera n estaciones sin problema
    elseif (N_est_cel.eq.N_est) then
	do i=1,N_est_cel
	    numeros(i)=i
	enddo
	do i=1,nceldas
	    idw_stations(:,i)=numeros
	enddo
	resultado=1 !se han asignado todas las estaciones a todas las celdas
    else
	idw_stations=-999
	resultado=2 !no ha sido posible asignar estaciones.
    endif
end subroutine
subroutine pre_tin_points_select(x,y,resultado,tin_stations,nceldas) !Arroja una matriz donde se indica a que triangulo pertenece cada celda 
    !Variables de entrada
    integer, intent(in) :: nceldas
    real, intent(in) :: x(nceldas),y(nceldas)
    !Variables de salida
    integer, intent(out) :: tin_stations(1,nceldas),resultado
    !f2py intent(in) :: nceldas,x,y
    !f2py intent(out) :: tin_stations,resultado
    !Variables locales
    integer i,j,flag,flag2,cont,N_tin,n(2)
    real ori1,ori2,ori3,oriP,xp,yp,x1,x2,x3,y1,y2,y3
    !Comienza la busqueda de pertenencia para cada entrada de la tabla
    flag2=0
    cont=0
    resultado=0
    n=shape(tin)
    N_tin=n(2)
    tin_stations=noData
    do i=1,nceldas
	!Inicializa una bandera que revisa si el punto se encuentra dentro de algún triangulo
	flag=0
	!Comienza a buscar el triangulo de pertenencia
	j=1
	do while (j.le.N_tin .and. flag.eq.0)
	    !Cambio de variables para programar facil
	    xp=x(i); yp=y(i)
	    x1=punt_coord(1,TIN(1,j)); x2=punt_coord(1,TIN(2,j)); x3=punt_coord(1,TIN(3,j))
	    y1=punt_coord(2,TIN(1,j)); y2=punt_coord(2,TIN(2,j)); y3=punt_coord(2,TIN(3,j))
	    !Calcula la orientación del triangulo    
	    oriP=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)
	    ori1=(x1-xp)*(y2-yp)-(y1-yp)*(x2-xp)
	    ori2=(xp-x3)*(y2-y3)-(yp-y3)*(x2-x3)
	    ori3=(x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)
	    !Si las cuatro orientaciones son o positivas o negativas el punto se encuentra dentro
	    if ((oriP.gt.0 .and. ori1.gt.0 .and. ori2.gt.0 .and. ori3.gt.0) .or. (oriP.lt.0 .and. ori1.lt.0 &
	    &.and. ori2.lt.0 .and. ori3.lt.0)) then
		tin_stations(1,i)=j
		flag=1
	    else
		j=j+1
	    end if
	end do
    end do
end subroutine

	
end module
