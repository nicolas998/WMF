!cuencas: Modulo para el trazado de cuencas y corrientes, y posterior analisis
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
!Cuencas es un modulo escrito en fortran y diseñado para ser operado
!tanto desde fortran como desde python, para usar desde fortran 
!compilar: "gfortran -c cuencas.f90", para usar desde python debe 
!ser compilado con: "f2py -c -m nombre_modulo cuencas.f90". para usar 
!esta funcionalidad se debe tener instalado el modulo "numpy" de python
!y un compilador de c y fortran instalado, si se trabaja en linux 
!basta con tener instalado gcc y gfortran junto con numpy.
!
!cuencas obtiene la red de drenaje y la cuenca a partir de un par coordenado
!de puntos, y obtiene parametros de los mismos, para su funcionamiento 
!requiere un mapa de direcciones y su correspondiente mapa de elevacion 
!digital.  El mapa de direcciones debe estar en formato "teclas de teclado numerico"
!, en caso de que estas no se encuentren en dicho formato, el modulo incluye 
!conversores de: GRASS, ArcGIS.
!
!Gran parte de los resultados son arrojados como arreglos de fortran o
!numpy segun el caso.
!
!Gran parte de los resultados del modulo pueden ser asociados a los modulos:
!lluvia, modelacion. Estos pueden ser compilados en conjunto.
!-----------------------------------------------------------------------
!Esrito por: Nicolas Velasquez Giron
!email: nicolas.velasquezgiron@gmail.com

module cu

!-----------------------------------------------------------------------
!Variableses globales del tipo de mapa leido
!-----------------------------------------------------------------------
!Globales
real xll,yll !coordenadas de la esquina inferior izquierda
real noData !Valor que representa los no datos
real dx !Largo del mapa leido
real dxP !largo proyectado del mapa, debe ser indicado si no se cononce
integer ncols,nrows !cantidad de columnas y filas del mapa
real, allocatable :: stream_temp(:,:) !Vector temporal para el trazado de la corriente
integer, allocatable :: basin_temp(:,:) !Vector temporal para el trazado de la cuenca
real, allocatable :: perim_temp(:,:) !Vector con el perimetro de la cuenca
integer, allocatable :: sub_basins_temp(:,:) !Vector con las sub-cuencas
real, allocatable :: ppal_stream_temp(:,:) !Vector con el cauce principal 
real, allocatable :: netxy_temp(:,:) !Vector con las coordenadas de la red hidrica
real area,perimetro,pend_media !Escalares: area, perimetro y pendiente media de la ultima cuenca trazada
real elevacion !Escalar: elevacion media de la cuenca
real centroX,centroY !Escalares: coordenadas del centroide X y Y de la cuenca trazada
!Para hacer sorting de cosas
public :: QsortC
!Para operaciones con matrices 
real, allocatable :: col_fil_temp(:,:)
!private :: Partition

!-----------------------------------------------------------------------
!Punto de inicio de funciones del modulo
contains

!-----------------------------------------------------------------------
!Herramientas varias
!-----------------------------------------------------------------------
subroutine coord2fil_col(x,y,col,fil) !Obtiene la fila columna a partir de un XY
    !Variables de entrada
    real, intent(in) :: x,y
    !variables de salida
    integer, intent(out) :: fil,col
    !f2py intent(in) :: x,y
    !f2py intent(out) :: fil,col
    !Pasa las coordenadas
    col=abs(ceiling((x-xll)/dx))
    fil=nrows-floor((y-yll)/dx)
    !si las coordenadas estan fuera del rango reporta como -999
    if (col.gt.ncols.or.col.le.0) col=-999
    if (fil.gt.nrows.or.fil.le.0) fil=-999
end subroutine


!-----------------------------------------------------------------------
!Herramientas lectura de mapas
!-----------------------------------------------------------------------
!Lectura de mapas acii
subroutine read_ASCII_hdr(ruta) !Lee el encabezado de un mapa 
    !variables de entrada
    character*255, intent(in) :: ruta
    !f2py intent(in) :: ruta
    !variables propias
    character*20 joder
    open(10,file=ruta,status='old')
	read(10,*) joder,ncols
	read(10,*) joder,nrows
	read(10,*) joder,xll
	read(10,*) joder,yll
	read(10,*) joder,dx
	read(10,*) joder,noData
    close(10)
end subroutine
subroutine read_float_ASCII_dat(ruta,nc,nr,Mat) !Lee mapa flotante
    !Variables de entrada
    character*255,intent(in) :: ruta
    integer, intent(in) :: nc,nr 
    !Variables de salida
    real, intent(out) :: Mat(nc,nr)
    !f2py intent(in) :: ruta
    !f2py intent(in) :: nc,nr
    !f2py intent(out) :: Mat
    !Variables locales
    integer i,j
    character*20 joder
    !Lectura del mapa ascii
    open(10,file=ruta,status='old')
	do i=1,6
	    read(10,*) joder
	enddo
	read(10,*) ((Mat(i,j),i=1,ncols),j=1,nrows)
    close(10)
end subroutine
subroutine read_int_ASCII_dat(ruta,nc,nr,Mat) !Lee mapa entero
    !Variables de entrada
    character*255,intent(in) :: ruta
    integer, intent(in) :: nc,nr 
    !Variables de salida
    integer, intent(out) :: Mat(nc,nr)
    !f2py intent(in) :: ruta
    !f2py intent(in) :: nc,nr
    !f2py intent(out) :: Mat
    !Variables locales
    integer i,j
    character*20 joder
    !Lectura del mapa ascii
    open(10,file=ruta,status='old')
	do i=1,6
	    read(10,*) joder
	enddo
	read(10,*) ((Mat(i,j),i=1,ncols),j=1,nrows)
    close(10)
end subroutine
subroutine read_ASCII_hdr_ng(ruta,nc,nr,xlln,ylln,dxn,noDataN) !Lee el encabezado de un mapa no global
    !variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(out) :: nc,nr
    real, intent(out) :: xlln,ylln,dxn,noDataN
    !f2py intent(in) :: ruta
    !f2py intent(out) :: nc,nr,xlln,ylln,dxn,noDataN
    !variables propias
    character*20 joder
    open(10,file=ruta,status='old')
	read(10,*) joder,nc
	read(10,*) joder,nr
	read(10,*) joder,xlln
	read(10,*) joder,ylln
	read(10,*) joder,dxn
	read(10,*) joder,noDataN
    close(10)
end subroutine
subroutine read_float_ASCII_ng(ruta,nc,nr,Mat) !Lee mapa flotante no global 
!Variables de entrada
    character*255,intent(in) :: ruta
    integer, intent(in) :: nc,nr 
    !Variables de salida
    real, intent(out) :: Mat(nc,nr)
    !f2py intent(in) :: ruta
    !f2py intent(in) :: nc,nr
    !f2py intent(out) :: Mat
    !Variables locales
    integer i,j
    character*20 joder
    !Lectura del mapa ascii
    open(10,file=ruta,status='old')
	do i=1,6
	    read(10,*) joder
	enddo
	read(10,*) ((Mat(i,j),i=1,nc),j=1,nr)
    close(10)
end subroutine
subroutine read_int_ASCII_ng(ruta,nc,nr,Mat) !Lee mapa entero no global 
!Variables de entrada
    character*255,intent(in) :: ruta
    integer, intent(in) :: nc,nr 
    !Variables de salida
    integer, intent(out) :: Mat(nc,nr)
    !f2py intent(in) :: ruta
    !f2py intent(in) :: nc,nr
    !f2py intent(out) :: Mat
    !Variables locales
    integer i,j
    character*20 joder
    !Lectura del mapa ascii
    open(10,file=ruta,status='old')
	do i=1,6
	    read(10,*) joder
	enddo
	read(10,*) ((Mat(i,j),i=1,nc),j=1,nr)
    close(10)
end subroutine
!Lectura de binarios de cuenca
subroutine read_int_basin(ruta,records,vect,nrecords,nceldas) !Lee los datos enteros de un binario de cuenca en los records ordenados
    !Variables de entrada
    integer, intent(in) :: nrecords,nceldas
    character*255, intent(in) :: ruta
    integer, intent(in) :: records(nrecords)
    !Variables de salida
    integer, intent(out) :: vect(nrecords,nceldas)
    !f2py intent(in) :: nrecords,nceldas,ruta,records
    !f2py intent(out) :: vect
    !Variables locales
    integer i 
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	do i=1,nrecords
	    read(10,rec=records(i)) vect(i,:)
	enddo
    close(10)
end subroutine
subroutine read_float_basin(ruta,records,vect,nrecords,nceldas) !Lee los datos flotantes de un binario de cuenca en los records ordenados
    !Variables de entrada
    integer, intent(in) :: nrecords,nceldas
    character*255, intent(in) :: ruta
    integer, intent(in) :: records(nrecords)
    !Variables de salida
    real, intent(out) :: vect(nrecords,nceldas)
    !f2py intent(in) :: nrecords,nceldas,ruta,records
    !f2py intent(out) :: vect
    !Variables locales
    integer i 
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	do i=1,nrecords
	    read(10,rec=records(i)) vect(i,:)
	enddo
    close(10)
end subroutine

!-----------------------------------------------------------------------
!Herramientas de Escritura y conversion de mapas de mapas
!-----------------------------------------------------------------------
!Escritura de mapas
subroutine write_float_ascii(ruta,MAPA) !Escribe mapas flotantes
    !variables de entrada
    !integer, intent(in) :: nc,nr
    character*255, intent(in) :: ruta
    real, intent(in) :: MAPA(:,:)
    !f2py intent(in) ruta,MAPA
    character Tcolumnas*5, formatoR*20
    integer i
    !Pasa el numero de columnas a valor textual
    write(Tcolumnas,'(I5)') ncols
    !Crea el formato
    formatoR='('//Tcolumnas//'F10.3)'
    open(10,file=ruta,status='replace')
	write(10,'(A14,I5)') 'ncols         ',ncols
	write(10,'(A14,I5)') 'nrows         ',nrows
	write(10,'(A14,f16.11)') 'xllcorner     ', xll
	write(10,'(A14,f15.11)') 'yllcorner     ',      yll
	write(10,'(A14,f12.8)') 'cellsize      ',dx
	write(10,'(A14,f8.2)') 'NODATA_value  ', noData
	do i=1,nrows; write(10,formatoR) MAPA(:,i); enddo
    close(10)
end subroutine
subroutine write_int_ascii(ruta,MAPA) !Escribe mapas enteros
    !variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(in) :: MAPA(:,:)
    !f2py intent(in) ruta,MAPA
    character Tcolumnas*5, formatoR*20
    integer i
    !Pasa el numero de columnas a valor textual
    write(Tcolumnas,'(I5)') ncols
    !Crea el formato
    formatoR='('//Tcolumnas//'I7)'
    open(10,file=ruta,status='replace')
	write(10,'(A14,I5)') 'ncols         ',ncols
	write(10,'(A14,I5)') 'nrows         ',nrows
	write(10,'(A14,f16.11)') 'xllcorner     ', xll
	write(10,'(A14,f15.11)') 'yllcorner     ',      yll
	write(10,'(A14,f12.8)') 'cellsize      ',dx
	write(10,'(A14,f8.1)') 'NODATA_value  ', noData
	do i=1,nrows; write(10,formatoR) MAPA(:,i); enddo
    close(10)
end subroutine
subroutine write_float_ascii_ng(ruta,MAPA,cols,rows,xll_ng,yll_ng,dx_ng) !Escribe mapas flotantes
    !variables de entrada
    integer, intent(in) :: cols,rows
    real, intent(in) :: xll_ng,yll_ng,dx_ng
    character*255, intent(in) :: ruta
    real, intent(in) :: MAPA(:,:)
    !f2py intent(in) ruta,MAPA,cols,rows,xll_ng,yll_ng,dx_ng
    character Tcolumnas*5, formatoR*20
    integer i
    !Pasa el numero de columnas a valor textual
    write(Tcolumnas,'(I5)') cols
    !Crea el formato
    formatoR='('//Tcolumnas//'F10.3)'
    open(10,file=ruta,status='replace')
	write(10,'(A14,I5)') 'ncols         ',cols
	write(10,'(A14,I5)') 'nrows         ',rows
	write(10,'(A14,f16.11)') 'xllcorner     ', xll_ng
	write(10,'(A14,f15.11)') 'yllcorner     ',      yll_ng
	write(10,'(A14,f12.8)') 'cellsize      ',dx_ng
	write(10,'(A14,f8.2)') 'NODATA_value  ', noData
	do i=1,rows; write(10,formatoR) MAPA(:,i); enddo
    close(10)
end subroutine
subroutine write_int_ascii_ng(ruta,MAPA,cols,rows,xll_ng,yll_ng,dx_ng) !Escribe mapas flotantes
    !variables de entrada
    integer, intent(in) :: cols,rows
    real, intent(in) :: xll_ng,yll_ng,dx_ng
    character*255, intent(in) :: ruta
    integer, intent(in) :: MAPA(:,:)
    !f2py intent(in) ruta,MAPA,cols,rows,xll_ng,yll_ng,dx_ng
    character Tcolumnas*5, formatoR*20
    integer i
    !Pasa el numero de columnas a valor textual
    write(Tcolumnas,'(I5)') cols
    !Crea el formato
    formatoR='('//Tcolumnas//'I7)'
    open(10,file=ruta,status='replace')
	write(10,'(A14,I5)') 'ncols         ',cols
	write(10,'(A14,I5)') 'nrows         ',rows
	write(10,'(A14,f16.11)') 'xllcorner     ', xll_ng
	write(10,'(A14,f15.11)') 'yllcorner     ',      yll_ng
	write(10,'(A14,f12.8)') 'cellsize      ',dx_ng
	write(10,'(A14,f8.2)') 'NODATA_value  ', noData
	do i=1,rows; write(10,formatoR) MAPA(:,i); enddo
    close(10)
end subroutine
!Escritura archivos binarios en formato para modelo
subroutine write_int_basin(ruta,vect,record,nceldas,estado)
    !Variables de entrada
    integer, intent(in) :: nceldas
    character*255, intent(in) :: ruta
    character*7, intent(in) :: estado
    integer, intent(in) :: vect(nceldas)
    integer, intent(in) :: record
    !f2py intent(in) :: ruta,vect,record,nceldas,estado
    !Escribe el archivo
    open(10,file=ruta,status=estado,form='unformatted',access='direct',RECL=4*nceldas)
	write(10,rec=record) vect 
    close(10)
end subroutine
subroutine write_float_basin(ruta,vect,record,nceldas,estado)
    !Variables de entrada
    integer, intent(in) :: nceldas
    character*255, intent(in) :: ruta
    character*7, intent(in) :: estado
    real, intent(in) :: vect(nceldas)
    integer, intent(in) :: record
    !f2py intent(in) :: ruta,vect,record,nceldas,estado
    !Escribe el archivo
    open(10,file=ruta,status=estado,form='unformatted',access='direct',RECL=4*nceldas)
	write(10,rec=record) vect 
    close(10)
end subroutine

!-----------------------------------------------------------------------
!Trazador de corrientes
!-----------------------------------------------------------------------
subroutine stream_find(x,y,DEM,DIR,nc,nf,nceldas) !Encuentra la corriente
    !Variables de entrada
    integer, intent(in) :: nc,nf
    real, intent(in) :: x,y
    real, intent(in) :: DEM(nc,nf)
    integer, intent(in) :: DIR(nc,nf)
    !Variables de salida
    integer, intent(out) :: nceldas
    !f2py intent(in) :: x,y,DIR,DEM,nc,nf
    !f2py intent(out) :: nceldas
    !variables locales
    integer kc,kf,flag,envia,c,i,dire,col,fil,cont
    character*20 a,b
    !aloja el vector de cauce
    if (.not. allocated(stream_temp)) allocate(stream_temp(4,ncols*nrows))
    !vector de mascara
    distancia=0
    cont=0
    flag=1
    stream_temp=0
    !col=coli; fil=fili
    call coord2fil_col(x,y,col,fil)
    dire=DIR(col,fil)
    do while (flag.eq.1)
		do kf=1,3
		    do kc=1,3
			!evalua si esa es la direccion
				envia=9-3*kf+kc
				if (dire.eq.envia) then
				    cont=cont+1
				    !Guarda valores en el vector de la corriente
				    stream_temp(1,cont)=xll+dx*(col-0.5)
				    stream_temp(2,cont)=yll+dx*(nrows-fil+0.5)
				    stream_temp(3,cont)=DEM(col,fil)
				    !actualiza el valor de la fila columna a evaluar
				    col=col+kc-2
				    fil=fil+kf-2
				    dire=DIR(col,fil)
				    !calcula la distancia acumulada
				    if (mod(envia,2).eq.0) then
						distancia=distancia+dxP
				    else
						distancia=distancia+dxP*1.4
				    endif
				    stream_temp(4,cont)=distancia
				endif
		    enddo
		enddo
		!Evalua si puede seguir o n
		if (col.gt.ncols.or.fil.gt.nrows.or.col.le.0 &
			&.or.fil.le.0 .or.dire.eq.noData .or. dire .le. 0) then
		    flag=0
		endif
    end do
    !copia el vector de salida
    nceldas=cont
end subroutine
subroutine stream_cut(nceldas,stream_f) !Genera una copia corta de la corriente
    !Variables de entrada
    integer, intent(in) :: nceldas
    real, intent(out) :: stream_f(4,nceldas)
    !f2py intent(in) :: nceldas
    !f2py intent(out) :: stream_f
    stream_f=stream_temp(:,1:nceldas)
end subroutine
subroutine stream_kml(corr,ruta,nceldas) !escribe un kml de la corriente
    !Variables de entrada
    integer, intent(in) :: nceldas
    character*255, intent(in) :: ruta
    real, intent(in) :: corr(4,nceldas)
    !f2py intent(in) :: nceldas,corr,ruta
    !Variables internas
    integer i
    !Escribe el archivo en formato kml
    open(10,file=ruta,status='replace')
	write(10,'(A38)')'<?xml version="1.0" encoding="UTF-8"?>'
	write(10,*) '<kml xmlns="http://www.opengis.net/kml/2.2">'
	write(10,*) '<Document>'
	write(10,*) '<Placemark>'
	write(10,*) '<name>Corriente</name>'
	write(10,*) '<description> corriente </description>'
	write(10,*) '<Style>'
	write(10,*) '<LineStyle>'
	write(10,*) '<color>ff660000</color>'
	write(10,*) '<width>3</width>'
	write(10,*) '</LineStyle>'
	write(10,*) '</Style>'
	write(10,*) '<LineString>'
	write(10,*) '<coordinates>'	
	do i=1,nceldas-1
	    write(10,'(1F15.5,A1,1F15.5,A1,I1)') corr(1,i),',',corr(2,i),',',0
	enddo
	write(10,*) '</coordinates>'
	write(10,*) '</LineString>'
	write(10,*) '</Placemark>'
	write(10,*) '</Document>'
	write(10,*) '</kml>'
    close(10)
end subroutine
subroutine stream_reset !Se usa si se va a cambiar de mapa de referencia
    deallocate(stream_temp)
end subroutine

!-----------------------------------------------------------------------
!Trazador de cuencas
!-----------------------------------------------------------------------
!funciones de trazado
subroutine basin_find(x,y,DIR,nc,nr,nceldas) !encuentra toda la cuenca
    !Variables de entrada
    real, intent(in) :: x,y
    integer, intent(in) :: nc, nr
    integer, intent(in) :: DIR(nc,nr)
    !Variables de salida
    integer, intent(out) :: nceldas
    !f2py intent(in) :: x,y,DIR
    !f2py intent(out) :: nceldas
    !f2py intent(in) :: nc,nr
    !Variables locales
    integer kf,kc
    integer tenia,coli,fili,col2,row2,i,j,cont3,cont2,res(2,7),paso,longitud,celTot,c,f
    !aloja el vector de cauce
    if (.not. allocated(basin_temp)) allocate(basin_temp(3,ncols*nrows))
    !Encuentra la fila columna 
    call coord2fil_col(x,y,coli,fili)
    !Encuentra la cuenca que drena al punto
    celTot=ncols*nrows
    basin_temp=-1    
    basin_temp(1,celTot)=0 !Celda a la que drena
    basin_temp(2,celTot)=coli !Columna de la celda
    basin_temp(3,celTot)=fili !Fila de la celda
    tenia=1
    cont2=1
    col2=coli 
    row2=fili
    do while (col2>0) !Si agarra una entrada negativa del vector para
	!Encuentra las celdas que le drenan a la celda objetivo
	cont3=0
	!Busca alrededor
	do kf=1,3
	    do kc=1,3
		!posicion de la columna y la fila a evaluar
		c=col2+kc-2
		f=row2+kf-2
		!se fija si la columna y la fila a evaluar se encuentran en el mapa
		if (((c.le.ncols) .and. (c.gt.0)) .and. ((f.le.nrows) .and. (f.gt.0))) then
		    !se fija si la celda aevaluar le drena a la celda objetivo
		    if (DIR(c,f).eq.3*kf-kc+1) then
			!Incrementa la cantidad de celdas que le drenan
			cont3=cont3+1
			res(1,cont3)=c
			res(2,cont3)=f
		    end if
		end if
	    end do
	end do
	!llena el vector
	if (cont3.gt.0) then
	    do i=1,cont3
		!basin_temp(1,celTot-tenia+1-i)=tenia+i 
		basin_temp(1,celTot-tenia+1-i)=cont2
		basin_temp(2,celTot-tenia+1-i)=res(1,i) 
		basin_temp(3,celTot-tenia+1-i)=res(2,i)
	    end do
	end if
	!Actualiza "tenia" con el fin de saber el tamaño que lleva el vector final
	tenia=tenia+cont3 
	paso=celTot-cont2
	!Actualiza la fila y la columna, tomando la siguiente
	col2=basin_temp(2,paso)
	row2=basin_temp(3,paso)
	cont2=cont2+1
    end do
    nceldas=tenia
end subroutine
subroutine basin_cut(nceldas,basin_f) !Genera una matriz con la cuenca final
    !Variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(out) :: basin_f(3,nceldas)
    !f2py intent(in) :: nceldas
    !f2py intent(out) :: stream_f
    !Variables locales
    integer celTot
    !Prueba si si hay cuenca temporal
    if (allocated(basin_temp)) then
	celTot=ncols*nrows
	basin_f=basin_temp(:,celTot-nceldas+1:)
    else
	print *, 'Error: La variable basin_temp no se encuentra alojada'
    endif
end subroutine
subroutine basin_basics(basin_f,DEM,DIR,nc,nf,nceldas,acum,long,pend,elev) !calcula: acumulada, longitud y pendiente
    !variables de entrada
    integer, intent(in) :: nceldas,nc,nf
    integer, intent(in) :: basin_f(3,nceldas), DIR(nc,nf)
    real, intent(in) :: DEM(nc,nf)
    !variables de salida
    integer, intent(out) :: acum(nceldas)
    real, intent(out) :: long(nceldas),pend(nceldas),elev(nceldas)
    !f2py intent(in) :: nceldas,nc,nf,basin_f,DEM,DIR
    !f2py intent(out) :: acum,long,pend,elev
    integer i,drenaid,col_pos,fil_pos
    real X(nceldas),Y(nceldas)
    !Calcula la elevacion
    do i=1,nceldas 
	elev(i)=DEM(basin_f(2,i),basin_f(3,i))
    end do
    !Calcula Longitudes, acum y pendiente
    acum=1
    do i=1,nceldas
	!Determina la celda a la que se drena
	drenaid=nceldas-basin_f(1,i)+1
	!Obtiene la longitud de la celda
	prueba=mod(DIR(basin_f(2,i),basin_f(3,i)),2)
	if (prueba.eq.0.0) then
	    long(i)=dxp
	else
	    long(i)=dxp*sqrt(2.0)
	endif
	!Calcula el area acumulada
	if (basin_f(1,i).ne.0) then
	    acum(drenaid)=acum(drenaid)+acum(i)
	    pend(i)=abs(elev(i)-elev(drenaid))/long(i)
	else
	    call drain_colfil(DIR(basin_f(2,i),basin_f(3,i)),col_pos,fil_pos)
	    pend(i)=abs(elev(i)-DEM(basin_f(2,i)+col_pos,basin_f(3,i)+fil_pos))/long(i)
	endif
	!Si la pendiente es plana 0, le da un poco de pendiente
	if (pend(i).eq.0) pend(i)=0.001
    end do
    !Calcula escalares genericos de la cuenca
    area=nceldas*dxp**2/1e6 
    pend_media=sum(pend)/nceldas
    elevacion=sum(elev)/nceldas
    !Calcula las coordenadas del centroide de la cuenca
    call basin_coordXY(basin_f,X,Y,nceldas)
    call QsortC(X)
    call QsortC(Y)
    if (mod(nceldas,2).eq.0.0) then
	centroX=X(nceldas/2); centroY=Y(nceldas/2)
    else
	centroX=X((nceldas+1)/2); centroY=Y((nceldas+1)/2)
    endif
end subroutine
subroutine basin_acum(basin_f,nceldas,acum) !calcula: acumulada
    !variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    !variables de salida
    integer, intent(out) :: acum(nceldas)
    !f2py intent(in) :: nceldas,nc,nf,basin_f
    !f2py intent(out) :: acum
    integer i,drenaid
    !Calcula Longitudes, acum y pendiente
    acum=1
    do i=1,nceldas
		!Determina la celda a la que se drena
		drenaid=nceldas-basin_f(1,i)+1
		!Calcula el area acumulada
		if (basin_f(1,i).ne.0) then
		    acum(drenaid)=acum(drenaid)+acum(i)			    
		endif	
    end do
    !Calcula escalares genericos de la cuenca
    area=nceldas*dxp**2/1e6 
end subroutine
subroutine basin_findlong(basin_f,nceldas,longMax,punto) !Calcula la longitud de la cuenca
	!variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    !variables de salida
    real, intent(out) :: longMax
    integer, intent(out) :: punto
    !f2py intent(in) :: nceldas,basin_f
    !f2py intent(out) :: punto,longMax
    integer i,drenaid
    real long
    LongMax=0.0
    do i=1,nceldas-1
		long=sqrt( ((basin_f(2,i)-basin_f(2,nceldas))*dxP)**2 &
		& +((basin_f(3,i)-basin_f(3,nceldas))*dxP)**2)
		if (long.gt.LongMax) then 
			LongMax=long
			punto=i
		endif
    enddo
    LongMax=LongMax/1000.0
end subroutine
subroutine basin_time_to_out(basin_f,long,speed,time,nceldas)
    !variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: long(nceldas),speed(nceldas)
    !variables de salida
    real, intent(out) :: time(nceldas)
    !variables internas
    integer i,j,drenaid
    real prueba,tiempo_cel(nceldas),tiempo_acum
    logical flag
    !Definicion de f2p2
    !f2py intent(in) :: nceldas,basin_f,long,speed
    !f2py intent(out) :: time
    !Calcula el tiempo por celda
    tiempo_cel=long/speed
    !Calcula para cada el tiempo a la salida
    do i=1,nceldas    
	drenaid=nceldas-basin_f(1,i)+1	
	if (drenaid.lt.nceldas) then
	    flag=.true.
	    tiempo_acum=tiempo_cel(i)
	    do while (flag)
		tiempo_acum=tiempo_acum+tiempo_cel(drenaid)
		drenaid=nceldas-basin_f(1,drenaid)+1		    		
		if (drenaid.ge.nceldas) flag=.false.		    
	    enddo
	    time(i)=tiempo_acum
	else
	    time(i)=tiempo_cel(i)
	endif	
    enddo    
end subroutine
subroutine basin_arc_slope(basin,DEM,slope,nceldas,nc,nr) !Calcula la pendiente segun el algoritmo de arcgis
    !Variables de entrada
    integer, intent(in) :: nceldas, nc, nr
    integer, intent(in) :: basin(3,nceldas)
    real, intent(in) :: DEM(nc,nr)
    !Variables de salida
    real, intent(out) :: slope(nceldas)
    !f2py intent(in) :: nceldas, nc,nr, basin, DEM
    !f2py intent(out) :: slope
    !Variables locales 
    integer i,j,cel,col,fil
    real k(3,3),pen_dx,pen_dy
    !Calcula la pendiente para cada celda de la ladera
    do cel=1,nceldas
	col=basin(2,cel); fil=basin(3,cel)
	!Obtiene la elevacion en el punto 
	Elev=DEM(col,fil)
	!obtiene el kernel en el punto y posiciona la fila columna en el extremo sup izquierdo
	k=0.0
	col=col-1; fil=fil-1	
	do i=0,2
	    do j=0,2
			if (col+j .le. nc .and. col+j .ge. 1 .and. fil+i .le. nr .and. fil+i .ge. 1) then
			    if (DEM(col+j,fil+i).ne.noData) then
				k(j+1,i+1)=DEM(col+j,fil+i)
			    else
				k(j+1,i+1)=DEM(col+1,fil+1)
			    endif
			endif
	    enddo
	enddo
	!Calcula la pendiente con el kernel 
	pen_dx=((k(3,1)+2*k(3,2)+k(3,3))-(k(1,1)+2*k(1,2)+k(1,3)))/(8*dxp)
	pen_dy=((k(1,3)+2*k(2,3)+k(3,3))-(k(1,1)+2*k(2,1)+k(3,1)))/(8*dxp)
	slope(cel)=sqrt(pen_dx**2 + pen_dy**2)	
    enddo
end subroutine 
subroutine basin_reset !Se usa si se va a cambiar de mapa de referencia
    deallocate(basin_temp)
end subroutine
subroutine basin_coordXY(basin_f,X,Y,nceldas) !obtiene las coordenadas X,Y del centro de cada celda
    !Variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    !Variables de salida
    real, intent(out) :: X(nceldas),Y(nceldas)
    !f2py intent(in) :: nceldas, basin_f
    !f2py intent(out) :: X,Y
    !Calcula 
    X=xll+dx*(basin_f(2,:)-0.5)
    Y=yll+dx*((nrows-basin_f(3,:))+0.5)
end subroutine
!Funciones de delineacion de cuenca
subroutine basin_perim_find(basin_f,nperim,nceldas) !Encuentra los puntos X,Y del borde de la cuenca
    !Variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    !Variables de salida
    integer, intent(out) :: nperim
    !f2py intent(in) :: nceldas,basin_f
    !f2py intent(out) :: nperim
    !Variables locales
    integer flag1,flag2,kf,kc,f1,c1,f2,c2,cont,cont2,col2,fil2,i,longitud,col,fil
    integer OrtoCol(4),OrtoFil(4),mov1(4),mov2(4),DiagMov1(4),DiagMov2(4),DiagMov3(4),pos(4),posN(4)
    integer OrtoAdi(4)
    real yll_loc,xll_loc
    integer mascara(ncols,nrows),Nveces
    !Movimientos de la busqueda ortogonal
    OrtoCol=[-1,0,1,0]; OrtoFil=[0,-1,0,1]
    OrtoAdi=[0,1,1,0]
    !Movimientos de la busqueda diagonal
    DiagMov1=[-1,1,1,-1];DiagMov2=[-1,-1,1,1];DiagMov3=[0,1,0,-1]
    !Movimientos de busqueda de linderos
    mov1=[-1,0,0,-1]; mov2=[-1,-1,0,0]
    !Vector de posiciones 
    pos=[4,1,2,3];
    !Asigna tamano a la matrizs de vectores a usar y define las coordenadas iniciales    
    if (allocated(perim_temp).eqv. .false.) then
		allocate(perim_temp(2,ncols*nrows))
	endif
    perim_temp=0.0
    col2=basin_f(2,nceldas)
    fil2=basin_f(3,nceldas)
    xll_loc=xll-dx
    yll_loc=yll-dx
    !Genera la matriz de mascara de ceros y unos de la cuenca
    mascara=0
    do i=1,nceldas
		mascara(basin_f(2,i),basin_f(3,i))=1
    enddo
    !Encuentra el punto de partida
    flag1=1; flag2=1; cont=1; 
    !Evalua si hay lados ortogonales libres para iniciar a partir de la celda de salida
    do while (flag1.eq.1)	
		!Col Fil a evaluar		
		c1=col2+OrtoCol(cont)
		f1=fil2+OrtoFil(cont)
		!Evalua si ya encontro un lugar para iniciar
		if (mascara(c1,f1).eq.0) then
			!Guarda el punto para ser graficado
			perim_temp(1,1)=(col2-1+OrtoAdi(cont))*dx+xll_loc
			perim_temp(2,1)=(nrows+2-fil2+mov2(5-cont))*dx+yll_loc
			!Da la orden de finalizar 
			flag1=0
		endif
		!Evalua si ya debe terminar 
		if (cont.gt.4) then
			flag1=0
			flag2=0
		endif
		!Incrementa el contador		
		cont=cont+1		
    enddo    
    !print *, flag1,flag2,cont
    !Si no encontro lados ortogonales evalua si hay lados diagonales
    if (flag2.eq.0) then
		!Reinicia contador y bandera 1
		flag1=1; cont=1
		print *, 'Entro a diagonales'
		!Evalula si hay lados diagonales
		do while (flag1.eq.1)
		    !Puntos de las celdas para evaluar
		    c1=col2+DiagMov1(cont); f1=fil2+DiagMov2(cont)
		    c2=col2+DiagMov3(cont); f2=fil2+DiagMov3(5-cont)
		    !Evalua si el punto sirve como putno de arranque
		    if ((mascara(c1,f1).eq. 1 .and. mascara(c2,f2).eq.0).or.(mascara(c1,f1).eq.0.and. mascara(c2,f2).eq.1)) then
				!Finaliza la busqueda
				flag1=0
				!Guarda el punto para ser graficado
				print *, 'encontro por diagonales'
				perim_temp(1,1)=(col2-1+OrtoAdi(cont))*dx+xll_loc
				perim_temp(2,1)=(nrows+2-fil2+mov2(5-cont))*dx+yll_loc
		    endif
		    !se le suma un valor al contador
		    cont=cont+1
		enddo
    endif    
    !Una ves termina le resta al contador una unidad con el fin de identificar el punto de inicio
    cont=cont-1
    if (cont .eq. 3) then
		col2=col2+DiagMov3(cont)
		fil2=fil2+OrtoCol(cont)
    endif	
    !Busca los demas puntos
    flag1=1; cont2=2
    col=col2; fil=fil2
    do while (flag1.eq.1)
		!De acuerdo al valor de entrada gira los vectores de movimiento
		posN=cshift(pos,SHIFT=cont-1)
		!cambio=cshift(DiagMov3,SHIFT=cont-1)
		!Estados iniciales
		i=1; flag2=1	
		!Busca cual de las tres opciones es la que sirve
		do while (i.le.3 .and. flag2.eq.1)
		    !Puntos de las celdas para evaluar
		    c1=col2+mov1(posN(i)); f1=fil2+mov2(posN(i))
		    c2=col2+mov2(5-posN(i)); f2=fil2+mov1(posN(i))
		    !Evalua la condicion de existencia
		    if ((mascara(c1,f1).eq.1.and.mascara(c2,f2).eq.0).or.(mascara(c1,f1).eq.0.and.mascara(c2,f2).eq.1)) then
				!Toma la siguiente entrada
				cont=posN(i)
				!Guarda las coordenadas 
				perim_temp(1,cont2)=perim_temp(1,cont2-1)+dx*DiagMov3(cont)
				perim_temp(2,cont2)=perim_temp(2,cont2-1)-dx*OrtoCol(cont)
				cont2=cont2+1
				!Obtiene la nueva fila columna
				col2=col2+DiagMov3(cont)
				fil2=fil2+OrtoCol(cont)
				!indica la salida de la busqueda interna
				flag2=0					
		    endif
		    !si no cumple continua iterando
		    i=i+1
		enddo
		!observa si le nueva fil, col es igual a la inicial
		if (col2.eq.col.and.fil2.eq.fil) then
		    flag1=0
		endif
    enddo
    nperim=cont2-1
    perimetro=nperim*dxp/1000.0
end subroutine
subroutine basin_perim_cut(nperim,basin_perim) !Corta el perimetro de la cuenca
    !Variables de entrada
    integer, intent(in) :: nperim
    !Varaibles de salida
    real, intent(out) :: basin_perim(2,nperim)
    !f2py intent(in) :: nperim
    !f2py intent(out) :: basin_perim
    !copia y libera el vector
    basin_perim=perim_temp(:,1:nperim)
    deallocate(perim_temp)
end subroutine
subroutine basin_perim_kml(basin_p,ruta,nperim) !Escribe un kml con la cuenca
    !Variables de entrada
    integer, intent(in) :: nperim
    real, intent(in) :: basin_p(2,nperim)
    character*255, intent(in) :: ruta
    !Variables locales 
    character*10 Tarea, Tperim, Tpend, Tx, Ty, Telev
    character*50 TextCoord
    integer i
    !Convierte a texto el area el perimetro y la pendiente
    write(Tarea,'(F10.2)') area
    write(Tpend,'(F10.2)') pend_media*100
    write(Tperim,'(F10.2)') perimetro
    write(Telev,'(F10.2)') elevacion
    write(Tx,'(F10.4)') centroX
    write(Ty,'(F10.4)') centroY
    !Itera para escribir el kml de la cuenca     
    open(10,file=ruta,status='replace')
	write(10,'(A38)')'<?xml version="1.0" encoding="UTF-8"?>'
	write(10,*) '<kml xmlns="http://www.opengis.net/kml/2.2">'
	write(10,*) '<Document>'
	write(10,*) '		<Style>'
	write(10,*) '			<LineStyle>'
	write(10,*) '				<color>ff660000</color>'
	write(10,*) '				<width>3</width>'
	write(10,*) '			</LineStyle>'
	write(10,*) '			<PolyStyle>'
	write(10,*) '				<color>4C000000</color>'
	write(10,*) '				<fill>1</fill>'
	write(10,*) '				<outline>1</outline>'
	write(10,*) '			</PolyStyle>'
	write(10,*) '		</Style>'
	write(10,*) '	<Placemark>'
	write(10,*) '		<name>Cuenca Trazada</name>'
	write(10,*) '		<description>'
	write(10,*) '		<![CDATA['
	write(10,*) '			<html> <head> </head> <body> <div class=''popupEstaciones''>'
	write(10,*) '			<b>Propieades Geomorfol&oacute;gicas</b>'
	write(10,*) '			<hr>'
	write(10,*) '			<table border=''0'' width=''100%''>'
	write(10,*) '			<tr>'
	write(10,*) '				<td> <b>Area:</b> </td>'
	write(10,*) '				<td>'//Tarea//' Km<sup>2</sup> </td>'
	write(10,*) '			</tr>'
	write(10,*) '			<tr>'
	write(10,*) '				<td> <b>Pendiente Media:</b> </td>'
	write(10,*) '				<td>'//Tpend//' %</td>'
	write(10,*) '			</tr>'
	write(10,*) '			<tr>'
	write(10,*) '				<td> <b>Per&iacute;metro:</b> </td>'
	write(10,*) '				<td>'//Tperim//' Km</td>'
	write(10,*) '			</tr>'
	write(10,*) '			<tr>'
	write(10,*) '				<td> <b>Elevaci&oacute;n:</b> </td>'
	write(10,*) '				<td>'//Telev//' m.s.n.m</td>'
	write(10,*) '			</tr>'
	write(10,*) '			<tr>'
	write(10,*) '				<td> <b>Centroide:</b> </td>'
	write(10,*) '				<td>'//Tx//' Lon, '//Ty//' Lat</td>'
	write(10,*) '			</tr>'
	write(10,*) '		]]>'
	write(10,*) '		</description>'	
	write(10,*) '		<area>'//Tarea//'</area>'
	write(10,*) '		<elevacion>'//Telev//'</elevacion>'
	write(10,*) '		<centroX>'//TX//'</centroX>'
	write(10,*) '		<centroY>'//TY//'</centroY>'
	write(10,*) '		<ruta_bin>'//ruta//'</ruta_bin>'
	write(10,*) '		<Polygon>'
	write(10,*) '			<outerBoundaryIs>'
	write(10,*) '				<LinearRing>'
	write(10,*) '					<coordinates>'
	do i=1,nperim
	    write(TextCoord,'(1F15.5,A1,1F7.5,A1,I1)') basin_p(1,i),',',basin_p(2,i),',',0
	    TextCoord=ADJUSTL(TextCoord) 
	    TextCoord=trim(TextCoord)
	    write(10,*) TextCoord
	enddo
	write(10,*) '					</coordinates>'
	write(10,*) '				</LinearRing>'
	write(10,*) '			</outerBoundaryIs>'
	write(10,*) '		</Polygon>'
	write(10,*) '	</Placemark>'
	write(10,*) '</Document>'
	write(10,*) '</kml>'
    close(10)	
end subroutine 
!Funciones de interaccion con mapas
subroutine basin_float_var2map(basin_f,vect,MAPA,nc,nf,nceldas) !Convierte una variable de la cuenca en mapa
    !variables de entrada
    integer, intent(in) :: nc,nf,nceldas
    real, intent(in) :: vect(nceldas)
    integer, intent(in) :: basin_f(3,nceldas)
    !Variables de salida
    real, intent(out) :: MAPA(nc,nf)
    !f2py intent(in) :: basin_f,vect,nc,nf,nceldas
    !f2py intent(out) :: MAPA
    integer i
    MAPA=noData
    do i=1,nceldas
	MAPA(basin_f(2,i),basin_f(3,i))=vect(i)
    enddo
end subroutine
subroutine basin_int_var2map(basin_f,vect,MAPA,nc,nf,nceldas) !Convierte una variable de la cuenca en mapa
    !variables de entrada
    integer, intent(in) :: nc,nf,nceldas
    integer, intent(in) :: vect(nceldas)
    integer, intent(in) :: basin_f(3,nceldas)
    !Variables de salida
    integer, intent(out) :: MAPA(nc,nf)
    !f2py intent(in) :: basin_f,vect,nc,nf,nceldas
    !f2py intent(out) :: MAPA
    integer i
    MAPA=noData
    do i=1,nceldas
	MAPA(basin_f(2,i),basin_f(3,i))=vect(i)
    enddo
end subroutine
subroutine basin_int_map2var(basin_f,MAPA,vect,nc,nf,xll_m,yll_m,dx_m,noData_m,nceldas) !Convierte un mapa a una variable ordenada para la cuenca
    !Variables de entrada
    integer, intent(in) :: nc,nf,nceldas,noData_M
    integer, intent(in) :: basin_f(3,nceldas),Mapa(nc,nf)
    real, intent(in) :: xll_m,yll_m,dx_m
    !Variables de salida
    integer, intent(out) :: vect(nceldas)
    !f2py intent(in) :: nc,nf,nceldas,noData_M,basin_f,Mapa,xll_m,yll_m,dx_m
    !f2py intent(out) :: vect
    !Variables locales
    integer i,j,fila,columna,Media,cont
    real Xpos,Ypos,M
    !Calcula la media de los valores del mapa
    M=sum(Mapa,mask=Mapa.ne.noData_m)
    cont=count(Mapa.ne.noData_m)
    Media=ceiling(M/cont)
    !Asigna los valores del mapa al vector
    cont=0
    do i=1,nceldas
	!Para cada entrada de la tabla calcula el lat, long a partir de cual determina que valor le entra
	Xpos=xll+dx*(basin_f(2,i)-0.5)
	Ypos=yll+dx*((nrows-basin_f(3,i))+0.5)
	!Evalua si la posicion esta por dentro del mapa
	if (Xpos.gt.xll_m.and.Xpos.lt.(xll_m+dx_m*nc).and.Ypos.gt.yll_m.and.Ypos.lt.(yll_m+nf*dx_m)) then
	    !si esta por dentro le asigna la columna equivalente
	    columna=ceiling((Xpos-xll_m)/dx_m)
	    fila=ceiling((Ypos-yll_m)/dx_m)
	    fila=nf-fila+1
	    vect(i)=Mapa(columna,fila)
	    if (vect(i).eq.noDataM) then
		vect(i)=Media
	    endif
	else
	    !si esta por fuera del mapa le asigna la media de los valores del mapa
	    vect(i)=Media
	endif
    enddo
end subroutine
subroutine basin_float_map2var(basin_f,MAPA,vect,nc,nf,xll_m,yll_m,dx_m,noData_m,nceldas) !Convierte un mapa a una variable ordenada para la cuenca
    !Variables de entrada
    integer, intent(in) :: nc,nf,nceldas,noData_M
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: Mapa(nc,nf)
    real, intent(in) :: xll_m,yll_m,dx_m
    !Variables de salida
    real, intent(out) :: vect(nceldas)
    !f2py intent(in) :: nc,nf,nceldas,noData_M,basin_f,Mapa,xll_m,yll_m,dx_m
    !f2py intent(out) :: vect
    !Variables locales
    integer i,j,fila,columna,cont
    real Xpos,Ypos,Media
    !Calcula la media de los valores del mapa
    Media=sum(Mapa,mask=Mapa.ne.noData_m)
    cont=count(Mapa.ne.noData_m)
    Media=Media/cont
    !Asigna los valores del mapa al vector
    cont=0
    do i=1,nceldas
	!Para cada entrada de la tabla calcula el lat, long a partir de cual determina que valor le entra
	Xpos=xll+dx*(basin_f(2,i)-0.5)
	Ypos=yll+dx*((nrows-basin_f(3,i))+0.5)
	!Evalua si la posicion esta por dentro del mapa
	if (Xpos.gt.xll_m.and.Xpos.lt.(xll_m+dx_m*nc).and.Ypos.gt.yll_m.and.Ypos.lt.(yll_m+nf*dx_m)) then
	    !si esta por dentro le asigna la columna equivalente
	    columna=ceiling((Xpos-xll_m)/dx_m)
	    fila=ceiling((Ypos-yll_m)/dx_m)
	    fila=nf-fila+1
	    vect(i)=Mapa(columna,fila)
	    if (vect(i).eq.noData_m) then
		vect(i)=Media
	    endif
	else
	    !si esta por fuera del mapa le asigna la media de los valores del mapa
	    vect(i)=Media
	endif
    enddo
end subroutine
subroutine basin_map2basin(basin_f,nceldas,Mapa,xllM,yllM,ncolsM,nrowsM,dxM,noDataM,opcion,vec) !pasa un mapa a un vector con los valores correspondientes a la topologia de la cuenca
    !Variables de entrada
    integer, intent(in) :: ncolsM,nrowsM,nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: xllM,yllM,dxM,noDataM,Mapa(ncolsM,nrowsM)
    character*11, intent(in) :: opcion
    !Variables de salida
    real, intent(out) :: vec(nceldas)
    !f2py intent(in) :: ncolsM,nrowsM,nceldas,basin_f,xllM,yllM,dxM,noDataM,opcion,Mapa
    !f2py intent(out) :: vec
    !Variables internas
    integer i,j,cont,cantidad
    real Xpos,Ypos,MediaMapa
    integer fila,columna
    !Calcula la media del mapa
    cantidad=count(mapa.ne.noDataM)
    mediaMapa=sum(mapa,MASK=mapa.ne.noDataM)
    mediaMapa=mediaMapa/cantidad
    !Para cada entrada de la tabla determina cual es el valor que debe ir en el vector
    vec=noDataM
    cont=0
    do i=1,nceldas
		!Calcula la pos de la celda
		Xpos=xll+dx*(basin_f(2,i)-0.5)
		Ypos=yll+dx*((nrows-basin_f(3,i))+0.5)
		!Evalua si la posicion esta por dentro del mapa
		if (Xpos.gt.xllM.and.Xpos.lt.(xllM+dxM*ncolsM).and.Ypos.gt.yllM.and.Ypos.lt.(yllM+nrowsM*dxM)) then
		    !si esta por dentro le asigna la columna equivalente
		    columna=ceiling((Xpos-xllM)/dxM)
		    fila=ceiling((Ypos-yllM)/dxM)
		    fila=nrowsM-fila+1
		    vec(i)=Mapa(columna,fila)
		    !Si esta adentro pero es un valor no data hace una de las dos opciones
		    if (vec(i).eq.noDataM) then
				!if (opcion=='fill_mean') then !llena con la media de los valores del mapa
					vec(i)=mediaMapa
				!endif
		    endif
		endif
    enddo
end subroutine
subroutine basin_2map_find(basin,map_ncols,map_nrows,nceldas) !Determina los limites y la cantidad de filas y columnas de un mapa enmarcando la cuenca trazada
	!varaibles de entrada
	integer, intent(in) :: nceldas
	integer, intent(in) :: basin(3,nceldas)
	!Variables de salida
	integer, intent(out) :: map_ncols,map_nrows
	!f2py intent(in) :: nceldas,basin
	!f2py intent(out) :: map_ncols, map_nrows
	!variables locales
	integer col_min,col_max,fil_min,fil_max
	!Encuentra las columnas y filas minimas y maximas
	col_min=minval(basin(2,:)); col_max=maxval(basin(2,:))
	fil_min=minval(basin(3,:)); fil_max=maxval(basin(3,:))
	!Encuentra el numero de columnas y de filas nuevo
	map_ncols=col_max-col_min+1
	map_nrows=fil_max-fil_min+1
end subroutine
subroutine basin_2map(basin,var,mapa,map_ncols,map_nrows,map_xll,map_yll,nceldas) !Genera un mapa de la cuenca tomando los limites encontrados por basin_2map_find
	!Variables de entrada
	integer, intent(in) :: nceldas
	integer, intent(in) :: basin(3,nceldas)
	real, intent(in) :: var(nceldas)
	integer, intent(in) :: map_ncols,map_nrows
	!Variables de salida
	real, intent(out) :: mapa(map_ncols,map_nrows),map_xll,map_yll
	!f2py intent(in) :: nceldas, basin, var, map_ncols, map_nrows
	!f2py intent(out) :: mapa,map_xll,map_yll
	!Variables locales
	integer i,j
	integer col_min,col_max,fil_min,fil_max
	integer col_rel,fil_rel
	!Encuentra la fila columna maxima y minima
	col_min=minval(basin(2,:)); col_max=maxval(basin(2,:))
	fil_min=minval(basin(3,:)); fil_max=maxval(basin(3,:))
	!Encuentra el xll y el yll nuevos
	map_xll=xll+dx*(col_min-1)
	map_yll=yll+dx*(nrows-fil_max)
	!Aloja la matriz y comienza a llenarla de datos
	mapa=nodata
	do i=1,nceldas
		col_rel=basin(2,i)-col_min+1
		fil_rel=basin(3,i)-fil_min+1
		mapa(col_rel,fil_rel)=var(i)
	enddo
end subroutine
subroutine basin_point2var(basin_f,id_coord,xy_coord,res_coord,basin_pts,ncoord,nceldas) !Obtiene el vector basin_pts con los puntos de control 
    !Variables de entrada
    integer, intent(in) :: nceldas,ncoord
    integer, intent(in) :: basin_f(3,nceldas),id_coord(ncoord)
    real, intent(in) :: xy_coord(2,ncoord)
    !Variables de salida
    integer, intent(out) :: basin_pts(nceldas),res_coord(ncoord)    
    !f2py intent(in) :: nceldas,ncoord,basin_f,xy_coord
    !f2py intent(out) :: basin_pts,res_coord
    !Variables locales
    integer i,j,x_col,y_fil,esta,posit
    real x,y,m
    !Recorre todos los puntos de entrada
    !basin_pts=0
    xy_new=noData
    do i=1,ncoord
	!Calcula la fila columna del punto
	x=(xy_coord(1,i)-xll)/dx
	x_col=ceiling(x)
	y=nrows-(xy_coord(2,i)-yll)/dx
	y_fil=ceiling(y)
	!Entrega la posicion dentro del vector
	call find_xy_in_basin(basin_f,x_col,y_fil,posit,nceldas)
	!Evalua si esta en la cuenca	
	if (posit.gt.0) then
	    res_coord(i)=0 !El punto esta dentro de la cuenca
	    basin_pts(posit)=id_coord(i)
	else
	    res_coord(i)=1
	endif    	    
    enddo
end subroutine
subroutine basin_extract_var_by_point(basin_f,var,xy_coord,kernel,var_values,ncoord,nceldas) !Entrega el valor de una var de la cuenca a partir de puntos
	!Variables de entrada
    integer, intent(in) :: nceldas,ncoord,kernel
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: xy_coord(2,ncoord),var(nceldas)
    !Variables de salida
    real, intent(out) :: var_values(ncoord)    
    !f2py intent(in) :: nceldas,ncoord,basin_f,xy_coord,kernel,var
    !f2py intent(out) :: var_values
    !Variables locales
	integer i,x_col,y_fil,posit,c,f,cont,posit_temp,k
	real x,y,valor
	!Si se da un kernel malo se corrige y se avisa
	if (kernel .lt. 3) then
		k=3
		print *, 'Alerta: Kernel dado inferior a 3, kernel adoptado igual a 3'
	else
		k=kernel
	endif
	k=floor(k/2.0)
	!Busca en cada uno de los puntos 
	var_values=-9999.0	
	do i=1,ncoord
		x=(xy_coord(1,i)-xll)/dx
		x_col=ceiling(x)
		y=nrows-(xy_coord(2,i)-yll)/dx
		y_fil=ceiling(y)
		!Entrega la posicion dentro del vector
		call find_xy_in_basin(basin_f,x_col,y_fil,posit,nceldas)
		!solo evalua si el punto esta dentro de la cuenca
		if (posit .ne. 0) then
			cont=0
			valor=0
			do c=-k,k
				do f=-k,k
					call find_xy_in_basin(basin_f,x_col+c,y_fil+f,posit_temp,nceldas)
					if (posit .ne. 0) then
						cont=cont+1
						valor=valor+var(posit_temp)
					endif
				enddo
			enddo
			var_values(i)=valor/cont
		endif
	enddo
end subroutine
!subroutine basin_var2smooth(basin_f,var,varOut,nc,nf,nceldas,kernel=3)!Suaviza una variable tomando un tamano de kernel variable (kernel cuadrado)
!	!Variables de entrada
!	integer, intent(in) :: nf,nc,nceldas,kernel
!	integer, intent(in) :: basin_f(3,nceldas)
!	real, intent(in) :: var(nceldas)
!	!variables de salida
!	real, intent(out) :: varOut(nceldas)
!	!f2py intent(in) :: nf,nc,nceldas,basin_f,var
!	!f2py intent(out) :: varOut
!	!Variables locales 
!	integer i,j
!	real mapa1(nc,nf), mapa2(nc,nf)
!	!Convierte el vector a mapa
!	call basin_float_var2map(basin_f,var,mapa1,nc,nf,nceldas)
!	!itera sobre el mapa suavizando
!	do i=1,nc
!		do j=1,nf
!			if ()
!		enddo
!	enddo
!	!convierte el mapa suavizado en el vector de salida
!	call basin_float_map2var(basin_f,mapa2,varOut,cu.ncols,cu.nrows,cu.xll,cu.yll,cu.dx,cu.noData,nceldas)
!end	subroutine
!Funciones de propiedades de corrientes
subroutine basin_stream_nod(basin_f,acum,nceldas,umbral,cauce,nodos,trazado,n_nodos,n_cauce) !obtiene vectores: celdas cauce, nodos hidrologicos, puntos de trazado
    !variables de entrada
    integer, intent(in) :: nceldas !numero de celdas que componen la cuenca
    integer, intent(in) :: basin_f(3,nceldas),acum(nceldas) !vector de cuenca y de celdas acum
    integer, intent(in) :: umbral !umbral de area a partir del cual se considera cauce
    !Variables de salida
    integer, intent(out) :: cauce(nceldas),nodos(nceldas),trazado(nceldas),n_nodos,n_cauce
    !f2py intent(in) :: nceldas,basin_f,acum,umbral
    !f2py intent(out) :: cauce,nodos,trazado,n_nodos,n_cauce
    integer i,drenaid,contador
    !Encuentra las celdas que son cauce
    cauce=0    
    where(acum.ge.umbral) cauce=1
    n_cauce=count(cauce.eq.1)
    !Encuentra los nodos
    nodos=0
    n_nodos=0
    do i=1,nceldas
	!Determina la celda a la que se drena
	drenaid=nceldas-basin_f(1,i)+1
	!Calcula el area acumulada
	if (basin_f(1,i).ne.0 .and. cauce(i)==1) then
	    nodos(drenaid)=nodos(drenaid)+1
	endif
    enddo
    !Genera los nodos de nacimientos y elimina los demas
    where(nodos.eq.0 .and. cauce.eq.1) nodos=3
    nodos(nceldas)=2  !La salida tambien es un nodo
    where(nodos.lt.2) nodos=0 !Hace nodos todos los puntos intermedios    
    n_nodos=count(nodos.gt.0)
    !Encuentra los puntos de trazado
    trazado=0
    trazado(nceldas)=1
    do i=1,nceldas
	drenaid=nceldas-basin_f(1,i)+1
	if (basin_f(1,i).ne.0 .and. nodos(drenaid).ne.0 .and. cauce(i).eq.1) then
	    trazado(i)=1
	endif
    enddo
end subroutine
subroutine basin_stream_slope(basin_f,basin_elev,basin_long,nodos,n_cauce,stream_s,stream_l,nceldas) !Obtiene la pendiente por tramo de cauce
    !Variables de entrada
    integer, intent(in) :: nceldas,n_cauce
    integer, intent(in) :: basin_f(3,nceldas),nodos(nceldas)
    real, intent(in) :: basin_elev(nceldas),basin_long(nceldas)
    !Variables de salida
    real, intent(out) :: stream_s(nceldas),stream_l(nceldas)
    !f2py intent(in) :: nceldas,basin_f,nodos,basin_elev,basin_long,n_cauce
    !f2py intent(out) :: stream_s,stream_l
    !Variables locales
    integer i,j,cont,drenaid,celdas_tramo(n_cauce),flag
    real med_elev,med_long,S_tramo,X2,XY,elev_dist(2,n_cauce)
    !Itera por las celdas
    stream_s=noData
    stream_l=noData
    do i=1,nceldas
	!Si la celda es nodo y no es el nodo de salida evalua hasta la salida
	if (nodos(i).gt.0 .and. basin_f(1,i).ne.0) then
	    !Condiciones del nodo icial y bandera en estado encendido	    	    
	    elev_dist(1,1)=basin_elev(i)
	    elev_dist(2,1)=basin_long(i)
	    celdas_tramo(1)=i	    
	    flag=1; cont=1
	    drenaid=nceldas-basin_f(1,i)+1
	    !Recorre hasta llegar a otro nodo
	    do while (flag.eq.1) 		
		!Verifica si la siguiente celda es nodo o no
		if (nodos(drenaid).eq.0) then
		    !si no es acumula longitud y guarda elevacion y acumulado de long
		    cont=cont+1
		    elev_dist(1,cont)=basin_elev(drenaid)
		    elev_dist(2,cont)=elev_dist(2,cont-1)+basin_long(drenaid)		   
		    celdas_tramo(cont)=drenaid 
		else
		    !si es nodo para de acumular por lo que apaga la bandera
		    flag=0
		    !si fue nodo de una toma valores 
		    if (cont.eq.1) then
			cont=cont+1
			elev_dist(1,cont)=basin_elev(drenaid)
			elev_dist(2,cont)=elev_dist(2,cont-1)+basin_long(drenaid)
		    endif
		endif
		!Calcula la nueva celda destino
		drenaid=nceldas-basin_f(1,drenaid)+1
	    enddo
	    !Calcula la media de la elevacion y la media de la longitud total	    
	    S_tramo=abs((elev_dist(1,1)-elev_dist(1,cont))/(elev_dist(2,1)-elev_dist(2,cont)))
	    if (S_tramo.le.0) S_tramo=0.001
	    L_tramo=elev_dist(2,cont)
	    if (L_tramo.le.0) L_tramo=dxP
	    !Le asigna esa pendiente a todas las celdas del tramo
	    do j=1,cont
		stream_s(celdas_tramo(j))=S_tramo
		stream_l(celdas_tramo(j))=L_tramo
	    enddo	    
	endif
    enddo
end subroutine
subroutine basin_stream_type(basin_f,acum,umbrales,stream_types,numbrales,nceldas) !Obtiene un vector con los tipos de celdas que son la cuenca (usado para modelacion), umbrales en orden ascendente
    !Variables de entrada
    integer, intent(in) :: nceldas,numbrales
    integer, intent(in) :: basin_f(nceldas),acum(nceldas)
    integer, intent(in) :: umbrales(numbrales)
    !Variables de salida
    integer, intent(out) :: stream_types(nceldas)
    !f2py intent(in) :: nceldas,basin_f,acum,umbrales,numbrales
    !f2py intent(out) :: stream_types
    !Variables locales
    integer i,j
    !Obtiene los tipos de celdas
    !1: celdas tipo ladera
    !2: celdas tipo carcava o gully
    !3: celdas tipo cauce
    !Si se pone una cantidad de umbrales diferente a 2 no se hace caso a estas catgorias
    stream_types=1
    do i=1,numbrales
	do j=1,nceldas
	    if (acum(j).gt.umbrales(i)) then
		stream_types(j)=i+1
	    endif
	enddo
    enddo
end subroutine
subroutine basin_stream_point2stream(basin_f,cauce,id_coord,xy_coord,res_coord,basin_pts,xy_new,ncoord,nceldas) !Obtiene el vector basin_pts con los puntos de control ubicados
    !Variables de entrada
    integer, intent(in) :: nceldas,ncoord
    integer, intent(in) :: basin_f(3,nceldas),cauce(nceldas),id_coord(ncoord)
    real, intent(in) :: xy_coord(2,ncoord)
    !Variables de salida
    integer, intent(out) :: basin_pts(nceldas),res_coord(ncoord)
    real, intent(out) :: xy_new(2,ncoord)
    !f2py intent(in) :: nceldas,ncoord,basin_f,xy_coord,cauce
    !f2py intent(out) :: basin_pts,res_coord,xy_new
    !Variables locales
    integer i,j,x_col,y_fil,esta,posit
    real x,y,m
    !Recorre todos los puntos de entrada
    !basin_pts=0
    xy_new=noData
    do i=1,ncoord
	!Calcula la fila columna del punto
	x=(xy_coord(1,i)-xll)/dx
	x_col=ceiling(x)
	y=nrows-(xy_coord(2,i)-yll)/dx
	y_fil=ceiling(y)
	!Entrega la posicion dentro del vector
	call find_xy_in_basin(basin_f,x_col,y_fil,posit,nceldas)
		!Evalua si esta en la cuenca
	if (posit.gt.0) then
	    res_coord(i)=1 !El punto esta dentro de la cuenca
	    !Evalua si esta en el cauce
	    if (cauce(posit).gt.0) then
		!El punto esta en el cauce de la cuenca
		res_coord(i)=2 
		basin_pts(posit)=id_coord(i)
	    else
		!El punto no esta en el cauce, va a buscar el cauce mas proximo
		do while (cauce(posit).eq.0 .and. posit.lt.nceldas)
		    posit=nceldas+1-basin_f(1,posit)
		enddo		
			if (posit.lt.nceldas) then 
			    !si termino bien asigna el punto 
			    basin_pts(posit)=id_coord(i)
			    res_coord(i)=2
			endif
	    endif
	else
	    res_coord(i)=0 !El punto esta por fuera de la cuenca
	endif
	!Si la coordenada quedo dentro del cauce grava el punto xy
	if (res_coord(i).eq.2) then
		xy_new(1,i)=xll+dx*(basin_f(2,posit)-0.5)
		xy_new(2,i)=yll+dx*((nrows-basin_f(3,posit))+0.5)
	endif
    enddo
end subroutine
!funciones para trazar la red hidrica
subroutine basin_netxy_find(basin_f,nodos,cauceHort,nceldas,netsize) !Funcion para obtener las coordenadas de la red
	!Variables de entrada
	integer, intent(in) :: nceldas
	integer, intent(in) :: nodos(nceldas),basin_f(3,nceldas)
	real, intent(in) :: cauceHort(nceldas)
	!Variables de salida
	integer, intent(out) :: netsize
	!Variables locales
	integer i,drenaid,cont
	real X(nceldas), Y(nceldas)
	logical flag1,flag2
	!f2py intent(in) :: nceldas,nodHort,basin_f,cauceHort
	!f2py intent(out) :: netsize
	!Aloja la red temporal
	if (.not. allocated(netxy_temp)) allocate(netxy_temp(3,nceldas))
	netxy_temp=-999
	!Obtiene las coordenadas de todas las celdas
	call basin_coordXY(basin_f,X,Y,nceldas)
	!Busca las coordenadas de la red
	cont=0
	do i=1,nceldas
		if (nodos(i).ne.0 .and. i .lt. nceldas) then
			cont=cont+1
			drenaid=i
			flag1=.true.
			flag2=.false.
			para=0
			do while (flag1)
				netxy_temp(1,cont)=cauceHort(drenaid)
				netxy_temp(2,cont)=X(drenaid)
				netxy_temp(3,cont)=Y(drenaid)
				drenaid=nceldas-basin_f(1,drenaid)+1
				cont=cont+1
				if (flag2) flag1=.false.
				if (nodos(drenaid).ne.0 .or. drenaid .eq. nceldas) flag2=.true.
			enddo
		endif
	enddo
	netsize=cont
end subroutine
subroutine basin_netxy_cut(netsize,nceldas,net)
	!Variables de entrada
	integer, intent(in) :: nceldas, netsize
	!Variables de salida
	real, intent(out) :: net(3,netsize)
	!f2py intent(in) :: nceldas, netsize
	!f2py intent(out) :: net
	net=netxy_temp(:,:netsize)
end subroutine
!funciones de cauce principal
subroutine basin_ppalstream_find(basin_f,nodos,longCeldas,elev,nceldas,ppal_nceldas,punto)
	!Variables de entrada
	integer, intent(in) :: nceldas
	integer, intent(in) :: nodos(nceldas),basin_f(3,nceldas)
	real, intent(in) :: longCeldas(nceldas),elev(nceldas)
	!Variables de salida
	integer, intent(out) :: ppal_nceldas,punto
	!Variables locales
	integer i,cont,drenaid
	real LongMax,Long,X(nceldas),Y(nceldas)
	logical flag
	!f2py intent(in) :: nceldas,nodos,basin_f,longCeldas,elev
	!f2py intent(out) :: ppal_nceldas,punto
	!Aloja el vector temporal 
	if (.not. allocated(ppal_stream_temp)) allocate(ppal_stream_temp(4,nceldas))
	ppal_stream_temp=-999
	!Inicia la longitud maxima
	LongMax=0.0
	call basin_coordXY(basin_f,X,Y,nceldas)
	!Itera para todas las celdas y enuentra el cauce ppal
	do i=1,nceldas
		if (nodos(i).eq.3) then
			drenaid=nceldas-basin_f(1,i)+1
			Long=longCeldas(i)
			cont=1
			do while (drenaid .le. nceldas)
				Long=Long+longCeldas(drenaid)				
				drenaid=nceldas-basin_f(1,drenaid)+1
				cont=cont+1
			enddo
			if (Long .gt. LongMax) then
				LongMax=Long
				punto=i
				ppal_nceldas=cont
			endif
		endif
	enddo
	!Obtiene las propiedades del cuace principal
	drenaid=punto; cont=0; Long=0
	flag=.true.
	do while (flag)
		cont=cont+1
		Long=Long+longCeldas(drenaid)
		ppal_stream_temp(1,cont)=elev(drenaid)
		ppal_stream_temp(2,cont)=Long
		ppal_stream_temp(3,cont)=X(drenaid)
		ppal_stream_temp(4,cont)=Y(drenaid)
		drenaid=nceldas-basin_f(1,drenaid)+1
		if (drenaid .gt. nceldas) flag=.false.
	enddo
end subroutine
subroutine basin_ppalstream_cut(ppal_nceldas,nceldas,ppal_f)
	!Variables de entrada
	integer, intent(in) :: nceldas, ppal_nceldas
	!Variables de salida
	real, intent(out) :: ppal_f(4,ppal_nceldas)
	!f2py intent(in) :: nceldas, ppal_nceldas
	!f2py intent(out) :: ppal_f
	ppal_f=ppal_stream_temp(:,:ppal_nceldas)
end subroutine
subroutine basin_ppal_hipsometric(basin_f,elev,ppal_punto,intervalos,nceldas,ppal_nceldas,ppal_hipso,basin_hipso)
	!Variables de entrada
	integer, intent(in) :: nceldas,ppal_nceldas,ppal_punto,intervalos
	integer, intent(in) :: basin_f(3,nceldas)
	real, intent(in) :: elev(nceldas)
	!Variables de salida
	real, intent(out) :: ppal_hipso(2,ppal_nceldas),basin_hipso(2,intervalos)
	!Variables locales
	integer cont,drenaid,acum(nceldas)
	real TamIntervalo,maxElev
	!f2py intent(in) :: nceldas,ppal_nceldas, elev,basin_f, ppal_punto,intervalos
	!f2py intent(out) :: ppal_hipso
	call basin_acum(basin_f,nceldas,acum)
	!Calcula la curva por el cuace ppal
	drenaid=ppal_punto
	do cont=1,ppal_nceldas
		ppal_hipso(1,cont)=(acum(drenaid)*dxp**2)/1e6
		ppal_hipso(2,cont)=elev(drenaid)
		drenaid=nceldas-basin_f(1,drenaid)+1
	enddo
	!Calcula la curva para toda la cuenca
	maxElev=maxval(elev,dim=1)
	minElev=minval(elev,dim=1)
	TamIntervalo=(maxElev-minElev)/intervalos
	Cant=0
	do i=1,intervalos		
		Cant=Cant+count(elev .ge. minElev+TamIntervalo*i .and. elev .lt. minElev+TamIntervalo*(i+1),dim=1)
		basin_hipso(1,i)=nceldas-Cant
		basin_hipso(2,i)=(minElev+TamIntervalo*i+minElev+TamIntervalo*(i+1))/2.0				
	enddo
end subroutine
!Funciones de balance y regionalizacion
subroutine basin_qmed(basin_f,elev,precip,qmed,ETR,nceldas,etr_type &
	&, mu_choud) !Calcula el caudal medio de largo plazo por el metodo de turc
    !Variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas),elev(nceldas)
    real, intent(in) :: precip(nceldas)
    integer, intent(in) :: etr_type
    real, intent(in) :: mu_choud
    !Variables de salida
    real, intent(out) :: qmed(nceldas),ETR(nceldas)
    !Variables locales 
    real temp(nceldas),L(nceldas),razon(nceldas)
    real ETP(nceldas), Rn(nceldas),alpha,mu
    integer drenaid,i
    !f2py intent(in) :: nceldas, basin_f,precip,acum,elev
    !f2py intent(in), optional :: etr_type
    
    !f2py intent(out) :: qmed,etr    
    if (etr_type.eq. 1) then
	    !Calcula la ETR por turc
	    temp=27.72-0.0055*elev
	    L=300+25*temp+0.05*temp**3
	    razon=precip/L
	    where(razon.gt.0.316) ETR=precip/sqrt(0.9+((precip**2)/(L**2)))
	    where(razon.le.0.316) ETR=precip
	elseif (etr_type.eq. 2) then
		!Calcula por Cenicafe Budyko
		ETP=1700.17*exp(-0.0002*elev)
		ETR=(ETP*precip*tanh(precip/ETP)*(1-cosh(ETP/precip)+sinh(ETP/precip)))**0.5
	elseif (etr_type .eq. 3) then
		!Choudhury
		if (mu_choud.lt.0.85 .or. mu_choud.gt.1.9) then 
			mu=1.37 
		else
			mu=mu_choud
		endif
		Rn=precip/mu
		ETR=precip/(1+(precip/Rn)**1.91)**(1/1.91)	
	endif
    !Calula el caudal
    qmed = (dxp**2)*(precip-ETR)/31536000000.0
	do i=1,nceldas
		!Determina la celda a la que se drena
		drenaid=nceldas-basin_f(1,i)+1
		!Calcula el area acumulada
		if (basin_f(1,i).ne.0) then
		    qmed(drenaid)=qmed(drenaid)+qmed(i)			    
		endif	
    end do
    !qmed=qmed*(dxp**2)/31536000000.0
end subroutine
subroutine basin_qofer_qcap(basin_f,q_oferta,q_cap,qres,escazes,nceldas) !Resta el cadual de captacion al oferta
    !Variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas)
    real, intent(in) :: q_cap(nceldas),q_oferta(nceldas)
    !Variables de salida
    real, intent(out) :: qres(nceldas),escazes(nceldas)
    !Variables locales 
    integer i,drenaid
    real es
    !f2py intent(in) :: nceldas, basin_f,q_oferta,q_cap
    !f2py intent(out) :: qres,escazes
    !Acumula caudal una ves retirado el caudal captado
    qres=q_oferta
!    do i=1,nceldas
!		if (q_cap(i).gt.0) then
!		    drenaid=nceldas-basin_f(1,i)+1
!		    do while (drenaid.le.nceldas)
!				qres(drenaid)=qres(drenaid)-q_cap(i)
!				drenaid=nceldas-basin_f(1,drenaid)+1
!		    enddo
!		endif
!    enddo
    !Calcula el indice de escazes
    escazes=0
    do i=1,nceldas
		if (q_cap(i).gt.0) then
		    es = 100*(q_cap(i)/qres(i))
		    if (es .le. 100) escazes(i) = es		   
		    drenaid=nceldas-basin_f(1,i)+1
		    do while (drenaid.le.nceldas .and. q_cap(drenaid).eq.0.0)
				es=100*(q_cap(i)/qres(drenaid))
				if (es .le. 100) escazes(drenaid) = es
				drenaid=nceldas-basin_f(1,drenaid)+1
		    enddo
		endif
    enddo
end subroutine
!Funciones de acumulacion de flujo 
subroutine basin_propagate(basin_f,var,var_prop,nceldas)
	!Variables de entrada
	integer, intent(in) :: nceldas
	integer, intent(in) :: basin_f(3,nceldas)
	real, intent(in) :: var(nceldas)
	!Variables de salida 
	real, intent(out) :: var_prop(nceldas)
	!f2py intent(in) :: nceldas, basin_f, var
	!f2py intent(out) :: var_prop
	!Variables locales 
	integer i,drenaid
	!Transito de la variable aguas abajo
	var_prop = var
	do i=1,nceldas
		!propaga si hay algo que propagar
		if (var_prop(i).gt.0) then 
			drenaid=nceldas-basin_f(1,i)+1
			if (drenaid .le. nceldas) then 
				var_prop(drenaid) = var_prop(drenaid) + var_prop(i)
			endif
		endif
	enddo
end subroutine
!funciones de sub-cuencas (nodos y laderas)
subroutine basin_subbasin_nod(basin_f,acum,nceldas,umbral,cauce,nodos_fin,n_nodos) !obtiene vectores: celdas cauce, nodos hidrologicos, puntos de trazado
    !variables de entrada
    integer, intent(in) :: nceldas !numero de celdas que componen la cuenca
    integer, intent(in) :: basin_f(3,nceldas),acum(nceldas) !vector de cuenca y de celdas acum
    integer, intent(in) :: umbral !umbral de area a partir del cual se considera cauce
    !Variables de salida
    integer, intent(out) :: cauce(nceldas),nodos_fin(nceldas),n_nodos
    !f2py intent(in) :: nceldas,basin_f,acum,umbral
    !f2py intent(out) :: cauce,nodos,n_nodos
    integer i,drenaid,contador,nodos(nceldas),cont,cont2,posi
    logical flag1,flag2
    !Encuentra las celdas que son cauce
    cauce=0    
    where(acum.ge.umbral) cauce=1
    n_cauce=count(cauce.eq.1)
    !Encuentra los nodos
    nodos=0
    n_nodos=0
    do i=1,nceldas
		!Determina la celda a la que se drena
		drenaid=nceldas-basin_f(1,i)+1
		!Calcula el area acumulada por la red
		if (basin_f(1,i).ne.0 .and. cauce(i)==1) then
		    nodos(drenaid)=nodos(drenaid)+1
		endif
    enddo         
    !Determina quienes son los nodos
    where(nodos.lt.2) nodos=0 
    !Encuentra los nodos de verdad
    nodos_fin=0; nodos_fin(nceldas)=1    
    do i=1,nceldas
		!Verifica que la celda sea cauce
		if (cauce(i).eq.1) then
		    !Determina la celda a la que se drena
		    drenaid=nceldas-basin_f(1,i)+1
		    !Si la celda destino es nodo esta se hace nodo final
		    if (drenaid .le. nceldas) then !si la celda destino no es la ultima de la cuenca
				if (nodos(drenaid).ne.0) then
				    nodos_fin(i)=1
				endif
		    endif
		endif
    enddo    
    n_nodos=sum(nodos_fin)
    !Si no esta alojado el vector de sub-cuencas lo aloja, pone las condiciones del nodo de salida
    if (allocated(sub_basins_temp) .eqv. .true.) deallocate(sub_basins_temp)
    if (.not. allocated(sub_basins_temp)) allocate(sub_basins_temp(2,n_nodos))
    sub_basins_temp=0
    sub_basins_temp(1,n_nodos)=nceldas; sub_basins_temp(2,n_nodos)=0   
    !Comienza a iterar de abajo hacia arriba 
    cont=1;cont2=0; flag1=.true.
    do while (flag1)
	!Examina si hay nodo para evaluar o no
		if (sub_basins_temp(1,n_nodos-cont+1).ne.0) then	   
		    !Si hay nodo toma su posicion para buscar otros
		    posi=sub_basins_temp(1,n_nodos-cont+1)
		    do j=1,posi-1
				if (nodos_fin(posi-j).ne.0) then
				    !si es nodo averigua hasta abajo si este le drena al nodo anterior
				    drenaid=nceldas-basin_f(1,posi-j)+1 ; flag2=.true.
				    do while (flag2)
						if (nodos_fin(drenaid).eq.0) then			
						    if (drenaid.ne.posi) then
								drenaid=nceldas-basin_f(1,drenaid)+1					    
						    endif
						elseif (nodos_fin(drenaid).ne.0 .and. drenaid.eq.posi) then
						    cont2=cont2+1
						    sub_basins_temp(1,n_nodos-cont2)=posi-j
						    sub_basins_temp(2,n_nodos-cont2)=cont
						    nodos_fin(posi-j)=cont2+1
						    flag2=.false.			    
						else
						    flag2=.false.
						endif
				    enddo		    
				endif
		    enddo
		    !Cuando termina actualiza el contador
		    cont=cont+1
		    if (cont.gt.n_nodos) flag1=.false.
		else
		    !Si no hay mas nodos para evaluar termina la ejecucion
		    flag1=.false.
		endif
    enddo
end subroutine
subroutine basin_subbasin_cut(n_nodos,sub_basins) !Corta el vector de la topologia de las sub-cuencas
	!Variables de entrada
	integer, intent(in) :: n_nodos
	!Variables de salida
	integer, intent(out) :: sub_basins(2,n_nodos)
	!f2py intent(in) :: n_nodos
	!f2py intent(out) :: sub_basins
	sub_basins=sub_basins_temp
	deallocate(sub_basins_temp)
end subroutine
subroutine basin_subbasin_horton(sub_basins,sub_horton,nod_horton,n_nodos,nceldas) !Encuentra el orden de horton de cada nodo
	!Variables de entrada
	integer, intent(in) :: n_nodos,nceldas
	integer, intent(in) :: sub_basins(2,n_nodos)
	!Variables de salida
	integer, intent(out) :: sub_horton(n_nodos),nod_horton(nceldas)
	!f2py intent(in) :: n_nodos,nceldas,sub_basins
	!f2py intent(out) :: sub_horton,nod_horton
	!Variables locales
	integer i,j,posi,drenantes(10),cont,min_val,max_val
	!Inicio de variables de salida
	sub_horton=0
	!Identifica los nodos de orden 1
	do i=1,n_nodos
		posi=n_nodos-i+1
		conteo=count(sub_basins(2,:).eq.posi)
		if (conteo.eq.0) then
			sub_horton(i)=1
			nod_horton(sub_basins(1,i))=1
		endif
	enddo
	!Itera de arriba a abajo acumulando ordenes
	do i=1,n_nodos
	    !Evalua si no es un nodo orden 1
	    if (sub_horton(i).ne.1) then
		!Si esta vacio determina la posicion e itera hacia arriba buscando quienes le drenan
		posi=n_nodos-i+1
		cont=0; drenantes=0
		!Encuentra los ordenes de quienes le drenan
		do j=1,i-1				
		    if (sub_basins(2,j).eq.posi) then
			cont=cont+1
			drenantes(cont)=sub_horton(j)
		    endif
		enddo
		!Determina si su orden se incrementa o no
		max_val=maxval(drenantes(1:cont))
		min_val=minval(drenantes(1:cont))
		if (max_val.eq.min_val) then
		    sub_horton(i)=max_val+1
		else
		    sub_horton(i)=max_val
		endif
		!ASigna ese orden al nodo de la cuenca
		nod_horton(sub_basins(1,i))=sub_horton(i)
	    endif
	enddo
end subroutine
subroutine basin_subbasin_find(basin_f,nodos,sub_pert,sub_basin,n_nodos,nceldas) !Determina las laderas de cada uno de los nodos
    !Variables de entrada
    integer, intent(in) :: n_nodos,nceldas,basin_f(3,nceldas),nodos(nceldas)    
    !Variables de salida
    integer, intent(out) :: sub_pert(nceldas),sub_basin(n_nodos)
    !f2py intent(in) :: nceldas,basin_f,nodos
    !f2py intent(out) sub_pert,sub_ncel
    !Variables locales
    integer i,j,numero,drenaid,cont,celdas_grupo(nceldas)
    logical flag    
    !Itera de arriba hacia abajo y encuentra la pertenencia de las celdas a las sub-cuencas
    sub_pert=nodos; sub_ncel=0
    do i=1,nceldas	
	!Verifica si la celda ya pertenece a un nodo
	if (sub_pert(i).eq.0) then
	    !Si no pertenece a ningun nodo determina a donde drena
	    drenaid=nceldas-basin_f(1,i)+1
	    !Itera hasta encontrar un nodo
	    flag=.true.
	    do while(flag)
		if (nodos(drenaid).eq.0) then
		    !Si no es nodo esa celda es guardada en en el grupo y sigue buscando hacia abajo		
		    drenaid=nceldas-basin_f(1,drenaid)+1		    
		else 
		    !si es nodo: se asigna a las celdas recorridas el numero del nodo al cual se apunta
		    sub_pert(i)=nodos(drenaid)
		    flag=.false.		    
		endif
	    enddo
	endif
    enddo
end subroutine
subroutine basin_subbasin_long(sub_pert,cauce,long,sub_basin,sub_horton,sub_basin_long,max_long,nodo_max_long,n_nodos,nceldas) !Determina la longitud del cauce total en cada nodo
	!Variables de entrada
	integer, intent(in) :: n_nodos,nceldas
	integer, intent(in) :: sub_pert(nceldas),cauce(nceldas),long(nceldas)
	integer, intent(in) :: sub_basin(2,n_nodos),sub_horton(n_nodos)
	!Variables de salida
	real, intent(out) :: sub_basin_long(n_nodos),max_long
	integer, intent(out) :: nodo_max_long
	!f2py intent(in) :: n_nodos,nceldas,sub_pert ,cauce,long, sub_basin, sub_horton
	!f2py intent(out) :: sub_basin_long,max_long,nodo_max_long
	!Variables locales
	integer i,j, cauce_pert(nceldas),drenaid
	real cauce_long(nceldas),long_acum
	!Calcula el cauce con longitud y las pertenencias de cada cauce
	cauce_long=cauce*long
	cauce_pert=cauce*sub_pert
	!Itera para cada nodo
	do i=1,n_nodos
		posi=n_nodos-i+1
		sub_basin_long(i)=sum(cauce_long,MASK=cauce_pert.eq.posi)
	enddo
	!Calcula la maxima longitud
	max_long=0; nodo_max_long=0
	do i=1,n_nodos
		!SI es horton=1 busca hasta abajo la longitud
		if (sub_horton(i).eq.1) then
			!Calcula la longitud desde ese punto de salida
			long_acum=sub_basin_long(i)
			drenaid=n_nodos-sub_basin(2,i)+1
			do while (sub_basin(2,drenaid).ge.1) 
				long_acum=long_acum+sub_basin_long(drenaid)
				drenaid=n_nodos-sub_basin(2,drenaid)+1				
			enddo
			!Compara a ver si es mayor que la maxima
			if (long_acum.gt.max_long) then
				max_long=long_acum
				nodo_max_long=n_nodos-i+1
			endif
		endif
	enddo
end subroutine
subroutine basin_subbasin_map2subbasin(sub_pert,basin_var,subbasin_sum,&
	&n_nodos,nceldas,cauce,sum_mean) !Agrega una variable de la cuenca a laderas
	!Varialbes de entrada
	integer, intent(in) :: n_nodos,nceldas,sum_mean
	integer, intent(in) :: sub_pert(nceldas)
	real, intent(in) :: basin_var(nceldas)
	integer, intent(in), optional :: cauce(nceldas)
	!Variables de salida
	real, intent(out) :: subbasin_sum(n_nodos)
	!f2py intent(in) n_nodos,nceldas,sub_pert,basin_var,sum_mean
	!f2py intent(out) :: subbasin_var 
	!Variables locales
	integer i,posi,cont_valores 
	real suma_valores
	!Localiza 
	do i=1,n_nodos
		posi=n_nodos-i+1
		!Cambia las condiciones de acuerdo a si esta o no el cauce 
		if (present(cauce)) then 
			suma_valores=sum(basin_var,&
				&MASK=sub_pert.eq.posi.and.cauce.eq.1)
			cont_valores=count(sub_pert.eq.posi.and.cauce.eq.1)
		else
			suma_valores=sum(basin_var,MASK=sub_pert.eq.posi)
			cont_valores=count(sub_pert.eq.posi)
		endif
		if (cont_valores .gt. 0.0) then
			if (sum_mean .eq. 0) then
				subbasin_sum(i)=suma_valores/cont_valores
			elseif (sum_mean .eq. 1) then 
				subbasin_sum(i) = suma_valores
			endif
		else
			subbasin_sum(i)=0.0
		endif
	enddo
end subroutine	
subroutine basin_subbasin_stream_prop(sub_pert,cauce,long,slope,stream_slope &
	&,stream_long,nceldas,n_nodos)
	!Variables de entrada
	integer, intent(in) :: nceldas, n_nodos
	integer, intent(in) :: sub_pert(nceldas), cauce(nceldas)
	real, intent(in) :: long(nceldas), slope(nceldas)
	!Variables de salida
	real, intent(out) :: stream_slope(n_nodos), stream_long(n_nodos)
	!f2py intent(in) :: nceldas, n_nodos
	!f2py intent(in) :: sub_pert,cauce, long, slope
	!f2py intent(out) :: stream_slope, stream_long
 	!Variables locales 
	integer i
	real N	
	!itera dersde 1 hasta el maxvar
	do i = 1, n_nodos
		!Para cada una encuentra 
		stream_long(i) = sum(long, mask = (sub_pert .eq. i .and. cauce .eq. 1))
		N = count(sub_pert .eq. i .and. cauce .eq. 1)
		stream_slope(i) = sum(slope, mask = (sub_pert .eq. i &
			&.and. cauce .eq. 1))/N
	enddo	
end subroutine

!-----------------------------------------------------------------------
!Geomorfologia a partir de cuenca
!-----------------------------------------------------------------------
subroutine geo_hand(basin_f,basin_elev,basin_long,cauce,nceldas,&
	&hand_model,hdnd_model,a_quien) !Calcula: HAND: Height above the nearest drainage y HDND: Horizontal distance to the nearest drainage  
    !variables de entrada
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3,nceldas),cauce(nceldas)
    real, intent(in) :: basin_elev(nceldas),basin_long(nceldas)
    !Variables de salida
    real, intent(out) :: hand_model(nceldas),hdnd_model(nceldas)
    integer, intent(out) :: a_quien(nceldas)
    !f2py intent(in) :: nceldas,basin_f,basin_elev,basin_long
    !f2py intent(out) :: hand_model,hdnd_model,a_quien
    !Variables locales
    integer i,drenaid
    real L_sum
    !Todo lo que sea cauce le pone HAND iguala  cero
    where(cauce.eq.1) hand_model=0
    hdnd_model=basin_long
    where(cauce.eq.1) hdnd_model=0
    !Itera por cada celda 
    do i=1,nceldas
	!Prueba si no es cauce
	if (cauce(i).ne.1) then
	    !Itera hasta llegar al cauce
	    flag=1; i_temp=i; L_sum=basin_long(i)
	    do while (flag.eq.1)
		drenaid=nceldas-basin_f(1,i_temp)+1
		if (cauce(drenaid).eq.1) then !llego a un cauce
		    flag=0
		    hand_model(i)=basin_elev(i)-basin_elev(drenaid)
		    hdnd_model(i)=L_sum
		    a_quien(i)=drenaid
		elseif (basin_f(1,drenaid).eq.0) then !dreno a la salida de la cuenca
		    flag=0
		    hand_model(i)=0
		    hdnd_model(i)=0
		    a_quien(i)=0
		else !sigue drenando a celdas tipo ladera
		    i_temp=drenaid
		    L_sum=L_sum+basin_long(drenaid) 
		endif
	    enddo
	endif
    enddo
end subroutine
subroutine geo_hand_global(dem,dir,red,hand,nc,nr) !Genera el mapa de Hand para todo el mapa
	!Variables de entrada
	integer, intent(in) :: nc,nr
	real, intent(in) :: dem(nc,nr)
	integer, intent(in) :: red(nc,nr),dir(nc,nr)
	!variables de salida
	real, intent(out) :: hand(nc,nr)
	!f2py intent(in) :: dem,red,dir,nc,nr
	!f2py intent(out) :: hand
	!Variables locales 
	integer i,j,col,fil,flag,col_move,fil_move
	!Hace que todos los puntos done hay red hand=0
	hand=noData
	where(red.eq.1) hand=0
	!Busca por todas las celdas del mapa
	do i=1,nc
		do j=1,nr
			!Se fija si la entrada no es no data y que no sea parte de la red hidrica
			if (dir(i,j).ne.noData.and.red(i,j).eq.0) then
				!Itera aguas abajo hasta que la celda destino sea cauce
				col=i; fil=j; flag=1
				do while (flag.eq.1 .and. dir(col,fil) .ne. noData)
					!Calcula la columna y la fila de drenaje y se fija que este si sea un no
					call drain_colfil(dir(col,fil),col_move,fil_move)
					col=col+col_move; fil=fil+fil_move
					!Se fija que el punto de drenaje este en el mapa
					if (col.gt.0 .and. col .le. ncols .and. fil .gt. 0 .and. fil .le. nrows) then
						!se fija si el punto de drenaje es red, si lo es calcula hand y para, si no sigue, si es un punto noData para y hand=nodata
						if (red(col,fil).eq.1) then
							hand(i,j)=dem(i,j)-dem(col,fil)
							flag=0
						elseif (red(col,fil).eq.noData) then
							hand(i,j)=noData
							flag=0
						endif
					else
						!Si el punto de drenahe esta por fuera para y hand=nodata
						hand(i,j)=noData
						flag=0
					endif
				end do
			end if
		end do
	end do
end subroutine 
subroutine dir_reclass(Mat_in,Mat_out,nc,nr) !reclasifica los valores de direccion obtenidos a partir de r.watershed de GRASS
	!Variables de entrada
	integer, intent(in) :: nc,nr
	integer, intent(in) :: Mat_in(nc,nr)
	!Variables de salida
	integer, intent(out) :: Mat_out(nc,nr)
	Mat_out=noData
	where(Mat_in.eq.3) Mat_out=7
	where(Mat_in.eq.2) Mat_out=8
	where(Mat_in.eq.1) Mat_out=9
	where(Mat_in.eq.8) Mat_out=6
	where(Mat_in.eq.7) Mat_out=3
	where(Mat_in.eq.6) Mat_out=2
	where(Mat_in.eq.5) Mat_out=1
	where(Mat_in.eq.4) Mat_out=4
end subroutine

!-----------------------------------------------------------------------
!Funciones para hacer sort de algo
!Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
!Based on algorithm from Cormen et al., Introduction to Algorithms,
!-----------------------------------------------------------------------
recursive subroutine QsortC(A)
  real, intent(in out), dimension(:) :: A
  integer :: iq
    !f2py intent(inout) :: A
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine 
subroutine Partition(A, marker) !subrutina utilizada por QsortC para hacer sort
  real, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine

!-----------------------------------------------------------------------
!Herramientas para proceso de mapas (Experimental)
!-----------------------------------------------------------------------
subroutine find_colrow_inArg(mat,umbral1,umbral2,nc,nf,cont)
	!Variables de entrada
	integer, intent(in) :: nc,nf
	real, intent(in) :: mat(:,:),umbral1,umbral2
	!Variables de salida
	integer, intent(out) :: cont
	!f2py intent(in) :: nc,nf,umbral1,umbral2,mat
	!f2py intent(out) :: cont
	!variables locales
	integer i,j
	!Busca
	if (.not. allocated(col_fil_temp)) allocate(col_fil_temp(2,nc*nf))
	cont=0
	do i=1,nc
		do j=1,nf
			if (mat(i,j).gt. umbral1 .and. mat(i,j) .lt. umbral2) then
				cont=cont+1
				col_fil_temp(1,cont)=xll+dx*(i-0.5)
				col_fil_temp(2,cont)=yll+dx*(nrows-j+0.5)
			endif
		enddo
	enddo
end subroutine
subroutine cut_colrow_inArg(cont,col_fil)
	!variables de entrada
	integer, intent(in) :: cont
	!variables de salida 
	real, intent(out) :: col_fil(2,cont)
	!f2py intent(in) :: cont
	!f2py intent(out) :: col_fil
	!Opera	
	if (allocated(col_fil_temp)) then 
		col_fil=col_fil_temp(:,:cont)
		deallocate(col_fil_temp)
	endif
end subroutine

subroutine dem_detect_clouds(image,Grad,kerX,kerY,nc,nf) !detecta las nubes en un dem
	!Variables de entrada
	integer, intent(in) :: nc,nf
	real, intent(in) :: KerX(3,3),KerY(3,3)
	real, intent(in) :: image(nc,nf)
	!variables de salida 
	real, intent(out) :: Grad(nc,nf)
	!f2py intent(in) :: nc,nf,image
	!f2py intent(out) :: Grad
	!Variables locales
	integer i,j
	real DifX(nc,nf),DifY(nc,nf),imageTemp(nc,nf),sumaX,sumaY
	logical mascara(nc,nf)
	!itera por todo el mapa calculando bordes
	imageTemp=image
	where(imageTemp .eq. nodata) imageTemp=0.0
	mascara=.false.
	where(imageTemp .eq. nodata) mascara=.true.
	DifX=0; DifY=0
	do i=2,nc-1
		do j=2,nf-1
			sumaX=0
			sumaY=0
			do ki=1,3
				do kj=1,3
					sumaX=sumaX+imageTemp(i-2+ki,j-2+kj)*KerX(ki,kj)
					sumaY=sumaY+imageTemp(i-2+ki,j-2+kj)*KerY(ki,kj)
				enddo
			enddo
			if (any(mascara(i-1:i+1,j-1:j+1))) then				
				DifX(i,j)=0.0
				DifY(i,j)=0.0
			else	
				DifX(i,j)=sumaX
				DifY(i,j)=sumaY
			endif
		enddo
	enddo
	!Calcula el gradiente de diferencias 
	Grad=sqrt(DifX**2+DifY**2)
end subroutine
subroutine dem_correct_dem_w_dem(dem_in,dem_out,dem_w,mask,nc,nf,nc_m,nf_m,xll_m,yll_m,dx_m)
	!Variables de entrada
	integer, intent(in) :: nc,nf,nc_m,nf_m
	real, intent(in) ::	dem_w(nc_m,nf_m),xll_m,yll_m,dx_m
	real, intent(in) :: dem_in(nc,nf)
	integer, intent(in) :: mask(nc,nf)
	!Variables de salida	
	real, intent(out) :: dem_out(nc,nf)
	!Declaracion f2py
	!f2py intent(in) :: dem_w,nc,nf,xll_m,yll_m,dx_m,mask,dem_in
	!f2py intent(out) :: dem_out
	!Variables locales
	integer i,j,fila,columna
	!Itera sobre el mapa 
	dem_out=dem_in
	do i=1,ncols
		do j=1,nrows
			if (mask(i,j) .eq. 1.0) then
				!Calcula la posicion
				Xpos=xll+dx*(i-0.5)
				Ypos=yll+dx*((nrows-j)+0.5)
				!Evalua si la posicion esta por dentro del mapa
				if (Xpos.gt.xll_m.and.Xpos.lt.(xll_m+dx_m*nc).and.Ypos.gt.yll_m.and.Ypos.lt.(yll_m+nf*dx_m)) then				
				    !si esta por dentro le asigna el valor de la celda equivalente
				    columna=ceiling((Xpos-xll_m)/dx_m)
				    fila=ceiling((Ypos-yll_m)/dx_m)
				    fila=nf_m-fila+1
				    dem_out(i,j)=dem_w(columna,fila)
				endif
			endif
		enddo
	enddo
end subroutine 

subroutine DEM_Slope(DEM,nc,nf,Slope) !Calcula la pendiente
	!Variables de entrada
	integer, intent(in) :: nc,nf
	real, intent(in) :: DEM(nc,nf)
	!Variables de salida
	real, intent(out) :: Slope(nc,nf)
	!f2py intent(in) :: nc,nf,DEM
	!f2py intent(out) :: Slope
	!Variables locales
	integer i,j,ki,kj
	real x,y,s,st
	!Recorre toda la matriz
	Slope=nodata
	do i=2,nc-1
		do j=2,nf-1
			!Recorre todo el kernel de 3x3 para cada celda
			s=0
			do ki=-1,1
				do kj=-1,1
					y=DEM(i,j)-DEM(i+ki,j+kj)
					if (ki .ne. 0 .and. kj .ne. 0) x=1.42*dxp
					st=y/x
					if (st .gt. s) s=st
				enddo
			enddo
			!Toma la mayor pendiente como la pend de la celda
			Slope(i,j)=s
		enddo
	enddo
end subroutine
!subroutine DEM_Pitfill(DEM,nc,nf,DEMfill) !llena huecos en el DEM
!	!Variables de entrada
!	integer, intent(in) :: nc,nf
!	real, intent(in) :: DEM(nc,nf)
!	!Variables de salida
!	real, intent(out) :: DEMfill(nc,nf)
!	!f2py intent(in) :: nc,nf,DEM
!	!f2py intent(out) :: DEMfill
!	!Variables locales
!	integer i,j
!	real ValCel,MinCel,MinCelMask
!	real DEMtemp(nc+2,nf+2)
!	!Recorre toda la matriz
!	DEMtemp(2:nc+1,2:nf+1)=DEM
!	DEMfill=DEM
!	do i=2,nc+1
!		do j=2,nf+1
!			!En cada punto mira si este es el mas bajo o hay iguales 
!			!a el en ese punto
!			ValCel=DEMtemp(i,j)
!			MinCel=minval(DEMtemp(i-1:i+1,j-1:j+1))
!			MinCelMask=minval(DEMtemp(i-1:i+1,j-1:j+1),&
!			&mask=DEMtemp(i-1:i+1,j-1:j+1).ne.MinCel)
!			if (ValCel .eq. MinCel) then
!				DEMfill(i-1,j-1)=MinCelMask
!			endif
!		enddo
!	enddo
	
!end subroutine
!subroutine DEM_carve(DEM,nc,nf,DEMcarve,DIR,Pits) !Erosiona Soille, 2003
!	!Variables de entrada
!	integer, intent(in) :: nc,nf
!	real, intent(in) :: DEM(nc,nf)
!	!Variables de salida
!	real, intent(out) :: DEMcarve(nc,nf)
!	integer, intent(out) :: DIR(nc,nf),Pits(nc,nf)
!	!f2py intent(in) :: nc,nf,DEM
!	!f2py intent(out) :: DEMcarve,DIR,Pits
!	!Variables locales
!	integer i,j,z,ki,kj,cont,cont2,cont3,dirTemp,x,y,x2,y2
!	real ValCel,MinCel,MinCelMask,ValMinN,ValMinV
!	real DEMtemp(nc,nf),PitVal(nc*nf),PitListVal(nc)
!	integer PitX(nc),PitY(nc),PitListX(nc),PitListY(nc)
!	logical flag1
!	!Copia matrices de trabajo
!	DEMtemp=DEM
!	DEMcarve=DEM
!	Pits=0; PitX=-999; PitY=-999
!	!Paso 1 Identifica PITS
!	!Recorre toda la matriz buscando pits	
!	cont=0
!	do i=2,nc-1
!		do j=2,nf-1
!			!En cada punto mira si este es el mas bajo o hay iguales 
!			!a el en ese punto
!			ValCel=DEMtemp(i,j)
!			MinCel=minval(DEMtemp(i-1:i+1,j-1:j+1))
!			if (ValCel .eq. MinCel) then
!				cont=cont+1
!				Pits(i,j)=1
!				PitX(cont)=i
!				PitY(cont)=j
!				PitVal(cont)=MinCel
!			endif
!		enddo
!	enddo
!!	!Paso 2 en cada pit hace el procedimiento de soille et al 2003
!	PitListVal=0
!	do i=1,cont
!		!Copia para iniciar la lista de busqueda
!		PitListX(1)=PitX(cont)
!		PitListY(1)=PitY(cont)
!		PitListVal(1)=PitVal(cont)
!		flag1=.true.
!		cont2=1		
!		ValMinV=PitVal(cont)
!		!Itera hasta encontrar una salida
!		do z=1,10							
!			do j=1,cont2
!				x=PitListX(j); y=PitListY(j)
!				!Encuentra el siguiente valor minimo del pit y 
!				!cuantos hay con ese valor
!				ValMinN=minval(DEMtemp(x-1:x+1,y-1:y+1),&
!					&mask=DEMtemp(x-1:x+1,y-1:y+1).gt.ValMinV)								
!				!Encuentra exactamente sus posiciones y busca
!				! a ver si alguna de las nuevas drena fuera del pit	
!				cont3=cont2
!				do ki=-1,1
!					do kj=-1,1
!						if (DEMtemp(x+ki,y+kj).eq.ValMinN) then
!							cont3=cont3+1
!							PitListX(cont3)=x+ki
!							PitListY(cont3)=y+kj
!							PitListVal(cont3)=DEMtemp(x+ki,y+kj)
!							!Busca si le drena a alguien
!							x2=x+ki; y2=y+kj
!							call kernel_DIR(DEMtemp(x2-1:x2+1,y2-1:y2+1),&
!								&dirTemp)
!							!Si drena por fuera termina la inundacion
!							if (dirTemp.ne.0) then
!								flag1=.false.
!								Pits(x2,y2)=2
!							endif									
!						endif
!					enddo
!				enddo
!				DEMtemp(x,y)=ValMinN			
!				!finaliza la busqueda interna
!			enddo			
!			ValMinV=ValMinN
!			!imprimer para revisar
!			print *, '-------------------------'
!			print *, ValMin
!			print *, PitListVal(1:cont3)
!			print *, flag1
!			print *, cont2
!			print *, cont3
!			cont2=cont3
!		enddo
!		!Otro pit		
!	enddo
!end subroutine
!subroutine kernel_DIR(DEMkernel,DIR)
!	!variables de entrada
!	real, intent(in) :: DEMkernel(3,3)
!	!variables de salida
!	integer, intent(out) :: DIR
!	!Variables locales 
!	integer i,j, directions(9),cont
!	real y,x,s,st
!	!busca una direccion
!	s=0
!	directions=(/7,4,1,8,0,2,9,6,3/)
!	cont=1
!	DIR=0
!	do i=1,3
!		do j=1,3
!			y=DEMkernel(2,2)-DEMkernel(i,j)
!			if (cont.eq.1 .or. cont.eq.3 .or. cont.eq.7 .or. cont .eq. 9) then
!				x=1.41*x
!			endif
!			st=y/x
!			cont=cont+1
!			if (st .gt. s) then
!				s=st
!				DIR=directions(cont)
!			endif
!		enddo
!	enddo
!end subroutine


!subroutine DEM_Process(DEM,nc,nf,DIR,Slope,DEMProc) !Subrutina 1 de procesamiento
!	!Variables de entrada
!	integer, intent(in) :: nc,nf
!	real, intent(in) :: DEM(nc,nf)
!	!Variables de salida
!	real, intent(out) :: DEMProc(nc,nf),Slope(nc,nf)
!	integer, intent(out) :: DIR(nc,nf)
!	!f2py intent(in) :: nc,nf,DEM
!	!f2py intent(out) :: DEMProc,DIR,Slope
!	!Variables locales
!	integer i,j,ki,kj
!	real x,y,s,st
!	integer directions(9),cont
!	!Recorre toda la matriz
!	directions=(/7,4,1,8,0,2,9,6,3/)
!	DIR=nodata
!	do i=2,nc-1
!		do j=2,nf-1
!			!Recorre todo el kernel de 3x3 para cada celda
!			s=0; cont=0
!			do ki=-1,1
!				do kj=-1,1
!					y=DEM(i,j)-DEM(i+ki,j+kj)
!					x=dxp
!					cont=cont+1	
!					if (ki .ne. 0 .and. kj .ne. 0) x=1.41*dxp
!					st=y/x
!					if (st .gt. s) then						
!						s=st
!						DIR(i,j)=directions(cont)						
!					endif									
!				enddo
!			enddo			
!			Slope(i,j)=s			
!		enddo
!	enddo

!end subroutine


!--------------------------------------------------------------
!Esneider


end module

!-----------------------------------------------------------------------
!Subrutinas extras del moduo que son usadas internamente por cuencas
!-----------------------------------------------------------------------
subroutine find_xy_in_basin(basin_f,col,fil,posit,nceldas) !Encuentra la posicion de un par col,fil en un vector tipo cuenca 
    !Variables de entrada
    integer, intent(in) :: nceldas,col,fil
    integer, intent(in) :: basin_f(3,nceldas)
    integer, intent(out) :: posit
    !Variables locales
    integer i,flag
    !Busca la posicion
    posit=0; flag=1
    i=1
    do while(i.le.nceldas.and.flag.eq.1)
	if (basin_f(2,i).eq.col.and.basin_f(3,i).eq.fil) then
	    flag=0 !Encontro el punto de salida
	    posit=i
	else
	    i=i+1
	endif
    enddo
    !Si posit a la salida ==0 es porque el punto esta fuera de la cuenca
end subroutine
subroutine drain_colfil(dir,col_obj,fil_obj) !Encuentra lo que hay que sumar o restar a fil col para llegar a fil col donde drena
    !Variables de entrada
    integer, intent(in) :: dir
    integer, intent(out) :: col_obj,fil_obj
    !encuentra
    if (dir.eq.7) then
	col_obj=-1; fil_obj=-1
    elseif (dir.eq.8) then
	col_obj=0; fil_obj=-1
    elseif (dir.eq.9) then
	col_obj=1; fil_obj=-1
    elseif (dir.eq.4) then
	col_obj=-1; fil_obj=0
    elseif (dir.eq.6) then
	col_obj=1; fil_obj=0
    elseif (dir.eq.1) then
	col_obj=-1; fil_obj=1
    elseif (dir.eq.2) then
	col_obj=0; fil_obj=1
    elseif (dir.eq.3) then
	col_obj=1; fil_obj=1
	else
	col_obj=noData; fil_obj=noData
    endif
end subroutine

