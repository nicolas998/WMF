!Modelos: acoplado con cuencas presenta una serie de modelos hidrologicos distribuidos
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

module modelos

!La lluvia es controlada por el modulo lluvia.
!use lluvia

!Todas las variables tienen que ser definidas dentro del modulo
implicit none

!-----------------------------------------------------------------------
!Variableses globales para la modelacion
!-----------------------------------------------------------------------
!Variables de la cuenca 
real xll,yll !coordenadas de la esquina inferior izquierda
real noData !Valor que representa los no datos
real dx !Largo del mapa leido
real dxP !largo proyectado del mapa, debe ser indicado si no se cononce
integer ncols,nrows !cantidad de columnas y filas del mapa
integer nceldas !Cantidad de celdas que componen la cuenca

!Variables de rutas para leer y escribir informacion
character*500 rute_rain !Lectura de archivos de lluvia, que tengan forma de la cuenca
character*500 rute_speed !Escritura de mapas de velocidad [N_reg,4,Nceldas]
character*500 rute_storage !Escritura de mapas de almacenamiento [N_reg,5,Nceldas]

!Variables de propiedades geomorfologicas
integer, allocatable :: drena(:,:) !Indicador de la topologia de la cuenca 
integer, allocatable :: unit_type(:,:) !Tipo de celda: 1. ladera, 2. Carcava, 3. Cauce fijo
real, allocatable :: hill_long(:,:), slope_hill(:,:) !Long celda, pendiente [1,Nceldas]
real, allocatable :: stream_long(:,:), slope_stream(:,:) !Longitud del cauce pendiente [1,Nceldas]
real, allocatable :: stream_width(:,:) !Ancho de los rios.
real, allocatable :: elem_area(:,:) !Area de cada elemento de la cuenca, si es celdas es dx**2, si es laderas es otra cosa


!Variables de propiedades fisicas
real, allocatable :: v_coef(:,:) !Coeficientes de velocidades verticales [LT-1] [4,Nceldas]: 1. Evp, 2. Ks, 3. Kp, 4. kpp
real, allocatable :: v_exp(:,:)  !Exponentes de velocidades verticales no lineales [adim] [4,Nceldas]. (no implementado)
real, allocatable :: h_coef(:,:) !Coeficientes de velocidades horizontales [L] [4,Nceldas]: 1. runoff, 2. sub-sup, 3. subte, 4. cauce
real, allocatable :: h_exp(:,:)  !Exponentes de velocidades horizontales, aplica en casos no lineales
real, allocatable :: Max_capilar(:,:) !Maximo almacenamiento capilar [L] [1,Nceldas]
real, allocatable :: Max_gravita(:,:) !Maximo almacenamiento gravitacional [L] [1,Nceldas]

!Variables de control y configuracion del modelo 
real dt !Delta del tiempo de modelacion
integer rain_first_point !Primer punto de lectura de lluvia
integer sim_sediments !Simula (1) o no (0) sedimentos
integer sim_slides !Simula (1) o no (0) deslizamientos
integer save_storage !Guarda (1) o no (0) almacenamiento de los tanques en cada intervalo.
integer save_speed !Guarda (1) o no (0) las velocidades en cada intervalo
integer speed_type(3) !Tipo de velocidad para tanque 1, 2 y 3: 1: lineal, 2: Cinematica Potencial, de momento no hay mas implementadas 
integer, allocatable :: control(:,:) !Celdas de la cuenca que tienen puntos de control
integer, allocatable :: control_h(:,:) !Celdas de la cuenca que son puntos de control de humedad
integer, allocatable :: guarda_cond(:,:) !Intervalos de tiempo en que se hace guardado de condiciones

!Variables de resultados globales (siempre van a estar ahi)
real, allocatable :: Storage(:,:) !Almacenamiento de los 5 tanques del modelo
real, allocatable :: Speed_map(:,:) !Mapa de velocidades en los tanques, solo se activa si se da la opcion save_speed
real, allocatable :: Mean_Rain(:,:) !Serie de lluvia promedio en la cuenca [mm]

!variables de sedimentos Par aalojar volumens y demas
real sed_factor !factor para calibrar el modelo de sedimentos
real wi(3), Qskr,G,diametro(3)
real ERO(3),EROt(3),DEP(3),DEPt(3)
real, allocatable :: VolERO(:),VolDEPo(:) !Volumen erosionado por celda, Volumen depositado en la celda
real, allocatable :: Vs(:,:),Vd(:,:)
real, allocatable :: VSc(:,:),Vdc(:,:)
!Mapas para el calculo de sedimentos
real, allocatable :: Krus(:,:),Crus(:,:),Prus(:,:) !Factores de la RUSLE
real, allocatable :: PArLiAc(:,:) !Porcentaje de arenas, limos y arcillas

!Variables del modelo de deslizamientos
!Factores de seguridad, mapas de riesgos y de deslizamientos
real FS !Factor de seguridad, por defecto se deja en 1.
real, allocatable :: RiskVector(:,:) !Vector de riesgos: 1. No riesgo, 2. Evaluable, 3. Siempre en riesgo.
real, allocatable :: SlideOcurrence(:,:) !Vector con las ocurrencias de deslizamientos de acuerdo al modelo.
real, allocatable :: Zcrit(:,:), Zmin(:,:), Zmax(:,:), Bo(:,:) !Profundiad critica, minima, maxima [mts] y angulo critico [rad]
!Propiedades fisicas del suelo 
real, allocatable :: Zs(:,:) !Profundidad del suelo [mm]
real GammaW !Densidad del agua, es un solo valor
real, allocatable :: GammaS(:,:) !Densidad del suelo 
real, allocatable :: Cohesion(:,:) !Cohesion del suelo [Kpa]
real, allocatable :: FrictionAngle(:,:) !Angulo de friccion [rad]

contains


!-----------------------------------------------------------------------
!Modelo
!-----------------------------------------------------------------------

subroutine shia_v1(calib,N_cel,N_cont,N_contH,N_reg,Q,&
	& Qsed, Hum, balance, StoOut)
    
	!Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_cont,N_contH
    real, intent(in) :: calib(10)
    
	!Variables de salia
    real, intent(out) :: Hum(N_contH,N_reg),Q(N_cont,N_reg),Qsed(3,N_cont,N_reg) !Control humedad en el suelo, Control caudales 
    real, intent(out) :: StoOut(5,N_cel),balance !Almacenamiento en tanques, balance total de la modelacion
     
	!Variables de la lluvia
	real Rain(N_cel) !Lluvia leida en el intervalo de tiempo [mm] [1,N_cel]
	integer Res !Dice si se leyo bien o no la lluvia 
	real rain_sum
	!Variables de iteracion
	integer celda,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    integer control_cont, controlh_cont, i
    !Variables para el balance
    real entradas, salidas, inicial
	!Variables de conversion
	real m3_mmHill(N_cel) !Convierte mm a m3 o viceversa se basa en el area del elemento.
	real m3_mmRivers(N_cel)
	!Variables de almacenamiento en tanques, flujo vertical y flujo horizontal
	real vflux(4) !Flujo vertical [mm]
	real hflux(4) !Flujo horizontal [mm]
	real Evp_loss !Salida por evaporcion del sistema [mm]
	!Variables de velocidad vertical y horizontal
	real vspeed(4,N_cel) !Velocidad vertical [cm/h]
	real hspeed(4,N_cel) !Velocidad horizontal [cm/h] o [m/s]
	real section_area !Area de la seccion resuleta en ladera o en el canal 
	!Variables Max storage en tanques 1 y 3
	real H(2,N_cel)
	!Variables sub-modelo de sedimentos
	real Area_coef(nceldas) !Coeficiente para el calculo del lateral en cada celda del tanque 2 para calcuo de sedimentos
    real Vsal_sed(3) !Volumen de salida de cada fraccion de sedimentos [m3/seg]
	
	!Inicia la variable global de lluvia promedio sobre la cuenca
	if (allocated(Mean_Rain)) deallocate(Mean_Rain)
	allocate(Mean_Rain(1,N_reg))
	!Establece variable de conversion
	m3_mmHill=elem_area(1,:)/1000.0
	m3_mmRivers=(stream_long(1,:)*stream_width(1,:))/1000.0
	Q=0.0
	control_cont=1
	controlh_cont=1
	!Inicia variables para realizar el balance en la cuenca
	StoOut=Storage
	inicial=sum(Storage)
	entradas=0
	salidas=0
	!Calcula parametros que cambien con el tiempo
	! Calib 1 2 3 y 4: Velocidad vertical, 1. vel evp, los demas vel prof.
	! Calib 5 6 7 y 8: Velocidad Hztal tanques 2 3 4 y 5
	! Calib 9 y 10: Alm maximo capilar y gravitacional
	do i=1,4
		vspeed(i,:)=v_coef(i,:)*Calib(i)*dt ! Cantidad [mm] que baja de tanque a tanque
		hspeed(i,:)=h_coef(i,:)*Calib(i+4) ! Velocidad [mm/seg] que se mueve
	enddo
	!Calcula parametros estaticos en el tiempo
	H(1,:)=Max_capilar(1,:)*Calib(9)
	H(2,:)=Max_gravita(1,:)*Calib(10)
	!Si guarda velocidad aloja la variable 
	if (save_speed.eq.1) then 
		if (allocated(Speed_map)) deallocate(Speed_map)
		allocate(Speed_map(4,N_cel))
	endif
	!Si va a ejecutar sedimentos, aloja variables y arregla coeficientes
	if (sim_sediments.eq.1) then 
		call allocate_sed
		Qsed=0  
	endif
	!Si va a ejecutar deslizamientos, aloja variables y arregla coeficientes 
	if (sim_slides .eq. 1) then 
		call slide_allocate
	endif
	
	!Iter in the time
	do tiempo=1,N_reg
		
		!Lee la lluvia 
		call read_float_basin(rute_rain, rain_first_point+tiempo-1,&
			& N_cel, Rain, Res) 
		
		!Iter around the cells or hills
		do celda=1,N_cel
			
			!determina el elemento objetivo y realiza balance de lluvia
			drenaid=N_cel-drena(1,celda)+1
			entradas=entradas+Rain(celda)
			rain_sum=rain_sum+Rain(celda)
			
			!Determina el flujo vertical entre tanques
			vflux(1) = max(0.0, Rain(celda)-H(1,celda)+StoOut(1,celda))
			StoOut(1,celda)=StoOut(1,celda)+Rain(celda)-vflux(1)
			Evp_loss=min(vspeed(1,celda)*(StoOut(1,celda)/H(1,celda))**0.6,&
				&StoOut(1,celda))
			StoOut(1,celda)=StoOut(1,celda)-Evp_loss
			salidas=salidas+Evp_loss
			do i=1,3
				vflux(i+1)=min(vflux(i),vspeed(i+1,celda))
				StoOut(i+1,celda)=StoOut(i+1,celda)+vflux(i)-vflux(i+1)
			enddo
			StoOut(2,celda)=StoOut(2,celda)+max(0.0, StoOut(3,celda)&
				&+vflux(3)-vflux(4)-H(2,celda))
			salidas=salidas+vflux(4)
			
			!Calcula el flujo que sale de los tanques 2 a 4
			do i=1,3
				!calcula la velocidad de transferencia
				select case(speed_type(i))
					!Caso por defecto es lineal 
					case(1)
						hflux(i)=(1-hill_long(1,celda)/(hspeed(i,celda)*dt+&
							& hill_long(1,celda)))*StoOut(i+1,celda)
					!Caso no lineal potencial 
					case(2)	
						call calc_speed(StoOut(i+1,celda)*m3_mmHill(celda), h_coef(i,celda),&
							& h_exp(i,celda), hill_long(1,celda), hspeed(i,celda), section_area)
						hflux(i)=section_area*hspeed(i,celda)*dt/m3_mmHill(celda) ![mm]
				end select
				!Actualiza el almacenamiento
				StoOut(i+1,celda)=StoOut(i+1,celda)-hflux(i)	
			enddo	
				
			!Envia los flujos que salieron de acuerdo al tipo de celda 
			if (unit_type(1,celda).eq.1) then
				if (drena(1,celda).ne.0) then
					StoOut(2:4,drenaid)=StoOut(2:4,drenaid)+hflux(1:3)
				else
					Q(1,tiempo)=Q(1,tiempo)+sum(hflux(1:3))*m3_mmHill(celda) ![m3/s]
					salidas=salidas+sum(hflux(1:3))
				endif
			elseif (unit_type(1,celda).gt.1) then
				!Envia los flujos de acuerdo al tipo de celda
				StoOut(4,drenaid)=StoOut(4,drenaid)+hflux(3)*&
					&(3-unit_type(1,celda))
				StoOut(5,celda)=StoOut(5,celda)+sum(hflux(1:2))+&
					& hflux(3)*(unit_type(1,celda)-2)
				!Resuelve el transporte en el canal 
				call calc_speed(StoOut(5,celda)*m3_mmRivers(celda), h_coef(4,celda),&
					& h_exp(4,celda), stream_long(1,celda), hspeed(4,celda), section_area)
				hflux(4)=min(section_area*hspeed(4,celda)*dt*Calib(8)/m3_mmRivers(celda),StoOut(5,celda)) ![mm]
				!Envia el flujo en el cauce aguas abajo
				if (drena(1,celda).ne.0) then
					StoOut(5,celda)=StoOut(5,celda)-hflux(4)
				else
					Q(1,tiempo)=hflux(4)*m3_mmRivers(celda)/dt ![m3/s]      
					salidas=salidas+hflux(4) ![mm]
				endif
			endif
			
			!Si evalua tte de sedimentos calcula el tte en ladera y cauce
			if (sim_sediments .eq. 1) then
				Area_coef=m3_mmHill(celda)/(hill_long(1,celda)+hspeed(1,celda)*dt)
				section_area=StoOut(2,celda)*Area_coef(celda)
				!Calculo en ladera
				call sed_hillslope(sed_factor, StoOut(2,celda), hspeed(1,celda)&
					&, slope_hill(1,celda), section_area, celda, drenaid, unit_type(1,celda))
				!Calculo en cauce
				call sed_channel(StoOut(5,celda),hspeed(4,celda),hflux(4)*m3_mmRivers(celda)/dt,&
					& slope_stream(1,celda), section_area,celda,drenaid,Vsal_sed)
				!Si es la salida registra los sedimentos en la salida 
				if (drena(1,celda).eq.0) then
					do i=1,3
						Qsed(i,1,tiempo)=Vsal_sed(i)
					enddo  
				endif
			endif
			
			!Si evalua deslizamientos
			if (sim_slides.eq.1) call slide_ocurrence(celda, StoOut(3,celda)&
				&, H(2,celda))
			
			!Record de variables y resultados del modelo
			!Caudales en el punto de control
			if (control(1,celda).ne.0) then
				Q(control_cont,tiempo)=hflux(4)*m3_mmRivers(celda)/dt ![m3/s]      
				!Si se simularon sedimentos los guarda 
				if (sim_sediments.eq.1) then
					do i=1,3
						Qsed(i,control_cont,tiempo)=Vsal_sed(i)
				    enddo
				endif
				control_cont=control_cont+1
			endif
			!Humedad en puntos de control
			if (control_h(1,celda).ne.0) then 
				Hum(controlh_cont,tiempo)=sum((/ StoOut(1,celda), StoOut(3,celda)/))
				controlh_cont=controlh_cont+1
			endif
				
		enddo
		
		!Obtiene la lluvia promedio para el intervalo de tiempo
		Mean_Rain(1,tiempo)=rain_sum/N_cel
		
		!Mapas de velocidad del flujo, muestra la velocidad promedio 
		if (save_storage .eq. 1) then
			call write_float_basin(rute_storage,StoOut,tiempo,N_cel,5)
		endif
		if (save_speed .eq. 1) then
			call write_float_basin(rute_storage,hspeed,tiempo,N_cel,4)
		endif
		
	enddo
	
	!Realiza el balance, para ver si no se perdio mucha agua 
	balance=inicial+entradas-salidas-sum(StoOut)
	
end subroutine


!-----------------------------------------------------------------------
!Subrutinas Lectura y escritura de mapas
!-----------------------------------------------------------------------
!Lee los datos flotantes de un binario de cuenca en los records ordenados
subroutine read_float_basin(ruta, record, N_cel, vect, Res) 
    !Variables de entrada
    integer, intent(in) :: record, N_cel
    character*500, intent(in) :: ruta
    !Variables de salida
    real, intent(out) :: vect(N_cel)
    integer, intent(out) :: Res
    !f2py intent(in) :: record, N_cel, ruta
    !f2py intent(out) :: vect, Res    
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',&
		& RECL=4*N_cel)
	    read(10,rec=record,iostat=Res) vect
	    if (Res.ne.0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
	close(10)
end subroutine
!Lee los datos flotantes de un binario de cuenca en los records ordenados
subroutine write_float_basin(ruta,vect,record,N_cel,N_col) 
    !Variables de entrada
    integer, intent(in) :: record, N_cel, N_col
    character*255, intent(in) :: ruta
    real, intent(in) :: vect(N_col,N_cel)
    !Variables internas
    character*10 estado
    !Escritura     
    estado='old'
    if (record.gt.1) estado='replace'
    open(10,file=ruta,form='unformatted',status=estado,access='direct',RECL=4*N_col*N_cel)
	write(10,rec=record) vect
    close(10)
end subroutine

!-----------------------------------------------------------------------
!Subrutinas de solucion
!-----------------------------------------------------------------------
!Solucion de la onda cinematica generica, aplica para ecuaciones 
!del tipo v=c*Area**Exp.
subroutine calc_speed(sm, coef, expo, elem_long, speed, area)
	!Variables de entrada
	real, intent(in) :: sm,coef, expo, elem_long
	!Variables de salidqa
	real, intent(out) :: Area
	real, intent(inout) :: speed
	!Variables locales
	real new_speed
	integer i 
	!Itera la cantidad de veces niter para solucionar la ecuacion
	do i=1,4
	    Area=sm/(elem_long+speed*dt) ![m2] Calcula el area de la seccion
	    new_speed=coef*(Area**expo) ![m/seg] Calcula la velocidad nueva
	    speed=(2*new_speed+speed)/3 ![m/seg] Promedia la velocidad
	enddo		
end subroutine 

!-----------------------------------------------------------------------
!Subrutinas de sedimentos
!-----------------------------------------------------------------------
subroutine sed_allocate !Funcion para alojar variables si se van a calcular sed
    !Tamaño a los vectores de sedimentos en suspención y depositads
    if (allocated(VS) .eqv. .false.) allocate(VS(3,nceldas))
    if (allocated(VD) .eqv. .false.) allocate(VD(3,nceldas))
    if (allocated(VolERO) .eqv. .false.) allocate(VolERO(nceldas))
    if (allocated(VolDEPo) .eqv. .false.) allocate(VolDEPo(nceldas))
    if (allocated(VSc) .eqv. .false.) allocate(VSc(3,nceldas))
    if (allocated(VDc) .eqv. .false.) allocate(VDc(3,nceldas))    
    !estado inicial del almacenamiento de sedimentos
    VS=0 ![m3]
    VD=0![m3]
    VSc=0 ![m3]
    VDc=0 ![m3]
    VolERO=0; VolDEPo=0 ![m3]
    EROt=0; DEPt=0
end subroutine
subroutine sed_hillslope(alfa,S2,v2,So,area_sec,celda,drena_id,tipo) !Subrutina para calcular los sedimentos en ladera
    !Variables de entrada
    real, intent(in) :: S2,v2,So,area_sec,alfa
    integer, intent(in) :: celda,drena_id,tipo
    !Variables locales 
    integer i
    real qsSUStot,qsBMtot,qsEROtot,qsSUS(3),qsBM(3),qsERO(3) !deposit y transporte
    real SUStot,DEPtot
    real totXSScap, REScap !Capacidades de transporte
    real q,Qskr,cap,Adv,Te(3) !Caudal lateral y Capacidad de arrastre
    !Inicia Variables propias    
    DEP=0
    qsSUStot=0;qsBMtot=0;qsEROtot=0
    qsSUS=0;qsBM=0;qsERO=0                
    SUStot=0;DEPtot=0    
    !Calcula capacidad de arrastre
    q=area_sec*v2/dxp
    Qskr=58390*alfa*(So**1.664)*(q**2.035)*Krus(1,celda)*Crus(1,celda)*Prus(1,celda)*dxp*dt
    !Calcula Depositacion y atrapamiento      
    do i=1,3
	if (S2/1000>wi(i)*dt) then
	    Te(i)=wi(i)*dt*1000/S2
	else
	    Te(i)=1
	endif
	!Calcula los sedimentos depositados en la celda
	Vd(i,celda)=Vd(i,celda)+Te(i)*Vs(i,celda)
	Vs(i,celda)=Vs(i,celda)-Te(i)*Vs(i,celda) 
	DEP(i)=Te(i)*Vs(i,celda)
	!Calcula totales de suspendidos y depositados
	SUStot=SUStot+Vs(i,celda)
	DEPtot=DEPtot+Vd(i,celda)	
    enddo    
    !Transporta los sedimentos suspendidos
    do i=1,3	
	if (Vs(i,celda).gt.0.0) then
	    if (Qskr .lt. SUStot) then
		cap=Qskr*Vs(i,celda)/SUStot  ![m3]                         
		!Volumen que se puede llevar por advección
		Adv=Vs(i,celda)*v2*dt/(dxp+v2*dt) ![m3]    
		!Transporta por suspención lo mayor entre: lo que se puede llevar la capacidad
		!y lo que se puede llevar por advección
		qsSUS(i)=max(adv,cap) ![m3]
	    else
		qsSUS(i)=Vs(i,celda)
	    endif
	endif
	qsSUStot=qsSUStot+qsSUS(i)
	Vs(i,celda)=Vs(i,celda)-qsSUS(i)
	if (tipo.eq.1) then
	    Vs(i,drena_id)=Vs(i,drena_id)+qsSUS(i)
	else
	    Vsc(i,celda)=Vsc(i,celda)+qsSUS(i)
	endif
    enddo	
    !Capacidad de transporte de exceso
    totXSScap=max(0.0,Qskr-qsSUStot)
    !Transporta los sedimentos depositados	
    if (totXSScap>0.0 .and. DEPtot>0.0) then
	do i=1,3			
	    if (Vd(i,celda)>0.0) then    
		!Observa si la cap excedente es mayor a la cant total depositada
		if (totXSScap<DEPtot) then
		    !tta porcentaje del volumen depositado de la fracción
		    qsBM(i)=totXSScap*Vd(i,celda)/DEPtot ![m3]
		else
		    !tta todo el volumen depositado de la fracción
		    qsBM(i)=Vd(i,celda) ![m3]
		endif
	    endif
	    !Acumula lo que se tta de cada fracción en el tot de transportado por cama
	    qsBMtot=qsBMtot+qsBM(i) ![m3]
	    !Actualiza lo que se fue del vol depositado de cada fracción
	    DEP(i)=DEP(i)-qsBM(i) ![m3]
	    if (DEP(i).lt.0) DEP(i)=0	    
	    Vd(i,celda)=Vd(i,celda)-qsBM(i)
	    if (tipo.eq.1) then
		Vs(i,drena_id)=Vs(i,drena_id)+qsBM(i)
	    else
		Vsc(i,celda)=Vsc(i,celda)+qsBM(i)
	    endif
	enddo
    endif
    DEPt=DEPt+DEP
    !Capacidad residual 
    REScap=max(0.0, totXSScap-qsBMtot) ![m3]
    !Si la capacidad residual es mayor a cero, se presenta erosión
    if (REScap>0.0) then
	!Evalua para cada fracción de tamaño
	do i=1,3
	    !Eroda y agrega a esa fracción de acuerdo a la porción presente en el suelo
	    qsERO(i)=PArLiAc(i,celda)*REScap/100.0 ![m3]
	    qsEROtot=qsEROtot+qsERO(i) ![m3]
	    if (tipo.eq.1) then
		Vs(i,drena_id)=Vs(i,drena_id)+qsERO(i)
	    else
		Vsc(i,celda)=Vsc(i,celda)+qsERO(i)
	    endif
	end do
    endif
    EROt=EROt+qsERO    
    !Actualiza el mapa total de erosion y depositacion
    VolERO(celda)=VolERO(celda)+sum(qsERO)
    VolDEPo(celda)=VolDEPo(celda)+sum(DEP)
end subroutine
subroutine sed_channel(S5,v5,Q5,So,area_sec,celda,drena_id,VolSal) !subrutina para calcular seduimentos en el canals
    !Variables de entrada
    real, intent(in) :: S5,v5,So,area_sec,Q5
    integer, intent(in) :: celda,drena_id
    !Variables de salida
    real, intent(out) :: VolSal(3)
    !Variables locales 
    integer i
    real, parameter :: grav=9.8
    real, parameter :: Gsed=2.65
    real Cw(3),Te(3) !Coef de Engelgud y Hansen 
    real qsSUS(3),qsBM(3) !deposit y transporte
    real XSScap,VolEH !Capacidades de transporte
    real AdvF,supply,L,Rh !Caudal lateral y Capacidad de arrastre
    !Inicia Variables propias
    DEP=0
    qsSUS=0;qsBM=0                
    !Longitud de la seccion y Radio Hidraulico
    if (S5.gt.0.0) then 
	L=area_sec*1000/S5
	Rh=L*S5/(1000*L+2*S5)
    else
	Rh=0.0
    endif
    !Calcula Depositacion y atrapamiento      
    do i=1,3
	!Concentracion de sedimentos por Engelund y Hansen        
	Cw(i)=0.05*(Gsed/(Gsed-1))*(v5*So)/(sqrt((Gsed-1)*grav*diametro(i)/1000.0))*sqrt((Rh*So)/((Gsed-1)*(diametro(i)/1000.0)))
	!Tasa de atrapamiento
	if (S5/1000>wi(i)*dt) then
	    Te(i)=wi(i)*dt*1000/S5
	else
	    Te(i)=1
	endif
	!Calcula los sedimentos depositados en la celda
	VDc(i,celda)=VDc(i,celda)+Te(i)*Vsc(i,celda)
	Vsc(i,celda)=Vsc(i,celda)-Te(i)*Vsc(i,celda)
	DEP(i)=DEP(i)+Te(i)*Vsc(i,celda)
    enddo	    
    !Mueve los sedimentos suspendidos
    do i=1,3
    !Calcula lo que hay de seimdneots de la fracción
	supply=Vsc(i,celda)+Vdc(i,celda) ![m3] 
	!inicializa el volumen suspendido de la fracción en cero
	qsSUS(i)=0.0   
	!Si lo que hay es mayor que cero, lo transporta
	if (supply.gt.0.0) then
	    !Calcula el vol que puede ser transportado por EH
	    VolEH=Q5*Cw(i)*dt/2.65 ![m3]   !Ojo: Creemos que falta multiplicar por dt
	    !Calcula factor de tte por advección
	    AdvF=min(1.0,v5*dt/(dxp+v5*dt)) ![adim]
	    !Realiza el transporte por advección
	    qsSUS(i)=Vsc(i,celda)*AdvF ![m3]
	    !Calcula la capacidad de Exceso sobrante
	    XSScap=max(0.0,VolEH-qsSUS(i))
	    !Calcula lo que se va de lo depositado
	    qsBM(i)=Vdc(i,celda)*Advf ![m3]
	    qsBM(i)=min(XSScap,qsBM(i)) ![m3] 	    
	    DEP(i)=DEP(i)-qsBM(i)
	    if (DEP(i).lt.0.0) DEP(i)=0.0           
	    !Actulia los almacenamientos 
	    Vsc(i,celda)=Vsc(i,celda)-qsSUS(i)
	    Vdc(i,celda)=Vdc(i,celda)-qsBM(i)
	    !Envia aguas abajo lo que se ha transportado
	    if (drena_id.ne.0.0) Vsc(i,drena_id)=Vsc(i,drena_id)+qsSUS(i)+qsBM(i)
	    !Reporta el volumen saliente en el intervalo de tiempo
	    VolSal(i)=(qsSUS(i)+qsBM(i))/dt
	endif
    enddo
    DEPt=DEPt+DEP
    VolDEPo(celda)=VolDEPo(celda)+sum(DEP)
end subroutine

!-----------------------------------------------------------------------
!Subrutinas de deslizamientos
!-----------------------------------------------------------------------
subroutine slide_allocate !Funcion para alojar variables de deslizamientos
	!Aloja variables de umbrales para deslizamientos
	if (allocated(Zmin) .eqv. .false.) allocate(Zmin(1,nceldas))
	if (allocated(Zmax) .eqv. .false.) allocate(Zmax(1,nceldas))
	if (allocated(Zcrit) .eqv. .false.) allocate(Zcrit(1,nceldas))
	if (allocated(Bo) .eqv. .false.) allocate(Bo(1,nceldas))
	!Aloja variables de resultados
	if (allocated(SlideOcurrence) .eqv. .false.) allocate(SlideOcurrence(1,nceldas))
	if (allocated(RiskVector) .eqv. .false.) allocate(RiskVector(1,nceldas))
	!Calcula variables de acuerdo a las propiedades fisicas del suelo 
	!Profundidad critica de inmunidad
	Zmin = Cohesion/((GammaW*(COS(slope_hill))**2*TAN(FrictionAngle))&
		&+(GammaS*(COS(slope_hill))**2*(TAN(slope_hill)-TAN(FrictionAngle))))	
	!Calculating the minimum value of landslide-triggering saturated depth
	Zcrit = (GammaS/GammaW)*Zs*(1.0-(TAN(slope_hill)/TAN(FrictionAngle)))&
		&+(Cohesion/(GammaW*(COS(slope_hill))**2*TAN(FrictionAngle)))
	!Calculating soil thickness for unstable conditions
	Zmax = Cohesion/((GammaS*(COS(slope_hill))**2)*(TAN(slope_hill)-TAN(FrictionAngle)))
	!Calculating slope angle for unconditional stable conditions
	Bo = ATAN(-TAN(FrictionAngle*(GammaW-GammaS)/GammaS))
	!Calcula el mapa de suceptibilidad 
	RiskVector=0 !Se asume todo el vector estable
	where(slope_hill.gt.Bo .and. Zs .gt. Zmin) RiskVector=1 ! Condicionado
	where(slope_hill.gt.Bo .and. Zs .gt. Zmax) RiskVector=2 ! Inestable
	!Inicia en cero el vector de deslizamientos, como si no ocurrieran
	SlideOcurrence=0
end subroutine 
subroutine slide_ocurrence(cell,StorageT3,MaxStoT3) !Evalua la ocurrencia o no de deslizamientos
	!Variables de entrada
	real StorageT3, MaxStoT3
	integer, intent(in) :: cell
	!Variables locales
	real Zw,Num,Den !Profunidad perched
	!Evalua si el suelo es suceptible 
	if (RiskVector(1,cell) .eq. 1 .and. SlideOcurrence(1,cell) .ne. 1) then
		!Calcula la profundidad encharcada en el tanque gravitacional
		Zw=Zs(1,cell)*(StorageT3/MaxStoT3)
		!Evalua si la profundida emparamada es mayor o igual a la critica
		if (Zw .ge. Zcrit(1,cell)) then
			!En caso afirmativo hay falla en el suelo 
			SlideOcurrence(1,cell)=1
			!Actualiza el tipo de celda y aguas abajo
			call slide_hill2gullie(cell)
		else 
			Num=Cohesion(1,cell)+(GammaS(1,cell)*Zs(1,cell)-Zw*GammaW)&
				&*(cos(slope_hill(1,cell)))**2*TAN(FrictionAngle(1,cell))
			Den=GammaS(1,cell)*Zs(1,cell)*sin(slope_hill(1,cell))*cos(slope_hill(1,cell))
			!Prueba si es menor al factorde seguridad 
			if (Num/Den .le. FS) then 
				!Si esta vaiana es mas baja que el factor de seguridad desliza 
				SlideOcurrence(1,cell)=1
				!Actualiza el tipo de celdas aguasd abajo
				call slide_hill2gullie(cell)
			endif
		endif
	endif
end subroutine
subroutine slide_hill2gullie(cell) !Cuando una celda se vuelve en carcava asegura que aguas abajo tambien lo sea
	!Variables de entrada
	integer, intent(in) :: cell	
	!Variables locales 
	logical Flag
	integer localCell,direction
	!Ejecuta hasta encontrar otra carcava aguas abajo
	if (unit_type(1,cell).le.2) then 
		unit_type(1,cell)=2
		flag=.True.
		localCell=cell
		do while (flag)
			direction=nceldas-drena(1,localCell)+1
			if (unit_type(1,direction) .ge. unit_type(1,localCell)) then
				flag=.False.
			else
				unit_type(1,direction)=unit_type(1,localCell)
				localCell=direction
			endif
		enddo
	endif
end subroutine 

end module
