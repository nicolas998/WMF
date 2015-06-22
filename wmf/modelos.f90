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
real dt !Delta del tiempo de modelacion
integer ncols,nrows !cantidad de columnas y filas del mapa
integer nceldas !Cantidad de celdas que componen la cuenca
!Variables de mapas que pueden ser usados por modelos
integer, allocatable :: tipo_celda(:,:) !Tipo de celda
real, allocatable :: acum_cel(:,:) !Area acumlada km2
integer, allocatable :: drena(:,:) !Indicador de la topologia de la cuenca 
real, allocatable :: L_celda(:,:),pend_celda(:,:) !Long celda, pend de celda
real, allocatable :: pend_cauce(:,:),Elev(:,:),X(:,:),Y(:,:) !pend de cauce, elevacion, coordenadas
real, allocatable :: Hu(:,:),H3(:,:),Man(:,:) !Alm capilar, Alm maximo sub-sup, Manning
real, allocatable :: Ks(:,:),Kp(:,:),EVP(:,:) !Cond sub-sup, Cond supterranea, EV potencial
!conversores de las variables ks y kp
real conver_ks,conver_kp !Conversores de unidades de ks y kp
!Variables de control 
integer, allocatable :: control(:,:) !Celdas de la cuenca que tienen puntos de control
integer, allocatable :: control_h(:,:) !Celdas de la cuenca que son puntos de control de humedad
integer, allocatable :: guarda_cond(:,:) !Intervalos de tiempo en que se hace guardado de condiciones
!Variables de mapas modelo shia
real, allocatable :: S(:,:) !Almacenamiento de los 5 tanques del modelo
!Variables de evaluacion
integer eval !Variable de evaluacion de que las variables necesarias para el modelo estan asignadas
!Variables de lluvia
real, allocatable :: coord(:,:) !Variable con las coordenadas de las estaciones de lluvia
integer, allocatable :: tin(:,:) !Variable que contiene los triangulos conformados a partir de las coordenadas
integer, allocatable :: tin_perte(:,:) !Variable con la pertenencia de los triangulos
real, allocatable :: lluvia(:,:) !Registros de lluvia
real, allocatable :: evp_p(:,:) !Registros de evp tomados de series
real, allocatable :: campo(:,:) !Campo en el ultimo intervalo
real, allocatable :: campo_medio(:,:) !Promedio de los campos de lluvia sobre el modelo en el evento 
!Variables de regionalizacion del flujo
real epsilo,exp1 !coeficiente manning y exponente manning en carcavas [0.5, 0.64]
real alfa(2),sigma(3),Omega,c1,k,fhi !Parametros de la OCG (Velez, 2001)
real w1,w2,w3,B !resultados regionales de la OCG
!Variables de la velocidad calculada en laderas y cauce
real, allocatable :: vel_ladera(:,:), vel_cauce(:,:)
!variables de sedimentos Par aalojar volumens y demas
real wi(3), Qskr,G,diametro(3)
real ERO(3),EROt(3),DEP(3),DEPt(3)
real, allocatable :: VolERO(:),VolDEPo(:) !Volumen erosionado por celda, Volumen depositado en la celda
real, allocatable :: Vs(:,:),Vd(:,:)
real, allocatable :: VSc(:,:),Vdc(:,:)
!Mapas para el calculo de sedimentos
real, allocatable :: Krus(:,:),Crus(:,:),Prus(:,:) !Factores de la RUSLE
real, allocatable :: PArLiAc(:,:) !Porcentaje de arenas, limos y arcillas

contains

!-----------------------------------------------------------------------
!Funciones Varias
!-----------------------------------------------------------------------
!Funciones para estimar variables del modelo
subroutine OCG_params
    real eB
    B=Omega*(c1*k**(alfa(1)-alfa(2)))**((2/3.0)-alfa(2))
    eB=1.0/(1+alfa(2)*((2/3.0)-sigma(2)))
    w1=((2/3.0)-sigma(2))*(1.0-alfa(2))*eB
    w2=(1+alfa(2)*((2/3.0)-sigma(2)))/(fhi*((2/3.0)-sigma(2))*(alfa(1)-alfa(2))+sigma(1))
    w2=(fhi*(0.667-sigma(2))*(alfa(2)-alfa(1))+sigma(1))*eB
    w3=(0.5-sigma(3))*eB
    B=B**(-eB)
    !Parametros de onda cinematica en ladera
    if (epsilo .eq. 0.0) epsilo=0.50 !Valores por defecto de foster y lane 1981
    if (exp1 .eq. 0.0) exp1=0.64
end subroutine
!verificador del modelo shia
subroutine shia_verify
    eval=0
    if (allocated(tipo_celda).and.allocated(drena).and.allocated(L_celda).and.allocated(pend_celda)&
    &.and.allocated(pend_cauce).and.allocated(X).and.allocated(Y).and.allocated(Hu)&
    &.and.allocated(Man).and.allocated(ks).and.allocated(kp).and.allocated(evp).and.allocated(S)) then
	eval=1
    endif
end subroutine
!Arroja una matriz donde se indica a que triangulo pertenece cada celda 
subroutine pre_tin_points_select(resultado,nceldas) 
    !Variables de entrada
    integer, intent(in) :: nceldas
    !Variables de salida
    integer, intent(out) :: resultado
    !f2py intent(in) :: nceldas
    !f2py intent(out) :: resultado
    !Variables locales
    integer i,j,flag,flag2,cont,N_tin,n(2)
    real ori1,ori2,ori3,oriP,xp,yp,x1,x2,x3,y1,y2,y3
    !Comienza la busqueda de pertenencia para cada entrada de la tabla
    if (allocated(tin_perte).eqv. .false.) then
	allocate(tin_perte(1,nceldas))
    endif
    flag2=0
    cont=0
    resultado=0
    n=shape(tin)
    N_tin=n(2)
    tin_perte=noData
    do i=1,nceldas
	!Inicializa una bandera que revisa si el punto se encuentra dentro de algún triangulo
	flag=0
	!Comienza a buscar el triangulo de pertenencia
	j=1
	do while (j.le.N_tin .and. flag.eq.0)
	    !Cambio de variables para programar facil
	    xp=x(1,i); yp=y(1,i)
	    x1=coord(1,TIN(1,j)); x2=coord(1,TIN(2,j)); x3=coord(1,TIN(3,j))
	    y1=coord(2,TIN(1,j)); y2=coord(2,TIN(2,j)); y3=coord(2,TIN(3,j))
	    !Calcula la orientación del triangulo    
	    oriP=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)
	    ori1=(x1-xp)*(y2-yp)-(y1-yp)*(x2-xp)
	    ori2=(xp-x3)*(y2-y3)-(yp-y3)*(x2-x3)
	    ori3=(x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)
	    !Si las cuatro orientaciones son o positivas o negativas el punto se encuentra dentro
	    if ((oriP.gt.0 .and. ori1.gt.0 .and. ori2.gt.0 .and. ori3.gt.0) .or. (oriP.lt.0 .and. ori1.lt.0 &
	    &.and. ori2.lt.0 .and. ori3.lt.0)) then
		tin_perte(1,i)=j
		flag=1
	    else
		j=j+1
	    end if
	end do
    end do
end subroutine
!Lee los datos flotantes de un binario de cuenca en los records ordenados
subroutine read_float_basin(ruta,records,vect,Res,nrecords,nceldas) 
    !Variables de entrada
    integer, intent(in) :: nrecords,nceldas
    character*255, intent(in) :: ruta
    integer, intent(in) :: records(nrecords)
    !Variables de salida
    real, intent(out) :: vect(nrecords,nceldas)
    integer, intent(out) :: Res
    !f2py intent(in) :: nrecords,nceldas,ruta,records
    !f2py intent(out) :: vect
    !Variables locales
    integer i 
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	do i=1,nrecords
	    read(10,rec=records(i),iostat=Res) vect(i,:)
	    if (Res.ne.0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
	enddo
    close(10)
end subroutine
!Lee los datos flotantes de un binario de cuenca en los records ordenados
subroutine write_float_basin(ruta,vect,nrecords,nceldas) 
    !Variables de entrada
    integer, intent(in) :: nrecords,nceldas
    character*255, intent(in) :: ruta
    real, intent(in) :: vect(nrecords,nceldas)
    !f2py intent(in) :: nrecords,nceldas,vect,ruta
    !Variables locales
    integer i 
    !Escritura     
    open(10,file=ruta,form='unformatted',status='replace',access='direct',RECL=4*nceldas)
		do i=1,nrecords
		    write(10,rec=i) vect(i,:)	    
		enddo
    close(10)
end subroutine

!-----------------------------------------------------------------------
!Subrutinas de transporte de sedimentos
subroutine allocate_sed !Funcion para alojar variables si se van a calcular sed
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
subroutine hillslope_sed(alfa,S2,v2,So,area_sec,celda,drena_id,tipo) !Subrutina para calcular los sedimentos en ladera
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
subroutine channel_sed(S5,v5,Q5,So,area_sec,celda,drena_id,VolSal) !subrutina para calcular seduimentos en el canals
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
!Versiones modelo shia
!-----------------------------------------------------------------------

!model_tin_lk_sed: Modelo shia completo, contiene: velocidad lineal en ladera 
!y velocidad cinematica en canal, no se usa la OCG. la lluvia es interpolada
!por idw o tin de acuerdo a lo indicado.
subroutine shia_tin_lk_sed(calib,N_cel,N_cont,N_reg,si,Q,Qsed,balance,mean_rain,sp)
    !Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_cont
    real, intent(in) :: calib(10), si(5,N_cel)
    !Variables de salia
    real, intent(out) :: Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance(5,N_reg)
    real, intent(out) :: Qsed(3,N_cont,N_reg)
    !Variables locales 
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid,Res !Vector con el id a donde drena cada celda
    real Rain,rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
    character*3 corr !local de si se va a hacer correccion o no.
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real ax(nceldas),ay(nceldas),bx(nceldas),by(nceldas),cx(nceldas),cy(nceldas),dx(nceldas),dy(nceldas)
    real det1,det2,det3,det4,az,bz,cz,Cel_x,Cel_y,coef1(nceldas),coef2 !Variables de interpolacion tin
    integer Cel_pert !Variables de interpolacion tin
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas,Satras,Sadelante !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real Area_coef(nceldas) !Coeficiente para el calculo del lateral en cada celda del tanque 2 para calcuo de sedimentos
    real Vsal_sed(3) !Volumen de salida de cada fraccion de sedimentos [m3/seg]
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante    
    !Definiciones para f2py
    !f2py intent(in) :: N_cel,N_reg,N_cont,calib,si
    !f2py intent(out) :: Q,Qsed,sp,mean_rain,balance
    !Preambulo a la ejecucion del modelo
    !Calcula la cantidad de agua que es retenida en el capilar
    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
    !Calcula las velocidades horizontales iniciales
    v_ladera=Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
    v_cauce=0.5 ![m/s] Vel ladera, Estado constante
    pend_man=B*(acum_cel(1,:)**w2)*(pend_celda(1,:)**w3)!(pend_celda(1,:)**(0.5))/man(1,:)
    ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
    kph=kp(1,:)*Calib(8)*(conver_kp/1000)!*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
    !Calcula el porcentaje de flujo horizontal saliente de cada celda
    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
    !Calcula las constantes para interpolacion por TIN
    ax=coord(1,tin(1,tin_perte(1,:))) ; ay=coord(2,tin(1,tin_perte(1,:)))
    bx=coord(1,tin(2,tin_perte(1,:))) ; by=coord(2,tin(2,tin_perte(1,:)))
    cx=coord(1,tin(3,tin_perte(1,:))) ; cy=coord(2,tin(3,tin_perte(1,:)))
    coef1=(bx-ax)*(cy-ay)-(cx-ax)*(by-ay)
    !Convertor de mm a mt3 para obtener caudal a la salida
    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
    !Constante para el calculo del area lateral para el calculo de sedimentos 
    Area_coef=m3_mm/(dxp+v_ladera*dt)
    Qsed=0 !Inicia el caudal de sedimentos en cero
    !Establece el almacenamiento inicial 
    Sp=Si
    !ASigna tamanos a las variables de las velocidades
    if (allocated(vel_ladera) .eqv. .false.) allocate(vel_ladera(1,nceldas))
    vel_ladera(1,:)=v_ladera
    if (allocated(vel_cauce) .eqv. .false.) allocate(vel_cauce(1,nceldas))
    vel_cauce(1,:)=0
    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
    inicial=sum(Sp(:,:))
    entradas=0
    salidas=0
    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
    mean_rain=0; balance=0
    if (eval.eq.1) then
	    !Calcula el tiempo inicial
	    call etime(tiempo_ejec,tiempo_result)
	    !Inicia variables de sedimentos
	    call allocate_sed
	    !Itera para cada intervalo de tiempo
	    do tiempo=1,N_reg
		!Reinicia el conteo de la lluvia y el contador de los puntos de control
		rain_sum=0
		control_cont=2	    	    
		!Itera para todas las celdas
		do celdas=1,N_cel
			!Interpola la lluvia
			Cel_pert=tin_perte(1,celdas)
			Cel_x=x(1,celdas) ; Cel_y=y(1,celdas)
			az=max(lluvia(tin(1,Cel_pert),tiempo),0.0)
			bz=max(lluvia(tin(2,Cel_pert),tiempo),0.0)
			cz=max(lluvia(tin(3,Cel_pert),tiempo),0.0)
			det1=(cx(celdas)-ax(celdas))*(Cel_y-ay(celdas)) 
			det2=(by(celdas)-ay(celdas))*(Cel_x-ax(celdas)) 
			det3=(cy(celdas)-ay(celdas))*(Cel_x-ax(celdas))
			det4=(Cel_y-ay(celdas))*(bx(celdas)-ax(celdas))
			coef2=det1*(bz-az)+det2*(cz-az)-det3*(bz-az)-det4*(cz-az)
			!Calcula el flujo vertical por la celda
			R1=max(az-coef2/coef1(celdas),0.0)		
			R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]		
			R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 				
			R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
			R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero				
			!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
			drenaid=N_cel-drena(1,celdas)+1	    		
			!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
			Satras=Sp(1,celdas)
			Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
			E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]		
			Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]		
			balance(1,tiempo)=balance(1,tiempo)+Sp(1,celdas)-Satras-R1+R2+E1
			salidas=salidas+E1		 
			!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
			Satras=Sp(2,celdas)
			Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
			Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
			Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         		
			balance(2,tiempo)=balance(2,tiempo)+Sp(2,celdas)-Satras-R2+R3+Es2
			!Calculo de sedimentos en ladera
			Area=Sp(2,celdas)*Area_coef(celdas)
			call hillslope_sed(calib(10),Sp(2,celdas),v_ladera(celdas),pend_celda(1,celdas),&
			&Area,celdas,drenaid,tipo_celda(1,celdas))
			!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
			Satras=Sp(3,celdas)
			Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm] 
			Es3=E3(celdas)*Sp(3,celdas) ![mm]
			Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
			balance(3,tiempo)=balance(3,tiempo)+Sp(3,celdas)-Satras-R3+R4+Es3
			!Tanque 4: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
			Satras=Sp(4,celdas)
			Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
			Es4=E4(celdas)*Sp(4,celdas) ![mm]
			Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    				
			balance(4,tiempo)=balance(4,tiempo)+Sp(4,celdas)-Satras-R4+R5+Es4
			!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
			select case (tipo_celda(1,celdas))
			    case(1)
				    !Sigue drenando a su tanque igual aguas abajo
				    if (drena(1,celdas).ne.0) then
					Sp(2,drenaid)=Sp(2,drenaid)+Es2
					Sp(3,drenaid)=Sp(3,drenaid)+Es3
					Sp(4,drenaid)=Sp(4,drenaid)+Es4
				    else
					Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
					salidas=salidas+Es2+Es3+Es4
				    endif
			    case(2,3)			
				!El tanque 5 recibe agua de los tanques 2 y 3 fijo
				Satras=Sp(5,celdas)
				Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
				Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
				!Copia variables de vectores a escalares par aahorrar tiempo
				vel=v_cauce(celdas) ![m/seg]
				S5=Sp(5,celdas) ![mm]			
				L=L_celda(1,celdas) ![mts]
				pm=pend_man(celdas)
				!Calcula mediante onda cinematica el flujo a la salida
				do i=1,4
				    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
				    vn=(Area**(2/3))*pm ![m/seg]
				    vel=(2*vn+vel)/3 ![m2/seg]
				enddo
				v_cauce(celdas)=vel		            			
				vel_cauce(1,celdas)=vel_cauce(1,celdas)+vel
				Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
				Sp(5,celdas)=Sp(5,celdas)-Es5			
				balance(5,tiempo)=balance(5,tiempo)+Sp(5,celdas)-Satras-Es2-Es3-Es4*(tipo_celda(1,celdas)-2)+Es5
				!Sedimentos en el cauce
				call channel_sed(Sp(5,celdas),vel,Es5*m3_mm/dt,pend_celda(1,celdas),Area,celdas,drenaid,Vsal_sed)
				if (drena(1,celdas).ne.0) then 
					Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]							
				endif			
			end select
			!si es la salida de la cuenca guarda
			if (drena(1,celdas).eq.0) then 		   
			    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]    
			    do i=1,3
					Qsed(i,1,tiempo)=Vsal_sed(i)
				enddo  
			    salidas=salidas+Es5 ![mm]
			endif
			!si es un punto de control guarda igualmente
			if (control(1,celdas).ne.0) then
			    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
			    do i=1,3
					Qsed(i,control_cont,tiempo)=Vsal_sed(i)
			    enddo
			    control_cont=control_cont+1
			endif
			!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
			rain_sum=rain_sum+R1
		enddo
		balance(1,tiempo)=balance(1,tiempo)/N_cel
		balance(2,tiempo)=balance(2,tiempo)/N_cel
		!calcula cantidad total de lluvia que entra al sistema
		entradas=entradas+rain_sum
		!Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
		mean_rain(tiempo)=rain_sum/N_cel
	    enddo
	    !Calcula y muestra el tiempo de ejecucion
	    vel_cauce=vel_cauce/N_reg
	    call etime(tiempo_ejec,tiempo_result)
	    print *, 'Tiempo ejecucion: ', tiempo_result
	else
	    print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
    endif
end subroutine

!model_tin_lk: Modelo shia completo, contiene: velocidad lineal en ladera 
!y velocidad cinematica en canal, no se usa la OCG. la lluvia es interpolada
!por idw o tin de acuerdo a lo indicado.
subroutine shia_tin_lk(calib,N_cel,N_cont,N_reg,N_est,si,Q,balance,mean_rain,campo_precip,sp &
	&,lluv_param,ruta_store,kinematic_hill)
    !Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_est,N_cont
    real, intent(in) :: calib(9),si(5,N_cel)
    character, intent(in), optional :: ruta_store*255
    integer, intent(in), optional :: kinematic_hill
    real, intent(in) :: lluv_param(N_est,N_reg)
    !Variables de salia
    real, intent(out) :: Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance(5,N_reg),campo_precip(1,N_cel)
    !Variables locales 
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid,Res !Vector con el id a donde drena cada celda
    real Rain,rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
    character*3 corr !local de si se va a hacer correccion o no.
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas),FactorVlad(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real ax(nceldas),ay(nceldas),bx(nceldas),by(nceldas),cx(nceldas),cy(nceldas),dx(nceldas),dy(nceldas)
    real det1,det2,det3,det4,az,bz,cz,Cel_x,Cel_y,coef1(nceldas),coef2 !Variables de interpolacion tin
    real Precip(N_est,N_reg)
    integer Cel_pert !Variables de interpolacion tin
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas,Satras,Sadelante !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante    
    character ruta_temp*270,tiempo_texto*10 !Ruta temporal para el guardado de condiciones intermedias
    !Definiciones para f2py
    !f2py intent(in) :: N_cel,N_reg,N-est,N_cont,calib,si,lluv_param
    !f2py intent(in),optional :: ruta_store,kinematic_hill
    !f2py intent(out) :: Q,sp,mean_rain,balance,campo_precip    
    !Preambulo a la ejecucion del modelo    
    !Calcula la cantidad de agua que es retenida en el capilar
    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua    
    !Velocidad inicial en ladera
    if (allocated(vel_ladera) .eqv. .false.) allocate(vel_ladera(1,nceldas))    
    if (kinematic_hill .ne. 1) then		
		v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
		E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
		vel_ladera(1,:)=v_ladera
    else
		FactorVlad=(epsilo*pend_celda(1,:)**0.5)/man(1,:) !Adim, la velocidad se calcula con una onda cinematica
		v_ladera=1
		vel_ladera(1,:)=0.0
	endif
    !Velocidad inicial en cauce
    v_cauce=1 ![m/s] Punto inicial de la velocidad en el cauce
    pend_man=B*(acum_cel(1,:)**w2)*(pend_celda(1,:)**w3) !Onda cinematica geomorfoliogica    
    if (allocated(vel_cauce) .eqv. .false.) allocate(vel_cauce(1,nceldas))
    vel_cauce(1,:)=0.0 !Variable global de velocidad en cauce [es un control]    
    !Velocidades constantes sub-superficial y subterranea
    ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
    kph=kp(1,:)*Calib(8)*(conver_kp/1000)!*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
    !Calcula el porcentaje de flujo horizontal saliente de cada celda  en tanques lineales
    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))    
    !Calcula las constantes para interpolacion por TIN
    ax=coord(1,tin(1,tin_perte(1,:))) ; ay=coord(2,tin(1,tin_perte(1,:)))
    bx=coord(1,tin(2,tin_perte(1,:))) ; by=coord(2,tin(2,tin_perte(1,:)))
    cx=coord(1,tin(3,tin_perte(1,:))) ; cy=coord(2,tin(3,tin_perte(1,:)))
    coef1=(bx-ax)*(cy-ay)-(cx-ax)*(by-ay)
    !Convertor de mm a mt3 para obtener caudal a la salida
    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
    !Establece el almacenamiento inicial 
    Sp=Si    
    !Asigna tamanos de las variables de la lluvia
    if (allocated(campo) .eqv. .false.) allocate(campo(1,nceldas))    
    if (allocated(campo_medio) .eqv. .false.) allocate(campo_medio(1,nceldas))
    campo_medio(1,:)=0
    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
    inicial=sum(Sp(:,:))
    entradas=0
    salidas=0    
    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
    mean_rain=0; balance=0
    if (eval.eq.1) then
	!Calcula el tiempo inicial
	call etime(tiempo_ejec,tiempo_result)
	!Itera para cada intervalo de tiempo
	do tiempo=1,N_reg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2	    	    
	    campo(1,:)=0
	    !Itera para todas las celdas
	    do celdas=1,N_cel
		!Interpola la lluvia
		Cel_pert=tin_perte(1,celdas)
		Cel_x=x(1,celdas) ; Cel_y=y(1,celdas)
		az=max(lluv_param(tin(1,Cel_pert),tiempo),0.0)
		bz=max(lluv_param(tin(2,Cel_pert),tiempo),0.0)
		cz=max(lluv_param(tin(3,Cel_pert),tiempo),0.0)
		det1=(cx(celdas)-ax(celdas))*(Cel_y-ay(celdas)) 
		det2=(by(celdas)-ay(celdas))*(Cel_x-ax(celdas)) 
		det3=(cy(celdas)-ay(celdas))*(Cel_x-ax(celdas))
		det4=(Cel_y-ay(celdas))*(bx(celdas)-ax(celdas))
		coef2=det1*(bz-az)+det2*(cz-az)-det3*(bz-az)-det4*(cz-az)
		!Calcula el flujo vertical por la celda
		R1=max(az-coef2/coef1(celdas),0.0)		
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]		
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 				
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero				
		!Calcula la lluvia en el intervalo
		campo(1,celdas)=R1
		campo_medio(1,celdas)=campo_medio(1,celdas)+R1
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=N_cel-drena(1,celdas)+1	    				
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		Satras=Sp(1,celdas)
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]		
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]		
		balance(1,tiempo)=balance(1,tiempo)+Sp(1,celdas)-Satras-R1+R2+E1
		salidas=salidas+E1		 
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		Satras=Sp(2,celdas)
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		if (kinematic_hill .ne. 1) then
			Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		else 
			pm=FactorVlad(celdas); vel=v_ladera(celdas)
			do i=1,3
			    Area=Sp(2,celdas)*m3_mm/(dxp+vel*dt) ![m2] Calcula el area de la seccion
			    vn=pm*(Area**(0.6666*exp1)) ![m/seg] Calcula la velocidad nueva
			    vel=(2*vn+vel)/3.0 ![m/seg] Promedia la velocidad
			enddo
			v_ladera(celdas)=vel ![m/seg]
			vel_ladera(1,celdas)=vel_ladera(1,celdas)+vel
			Es2=min(Area*vel*Calib(4)*dt/m3_mm,Sp(2,celdas)) ![mm]		
		endif
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         		
		balance(2,tiempo)=balance(2,tiempo)+Sp(2,celdas)-Satras-R2+R3+Es2
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
		Satras=Sp(3,celdas)
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm] 
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		balance(3,tiempo)=balance(3,tiempo)+Sp(3,celdas)-Satras-R3+R4+Es3
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Satras=Sp(4,celdas)
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    				
		balance(4,tiempo)=balance(4,tiempo)+Sp(4,celdas)-Satras-R4+R5+Es4		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(1,celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2+Es3+Es4
			endif
		    case(2,3)			
			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
			Satras=Sp(5,celdas)
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante onda cinematica el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(w1))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			vel_cauce(1,celdas)=vel_cauce(1,celdas)+vel
			Es5=min(Area*v_cauce(celdas)*Calib(9)*dt/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5			
			balance(5,tiempo)=balance(5,tiempo)+Sp(5,celdas)-Satras-Es2-Es3-Es4*(tipo_celda(1,celdas)-2)+Es5
			if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select
		!si es la salida de la cuenca guarda		
		if (drena(1,celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(1,celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    balance(1,tiempo)=balance(1,tiempo)/N_cel
	    balance(2,tiempo)=balance(2,tiempo)/N_cel
	    !calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/N_cel
	    !Guardado de condiciones intermedias de simulacion
	    if (allocated(guarda_cond)) then
			if (guarda_cond(1,tiempo).ne.0.0) then				
				write(tiempo_texto,'(I10)') guarda_cond(1,tiempo) 
				tiempo_texto=adjustl(tiempo_texto)
				ruta_temp=trim(ruta_store) // 'Store_' // trim(tiempo_texto) // '.sto'
				ruta_temp=adjustl(ruta_temp); ruta_temp=trim(ruta_temp)				
				call write_float_basin(ruta_temp,Sp,5,nceldas) 
			endif
	    endif
	enddo
	!Calcula y muestra el tiempo de ejecucion
	vel_cauce=vel_cauce/N_reg
	if (present(kinematic_hill)) vel_ladera=vel_ladera/N_reg	
	campo_medio=campo_medio/N_reg
	campo_precip=campo_medio
	call etime(tiempo_ejec,tiempo_result)
	!print *, 'Tiempo ejecucion: ', tiempo_result
    else
	print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
    endif
end subroutine

!model_interp_lk: Modelo shia completo, contiene: velocidad lineal en ladera 
!y velocidad cinematica en canal, no se usa la OCG. la lluvia es interpolada
!por idw. 
subroutine shia_idw_lk(calib,N_cel,N_cont,N_reg,Num_est,pp,si,Q,balance,mean_rain,sp)
    !Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_cont,Num_est
    real, intent(in) :: calib(9),si(5,N_cel)
    real, intent(in) :: pp
    !Variables de salia
    real, intent(out) :: Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance(5,N_reg)
    !Variables locales 
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid,Res !Vector con el id a donde drena cada celda
    real Rain,rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
    character*3 corr !local de si se va a hacer correccion o no.
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real W(Num_est,nceldas),Cel_x,Cel_y,Wr
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante    
    real Satras !Guarda el almacenamiento anterior en el tanque [mm]
    !Definiciones para f2py
    !f2py intent(in) :: N_cel,N_reg,N_cont,calib,Num_est,pp,si
    !f2py intent(out) :: Q,sp,mean_rain,balance
    !Preambulo a la ejecucion del modelo
    !Calcula la cantidad de agua que es retenida en el capilar
    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
    !Calcula las velocidades horizontales iniciales
    v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
    v_cauce=1 ![m/s] Vel ladera, Estado constante
    pend_man=(pend_celda(1,:)**(0.5))/man(1,:)
    ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
    kph=kp(1,:)*Calib(8)*(conver_kp/1000)!*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
    !Calcula el porcentaje de flujo horizontal saliente de cada celda
    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))    
    !E4=0.08
    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
    !Convertor de mm a mt3 para obtener caudal a la salida
    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
    !Establece el almacenamiento inicial 
    Sp=Si
    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
    inicial=sum(Sp(:,:))
    entradas=0
    salidas=0
    !Calcula la ponderacion de cada celda a cada estacion
    do i=1,Num_est
	W(i,:)=1.0/(sqrt(((coord(1,:)-x(1,:))**2+(coord(2,:)-y(1,:))**2)))**pp
    end do
    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
    mean_rain=0
    if (eval.eq.1) then
	!Calcula el tiempo inicial
	call etime(tiempo_ejec,tiempo_result)
	!Itera para cada intervalo de tiempo
	do tiempo=1,N_reg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2	    	    
	    !Itera para todas las celdas
	    do celdas=1,N_cel
		!Interpola la lluvia
		Wr=sum(W(:,celdas)*lluvia(:,tiempo),mask=lluvia(:,tiempo).gt.0.0)
		R1=max(Wr/sum(W(:,celdas),mask=lluvia(:,tiempo).ge.0.0),0.0)		
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]		
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 				
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero				
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=N_cel-drena(1,celdas)+1	    		
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		Satras=Sp(1,celdas)
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]		
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]		
		balance(1,tiempo)=balance(1,tiempo)+Sp(1,celdas)-Satras-R1+R2+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		Satras=Sp(2,celdas)
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         		
		balance(2,tiempo)=balance(2,tiempo)+Sp(2,celdas)-Satras-R2+R3+Es2
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
		Satras=Sp(3,celdas)
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm] 
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		balance(3,tiempo)=balance(3,tiempo)+Sp(3,celdas)-Satras-R3+R4+Es3
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Satras=Sp(4,celdas)
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    				
		balance(4,tiempo)=balance(4,tiempo)+Sp(4,celdas)-Satras-R4+R5+Es4
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(1,celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]			    
			endif
		    case(2,3)			
			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
			Satras=Sp(5,celdas)
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante onda cinematica el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5			
			
			balance(5,tiempo)=balance(5,tiempo)+Sp(5,celdas)-Satras-Es2-Es3-Es4*(tipo_celda(1,celdas)-2)+Es5
			if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select
		!si es la salida de la cuenca guarda
		if (drena(1,celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(1,celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !balance(1,tiempo)=balance(1:4,tiempo)/N_cel	    
	    !calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/N_cel
	enddo
	!Calcula y muestra el tiempo de ejecucion
	call etime(tiempo_ejec,tiempo_result)
	print *, 'Tiempo ejecucion: ', tiempo_result
    else
	print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
    endif    
end subroutine

!model_idw_klk: Modelo shia completo, contiene: velocidad lineal en ladera,
!velociad no lineal sub-superficial, velocidad cinematica en canal, no se usa la OCG. 
!la lluvia es interpolada por idw. 
subroutine shia_idw_klk(calib,N_cel,N_cont,N_contH,N_reg,Num_est,pp,Q,Hum,Etr,Infiltra,mean_rain,sp)
    !Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_cont,Num_est,N_contH
    real, intent(in) :: calib(10)
    real, intent(in) :: pp
    !Variables de salia
    real, intent(out) :: Hum(N_contH,N_reg),Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),Etr(N_contH,N_reg)
    real, intent(out) :: Infiltra(N_contH,N_reg)
    !Variables locales 
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid,Res !Vector con el id a donde drena cada celda
    real Rain,rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
    character*3 corr !local de si se va a hacer correccion o no.
    real Hu_loc(nceldas),H3_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5,Z3 !lluvia (R1) y flujo que drena por los tanques, y retorno del tanque 3
    real v_ladera(nceldas),v_cauce(nceldas),v_sub(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real W(Num_est),Cel_x,Cel_y,Wr
    integer Cel_pert !Variables de interpolacion tin
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    integer i,control_cont,cont_h,cont !iterador,contador de puntos de control, contador de lluvia faltante    
    !Definiciones para f2py
    !f2py intent(in) :: N_cel,N_reg,N_cont,calib,Num_est,pp,N_contH
    !f2py intent(out) :: Q,Hum,sp,mean_rain,Etr,Infiltra
    !Preambulo a la ejecucion del modelo
    !Convertor de mm a mt3 para obtener caudal a la salida    
    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
    !Calcula la cantidad de agua que es retenida en el capilar
    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
    H3_loc=H3(1,:)*Calib(10) !Cantidad de agua que es almacenada gravitacionalmente
    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
    !Calcula las velocidades horizontales iniciales
    v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
    v_cauce=1; v_sub=0.1 ![m/s] Vel inicial cauce, Vel inicial sub-superficial
    pend_man=(pend_celda(1,:)**(0.5))/man(1,:)
    ksh=ks(1,:)*(conver_ks/1000) ![m/s] Convierte la velosidad de unidades		    
    ksh=(Calib(6)*ksh*pend_celda(1,:)*dxP**2)/(3*(H3_loc*m3_mm)**2) !Expresion constante de kubota y Sivapalan, 1995        
    kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
    !Calcula el porcentaje de flujo horizontal saliente de cada celda
    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))    
    !Establece el almacenamiento inicial 
    Sp=S
    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
    inicial=sum(Sp(:,:))
    entradas=0
    salidas=0
    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
    mean_rain=0 
    if (eval.eq.1) then
	!Calcula el tiempo inicial
	call etime(tiempo_ejec,tiempo_result)
	!Itera para cada intervalo de tiempo
	do tiempo=1,N_reg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2	    	    
	    cont_h=0
	    !Itera para todas las celdas
	    do celdas=1,N_cel
		!Interpola la lluvia		
		Cel_x=x(1,celdas) ; Cel_y=y(1,celdas)		
		do i=1,Num_est
		    W(i)=1.0/(sqrt(((coord(1,i)-Cel_x)**2+(coord(2,i)-Cel_y)**2)))**pp
		end do				
		Wr=sum(W*lluvia(:,tiempo),mask=lluvia(:,tiempo).gt.0.0)
		R1=Wr/sum(W,mask=lluvia(:,tiempo).ge.0.0)		
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=N_cel-drena(1,celdas)+1	    
		!Al sistema le entra la lluvia
		entradas=entradas+R1						
		!Flujo vertical entre tanques		
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		R5=min(R4,mm_kpp(celdas)) !pierde el minimo entre lo que hay y la conductividad del acuifero						
		!Determina si hay o no retorno del tanque 3 
		Z3=max(0.0,Sp(3,celdas)+R3-R4-H3_loc(celdas))	    		   
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial				
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(evp_p(1,tiempo)*Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 			    
		Sp(2,celdas)=Sp(2,celdas)+R2+Z3-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         				
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial				
		Sp(3,celdas)=Sp(3,celdas)+R3-Z3-R4 ![mm] 
		pm=ksh(celdas); vel=v_sub(celdas) !Copia variables por velocidad
		do i=1,4
		    Area=Sp(3,celdas)*m3_mm/(dxp+vel*dt) ![m2] Calcula el area de la seccion
		    vn=pm*(Area**2) ![m/seg] Calcula la velocidad nueva
		    vel=(2*vn+vel)/3.0 ![m/seg] Promedia la velocidad
		enddo		
		v_sub(celdas)=vel ![m/seg]
		Es3=Area*vel*dt/m3_mm ![mm]		
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		if (control_h(1,celdas).ne.0) then
		    cont_h=cont_h+1
		    Hum(cont_h,tiempo)=(Sp(3,celdas)+Sp(1,celdas))/(Hu_loc(celdas)+H3_loc(celdas))
		    Etr(cont_h,tiempo)=E1
		    Infiltra(cont_h,tiempo)=R3
		endif
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo subterraneo
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(1,celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2+Es3+Es4
			endif
		    case(2,3)			
			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
			!Copia variables de vectores a escalares par ahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante onda cinematica el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5
			if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select		
		!si es la salida de la cuenca guarda
		if (drena(1,celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif		
		!si es un punto de control guarda igualmente		
		if (control(1,celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif		
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula cantidad total de lluvia que entra al sistema	    
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/N_cel	   
	enddo
	!Calcula y muestra el tiempo de ejecucion
	call etime(tiempo_ejec,tiempo_result)
	print *, 'Tiempo ejecucion: ', tiempo_result
    else
	print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
    endif
end subroutine

!model_idw_kkk: Modelo shia completo, contiene: velocidad no lineal en ladera,
!velociad no lineal sub-superficial, velocidad cinematica en canal, no se usa la OCG. 
!la lluvia es interpolada por idw. 
subroutine shia_idw_kkk(calib,N_cel,N_cont,N_contH,N_reg,Num_est,pp,Q,Hum,balance,mean_rain,sp)
    !Variables de entrada
    integer, intent(in) :: N_cel,N_reg,N_cont,Num_est,N_contH
    real, intent(in) :: calib(10)
    real, intent(in) :: pp
    !Variables de salia
    real, intent(out) :: Hum(N_contH,N_reg),Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance
    !Variables locales 
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid,Res !Vector con el id a donde drena cada celda
    real Rain,rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
    character*3 corr !local de si se va a hacer correccion o no.
    real Hu_loc(nceldas),H3_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5,Z3 !lluvia (R1) y flujo que drena por los tanques, y retorno del tanque 3
    real v_ladera(nceldas),v_cauce(nceldas),v_sub(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real pend_man_lad(nceldas) !combinacion pendiente manning y coeficientes para flujo en laderas
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real W(Num_est),Cel_x,Cel_y,Wr
    integer Cel_pert !Variables de interpolacion tin
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    integer i,control_cont,cont_h,cont !iterador,contador de puntos de control, contador de lluvia faltante    
    !Definiciones para f2py
    !f2py intent(in) :: N_cel,N_reg,N_cont,calib,Num_est,pp,N_contH
    !f2py intent(out) :: Q,Hum,sp,mean_rain,balance
    !Preambulo a la ejecucion del modelo
    !Convertor de mm a mt3 para obtener caudal a la salida    
    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
    !Calcula la cantidad de agua que es retenida en el capilar
    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
    H3_loc=H3(1,:)*Calib(10) !Cantidad de agua que es almacenada gravitacionalmente
    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
    !Calcula las velocidades horizontales iniciales
    v_ladera=0.5 ![m/s] Vel ladera, Estado constante
    v_cauce=1; v_sub=0.1 ![m/s] Vel inicial cauce, Vel inicial sub-superficial
    pend_man=(pend_celda(1,:)**(0.5))/man(1,:) !combinacion pend y manning para flujo en cauce
    pend_man_lad=Calib(4)*(epsilo/man(1,:))*(pend_celda(1,:)**0.5)
    ksh=ks(1,:)*(conver_ks/1000) ![m/s] Convierte la velosidad de unidades		    
    ksh=(Calib(6)*ksh*pend_celda(1,:)*dxP**2)/(3*(H3_loc*m3_mm)**2) !Expresion constante de kubota y Sivapalan, 1995        
    kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
    !Calcula el porcentaje de flujo horizontal saliente de cada celda
    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))    
    !Establece el almacenamiento inicial 
    Sp=S
    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
    inicial=sum(Sp(:,:))
    entradas=0
    salidas=0
    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
    mean_rain=0 
    if (eval.eq.1) then
	!Calcula el tiempo inicial
	call etime(tiempo_ejec,tiempo_result)
	!Itera para cada intervalo de tiempo
	do tiempo=1,N_reg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2	    	    
	    cont_h=0
	    !Itera para todas las celdas
	    do celdas=1,N_cel
		!Interpola la lluvia		
		Cel_x=x(1,celdas) ; Cel_y=y(1,celdas)		
		do i=1,Num_est
		    W(i)=1.0/(sqrt(((coord(1,i)-Cel_x)**2+(coord(2,i)-Cel_y)**2)))**pp
		end do		
		Wr=sum(W*lluvia(:,tiempo),mask=lluvia(:,tiempo).gt.0.0)
		R1=Wr/sum(W,mask=lluvia(:,tiempo).ge.0.0)		
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=N_cel-drena(1,celdas)+1	    
		!Al sistema le entra la lluvia
		entradas=entradas+R1				
		!Flujo vertical entre tanques		
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		R5=min(R4,mm_kpp(celdas)) !pierde el minimo entre lo que hay y la conductividad del acuifero				
		!Determina si hay o no retorno del tanque 3 
		Z3=max(0.0,Sp(3,celdas)+R3-R4-H3_loc(celdas))	    		   
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial				
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(evp_p(1,tiempo)*Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1		
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 			    
		Sp(2,celdas)=Sp(2,celdas)+R2+Z3-R3 ![mm] !Actualiza alm con la lluvia nueva
		pm=pend_man_lad(celdas); vel=v_ladera(celdas)
		do i=1,4
		    Area=Sp(2,celdas)*m3_mm/(dxp+vel*dt) ![m2] Calcula el area de la seccion
		    vn=pm*(Area**2) ![m/seg] Calcula la velocidad nueva
		    vel=(2*vn+vel)/3.0 ![m/seg] Promedia la velocidad
		enddo
		v_ladera(celdas)=vel ![m/seg]
		Es2=Area*vel*dt/m3_mm ![mm]		
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         		
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial		
		Sp(3,celdas)=Sp(3,celdas)+R3-Z3-R4 ![mm] 
		pm=ksh(celdas); vel=v_sub(celdas) !Copia variables por velocidad
		do i=1,4
		    Area=Sp(3,celdas)*m3_mm/(dxp+vel*dt) ![m2] Calcula el area de la seccion
		    vn=pm*(Area**2) ![m/seg] Calcula la velocidad nueva
		    vel=(2*vn+vel)/3.0 ![m/seg] Promedia la velocidad
		enddo		
		v_sub(celdas)=vel ![m/seg]
		Es3=Area*vel*dt/m3_mm ![mm]		
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		if (control_h(1,celdas).ne.0) then
		    cont_h=cont_h+1
		    Hum(cont_h,tiempo)=Sp(3,celdas)+Sp(1,celdas)
		endif		
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo subterraneo
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(1,celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2+Es3+Es4
			endif
		    case(2,3)			
			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
			!Copia variables de vectores a escalares par ahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante onda cinematica el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5
			if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select		
		!si es la salida de la cuenca guarda
		if (drena(1,celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif		
		!si es un punto de control guarda igualmente		
		if (control(1,celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif		
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula cantidad total de lluvia que entra al sistema	    
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/N_cel	   
	enddo
	!Calcula y muestra el tiempo de ejecucion
	call etime(tiempo_ejec,tiempo_result)
	print *, 'Tiempo ejecucion: ', tiempo_result
    else
	print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
    endif
end subroutine

!model_interp_lk: Modelo shia completo, contiene: velocidad lineal en ladera 
!y velocidad cinematica en canal, no se usa la OCG. la lluvia es interpolada
!por idw o tin de acuerdo a lo indicado.
!subroutine shia_interp_lk(calib,tipo,N_cel,N_cont,N_reg,corr_sino,pot,Q,balance,mean_rain,sp)
!    !Variables de entrada
!    integer, intent(in) :: N_cel,N_reg,N_cont
!    real, intent(in) :: calib(9)
!    real, intent(in), optional :: pot
!    character*3, intent(in) :: tipo
!    character*3, intent(in), optional :: corr_sino
!    !Variables de salia
!    real, intent(out) :: Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance
!    !Variables locales 
!    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
!    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
!    integer drenaid,Res !Vector con el id a donde drena cada celda
!    real Rain(nceldas),rain_sum,pot_local !campo de lluvia  para el intervalo de tiempo, potencial para idw local,
!    character*3 corr !local de si se va a hacer correccion o no.
!    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
!    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
!    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
!    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
!    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
!    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
!    real m3_mm !Variable para convertir m3 a mm2
!    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
!    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
!    !Definiciones para f2py
!    !f2py intent(in) :: N_cel,N_reg,N_cont,calib,tipo
!    !f2py intent(in), optional :: pot
!    !f2py intent(out) :: Q,sp,mean_rain,balance
!    !Preambulo a la ejecucion del modelo
!    !Calcula la cantidad de agua que es retenida en el capilar
!    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
!    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
!    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
!    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
!    !Calcula las velocidades horizontales iniciales
!    v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
!    v_cauce=1 ![m/s] Vel ladera, Estado constante
!    pend_man=(pend_celda(1,:)**(0.5))/man(1,:)
!    ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
!    kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
!    !Calcula el porcentaje de flujo horizontal saliente de cada celda
!    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
!    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
!    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
!    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
!    !Convertor de mm a mt3 para obtener caudal a la salida
!    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
!    !Establece el almacenamiento inicial 
!    Sp=S
!    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
!    inicial=sum(Sp(:,:))
!    entradas=0
!    salidas=0
!    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
!    mean_rain=0
!    !Si se va a interpolar con idw verifica el potencial a usar
!    if (tipo.eq.'idw') then
!	if (present(pot)) then
!	    if (pot.lt.1.0) then
!		pot_local=1.0
!	    else
!		pot_local=pot
!	    endif
!	else
!	    pot_local=1.0
!	endif
!    elseif (tipo.eq.'tin') then
!	pot_local=1.0
!    else
!	eval=0
!    endif
!    !Verifica si van a hacer correcciones o no
!	if (present(corr_sino)) then
!	    if (len_trim(corr_sino).gt.0) then
!			corr=corr_sino
!	    else
!			corr='no'
!	    endif
!	endif
!    if (eval.eq.1) then
!		!Itera para cada intervalo de tiempo
!		do tiempo=1,N_reg
!		    !Reinicia el conteo de la lluvia y el contador de los puntos de control
!		    rain_sum=0
!		    control_cont=2
!		    !interpola la lluvia para todas las celdas en el intervalo de tiempo
!		    if (tipo.eq.'tin') then
!				call interpolation_tin_one(Rain,x,y,tiempo,nceldas,correccion=corr_sino)
!		    elseif (tipo.eq.'idw') then
!				call interpolation_idw_one(Rain,x,y,tiempo,nceldas,pp=pot_local,correccion=corr_sino)
!		    endif		    
!		    !Itera para todas las celdas
!		    do celdas=1,N_cel
!				!Verifica la lluvia
!				if (Rain(celdas).lt.0 .or. Rain(celdas).eq.0) then
!				    R1=0.0
!				else
!				    R1=Rain(celdas)
!				endif
!				!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
!				drenaid=N_cel-drena(1,celdas)+1	    
!				!Al sistema le entra la lluvia
!				entradas=entradas+R1
!				!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
!				R2=max(0.0,Rain(celdas)-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
!				Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
!				E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
!				Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
!				salidas=salidas+E1
!				!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
!				R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
!				Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
!				Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
!				Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
!				!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
!				R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
!				Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm] servicioalcliente@flybox.com.co
!				Es3=E3(celdas)*Sp(3,celdas) ![mm]
!				Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
!				!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
!				R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
!				Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
!				Es4=E4(celdas)*Sp(4,celdas) ![mm]
!				Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
!				!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
!				select case (tipo_celda(1,celdas))
!				    case(1)
!					!Sigue drenando a su tanque igual aguas abajo
!					if (drena(1,celdas).ne.0) then
!					    Sp(2,drenaid)=Sp(2,drenaid)+Es2
!					    Sp(3,drenaid)=Sp(3,drenaid)+Es3
!					    Sp(4,drenaid)=Sp(4,drenaid)+Es4
!					else
!					    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
!					    salidas=salidas+Es2+Es3+Es4
!					endif
!				    case(2,3)			
!					!El tanque 5 recibe agua de los tanques 2 y 3 fijo
!					Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
!					Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
!					!Copia variables de vectores a escalares par aahorrar tiempo
!					vel=v_cauce(celdas) ![m/seg]
!					S5=Sp(5,celdas) ![mm]			
!					L=L_celda(1,celdas) ![mts]
!					pm=pend_man(celdas)
!					!Calcula mediante onda cinematica el flujo a la salida
!					do i=1,4
!					    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
!					    vn=(Area**(2/3))*pm ![m/seg]
!					    vel=(2*vn+vel)/3 ![m2/seg]
!					enddo
!					v_cauce(celdas)=vel		            			
!					Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
!					Sp(5,celdas)=Sp(5,celdas)-Es5
!					if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
!				end select
!				!si es la salida de la cuenca guarda
!				if (drena(1,celdas).eq.0) then 		   
!				    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
!				    salidas=salidas+Es5 ![mm]
!				endif
!				!si es un punto de control guarda igualmente
!				if (control(1,celdas).ne.0) then
!				    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
!				    control_cont=control_cont+1
!				endif
!				!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
!				rain_sum=rain_sum+R1
!		    enddo
!		    !calcula cantidad total de lluvia que entra al sistema
!		    entradas=entradas+rain_sum
!		    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
!		    mean_rain(tiempo)=rain_sum/N_cel
!		enddo
!    else
!		print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
!    endif
!end subroutine

!model_field_lk: Modelo shia completo, contiene: velocidad lineal en ladera 
!y velocidad cinematica en canal, no se usa la OCG. lluvia es leida por 
!archivo, generado por ll.interpolation_idw o ll.interpolation_tin, o radar
!con el mismo formato.
!subroutine shia_field_lk(calib,ruta_rain,N_cel,N_cont,N_reg,Q,balance,mean_rain,sp)
!    !Variables de entrada
!    integer, intent(in) :: N_cel,N_reg,N_cont
!    real, intent(in) :: calib(9)
!    character*255, intent(in) :: ruta_rain
!    !Variables de salia
!    real, intent(out) :: Q(N_cont,N_reg),sp(5,N_cel),mean_rain(N_reg),balance
!    !Variables locales 
!    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
!    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
!    integer drenaid,Res !Vector con el id a donde drena cada celda
!    real Rain(nceldas),rain_sum !campo de lluvia  para el intervalo de tiempo
!    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
!    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
!    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
!    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
!    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
!    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
!    real m3_mm !Variable para convertir m3 a mm2
!    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
!    real X,Y !Posiciones actuales de la celda
!    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
!    !Definiciones para f2py
!    !f2py intent(in) :: N_cel,N_reg,N_cont,calib
!    !f2py intent(out) :: Q,sp,mean_rain,balance
!    !Preambulo a la ejecucion del modelo
!    !Calcula el tiempo inicial
!    call cpu_time(tiempo_i)
!    !Calcula la cantidad de agua que es retenida en el capilar
!    Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
!    mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
!    mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que precola el agua
!    mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm] Velocidad a la que se dan perdidas de agua
!    !Calcula las velocidades horizontales iniciales
!    v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
!    v_cauce=1 ![m/s] Vel ladera, Estado constante
!    pend_man=(pend_celda(1,:)**(0.5))/man(1,:)
!    ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
!    kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
!    !Calcula el porcentaje de flujo horizontal saliente de cada celda
!    E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
!    E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
!    E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
!    E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
!    !Convertor de mm a mt3 para obtener caudal a la salida
!    m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
!    !Establece el almacenamiento inicial 
!    Sp=S
!    !Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
!    inicial=sum(Sp(:,:))
!    entradas=0
!    salidas=0
!    !Determina como lluvia anterior y la lluvia media para cada celda el valor cero
!    mean_rain=0
!    if (eval.eq.1) then
!	!Itera para cada intervalo de tiempo
!	do tiempo=1,N_reg
!	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
!	    rain_sum=0
!	    control_cont=2
!	    !lee la lluvia para el intervalo
!	    call read_float_basin(ruta_rain,tiempo,Rain,res,1,N_cel)   
!	    !Itera para todas las celdas
!	    do celdas=1,N_cel
!		!Verifica la lluvia
!		if (Rain(celdas).lt.0 .or. Rain(celdas).eq.0) then
!		    R1=0.0
!		else
!		    R1=Rain(celdas)
!		endif
!		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
!		drenaid=N_cel-drena(1,celdas)+1	    
!		!Al sistema le entra la lluvia
!		entradas=entradas+R1
!		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
!		R2=max(0.0,Rain(celdas)-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
!		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
!		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
!		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
!		salidas=salidas+E1
!		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
!		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
!		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
!		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
!		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
!		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
!		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
!		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm] 
!		Es3=E3(celdas)*Sp(3,celdas) ![mm]
!		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
!		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
!		R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
!		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
!		Es4=E4(celdas)*Sp(4,celdas) ![mm]
!		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
!		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
!		select case (tipo_celda(1,celdas))
!		    case(1)
!			!Sigue drenando a su tanque igual aguas abajo
!			if (drena(1,celdas).ne.0) then
!			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
!			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
!			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
!			else
!			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
!			    salidas=salidas+Es2+Es3+Es4
!			endif
!		    case(2,3)			
!			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
!			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
!			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
!			!Copia variables de vectores a escalares par aahorrar tiempo
!			vel=v_cauce(celdas) ![m/seg]
!			S5=Sp(5,celdas) ![mm]			
!			L=L_celda(1,celdas) ![mts]
!			pm=pend_man(celdas)
!			!Calcula mediante onda cinematica el flujo a la salida
!			do i=1,4
!			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
!			    vn=(Area**(2/3))*pm ![m/seg]
!			    vel=(2*vn+vel)/3 ![m2/seg]
!			enddo
!			v_cauce(celdas)=vel		            			
!			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
!			Sp(5,celdas)=Sp(5,celdas)-Es5
!			if (drena(1,celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
!		end select
!		!si es la salida de la cuenca guarda
!		if (drena(1,celdas).eq.0) then 		   
!		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
!		    salidas=salidas+Es5 ![mm]
!		endif
!		!si es un punto de control guarda igualmente
!		if (control(1,celdas).ne.0) then
!		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
!		    control_cont=control_cont+1
!		endif
!		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
!		rain_sum=rain_sum+R1
!	    enddo
!	    !calcula cantidad total de lluvia que entra al sistema
!	    entradas=entradas+rain_sum
!	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
!	    mean_rain(tiempo)=rain_sum/N_cel
!	enddo
!    else
!	print *, 'Aviso: El modelo no ha ejecutado faltan variables por asignar dimension!!!'
!    endif
!end subroutine


end module

!-----------------------------------------------------------------------
!Subrutinas extras del moduo que son usadas internamente por modelos
!-----------------------------------------------------------------------

