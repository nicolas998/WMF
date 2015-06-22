module modelos
implicit none

!-----------------------------------------------------------------------
!Variableses globales para la modelacion
!-----------------------------------------------------------------------
!El tiempo
real dt
!Magnitud espacial de las celdas
real dxP
!Variable constante de la ecuacion OCG
real, allocatable :: part_ocg(:,:)
!Variables enteras provenientes del terreno
integer, allocatable :: tipo_celda(:,:),acum_cel(:,:)
!Variables de mapas, variables correspondientes a los mapas de usos de suelo, tipos, etc.
real, allocatable :: Hu(:,:),H3(:,:),Man(:,:),Ks(:,:),Kp(:,:),EVP(:,:)
!conversores de las variables ks y kp
real conver_ks,conver_kp
!Variables propias del terreno 
real, allocatable :: L_celda(:,:),pend_celda(:,:),pend_cauce(:,:),Elev(:,:),XY_celdas(:,:)
!almacenamiento del modelo
real, allocatable :: S(:,:)
!Valor de los datos nodata
real noData
!Coordenadas de la lluvia y registros de lluvia
real, allocatable :: coord_rain(:,:), rain(:,:)
integer N_reg,N_est
!Cantidad de celdas en la cuenca
integer nceldas
!indicadores de chequeo para permitir operar el modelo
integer :: check_all=1,check_topology=1,check_terrain=1,check_map=1,check_store=1,check_rain=1

!-----------------------------------------------------------------------
!Punto de inicio de funciones del modulo
contains

!-----------------------------------------------------------------------
!Funciones Varias
!-----------------------------------------------------------------------
!funciones para lectura de variables para shia
subroutine shia_read_map(ruta,records,nceldas) !Lee las variables de mapas externos
    !Variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(in) :: nceldas
    integer, intent(in) :: records(5)
    !f2py intent(in) :: ruta,nceldas,nreg,registros
    !Inicia las variables de los mapas
    allocate(Hu(1,nceldas),Man(1,nceldas),Ks(1,nceldas),kp(1,nceldas),EVP(1,nceldas))
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	read(10,rec=records(1)) Hu(1,:)
	read(10,rec=records(2)) man(1,:)
	read(10,rec=records(3)) ks(1,:)
	read(10,rec=records(4)) kp(1,:)
	read(10,rec=records(5)) EVP(1,:)
    close(10)
end subroutine
subroutine shia_read_terrain(ruta,records,nceldas) !Lee las variables obtenidas del terreno
    !Variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(in) :: nceldas
    integer, intent(in) :: records(5)
    !f2py intent(in) :: ruta,nceldas,nreg,registros
    !Inicia las variables de los mapas
    allocate(L_celda(1,nceldas),pend_cauce(1,nceldas),pend_celda(1,nceldas),Elev(1,nceldas),XY_celdas(2,nceldas))
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	read(10,rec=records(1)) L_celda(1,:)
	read(10,rec=records(2)) pend_cauce(1,:)
	read(10,rec=records(3)) pend_celda(1,:)
	read(10,rec=records(4)) Elev(1,:)
	read(10,rec=records(5)) XY_celdas(1,:)
    close(10)
end subroutine
subroutine shia_read_topology(ruta,records,nceldas) !Lee las variables obtenidas de la topologia del terreno
    !Variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(in) :: nceldas
    integer, intent(in) :: records(2)
    !f2py intent(in) :: ruta,nceldas,nreg,registros
    !Inicia las variables de los mapas
    allocate(tipo_celda(1,nceldas),acum_cel(1,nceldas))
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	read(10,rec=records(1)) tipo_celda(1,:)
	read(10,rec=records(2)) acum_cel(1,:)
    close(10)
end subroutine
subroutine shia_read_store(ruta,records,nceldas) !Lee las variables de almacenamiento 
    !Variables de entrada
    character*255, intent(in) :: ruta
    integer, intent(in) :: nceldas
    integer, intent(in) :: records(5)
    !f2py intent(in) :: ruta,nceldas,nreg,registros
    !Inicia las variables de los mapas
    allocate(S(5,nceldas))
    !Lectura 
    open(10,file=ruta,form='unformatted',status='old',access='direct',RECL=4*nceldas)
	read(10,rec=records(1)) S(1,:)
	read(10,rec=records(2)) S(2,:)
	read(10,rec=records(3)) S(3,:)
	read(10,rec=records(4)) S(4,:)
	read(10,rec=records(5)) S(5,:)
    close(10)
end subroutine
!funciones para desalojar variables de shia
subroutine shia_free_map !libera las variables de mapas
    deallocate(Hu,Man,Ks,kp,EVP)
end subroutine
subroutine shia_free_terrain !libera las variables de mapas
    deallocate(L_celda,pend_cauce,pend_celda,Elev,XY_celdas)
end subroutine
subroutine shia_free_topology !libera las variables de mapas
    deallocate(tipo_celda,acum_cel)
end subroutine
subroutine shia_free_store !libera las variables de mapas
    deallocate(S)
end subroutine
subroutine shia_free_all !Libera todas
    deallocate(Hu,Man,Ks,kp,EVP)
    deallocate(L_celda,pend_cauce,pend_celda,Elev)
    deallocate(tipo_celda,acum_cel)
    deallocate(S)
end subroutine
!funciones para alojar variables sin leerlas
subroutine shia_set_map
    !Inicia las variables de los mapas
    allocate(Hu(1,nceldas),Man(1,nceldas),Ks(1,nceldas),kp(1,nceldas),EVP(1,nceldas))
end subroutine
subroutine shia_set_terrain
    !Inicia las variables de los mapas
    allocate(L_celda(1,nceldas),pend_cauce(1,nceldas),pend_celda(1,nceldas),Elev(1,nceldas),XY_celdas(2,nceldas))
end subroutine
subroutine shia_set_topology
    !Inicia las variables de los mapas
    allocate(tipo_celda(1,nceldas),acum_cel(1,nceldas))
end subroutine
subroutine shia_set_store
    !Inicia las variables de los mapas
    allocate(S(5,nceldas))
end subroutine
!funciones para leer variables de lluvia
subroutine rain_dimension_point
    !Dimensiona
    allocate(coord_rain(2,N_est),rain(N_est,N_reg))
end subroutine
!Funciones para estimar variables del modelo
subroutine OCG_params(params,exp_param,calculo_ocg)
    !Variables de entrada
    real, intent(in) :: params(9) ! c1, Omega, k, alfa1, alfa2, sigma1, sigma2, sigma3, fhi
    real, intent(out) :: exp_param(4) !los cuatro parametros para la onda, e_cau,e_pend,e_acum,B
    integer, intent(out) :: calculo_ocg !0 si si la calculo, 1 si no
    !f2py intent(in) :: params
    !f2py intent(out) :: exp_param,calculo_ocg 
    !Internos 
    real c1, Omega, k, alfa1, alfa2, sigma1, sigma2, sigma3, fhi 
    real e_cau,e_pend,e_acum,B,e_B
    !Copia params
    sigma1=params(1); sigma2=params(2); sigma3=params(3)
    Omega=params(4); alfa1=params(5); alfa2=params(6)
    k=params(7); fhi=params(8); c1=params(9)
    !Calcula los exponentes y el coeficiente de la OGC
    B=Omega*(c1*k**(alfa1-alfa2))**((2/3)-sigma2)
    e_B=1/(1+alfa2*(0.667-sigma2))
    e_cau=(0.667-sigma2)*(1-alfa2)*e_B
    e_pend=(0.5-sigma3)*e_B
    e_acum=(fhi*(0.667-sigma2)*(alfa2-alfa1)+sigma1)*e_B
    e_B=-1*e_B
    B=B**e_B
    exp_param(1)=e_cau; exp_param(2)=e_pend; exp_param(3)=e_acum; exp_param(4)=B 
    !si la variable de area acumulada y pendiente se encuentran alojadas calcula vn=(Area**OCG(1)) --> *(So**OCG(2))*(Ac**OCG(3))*OCG(4) ![m/seg]
    if (allocated(pend_cauce).and.allocated(pend_celda).and.allocated(acum_cel).and.nceldas.gt.0) then
	!se fija si la variable ya ha sido alojada o no
	if (allocated(part_ocg) .eqv. .false.) then
	    allocate(part_ocg(1,nceldas))
	endif
	part_ocg(1,:)=noData
	where(tipo_celda(1,:).gt.1) part_ocg(1,:)=(pend_celda(1,:)**e_pend)*(((acum_cel(1,:)*(dxP**2))/1e6)**e_acum)*B ![adim]
	!where(tipo_celda(1,:).eq.3) part_ocg(1,:)=(pend_cauce(1,:)**e_pend)*(((acum_cel(1,:)*(dxP**2))/1e6)**e_acum)*B ![adim]
	calculo_ocg=0
    else
	calculo_ocg=1
    endif
end subroutine
!Funciones para determinar si el modelo ya se puede poner a funcionar
subroutine shia_idw_check !Determina si se puede correr shia
    check_all=0
    !Chequea variables de lluvia
    if (allocated(rain).and.allocated(coord_rain)) then
	if (N_est.gt.0 .and. N_reg.gt.0) then
	    check_rain=0
	endif
    else
	check_all=check_all+check_rain
    endif
    !Chequea variables de mapa
    if (allocated(Hu).and.allocated(Man).and.allocated(Ks).and.allocated(Kp).and.allocated(EVP)) then
	check_map=0
    else
	check_all=check_all+check_map
    endif
    !chequeva variables de la topologia
    if (allocated(acum_cel).and.allocated(tipo_celda)) then
	check_topology=0
    else
	check_all=check_all+check_topology
    endif
    !Chequea variables del terreno
    if (allocated(L_celda).and.allocated(pend_celda).and.allocated(pend_cauce).and.allocated(Elev).and.allocated(XY_celdas)) then
	check_terrain=0
    else
	check_all=check_all+check_terrain
    endif
    !chequea almacenamiento
    if (allocated(S)) then
	check_store=0
    else
	check_all=check_all+check_store
    endif
end subroutine 

!-----------------------------------------------------------------------
!Versiones modelo shia
!-----------------------------------------------------------------------
!Codigos:
subroutine shia_idw_lllk(drena,control,Calib,a_ocg,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) !Modelo con los tanques 2,3 y 4 estaticos, y el tanque 5 cinematico
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(9),a_ocg,pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(5,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_ejec(2),tiempo_result !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !Velocidades del flujo horizontal [m/s]
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Tiempos de recidencia verticales [seg]
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real E1,E2(nceldas),E3(nceldas),E4(nceldas) !Porcentaje del flujo horizontal que fluye [%] 
    real Es2,Es3,Es4,E5 !Cantidad de flujo que finalmente sale de cada tanque en cada intervalo de tiempo  
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real Area,vel,vn,S5,OCG !Area acum,velocidad, almacenamiento tanque 5 y OCG; todos son temporales (copias de vector)
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call etime(tiempo_ejec,tiempo_result)
	!Calcula el valor modificado del almacenamiento maximo capilar
	Hu_loc=Hu(1,:)*Calib(1)
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Estado constante
	v_cauce=1.0!1.4*Calib(9)*(pend_cauce(1,:)**0.5)/man(1,:) ![m/s] Estado inicial, recalculada de manera dinamica
	ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Estado constante
	kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Estado constante
	!calcula los tiempos de residencia en los tanques 2 3 y 4 
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm]
	mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm]
	mm_kpp=kp(1,:)*Calib(7)*dt*conver_kp ![mm]
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
	!Calcula conversores de unidades
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Copia el almacenamiento inicial 
	Sp=S
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		!Tanque 2: Calculo de la cantidad de agua que escurre por ladera y de la cantidad que se infiltra
		R3=min(R2,mm_ks(celdas)) ![mm]
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] 
		Es2=E2(celdas)*Sp(2,celdas)
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm]               			    
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial y la que se precola		
		R4=min(R3,mm_kp(celdas)) ![mm]
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm]
		Es3=E3(celdas)*Sp(3,celdas)
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]		
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo subterraneo y la que pasa a perdidas		
		R5=min(R4,mm_kpp(celdas)) ![mm]
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]                              
		Es4=E4(celdas)*Sp(4,celdas)
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]		   
		!Distribucion del Agua y flujo en carcavas y canales			
		select case(tipo_celda(1,celdas))
		    case(1)
			!Cada tanque drena a su igual aguas abajo
			Sp(2,drenaid)=Sp(2,drenaid)+Es2
			Sp(3,drenaid)=Sp(3,drenaid)+Es3
			Sp(4,drenaid)=Sp(4,drenaid)+Es4
		    case(2,3)
			!El tanque 5 recibe agua de los tanques 2 y 3 fijo
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4*(tipo_celda(1,celdas)-2) ![mm]
			Sp(4,drenaid)=Sp(4,drenaid)+Es4*(3-tipo_celda(1,celdas)) ![mm]
			!Si hay agua en el tanque calcula si no no
			if (Sp(5,celdas).gt.0)	then		
			    !Copia variables de vectores a escalares par aahorrar tiempo
			    vel=v_cauce(celdas) ![m/seg]
			    S5=Sp(5,celdas) ![mm]			
			    OCG=part_ocg(1,celdas) ![adim]
			    L=L_celda(1,celdas) ![mts]
			    !Calcula mediante OCG el flujo a la salida
			    do i=1,4
				Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
				vn=(Area**a_ocg)*OCG ![m/seg]
				vel=(2*vn+vel)/3 ![m2/seg]
			    enddo		            
			    !Devuelve valores a los vectores
			    v_cauce(celdas)=vel ![m/seg]	    
			    !Calcula la salida y actualiza tanques
			    E5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas)) ![mm]			    			    
			    Sp(5,celdas)=Sp(5,celdas)-E5 ![mm]			
			    if (drena(celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+E5	![mm]
			endif			
		end select	    				
		!Si es la salida de la cuenca guarda en la primera entrada del caudal
		if (drena(celdas).eq.0) Q(1,tiempo)=E5*m3_mm/dt ![m3/s]      
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=E5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call etime(tiempo_ejec,tiempo_result)
	print *, 'Tiempo ejecucion: ', tiempo_result
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine

!-----------------------------------------------------------------------
!Versiones modelos sencillos de transporte interpolados con IDW
!-----------------------------------------------------------------------
!model_hlikp: modelo con velocidad constante en la ladera, con infiltración y precolacion, con retencion capilar, velocidad cinematica en el canal.
subroutine model_idw_hlikp(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(9),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(5,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
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
	kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
	E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
	print *, E4(nceldas)	
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=S
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
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
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante OCG el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5
			if (drena(celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select
		!si es la salida de la cuenca guarda
		if (drena(celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hlikp: modelo con velocidad constante en la ladera, con infiltración y precolacion, con retencion capilar, velocidad cinematica en el canal.
subroutine model_idw_hlikp_campo(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(9),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(5,nceldas),mean_rain(nceldas),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4,R5 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S5 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas),mm_kpp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),E5(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
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
	kph=kp(1,:)*Calib(8)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel acuifero, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:))
	E5=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=S
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+peso*rain(pert_rain(i,celdas),tiempo)
			W=W+peso			
		    else
			cont=cont-1
		    endif		
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial
		R4=min(R3,mm_kp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		R5=min(R4,mm_kpp(celdas)) !precola el minimo entre lo que hay y la conductividad del acuifero		
		Sp(4,celdas)=Sp(4,celdas)+R4-R5 ![mm]
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
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
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante OCG el flujo a la salida
			do i=1,4
			    Area=Sp(5,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(9)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5
			if (drena(celdas).ne.0) Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]			
		end select
		!si es la salida de la cuenca guarda
		if (drena(celdas).eq.0) then 		   
		    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es5 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
		mean_rain(celdas)=R1
	    enddo
	    !calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    !mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine


!model_l: modelo con velocidad constante en la ladera y constante para toda la cuenca, sin infiltración, sin retencion capilar.
subroutine model_idw_lc(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib,pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real R1,R2 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas) !velocidad en ladera
    real E2(nceldas),Es2 !Cantidad de flujo que viaja en cada celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial	
	call cpu_time(tiempo_i)
	!Calcula las velocidades horizontales iniciales
	v_ladera=1*Calib ![m/s] Estado constante
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)		    
		endif 
		!Tanque 1: todo lo que cae fluye eventualmente		
		Sp(celdas)=Sp(celdas)+R1 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(celdas) ! [mm] De acuerdo al almacenamiento calcula cuanto se va	
		Sp(celdas)=Sp(celdas)-Es2 ![mm] Actualiza almacenamiento              			    
		!Si es la salida de la cuenca guarda en la primera entrada del caudal
		if (drena(celdas).eq.0) then
		    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s] 
		    salidas=salidas+Es2 ![mm] Lo que sale del sistema en mm     
		else
		    Sp(drenaid)=Sp(drenaid)+Es2
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance de agua 
	balance=100*((entradas-salidas-sum(Sp(:)))/entradas)	
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_l: modelo con velocidad constante en la ladera, sin infiltración, sin retencion capilar.
subroutine model_idw_l(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib,pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real R1,R2 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas) !velocidad en ladera
    real E2(nceldas),Es2 !Cantidad de flujo que viaja en cada celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Estado constante
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: todo lo que cae fluye eventualmente
		Sp(celdas)=Sp(celdas)+R1 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(celdas) ! [mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(celdas)=Sp(celdas)-Es2 ![mm] Actualiza almacenamiento              			    
		!Si es la salida de la cuenca guarda en la primera entrada del caudal
		if (drena(celdas).eq.0) then
		    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]      
		else
		    Sp(drenaid)=Sp(drenaid)+Es2
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula la cantidad total de lluvia sobre el sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance de agua en la cuenca
	balance=100*((entradas-salidas-sum(Sp(:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hl: modelo con velocidad constante en la ladera, sin infiltración, con retencion capilar.
subroutine model_idw_hl(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(3),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(2,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas) !velocidad en ladera
    real E2(nceldas),Es2,E1 !Cantidad de flujo que viaja en cada celda, Cantidad de flujo evaporado de la celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib(3)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Estado constante
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		Sp(2,celdas)=Sp(2,celdas)+R2 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento              			    
		!Si es la salida de la cuenca guarda en la primera entrada del caudal
		if (drena(celdas).eq.0) then
		    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]      
		else
		    Sp(2,drenaid)=Sp(2,drenaid)+Es2
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calculo de la cantidad de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance de agua en la cuenca
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hli: modelo con velocidad constante en la ladera, con infiltración, con retencion capilar.
subroutine model_idw_hli(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(5),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(3,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),ksh(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),Es2,Es3,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Sp(3,celdas)=Sp(3,celdas)+R3 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			endif
		    case(2,3)
			!Drena de nuevo al flujo superficial
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2+Es3
			else
			    Q(1,tiempo)=(Es2+Es3)*m3_mm/dt ![m3/s]
			endif
		end select
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula la cantidad total de lluvia que entra
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance en la cuenca
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hlic: modelo con velocidad constante en la ladera, con infiltración, con retencion capilar, velocidad constante en el canal.
subroutine model_idw_hlic(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(6),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(4,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),Es2,Es3,Es4,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
	v_cauce=1.4*Calib(6)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	E4=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
	print *, E4(nceldas)	
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Sp(3,celdas)=Sp(3,celdas)+R3 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2
			endif
		    case(2,3)
			!Drena al tanque del cauce de la misma celda
			Sp(4,celdas)=Sp(4,celdas)+Es2+Es3
			Es4=Sp(4,celdas)*E4(celdas)
			Sp(4,celdas)=Sp(4,celdas)-Es4
			if (drena(celdas).ne.0) Sp(4,drenaid)=Sp(4,drenaid)+Es4	![mm]			
		end select
		!si es la salida de la cuenca guarda
		if (drena(celdas).eq.0) then 		   
		    Q(1,tiempo)=Es4*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es4 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es4*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calculo de la cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hlik: modelo con velocidad constante en la ladera, con infiltración, con retencion capilar, velocidad cinematica en el canal.
subroutine model_idw_hlik(drena,control,Calib,pert_rain,pot_rain, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(6),pot_rain !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(4,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas), pm, L, Area, vel, vn, S4 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),Es2,Es3,Es4,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1.4*Calib(4)*(pend_celda(1,:)**0.5)/man(1,:) ![m/s] Vel ladera, Estado constante
	v_cauce=1 ![m/s] Vel ladera, Estado constante
	pend_man=(pend_celda(1,:)**(0.5))/man(1,:)
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	E4=1-L_celda(1,:)/(v_cauce*dt+L_celda(1,:))
	print *, E4(nceldas)	
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		Es2=E2(celdas)*Sp(2,celdas) ![mm] De acuerdo al almacenamiento calcula cuanto se va
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Sp(3,celdas)=Sp(3,celdas)+R3 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2
			endif
		    case(2,3)
			!Drena al tanque del cauce de la misma celda
			Sp(4,celdas)=Sp(4,celdas)+Es2+Es3
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S4=Sp(4,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante OCG el flujo a la salida
			do i=1,4
			    Area=Sp(4,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es4=min(Area*v_cauce(celdas)*dt*Calib(6)/m3_mm,Sp(4,celdas))
			Sp(4,celdas)=Sp(4,celdas)-Es4
			if (drena(celdas).ne.0) Sp(4,drenaid)=Sp(4,drenaid)+Es4	![mm]			
		end select
		!si es la salida de la cuenca guarda
		if (drena(celdas).eq.0) then 		   
		    Q(1,tiempo)=Es4*m3_mm/dt ![m3/s]      
		    salidas=salidas+Es4 ![mm]
		endif
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es4*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula cantidad total de lluvia que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hkik: modelo con velocidad cinematica en la ladera, con infiltración, con retencion capilar, velocidad cinematica en el canal.
subroutine model_idw_hkik(drena,control,Calib,pert_rain,pot_rain,Q,balance,mean_rain,Sp,corrio,npert,nceldas,ncontrol,nreg)
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(8),pot_rain!surco_param(2) !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(4,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas),surcos(nceldas), pm, L, Area, vel, vn, S4, S2 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),Es2,Es3,Es4,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    real surco_param
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1 ![m/s] Vel ladera, Estado constante
	v_cauce=1 ![m/s] Vel ladera, Estado constante
	pend_man=(pend_celda(1,:)**(0.5))/man(1,:) !Constante por celda para calcular velocidad cinematica en cauce
	surcos=(pend_celda(1,:)**0.5)*Calib(7)/man(1,:) !constante por celda para calcular velocidad cinematica en ladera
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**pot_rain
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 		
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		!Copia variables de vectores a escalares par aahorrar tiempo
		vel=v_ladera(celdas) ![m/seg]
		S2=Sp(2,celdas) ![mm]			
		L=L_celda(1,celdas) ![mts]
		pm=surcos(celdas)
		!Calcula mediante OCG el flujo a la salida
		do i=1,4
		    Area=Sp(2,celdas)*m3_mm/(L+vel*dt) ![m2]
		    vn=(Area**((2/3)*Calib(8)))*pm ![m/seg]
		    vel=(2*vn+vel)/3 ![m2/seg]
		enddo
		v_ladera(celdas)=vel		            			
		Es2=min(Area*v_ladera(celdas)*dt*Calib(4)/m3_mm,Sp(2,celdas)) ! [mm]		
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Sp(3,celdas)=Sp(3,celdas)+R3 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2
			endif
		    case(2,3)
			!Drena al tanque del cauce de la misma celda
			Sp(4,celdas)=Sp(4,celdas)+Es2+Es3
			!Copia variables de vectores a escalares par aahorrar tiempo
			vel=v_cauce(celdas) ![m/seg]
			S4=Sp(4,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante OCG el flujo a la salida
			do i=1,4
			    Area=Sp(4,celdas)*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es4=min(Area*v_cauce(celdas)*dt*Calib(6)/m3_mm,Sp(4,celdas))
			Sp(4,celdas)=Sp(4,celdas)-Es4
			if (drena(celdas).ne.0) then 
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4	![mm]
			else
			    Q(1,tiempo)=Es4*m3_mm/dt ![m3/s]      
			    salidas=salidas+Es4 ![mm]
			endif
		end select
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es4*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calcula la cantidad de lluvia total que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hkika: modelo con velocidad cinematica en la ladera, con infiltración, con retencion capilar, velocidad cinematica en el canal y almacenamiento en el acuifero.
subroutine model_idw_hkika(drena,control,Calib,pert_rain,params,Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(8),params(3) !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(5,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R4 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),v_cauce(nceldas),ksh(nceldas),kph(nceldas) !velocidad en ladera, Velocidad del flujo sub-superficial horizontal
    real pend_man(nceldas),surcos(nceldas), pm, L, Area, vel, vn, S5, S2 !combinado pendiente manning para calcular, Long, Area y velocidad para calcular flujo cinematico 
    real mm_ks(nceldas),mm_kp(nceldas) !Velocidad de la conductividad vertical
    real E2(nceldas),E3(nceldas),E4(nceldas),Es2,Es3,Es4,Es5,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm !Variable para convertir m3 a mm2
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	mm_kp=kp(1,:)*Calib(5)*dt*conver_kp ![mm] Velocidad a la que se precola el flujo
	!Calcula las velocidades horizontales iniciales
	v_ladera=1 ![m/s] Vel ladera, Estado constante
	v_cauce=1 ![m/s] Vel ladera, Estado constante
	pend_man=(pend_celda(1,:)**(0.5))/man(1,:) !Constante por celda para calcular velocidad cinematica en cauce
	surcos=(pend_celda(1,:)**0.5)*params(2)/man(1,:) !constante por celda para calcular velocidad cinematica en ladera
	ksh=ks(1,:)*Calib(6)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	kph=kp(1,:)*Calib(7)*(conver_kp/1000)*pend_celda(1,:) ![m/s] Vel flujo base, Estado constante
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:)) !Porcentaje constante de flujo sub-superficial saliente
	E4=1-L_celda(1,:)/(kph*dt+L_celda(1,:)) !Porcentaje constante de flujo subterraneo saliente	
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**params(2)
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 		
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		salidas=salidas+E1
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		!Copia variables de vectores a escalares par aahorrar tiempo
		vel=v_ladera(celdas) ![m/seg]
		S2=Sp(2,celdas) ![mm]			
		L=L_celda(1,celdas) ![mts]
		pm=surcos(celdas)
		!Calcula mediante OCG el flujo a la salida
		do i=1,4
		    Area=Sp(2,celdas)*m3_mm/(L+vel*dt) ![m2]
		    vn=(Area**((2/3)*params(3)))*pm ![m/seg]
		    vel=(2*vn+vel)/3 ![m2/seg]
		enddo
		v_ladera(celdas)=vel		            			
		Es2=min(Area*v_ladera(celdas)*dt*Calib(4)/m3_mm,Sp(2,celdas)) ! [mm]		
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		R4=min(R3,mm_kp(celdas)) !Precola el minimo entre la cantidad de flujo y la conductividad
		Sp(3,celdas)=Sp(3,celdas)+R3-R4 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Tanque 4: Calculo de la cantidad de agua que escurre como flujo subterraneo
		Sp(4,nceldas)=Sp(4,nceldas)+R4
		Es4=E4(celdas)*Sp(4,celdas) ![mm]
		Sp(4,celdas)=Sp(4,celdas)-Es4 ![mm]	     			    		
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			else
			    Q(1,tiempo)=(Es2+Es3+Es4)*m3_mm/dt ![m3/s]
			    salidas=salidas+Es2+Es3+Es4 ![mm]
			endif
		    case(2,3)
			!Drena al tanque del cauce de la misma celda de los tanques 2,3 y 4 dependiendo del tipo de celda
			if (tipo_celda(1,celdas).eq.2) then
			    Sp(5,celdas)=Sp(5,celdas)+Es2+Es3
			    Sp(4,drenaid)=Sp(4,drenaid)+Es4
			elseif (tipo_celda(1,celdas).eq.3) then
			    Sp(5,celdas)=Sp(5,celdas)+Es2+Es3+Es4
			endif
			!Copia variables de vectores a escalares par ahorrar tiempo
			Sp(5,celdas)=Sp(5,celdas)+Es2+Es3
			vel=v_cauce(celdas) ![m/seg]
			S5=Sp(5,celdas) ![mm]			
			L=L_celda(1,celdas) ![mts]
			pm=pend_man(celdas)
			!Calcula mediante OCG el flujo a la salida
			do i=1,4
			    Area=S5*m3_mm/(L+vel*dt) ![m2]
			    vn=(Area**(2/3))*pm ![m/seg]
			    vel=(2*vn+vel)/3 ![m2/seg]
			enddo
			v_cauce(celdas)=vel		            			
			Es5=min(Area*v_cauce(celdas)*dt*Calib(8)/m3_mm,Sp(5,celdas))
			Sp(5,celdas)=Sp(5,celdas)-Es5
			if (drena(celdas).ne.0) then 
			    Sp(5,drenaid)=Sp(5,drenaid)+Es5	![mm]
			else
			    Q(1,tiempo)=Es5*m3_mm/dt ![m3/s]      
			    salidas=salidas+Es5 ![mm]
			endif
		end select
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es5*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !Calcula la cantidad de lluvia total que entra al sistema
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hkir: modelo con velocidad cinematica en la ladera, con infiltración, con retencion capilar y retorno sub-superficial.
subroutine model_idw_hkir(drena,control,Calib,pert_rain,params, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(6),params(3) !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(3,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas),H3_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3,R3r !lluvia (R1) y flujo que drena por los tanques, y R3r retorno del tanque 3 al 2 por exceso
    real v_ladera(nceldas),surcos(nceldas),ksh(nceldas),vn,vel,Area,pm,S2 !velocidad en ladera, parametro de surcos, Velocidad del flujo sub-superficial horizontal
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E3(nceldas),Es2,Es3,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	H3_loc=H3(1,:)*Calib(6) !Cantidad de agua almacenada en el sub-suelo antes de volver a la celda superficial
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1 ![m/s] Vel ladera, Estado constante
	surcos=(pend_celda(1,:)**0.5)*params(2)/man(1,:) !constante por celda para calcular velocidad cinematica en ladera
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**params(3)
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]		
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]		
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		R3r=max(0.0,R3-H3_loc(celdas)+Sp(3,celdas))
		!if (celdas==100000) print *, R3,Sp(3,celdas),H3_loc(celdas),R3-H3_loc(celdas)+Sp(3,celdas)
		Sp(2,celdas)=Sp(2,celdas)+R3r
		Sp(3,celdas)=Sp(3,celdas)+R3-R3r ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Calcula velocidad superficial el ladera
		vel=v_ladera(celdas) ![m/seg]
		S2=Sp(2,celdas) ![mm]			
		L=L_celda(1,celdas) ![mts]
		pm=surcos(celdas)
		!Calcula mediante OCG el flujo a la salida
		do i=1,4
		    Area=Sp(2,celdas)*m3_mm/(L+vel*dt) ![m2]
		    vn=(Area**((2/3)*params(3)))*pm ![m/seg]
		    vel=(2*vn+vel)/3 ![m2/seg]
		enddo
		v_ladera(celdas)=vel		            			
		Es2=min(Area*v_ladera(celdas)*dt*Calib(4)/m3_mm,Sp(2,celdas)) ! [mm]		
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			endif
		    case(2,3)
			!Drena de nuevo al flujo superficial
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2+Es3
			else
			    Q(1,tiempo)=(Es2+Es3)*m3_mm/dt ![m3/s]
			endif
		end select
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula la cantidad total de lluvia que entra
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance en la cuenca
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine
!model_hki: modelo con velocidad cinematica en la ladera, con infiltración, con retencion capilar.
subroutine model_idw_hki(drena,control,Calib,pert_rain,params, Q,balance,mean_rain,Sp,corrio, npert,nceldas,ncontrol,nreg) 
    !variables de entrada
    integer, intent(in) :: nceldas,npert,ncontrol,nreg
    integer, intent(in) :: drena(nceldas),control(nceldas),pert_rain(npert,nceldas)
    real, intent(in) :: Calib(6),params(3) !La entrada 10 de Calib es el exponente de la seccion del area del flujo en la OCG
    !Variables de salida
    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(3,nceldas),mean_rain(nreg),corrio
    !f2py intent(in) :: nceldas,npert,ncontrol,drena,tipo_celda,control,pert_rain,Calib,pot_rain
    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
    !f2py note(<Esto es un ensayo>)
    !Variables locales
    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
    integer drenaid !Vector con el id a donde drena cada celda
    real Rain_past(nceldas),peso,W,Wr,rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
    real Hu_loc(nceldas) !Hu Local es una copia de Hu multiplicada por la claibracion [mm]
    real R1,R2,R3 !lluvia (R1) y flujo que drena por los tanques
    real v_ladera(nceldas),surcos(nceldas),ksh(nceldas),vn,vel,Area,pm,S2 !velocidad en ladera, parametro de surcos, Velocidad del flujo sub-superficial horizontal
    real mm_ks(nceldas) !Velocidad de la conductividad vertical
    real E3(nceldas),Es2,Es3,E1 !Cantidad de flujo que viaja en cada celda, cantidad de flujo que sale sub-superficial, Cantidad de flujo evaporado de la celda
    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
    real X,Y !Posiciones actuales de la celda
    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
    if (check_all.eq.0) then
	!Calcula el tiempo inicial
	call cpu_time(tiempo_i)
	!Calcula la cantidad de agua que es retenida en el capilar
	Hu_loc=Hu(1,:)*Calib(1) !Cantidad de agua almacenada en el suelo antes de escurrir
	mm_ks=ks(1,:)*Calib(3)*dt*conver_ks ![mm] Velocidad a la que infiltra el agua
	!Calcula las velocidades horizontales iniciales
	v_ladera=1 ![m/s] Vel ladera, Estado constante
	surcos=(pend_celda(1,:)**0.5)*params(2)/man(1,:) !constante por celda para calcular velocidad cinematica en ladera
	ksh=ks(1,:)*Calib(5)*(conver_ks/1000)*pend_celda(1,:) ![m/s] Vel sub-superficial, Estado constante		
	!Calcula el porcentaje de flujo horizontal saliente de cada celda
	E3=1-L_celda(1,:)/(ksh*dt+L_celda(1,:))
	!Convertor de mm a mt3 para obtener caudal a la salida
	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
	!Establece el almacenamiento inicial 
	Sp=0
	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
	inicial=sum(Sp(:,:))
	entradas=0
	salidas=0
	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
	Rain_past=0
	mean_rain=0
	!Itera para cada intervalo de tiempo
	do tiempo=1,nreg
	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
	    rain_sum=0
	    control_cont=2
	    !Itera para todas las celdas
	    do celdas=1,nceldas
		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
		drenaid=nceldas-drena(celdas)+1
		X=XY_celdas(1,celdas)
		Y=XY_celdas(2,celdas)
		!Itera por todas calculando el peso total y el peso total ponderado
		Wr=0; W=0
		cont=npert
		do i=1,npert		    
		    if (rain(pert_rain(i,celdas),tiempo).ne.noData) then
			peso=1/(sqrt((X-coord_rain(1,pert_rain(i,celdas)))**2+(Y-coord_rain(1,pert_rain(i,celdas)))**2))**params(3)
			Wr=Wr+1*rain(pert_rain(i,celdas),tiempo)
			W=W+1
		    else
			cont=cont-1
		    endif
		enddo
		!calcula la lluvia sobre la celda
		if (cont.gt.0) then 
		    !si una o mas estaciones tienen dato calcula
		    R1=Wr/W
		    Rain_past(celdas)=R1	
		else
		    !si todas las estaciones no reportan dato toma el valor anterior para la celda
		    R1=Rain_past(celdas)
		endif 
		!Al sistema le entra la lluvia
		entradas=entradas+R1
		!Tanque 1: Calculo de la cantidad de agua que es evaporada y la sobrante que pasa al flujo superficial
		R2=max(0.0,R1-Hu_loc(celdas)+Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)+R1-R2 ![mm] 
		E1=min(Evp(1,celdas)*Calib(2)*(Sp(1,celdas)/Hu_loc(celdas))**0.6,Sp(1,celdas)) ![mm]
		Sp(1,celdas)=Sp(1,celdas)-E1 ![mm]
		!Tanque 2: Flujo superficial, lo que no queda atrapado en capilar fluye 
		R3=min(R2,mm_ks(celdas)) !Infiltra el minimo entre lo que hay y la conductividad 		
		Sp(2,celdas)=Sp(2,celdas)+R2-R3 ![mm] !Actualiza alm con la lluvia nueva
		!Copia variables de vectores a escalares par ahorrar tiempo
		vel=v_ladera(celdas) ![m/seg]
		S2=Sp(2,celdas) ![mm]			
		L=L_celda(1,celdas) ![mts]
		pm=surcos(celdas)
		!Calcula mediante OCG el flujo a la salida
		do i=1,4
		    Area=Sp(2,celdas)*m3_mm/(L+vel*dt) ![m2]
		    vn=(Area**((2/3)*params(3)))*pm ![m/seg]
		    vel=(2*vn+vel)/3 ![m2/seg]
		enddo
		v_ladera(celdas)=vel		            			
		Es2=min(Area*v_ladera(celdas)*dt*Calib(4)/m3_mm,Sp(2,celdas)) ! [mm]		
		Sp(2,celdas)=Sp(2,celdas)-Es2 ![mm] Actualiza almacenamiento         
		!Tanque 3: Calculo de la cantidad de agua que escurre como flujo sub-superficial, nada precola		
		Sp(3,celdas)=Sp(3,celdas)+R3 ![mm]
		Es3=E3(celdas)*Sp(3,celdas) ![mm]
		Sp(3,celdas)=Sp(3,celdas)-Es3 ![mm]	     			    
		!Selecciona a donde enviar el flujo de acuerdo al tipo de celda
		select case (tipo_celda(1,celdas))
		    case(1)
			!Sigue drenando a su tanque igual aguas abajo
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2
			    Sp(3,drenaid)=Sp(3,drenaid)+Es3
			else
			    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s]
			endif
		    case(2,3)
			!Drena de nuevo al flujo superficial
			if (drena(celdas).ne.0) then
			    Sp(2,drenaid)=Sp(2,drenaid)+Es2+Es3
			else
			    Q(1,tiempo)=(Es2+Es3)*m3_mm/dt ![m3/s]
			endif
		end select
		!si es un punto de control guarda igualmente
		if (control(celdas).ne.0) then
		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
		    control_cont=control_cont+1
		endif
		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
		rain_sum=rain_sum+R1
	    enddo
	    !calcula la cantidad total de lluvia que entra
	    entradas=entradas+rain_sum
	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
	    mean_rain(tiempo)=rain_sum/nceldas
	enddo
	!Variable que al terminar en cero significa que el modelo se ejecuto
	corrio=0
	!Calcula el tiempo de ejecucion en segundos
	call cpu_time(tiempo_f)
	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
	!calcula el balance en la cuenca
	balance=100*((entradas-salidas-sum(Sp(:,:)))/entradas)
    else
	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
	corrio=1
	return
    endif    
end subroutine


!-----------------------------------------------------------------------
!Versiones modelos sencillos de transporte interpolados con TIN
!-----------------------------------------------------------------------
!subroutine model_tin_lc(drena,control,Calib,pert_tin,tin, Q,balance,mean_rain,Sp,corrio, ntin,nceldas,ncontrol,nreg) 
!    !variables de entrada
!    integer, intent(in) :: nceldas,ntin,ncontrol,nreg
!    integer, intent(in) :: drena(nceldas),control(nceldas),tin(3,ntin),pert_tin(nceldas)
!    real, intent(in) :: Calib 
!    !Variables de salida
!    real, intent(out) :: Q(ncontrol,nreg),balance,Sp(nceldas),mean_rain(nreg),corrio
!    !f2py intent(in) :: nceldas,ntin,ncontrol,drena,tipo_celda,control,tin,Calib,pot_rain,pert_tin
!    !f2py intent(out) :: Q, balance,S_private,mean_rain,corrio
!    !f2py note(<Esto es un ensayo>)
!    !Variables locales
!    real tiempo_i,tiempo_f !variables de tiempo de ejecucion del programa
!    integer celdas,tiempo !Iteradores para la cantidad de celdas y los intervalos de tiempo
!    integer drenaid !Vector con el id a donde drena cada celda
!    real Rain_past(nceldas),rain_sum !lluvia anterior, peso actual, suma pesos, suma pesos con lluvia, suma de la lluvia en el intervalo de tiempo
!    real R1,R2 !lluvia (R1) y flujo que drena por los tanques
!    real a(2,nceldas),b(2,nceldas),c(2,nceldas),d(2,nceldas),coef1(nceldas),coef2(3) !Coeficientes del determinante del interpolador de lluvia
!    real v_ladera(nceldas) !velocidad en ladera
!    real E2(nceldas),Es2 !Cantidad de flujo que viaja en cada celda
!    real m3_mm, L !Variable para convertir m3 a mm2, Longitud copiada 
!    real inicial,entradas,salidas !Cantidad de agua inicial, que entra al sistema y que sale [mm]
!    real X,Y !Posiciones actuales de la celda
!    integer i,control_cont,cont !iterador,contador de puntos de control, contador de lluvia faltante
!    !Evalua si permite o no ejecutar el modelo de acuerdo a los mapas dimensionados
!    if (check_all.eq.0) then
!	!Calcula el tiempo inicial	
!	call cpu_time(tiempo_i)
!	!Calcula las velocidades horizontales iniciales
!	v_ladera=1*Calib ![m/s] Estado constante
!	!Calcula el porcentaje de flujo horizontal saliente de cada celda
!	E2=1-L_celda(1,:)/(v_ladera*dt+L_celda(1,:))
!	m3_mm=(dxp**2)/1000 ![m2 * 1m/1000mm] --> [m3/mm] convierte tanque S5[mm] a [m3]
!	!Establece el almacenamiento inicial 
!	Sp=0
!	!Establece la cantidad de agua inicial dentro del sistema y la cantidad de salida y de entrada
!	inicial=sum(Sp(:))
!	entradas=0
!	salidas=0
!	!Determina como lluvia anterior y la lluvia media para cada celda el valor cero
!	Rain_past=0
!	mean_rain=0
!	!Determina los coeficientes para el calculo de la lluvia por tines
!	a(1,:)=coord_rain(1,TIN(1,pert_tin(:))); a(2,:)=coord_rain(2,TIN(1,pert_tin(:)))
!	b(1,:)=coord_rain(1,TIN(2,pert_tin(:))); b(2,:)=coord_rain(2,TIN(2,pert_tin(:)))
!	c(1,:)=coord_rain(1,TIN(3,pert_tin(:))); c(2,:)=coord_rain(2,TIN(3,pert_tin(:)))
!	d(1,:)=xll+ColFil(1,:)*dx-0.5*dx; d(2,:)=yll+(nrows-ColFil(2,:))*dx+0.5*dx
!	coef(:)=(b(1,:)-a(1,:))*(c(2,:)-a(2,:))-(c(1,:)-a(1,:))*(b(2,:)-a(2,:))			
!	!Itera para cada intervalo de tiempo
!	do tiempo=1,nreg
!	    !Reinicia el conteo de la lluvia y el contador de los puntos de control
!	    rain_sum=0
!	    control_cont=2
!	    !Itera para todas las celdas
!	    do celdas=1,nceldas
!		!Calcula la posicion de la celda objetivo y copia la posicion X,Y de la celda actual
!		drenaid=nceldas-drena(celdas)+1
!		X=XY_celdas(1,celdas)
!		Y=XY_celdas(2,celdas)	    
!		!resuelve el determinante y estima la lluvia sobre la celda en el intervalo de tiempo
!		coef2(1)=rain(TIN(1,Tperte(celdas)),tiempo); coef2(2)=rain(TIN(2,Tperte(celdas)),tiempo); coef2(3)=rain(TIN(3,Tperte(celdas)),tiempo)
!		do i=1,3
!		    if (coef2(i).eq.noData) then
!			coef2(i)=0
!		    endif
!		enddo
!		coef1=(c(1,celdas)-a(1,celdas))*(d(2,celdas)-a(2,celdas))*(coef2(2)-coef2(1))+(b(2,celdas)-a(2,celdas))*(coef2(3)-coef2(1))*(d(1,celdas)-a(1,celdas))-&
!		&(coef2(2)-coef2(1))*(c(2,celdas)-a(2,celdas))*(d(1,celdas)-a(1,celdas))-(d(2,celdas)-a(2,celdas))*(coef2(3)-coef2(1))*(b(1,celdas)-a(1,celdas))
!		R1=coef2(1)-coef1/coef(celdas) !lluvia sobre el pixel [mm]
!		if (R1.lt.0) R1=0		
!		!Tanque 1: todo lo que cae fluye eventualmente		
!		Sp(celdas)=Sp(celdas)+R1 ![mm] !Actualiza alm con la lluvia nueva
!		Es2=E2(celdas)*Sp(celdas) ! [mm] De acuerdo al almacenamiento calcula cuanto se va	
!		Sp(celdas)=Sp(celdas)-Es2 ![mm] Actualiza almacenamiento              			    
!		!Si es la salida de la cuenca guarda en la primera entrada del caudal
!		if (drena(celdas).eq.0) then
!		    Q(1,tiempo)=Es2*m3_mm/dt ![m3/s] 
!		    salidas=salidas+Es2 ![mm] Lo que sale del sistema en mm     
!		else
!		    Sp(drenaid)=Sp(drenaid)+Es2
!		endif
!		!si es un punto de control guarda igualmente
!		if (control(celdas).ne.0) then
!		    Q(control_cont,tiempo)=Es2*m3_mm/dt ![m3/s]      
!		    control_cont=control_cont+1
!		endif
!		!Lluvia media: suma la cantidad de lluvia que cae sobre la cuenca 
!		rain_sum=rain_sum+R1
!	    enddo
!	    !Calcula cantidad total de lluvia que entra al sistema
!	    entradas=entradas+rain_sum
!	    !Calculo de la lluvia media sobre la cuenca en el intervalo de tiempo
!	    mean_rain(tiempo)=rain_sum/nceldas
!	enddo
!	!Variable que al terminar en cero significa que el modelo se ejecuto
!	corrio=0
!	!Calcula el tiempo de ejecucion en segundos
!	call cpu_time(tiempo_f)
!	print *, 'Tiempo ejecucion: ', tiempo_f-tiempo_i
!	!calcula el balance de agua 
!	balance=100*((entradas-salidas-sum(Sp(:)))/entradas)
!	print *, 'inicial:', inicial
!	print *, 'balance:', balance
!	print *, 'entradas:', entradas
!	print *, 'salidas:' , salidas
!	print *, 'almacenado:', sum(Sp(:))
!    else
!	!Si falta alguna variable por asignar se retira sin ejecutar el modelo
!	corrio=1
!	return
!    endif    
!end subroutine

end module
