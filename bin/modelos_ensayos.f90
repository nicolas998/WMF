module modelos_ensayos


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




end module
