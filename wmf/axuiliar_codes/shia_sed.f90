subroutine shia_est_ocin_sed(lectura1,lectura2,lectura3,lectura4,lectura5,lectura6,lectura7,lectura8,lectura9,escribe,escribe2,escribe3,escribe4)
    use lecturas
    use sedimentos

!-------------------------------------------------------------------------------------------------------------------!
!Descripción

    !shia_est_ocin_sed: versión estatica del shia, en ésta versión se calcula la velocidad en el cauce de manera no lienal
    !empleando la Onda Cinemática Geomorfológica (Vélez 2001), además se calcula la velocidad en ladera mediante relación
    !no lineal entre la velocidad y el almacenamiento en el tanque.    
    !Se calculan sedimentos

!-------------------------------------------------------------------------------------------------------------------!
!Definición de variables
    
    !-----------------------------------------------------------------------!
    !Lectura y escritura
    character*255 lectura1,lectura2,lectura3,lectura4,lectura5,lectura6,lectura7,lectura8,escribe,escribe2,lectura9,escribe3,escribe4
    
    !-----------------------------------------------------------------------!
    !Variables Para la distribución de la lluvia
    real*4, allocatable :: a(:,:),b(:,:),c(:,:),d(:,:),coef(:)
    real*4 a3,b3,c3,coef1,R0        
    
    !-----------------------------------------------------------------------!
    !Variables empleadas en los tanques        
    !Variables de flujo
    real*4 R(5)           
    real*8 E,E2,E3,E4,E5
    !Almacenamiento en cada tanque
    real*4, allocatable :: S1(:),S2(:),S3(:),S4(:),S5(:)   
    !Parámetros
    real conver,prof,per                
    
    !-----------------------------------------------------------------------!
    !Variables empleadas por ahí
    integer i,j,z,k,cont,drenaid
    
    !-----------------------------------------------------------------------!
    !Variables de Onda Cinemática Geomorfológica
    real Bet1,eacum1,epend1,easec1,Ac,pen,vn !En el cauce
    real Bet2,eacum2,epend2,easec2 !En las carcavas
    real Bet,eacum,epend,easec !Las que se emplean dependiendo de si es car o cau   
    
    !-----------------------------------------------------------------------!
    !Variables para sedimentos
    real*4 Te(3),wi(3), alfaSED, Qskr, VolSed(3), VolDep(3), VolSal(3),Cw(3),Gesp
    real*4 dSED(3),ERO(3),EROt(3),DEP(3),DEPt(3)
    real*4, allocatable :: VolERO(:),VolDEPo(:) !Volumen erosionado por celda, Volumen depositado en la celda
    real*4, allocatable :: VS1(:),Vs2(:),VS3(:),VD1(:),VD2(:),VD3(:),PAre(:),PArc(:),PGra(:),Krus(:),Crus(:),Prus(:)
    real*4, allocatable :: VSc1(:),VSc2(:),VSc3(:),VDc1(:),VDc2(:),VDc3(:),tablaSED(:,:)
    real*4, allocatable :: Qsaux(:,:)
    
    !-----------------------------------------------------------------------!
    !Guarda resultados (por el momento sólo caudal)
    real*4, allocatable :: Q(:,:)
    
    !-----------------------------------------------------------------------!
    !Para realizar balance
    real*8 entradas, salidas, diferencia
    
    !-----------------------------------------------------------------------!
    !Variables para leer tabla de cuenca
    integer*4 filas,col,ncols,nrows,CantPtos
    integer ipen,iarea,dt
    real iv1,iv4
    real*8 xll,yll
    real*4 dx,nodata
    real, allocatable :: tabla(:,:)
    
    !-----------------------------------------------------------------------!
    !Variables de los triangulos, y distribución de lluvia    
    integer*4 numTri
    integer*4, allocatable :: tri(:,:)      
    integer*4 numReg,numEst    
    real*4, allocatable :: lluvia(:,:)
    real*8, allocatable :: coordEst(:,:)
    
    !-----------------------------------------------------------------------!
    !Variables par leer parámetros de calibración
    real R1,R2,R3,R4,R5,R6,R7,R8,R9,fc
    
    !-----------------------------------------------------------------------!
    !Variables par leer estado inicial de la cuenca    
    real si1,si2,si3,si4,si5

    !-----------------------------------------------------------------------!    
    !Variables de asignación de valores tabla
	integer, allocatable :: IDdren(:,:), ColFil(:,:), Ptos_cont(:), T_pert(:), nivel(:)
    real, allocatable :: ddx(:), Evp(:), Hu(:),v1(:),v2(:),v3(:),v4(:), Ks(:), Kp(:), Kpp(:), pend(:), man(:) 
    real*8, allocatable :: acum(:)
    
!-------------------------------------------------------------------------------------------------------------------!
!Lectura de datos    

    !Lectura de tabla
    call LeerTabla(lectura1,filas,col,xll,yll,CantPtos,dx,ncols,nrows,iv1,iv4,ipen,iarea,dt,tabla)
    !Lectura de Triangulos
    call LeerTriangulos(lectura2,numTri,Tri)
    !Lectura de datos y estaciones de precipitacion
    call leerLluvia(lectura3,numReg,lluvia,numEst,coordEst)    
    !Lectura de Parámetros de calibración
    call leeCalib(lectura4,R1,R2,R3,R4,R5,R6,R7,R8,R9,fc)
    !Lectura de estado inicial de la cuenca
    call leeInicial(lectura5,si1,si2,si3,si4,si5) 
    !Lectura de parametros OCG en cauce
    call parametrosOndaCau(lectura6,Bet1,eacum1,epend1,easec1)
    !Lectura de parametros OCG en carcavas
    call parametrosOndaCar(lectura9,Bet2,eacum2,epend2,easec2)        
    !Lectura de parametros de sedimentos
    call paramSED(lectura7,wi,dSED,alfaSED,Gesp)
    !Lectura de tabla con variables para calculo de sedimentos
    call leeTablaSED(lectura8,tablaSED)
    
    
!------------------------------------------****ACA EMPIEZA EL MODELO POR FIN****--------------------------------------!
!Corre el modelo
    
    !----------------------------------------------------------------------------------------------------------------!
    !Condiciones iniciales

	    !Asignaciones de valores de la tabla	    
	    tabla(16,filas+1)=4
	    allocate(IDdren(2,filas+1), ColFil(2,filas+1), Ptos_cont(filas+1), T_pert(filas+1), nivel(filas+1))
	    allocate(ddx(filas+1), Evp(filas+1), Hu(filas+1), Ks(filas+1),Kp(filas+1),Kpp(filas+1),man(filas+1))
	    allocate(v1(filas+1),v2(filas+1),v3(filas+1),v4(filas+1),pend(filas+1),acum(filas+1))
	    
	    !Tamaño a los vectores de sedimentos en suspención y depositados
	    allocate(VS1(filas+1),VS2(filas+1),VS3(filas+1),VD1(filas+1),VD2(filas+1),VD3(filas+1),VolERO(filas+1),VolDEPo(filas+1))
	    allocate(VSc1(filas+1),VSc2(filas+1),VSc3(filas+1),VDc1(filas+1),VDc2(filas+1),VDc3(filas+1))
	    
	    !Tamaño a los vectores de los parámetros de los sedimentos
	    allocate(Krus(filas+1),Crus(filas+1),Prus(filas+1),PAre(filas+1),PArc(filas+1),PGra(filas+1))
	    
	    !Parámetros iniciales (Valores estan en el archivo de parámetros)      
        allocate(S1(filas+1),S2(filas+1),S3(filas+1),S4(filas+1))
        allocate(S5(filas+1))
	    
	    !Calcula antes de correr algunos de los parámetros para interpolar la lluvia
        allocate(a(2,filas),b(2,filas),c(2,filas),d(2,filas),coef(filas));
        a=1; b=1; c=1; d=1; coef=1;
        
        !inicializa balance
        entradas=0; salidas=0
        
        !Itera para realizar calculos y asignar valores de la tabla
	    do i=1,filas	        
	        !----------------------------------------------------------------!
	        !Asignación de valotres de la tabla
	        IDdren(1,i)=tabla(1,i); IDdren(2,i)=tabla(2,i)
	        ColFil(1,i)=tabla(3,i); ColFil(2,i)=tabla(4,i)
	        Ptos_cont(i)=tabla(5,i); ddx(i)=tabla(6,i)
            T_pert(i)=tabla(8,i); nivel(i)=tabla(17,i)

            !----------------------------------------------------------------!
	        !Valores con calibracion	        
            Hu(i)=tabla(9,i)*R1 ![mm]
            EVP(i)=tabla(12,i)*R2 ![mm]
            Ks(i)=tabla(10,i)*R3*dt*fc ![mm]
            v1(i)=tabla(13,i)*R4 ![m/seg]
            Kp(i)=tabla(11,i)*R5*dt*fc ![mm]
            v2(i)=tabla(10,i)*R6*fc/1000 ![m/seg]
            Kpp(i)=tabla(11,i)*R7*dt*fc ![mm]
            v3(i)=tabla(11,i)*R8*fc/1000 ![m/seg]            
            pend(i)=tabla(17,i)/100
            if (pend(i)<0.001) then; pend(i)=0.01; endif    
            Acum(i)=tabla(7,i)*dx**2/1e6
            if (Acum(i)==0) then; Acum(i)=(dx**2)/1e6; endif 
            
            !----------------------------------------------------------------!
	        !Parámetros de los sedimentos
	        Krus(i)=tablaSED(1,i) ![adim] 
            Crus(i)=tablaSED(2,i) ![adim]
            Prus(i)=tablaSED(3,i) ![adim]
            PAre(i)=tablaSED(4,i) ![%]
            PGra(i)=tablaSED(5,i) ![%]
            PArc(i)=tablaSED(6,i) ![%]
            
            !----------------------------------------------------------------!
            !Estado inicial de los pixeles
            S1(i)=si1 ![mm]
            S2(i)=si2 ![mm]
            S3(i)=si3 ![mm]
            S4(i)=si4 ![mm]
            S5(i)=si5 ![m^3]        
	        
	        !----------------------------------------------------------------!
	        !Agrega datos a la entrada del balance
	        entradas = entradas + S1(i)+S2(i)+S3(i)+S4(i)+S5(i)	        
	        
	        !----------------------------------------------------------------!
	        !Calculo de valores para la iteración de la lluvia
	        a(1,i)=coordEst(1,Tri(1,T_pert(i))); a(2,i)=coordEst(2,Tri(1,T_pert(i)))
            b(1,i)=coordEst(1,Tri(2,T_pert(i))); b(2,i)=coordEst(2,Tri(2,T_pert(i)))
            c(1,i)=coordEst(1,Tri(3,T_pert(i))); c(2,i)=coordEst(2,Tri(3,T_pert(i)))
            d(1,i)=xll+ColFil(1,i)*dx-0.5*dx; d(2,i)=yll+(nrows-ColFil(2,i))*dx+0.5*dx
            coef(i)=(b(1,i)-a(1,i))*(c(2,i)-a(2,i))-(c(1,i)-a(1,i))*(b(2,i)-a(2,i))
	    end do               	           
        deallocate(tabla,tablaSED)        
            

    !----------------------------------------------------------------------------------------------------------------!
    !Comienza a correr       

        !Tamaño Vector Caudales simulados y de sedimentos simulados
        allocate(Q(CantPtos,numReg),Qsaux((CantPtos)*3,numReg))
        
        !Para convertir de mm a m^3
        conver=(dx**2)/1000 ![m2]*[m/mm]=[m3/mm]                
        
        !Estado inicial de las velocidades
        v1=0.1; v4=1; ![m/s]
        
        !estado inicial del almacenamiento de sedimentos
        VS1=0;VS2=0;VS3=0 ![m3]
        VD1=0;VD2=0;VD3=0 ![m3]
        VSc1=0;VSc2=0;VSc3=0 ![m3]
        VDc1=0;VDc2=0;VDc3=0 ![m3]
        VolERO=0; VolDEPo=0 ![m3]
        EROt=0; DEPt=0
        
        !Itera en cada uno de los intervalos de tiempo
        do i=1,numReg
            cont=0;    
            !itera sobre cada una de las celdas en el intervalo de tiempo
            do j=1,filas
            
            !Calcula el valor de a donde drena la cosa
            drenaid=IDdren(1,j)-IDdren(2,j)+j
            
            !--------------------------------------------------------------------------------------------------------!
            !Se calcula R1 para el pixel de acuerdo al triangulo de pertenencia
                
                !Obtiene los valores del determinante
                a3=lluvia(Tri(1,T_pert(j)),i); b3=lluvia(Tri(2,T_pert(j)),i); c3=lluvia(Tri(3,T_pert(j)),i)                
                !Obtiene los coeficientes para obtener z (o la lluvia)
                coef1=(c(1,j)-a(1,j))*(d(2,j)-a(2,j))*(b3-a3)+(b(2,j)-a(2,j))*(c3-a3)*(d(1,j)-a(1,j))-(b3-a3)*(c(2,j)-a(2,j))*(d(1,j)-a(1,j))-(d(2,j)-a(2,j))*(c3-a3)*(b(1,j)-a(1,j))                
                !Obtiene el valor de la lluvia
                R(1)=a3-coef1/coef(j) !lluvia sobre el pixel [mm]
                !Actualiza entradas del balance
                entradas=entradas+R(1)
                
            !--------------------------------------------------------------------------------------------------------!
            !Primer Tanque (Flujo Capilar)    
                
                !Calcula el excedente del almacenamiento capilar    (recordar: Sm13(1,j)=Sm1(j))
                !R(2)=R(1) - min(R(1)*(1-(S1(j)/Hu(j))**2),Hu(j)-S1(j)) ![mm]
                R(2)=max(0.0,R(1)-Hu(j)+S1(j))
                !Calcula el almacenamiento en el tanque
                S1(j)=min(S1(j)+R(1)-R(2),Hu(j)) ![mm] (Sm13(1,j)=Sm1)
                !Calcula lo que se evapora                
                E=min(Evp(j)*(S1(j)/Hu(j))**0.6,S1(j)) ![mm]
                !Actualiza el almacenamiento en el tanque
                S1(j)=S1(j)-E ![mm]
                !Actualiza entradas del balance
                salidas=salidas+E
                
            !--------------------------------------------------------------------------------------------------------!
            !Segundo Tanque (Flujo Superficial)
            
                !Agua que sigue gracias a la infiltración
                R(3)=min(R(2),Ks(j)) ![mm] (método del umbral) (recordar: Ksp(1,j)=Ks)
                !Calcula el almacenamiento
                S2(j)=S2(j)+R(2)-R(3) ![mm]                
                !Calcula la escorrentía superficial (Embalse lineal)                            
                E2=(1-ddx(j)/(v1(j)*dt+ddx(j)))*S2(j) ![mm]
                !Calcula el área transversal
                Area=S2(j)*conver/(ddx(j)+v1(j)*dt)
                !Actualiza el almacenamiento
                S2(j)=S2(j)-E2 ![mm] 
                !Para el balance
                if (j==filas) then
                    Salidas = Salidas + E2
                endif
                
                !******************************************************************
                !Sedimentos en Ladera                    
                    !Calcula depósito en el pixel sólo si hay suspendido                     
                    !Calcula Eficiencia de Atrapamiento (método 1)                        
                    call TrapEfic(Te,wi,S2(j),dt)  
                    !Calcula Eficiencia de Atrapamiento (método 2)
                    !call TrapEfic2(Te,wi,S2(j),ddx(j),v1(j))  
                    !Calcula los sedimentos depositados en el pixel                  
                    call CalcDep(VS1(j),VS2(j),VS3(J),VD1(j),VD2(j),VD3(j),Te,DEP)                                        
                    !Calcula capacidad 
                    call CalcQskr(Qskr,area,v1(j),Krus(j),Crus(j),Prus(j),ddx(j),dt,alfaSED,pend(j))
                    !Calcula lo que se entrega al siguiente pixel                                                           
                    call CalcSEDlad(VS1(j),VS2(j),VS3(J),VD1(j),VD2(j),VD3(j),VolSal,PAre(j),PGra(j),PArc(j),Te,Qskr,v1(j),dt,ddx(j),ERO,DEP)
                    !Suma la erosión total de cada fracción
                    EROt=EROt+ERO
                    !Suma la depositación total de cada fracción
                    DEPt=DEPt+DEP
                    !Actualiza la erosión en la celda
                    VolERO(j)=VolERO(j)+sum(ERO)
                    !Actualiza la depositación enla celda
                    VolDEPo(j)=VolDEPo(j)+sum(DEP)
                                                                             
            !--------------------------------------------------------------------------------------------------------!
            !Tercer Tanque (Flujo Subsuperficial)
            
                !Agua que precola al suelo
                R(4)=min(R(3),Kp(j)) ![mm] (método del umbral) (recordar: Ksp(2,j)=Kp)
                !Actualiza almacenamiento
                S3(j)=S3(j)+R(3)-R(4) ![mm] (Sm13(2,j)=Sm3)                
                !Calcula lo que sale del tanque                                               
                E3=(1-ddx(j)/(v2(j)*dt+ddx(j)))*S3(j) ![mm] Estilo juan Jose                                            
                !Actualiza almacenamiento
                S3(j)=S3(j)-E3 ![mm]
                !Para el balance
                if (j==filas) then             
                    Salidas = Salidas + E3
                endif 
                
            !--------------------------------------------------------------------------------------------------------!
            !Cuarto Tanque (Flujo Subterráneo)
            
                !Agua que entra al acuífero
                R(5)=min(R(4),Kpp(j)) ![mm] estilo juan jose
                !Actualiza el almacenamiento en el acuífero
                S4(j)=S4(j)+R(4)-R(5) ![mm]                                 
                !Calcula lo que sale del tanque                                                                
                E4=(1-ddx(j)/(v3(j)*dt+ddx(j)))*S4(j) ![mm] Estilo Juan jose                                                       
                !Actualiza el almacenamiento
                S4(j)=S4(j)-E4 ![mm]
                !Para el balance
                if (j==filas) then
                    Salidas = Salidas + E4
                endif                 
            
            !--------------------------------------------------------------------------------------------------------!
            !Distribución del agua
                
                select case(nivel(j))
                    case(1)
                        !Transporte Agua
                        S2(drenaid)=S2(drenaid)+E2 ![mm]
                        S3(drenaid)=S3(drenaid)+E3 ![mm]
                        S4(drenaid)=S4(drenaid)+E4 ![mm]
                        !Transporte Sedimentos
                        VS1(drenaid)=VS1(drenaid)+VolSal(1) !Vol en suspención fracción 1 [m3]
                        VS2(drenaid)=VS2(drenaid)+VolSal(2) !Vol en suspención fracción 2 [m3]
                        VS3(drenaid)=VS3(drenaid)+VolSal(3) !Vol en suspención fracción 3 [m3]
                    case(2)
                        !Transporte Agua
                        S5(j)=S5(j)+(E2+E3)*conver ![m3]
                        S4(drenaid)=S4(drenaid)+E4 ![mm]
                        !Transporte Sedimentos
                        VSc1(j)=VSc1(j)+VolSal(1) !Vol en suspención fracción 1 [m3]
                        VSc2(j)=VSc2(j)+VolSal(2) !Vol en suspención fracción 2 [m3]
                        VSc3(j)=VSc3(j)+VolSal(3) !Vol en suspención fracción 3 [m3]
                        !Selección de variables para carcava
                        easec=easec2; eapend=eapend2; eacum=eacum2; Bet=Bet2                       
                    case(3)
                        !Transporte Agua
                        S5(j)=S5(j)+(E2+E3+E4)*conver ![m3]
                        !Transporte Sedimentos
                        VSc1(j)=VSc1(j)+VolSal(1) !Vol en suspención fracción 1 [m3]
                        VSc2(j)=VSc2(j)+VolSal(2) !Vol en suspención fracción 2 [m3]
                        VSc3(j)=VSc3(j)+VolSal(3) !Vol en suspención fracción 3 [m3]
                        !Selección de variables para cauce
                        easec=easec2; eapend=eapend2; eacum=eacum2; Bet=Bet2
                end select    
            
            !--------------------------------------------------------------------------------------------------------!
            !Quinto Tanque (Cauce, no todas las celdas tinen cauce, y éste varia en 2 niveles)
                               
                !Si hay agua en el quinto tanque éste la pasa aguas abajo
                if (S5(j)>0) then                                             
		            !Asume una velocidad inicial		            		            
		            pen=pend(j) !Copia para ahorrar velocidad 
		            Ac=acum(j) !Copia para ahorrar velocidad 
		            do z=1,4
		                Area=S5(j)/(ddx(j)+(v4(j)*dt)) ![m2]
		                vn=(R9*(Area**easec)*(pen**epend)*(Ac**eacum))*Bet
		                v4(j)=(2*vn+v4(j))/3
		            enddo		            
		            !calcula lo que sale con un embalse lineal                                  
		            E5=Area*v4(j)*dt ![m3]  
		            !Sale del tanque lo que se drena al otro
                    S5(j)=S5(j)-E5 ![m3]           
                    if (j==filas) then; salidas=salidas+E5/conver; endif 
                else
                    E5=0; v4(j)=0; Area=0;
                endif   
                    
                !******************************************************************
                !Sedimentos en Cauce                                                
                    !Transforma la profundidad a mm
                    prof=S5(j)/conver ![mm]                        
                    !Calcula cuanto se deposita si hay sedimentos suspendido (método 1)
                    call TrapEfic(Te,wi,prof,dt)             
                    !Calcula lo que se queda atrapado
                    call CalcDep(VSc1(j),VSc2(j),VSc3(J),VDc1(j),VDc2(j),VDc3(j),Te,DEP)                                      
                    !Calcula la concentración usando Engelund y Hansen (1967)
                    call CalcEH(Cw,Gesp,dSED,v4(j),pen,Area,prof)
                    !Calcula los sedimentos que son transportados
                    call CalcSEDCorr(VSc1(j),VSc2(j),VSc3(j),VDc1(j),VDc2(j),VDc3(j),VolSal,Cw,E5,dt,v4(j),ddx(j),DEP)
                    !Actualiza depositación en la celda                        
                    VolDEPo(j)=VolDEPo(j)+sum(DEP)
                    !Suma la depositación total de cada fracción
                    DEPt=DEPt+DEP
                    
                !******************************************************************
                !Vuelve al agua normal                    
                                         
                !Entrega al tanque de corriente del pixel destino
                if (IDdren(2,j)/=0) then
                    !Transporta agua
                    S5(drenaid)=S5(drenaid)+E5  ![m3] 
                    !Transporte Sedimentos
                    VSc1(drenaid)=VSc1(drenaid)+VolSal(1) !Vol en suspención fracción 1 [m3]
                    VSc2(drenaid)=VSc2(drenaid)+VolSal(2) !Vol en suspención fracción 2 [m3]
                    VSc3(drenaid)=VSc3(drenaid)+VolSal(3) !Vol en suspención fracción 3 [m3]                     
                end if   
        
            !--------------------------------------------------------------------------------------------------------!
            !Guarda resultados en la variable Q
                if (Ptos_Cont(j)/=0) then 
                    cont=cont+1;
                    Q(Ptos_cont(j),i)=E5/dt ![m3/s]                          
                    Qsaux(1+3*(Ptos_Cont(j)-1),i)=VolSal(1)/dt ![m3/s]
                    Qsaux(2+3*(Ptos_Cont(j)-1),i)=VolSal(2)/dt ![m3/s]
                    Qsaux(3+3*(Ptos_Cont(j)-1),i)=VolSal(3)/dt ![m3/s]
                end if                                              
            end do
            write(*,*) i                   
        end do
        
        !Escribe balance
        do i=1,filas
            salidas=salidas+S1(i)+S2(i)+S3(i)+S4(i)+S5(i)/conver
        end do
        write(*,*) salidas
        write(*,*) entradas
        diferencia=((salidas-entradas)/salidas)*100
        write(*,*) 'Diferencia en el balance'
        write(*,*) diferencia
    !----------------------------------------------------------------------------------------------------------------1
    !Escribe resultados
    
        !Escribe el caudal simulado
        open(30,file=escribe,status='replace')
            write(30,*) 'Caudal Simulado'
            write(30,*) 'Balance:'
            write(30,*) 'Salidas:',salidas
            write(30,*) 'Entradas:',entradas
            write(30,*) 'Diferencia:',diferencia
            write(30,*) 'Puntos de control y salida de la cuenca'
            call mtrxwrt(Q,numReg,CantPtos,30)
        close(30)

        !Escribe los sedimentos simulados
        open(30,file=escribe2,status='replace')
            write(30,*) 'Sedimentos simulados'
            write(30,*) 'Erodado de cada fracción [m3]'
            write(30,*) 'Arenas          Limos          Arcillas'
            write(30,*) EROt(1),EROt(2),EROt(3)
            write(30,*) 'Depositado de cada fracción [m3]'
            write(30,*) 'Arenas          Limos          Arcillas'
            write(30,*) DEPt(1),DEPt(2),DEPt(3)      
            write(30,*) 'Producción de cada fracción [m3]'
            write(30,*) 'Arenas          Limos          Arcillas'
            write(30,*) EROt(1)-DEPt(1),EROt(2)-DEPt(2),EROt(3)-DEPt(3)
            write(30,*) 'Producción total de sedimentos [m3]'
            write(30,*) sum(EROt-DEPt)            
            write(30,*) 'Cada 3 columnas son un pt de contro, Arenas, Limos, Arcillas [m3s-1]'
            call mtrxwrtSED(Qsaux,numReg,(CantPtos)*3,30)
        close(30)
        
        !Escribe mapa de volumen erodado
        nodata=-9999
        call escriMapa(escribe3,VolERO,ColFil,ncols,nrows,xll,yll,dx,nodata,filas)

        !Escribe mapa de volumen depositado
        call escriMapa(escribe4,VolDEPo,ColFil,ncols,nrows,xll,yll,dx,nodata,filas)

end subroutine