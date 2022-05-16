!!    This program solves 2D Euler equations by using 
!!		HOCUS-BVD
! https://doi.org/10.1016/j.jcp.2020.110067
! https://arxiv.org/pdf/2012.09905
!**********************************************************************
! Written by Sainath Ch, Email: skywayman9@gmail.com
!**********************************************************************

	program euler

	implicit none


    integer 				:: i, j
	integer					:: N, rk_step
	integer					:: NTMAX=100000
	double precision		:: t_end, time , dt

	double precision		:: CFL =0.19d0, gamma =1.4d0
	!test case 1 Riemann problem
	! integer, parameter 		:: NX = 400, NY = 400, ghostp = 5, n_eqn = 4

	!test case 2 Shock-entropy
	! integer, parameter 		:: NX = 400, NY = 80, ghostp = 5, n_eqn = 4

	!test case 3 Richtmeyer - Meshkov
	integer, parameter 		:: NX = 320, NY = 80, ghostp = 4, n_eqn = 4

	!test case 4 Shock-Vortex
	! integer, parameter 		:: NX = 400, NY = 200, ghostp = 5, n_eqn = 4



	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2)
	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    integer, parameter		:: file_save=200

	integer 				:: time_ini,time_end,time_calc,test_case
	double precision 		:: start, finish

	common /grid/ dx, dy


	call cpu_time(start)
	write(*,*) 'Program start...'

	call system_clock(count = time_ini)


	test_case = 3


if(test_case .eq. 1) then
	xmin =  0.00d0
	xmax =  1.0d0

	ymin =  0.00d0
	ymax =  1.00d0
elseif(	test_case .eq. 2) then

	xmin =  -5.00d0
	xmax =  5.0d0

	ymin =  -1.00d0
	ymax =  1.00d0
elseif(	test_case .eq. 3) then
	xmin = 0.0d0
	xmax = 4.0d0

	ymin = 0.0d0
	ymax = 1.0d0

elseif(	test_case .eq. 4) then
	xmin =   0.00d0
	xmax =   2.0d0

	ymin =   0.00d0
	ymax =   1.00d0
endif

	
    dx = (xmax - xmin)/NX
	dy = (ymax - ymin)/NY

	do i = -ghostp, NX + ghostp

		  x(i) = xmin + (i-0.5d0)*dx
	enddo

	do j = -ghostp, NY + ghostp

	  	y(j) = ymin + (j-0.5d0)*dy
	
	enddo	



	call initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp,test_case,t_end)
	call timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)


	time = 0.0d0
	N=1
	
	  write(*,*)'*********************************************'
      write(*,*)'   time step N        time             '
      write(*,*)'*********************************************'

     ! Computations starts here

	do while(time.lt.t_end)
			call timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
			
			time = time + dt

			write(*,*) N ,time

			do rk_step = 0,2

				call boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp,test_case)


				do j = -ghostp, NY + ghostp
					do i = -ghostp, NX + ghostp

						cons(i,j,1,rk_step ) = density(i,j)
					    cons(i,j,2,rk_step ) = density(i,j)*u_vel(i,j)
					    cons(i,j,3,rk_step ) = density(i,j)*v_vel(i,j)
					    cons(i,j,4,rk_step ) = pressure(i,j)/(gamma-1.0d0) + 0.5d0*density(i,j)*(u_vel(i,j)**2.0d0 + v_vel(i,j)**2.0d0)

					enddo
				enddo


				call FX(density, u_vel, v_vel, pressure, residual, gamma,NX,NY,ghostp)

				call GY(density, u_vel, v_vel, pressure, residual, gamma,NX,NY,ghostp)
				
				
				call rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)

			enddo

			N=N+1
			

			if(MOD(N,file_save) .eq. 0) then

				call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

			endif	

			if (abs(time-t_end) .le. 1.0d-06) then

				write(*,*)'*********************************************'
           		write(*,*)'   Number of time steps = ',N
          	    write(*,*)'*********************************************'

          	    exit
          	endif

    enddo

    	
    		call output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)
    		call tecplot(N/file_save,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

			call system_clock(count = time_end)
    
   			time_calc = time_end-time_ini
    
    	write(*,'(A20,I10,A)')   'Calculation time ',time_calc,' [CPU ticks]'

    		call cpu_time(finish)
    
    	write(*,*) " Total CPU time to solution = ", finish-start, " seconds"

    	write(*,*) 'Program ends...'


	end program euler


	!***********************************************************************
	!*****                       Initial conditions                    *****
	!***********************************************************************


	subroutine initialconditions(x,y,density, u_vel, v_vel,pressure, sound, gamma,NX,NY,ghostp,test_case,t_end)

	implicit none

	integer 				:: i, j, test_case

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy, xmin, xmax, ymin, ymax, gamma, t_end

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: f,tau,temp

	double precision,parameter :: pi=acos(-1.0d0)

	double precision,parameter :: rc=0.05d0,epsilon=0.3d0,alpha=0.204


	if(test_case .eq. 1) then

	!Configuration-3 Riemann problem

	t_end =0.8d0

	do i = 1, NX

		do j =  1, NY

	    if ( x(i) .gt. 0.8d0 .and. y(j) .gt.  0.8d0) then
		density(i,j) = 1.5d0;
        u_vel(i,j)   = 0.0d0;
        v_vel(i,j)   = 0.0d0;
        pressure(i,j)= 1.5d0;
      
       elseif ( x(i) .le. 0.8d0 .and. y(j) .gt.  0.8d0) then
        
        density(i,j) = 33.0d0/62.0d0;
        u_vel(i,j)   = 4.0d0/11.0d0**0.5d0;
        v_vel(i,j)   = 0.00d0;
       	pressure(i,j)= 0.3d0;
       
       else if ( x(i) .le. 0.8d0 .and. y(j) .le.  0.8d0) then
        
        density(i,j) = 77.0d0/558.0d0;
        u_vel(i,j)   = 4.0d0/11.0d0**0.5d0;
        v_vel(i,j)   = 4.0d0/11.0d0**0.5d0;
        pressure(i,j)= 9.0d0/310.0d0;
       
       else if ( x(i) .gt. 0.8d0 .and. y(j) .le.  0.8d0) then
        
        density(i,j) = 33.0d0/62.0d0;
        u_vel(i,j)   = 0.0d0;
        v_vel(i,j)   = 4.0d0/11.0d0**0.5d0;
        pressure(i,j)= 0.3d0;
      	
      	endif

			
			sound(i,j)			= (gamma*pressure(i,j)/density(i,j))**0.5d0
		enddo	
	enddo

	elseif(test_case .eq. 2) then

	! Shock-entropy case
	t_end = 1.8d0

	do i = 1, NX
		do j =  1, NY

	   if ( x(i) .lt. -4.0d0) then
		density(i,j) = 27.0d0/7.0d0
        u_vel(i,j)   = 4.0d0*dsqrt(35.0d0)/9.0d0
        v_vel(i,j)   = 0.0d0;
        pressure(i,j)= 31.0d0/3.0d0
      
       
       else
        
        density(i,j) = 1.0d0+ (0.2d0*dsin(10.0d0*x(i)*acos(pi/6.0d0)+10.0d0*y(j)*asin(pi/6.0d0)))
        u_vel(i,j)   = 0.0d0;
        v_vel(i,j)   = 0.0d0;
        pressure(i,j)= 1.d0;
      	
      	endif

			
		sound(i,j)			= (gamma*pressure(i,j)/density(i,j))**0.5d0
		enddo	
	enddo

	

	elseif(test_case .eq. 3) then

	! Richtmeyer - Meshkov

	t_end = 9.0d0

	do i = 1, NX

		do j =  1, NY

			if(x(i) .lt. 2.9d0-0.1d0*sin(2.0d0*pi*(y(j)+0.25d0)) )then

				density(i,j) 	= 5.04d0
				u_vel(i,j)		= 0.0d0
				v_vel(i,j)		= 0.0d0
				pressure(i,j)	= 1.0d0				
				

			elseif(x(i) .lt. 3.2d0 )then

				density(i,j)	= 35517.0d0/35518.0d0 
				u_vel(i,j)		= 0.0d0
				v_vel(i,j)		= 0.0d0
				pressure(i,j)	= 1.0d0

			else
				density(i,j)	= 1.411d0
				u_vel(i,j)		=-665.0d0/1556.0d0
				v_vel(i,j)		= 0.0d0
				pressure(i,j)	= 1.628d0
				
			endif
			
			sound(i,j)			= (gamma*pressure(i,j)/density(i,j))**0.5d0
		enddo	
	enddo

	elseif(test_case .eq. 4) then

	!  Shock-vortex

	t_end =0.6d0

	do i = 1, NX

		do j =  1, NY


		if ( x(i) .lt. 0.5d0) then
		
		f = ( (x(i)-0.25d0)**2.0d0+(y(j)-0.5d0)**2.0d0 )**0.5d0

		tau=f/rc

		u_vel(i,j)   = dsqrt(1.4d0)+(epsilon*((y(j)-0.5d0)/rc)*dexp(alpha*(1.0d0-tau**2.0d0)));

        v_vel(i,j)   = 0.0d0- (epsilon*((x(i)-0.25d0)/rc)*dexp(alpha*(1.0d0-tau**2.0d0)));
		
		density(i,j) = (1.0d0-( ((gamma-1.0d0)*epsilon*epsilon*dexp(2.0d0*alpha*(1.0d0-tau**2.0d0)))&
       	               /(4.0d0*gamma*alpha)))**(1.0d0/(gamma-1.0d0));


        
        pressure(i,j)= density(i,j)**gamma;
      
       else


		v_vel(i,j)   = 0.0d0 

		temp = (gamma**0.5d0)+dsqrt(2.0d0)*((1.0d0-1.3d0)/dsqrt((gamma-1.0d0)+1.30d0*(gamma+1.0d0)))

		u_vel(i,j)   = temp 


       	density(i,j) = 1.0d0*(gamma-1.0d0+(gamma+1.0d0)*1.3d0)/(gamma+1.0d0+(gamma-1.0d0)*1.3d0)
                      
        pressure(i,j)= 1.3d0

		sound(i,j)	= (gamma*pressure(i,j)/density(i,j))**0.5d0
			endif



		enddo	
	enddo



	else

	write(*,*) 'More tests will be added later'

	endif


	end subroutine initialconditions

	!***********************************************************************
	!*****                       Output 			                   *****
	!***********************************************************************

	subroutine output(density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)

	integer 				:: i, j

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)


	open(unit=25, file="soln.txt",action="write",status="replace")
    do j = 1, NY
        do i = 1, NX
	 
	  write(25,'(6F25.8)') x(i),y(j),density(i,j),pressure(i,j),u_vel(i,j),v_vel(i,j)
	 enddo
	enddo
	close(25)


 	end subroutine output

 	subroutine tecplot(file,density,u_vel,v_vel,pressure,x,y,NX,NY,ghostp)
	implicit none


    integer 				:: i, j, file,l

	integer			 		:: NX, NY, ghostp

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	character(len=8) 		:: number*4, file_name

	write(number,'(i4.4)') file
	file_name="Rslt"//number
	open(unit=1,file=file_name//'.tec')
	
	write(1,*) 'TITLE="',file_name,'"'
    write(1,*) 'VARIABLES = "x","y","rho","vx","vy","Pre"'
	write(1,*) "ZONE I=",NX," J=",NY," F=POINT"


    do j = 1, NY
        do i = 1, NX

          write(1,'(6F25.8)') x(I), y(J), density(I,J), u_vel(I,J), v_vel(I,J), pressure(I,J)

	  enddo
    enddo
	
    close(1)



    end subroutine tecplot

	!***********************************************************************
	!*****                       Compute time step                	   *****
	!***********************************************************************

 	subroutine timestep(u_vel,v_vel,density,pressure,sound,CFL,time,t_end,dt,NX,NY,ghostp,gamma)
 	implicit none


 	integer 				:: i, j

	integer			 		:: NX, NY, ghostp
	double precision 		:: dx, dy,gamma

	double precision 		:: dt, time, t_end, dtnew, CFL

 	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
 	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: x_velocity, y_velocity

	common /grid/ dx, dy

	dt = 1.0d10

		do i = 0, NX

			do j = 0, NY
				sound(i,j) =  (gamma*pressure(i,j)/density(i,j))**0.5d0
				x_velocity =  ABS(u_vel(i,j)) + sound(i,j)
				y_velocity =  ABS(v_vel(i,j)) + sound(i,j)

				dtnew = min(dx/x_velocity, dy/y_velocity)

				if(dtnew .lt. dt) dt = dtnew
			enddo
		enddo

		dt = CFL*dt

		if ((time+dt) .gt. t_end ) then

			dt = t_end - time
		endif	


 	end subroutine timestep


 	!***********************************************************************
	!*****                       Boundary conditions                   *****
	!***********************************************************************


 	subroutine boundaryconditions(density, u_vel, v_vel, pressure, x, y, time, gamma,NX,NY,ghostp,test_case)
    implicit none

	integer 				:: i, j, NX, NY, ghostp,test_case

	double precision 		:: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
	double precision		:: dx, dy,gamma

	double precision		:: time

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)



	common /mesh/ dx, dy

	
	if(test_case .eq. 1) then

	! Set left and right boundary conditions.
	do i = 1,ghostp
		do j = 1, NY

		density(-i+1,j)	    = density(i,j)
		u_vel(-i+1,j)		= u_vel (i,j)
		v_vel(-i+1,j)		= v_vel (i,j)
		pressure(-i+1,j)	= pressure(i,j)

		
		density(NX+i,j)		= density(NX-i+1,j)
		u_vel(NX+i,j)		= u_vel(NX-i+1,j)
		v_vel(NX+i,j)		= v_vel(NX-i+1,j)
		pressure(NX+i,j)	= pressure(NX-i+1,j)

		enddo
	enddo

	! Set top and bottom boundary condition
	do j = 1,ghostp

		do i = 1, NX

			density(I,-J+1)  = density(i,j)
			u_vel(I,-J+1)    = u_vel(i,j)
			v_vel(I,-J+1)    = v_vel(i,j)
			pressure(I,-J+1) = pressure(i,j)

	    
	      	density(i,NY+j)		=  density(i,NY-j+1)
	      	u_vel(i,NY+j)		=  u_vel(i,NY-j+1)
	      	v_vel(i,NY+j)		=  v_vel(i,NY-j+1)
	      	pressure(i,NY+j)	=  pressure(i,NY-j+1)

	  enddo
	enddo

	elseif(test_case .eq. 2) then

	do i = 1,ghostp
		do j = 1, NY

		density(-i+1,j) = 27.0d0/7.0d0
        u_vel(-i+1,j)   = 4.0d0*dsqrt(35.0d0)/9.0d0
        v_vel(-i+1,j)   = 0.0d0;
        pressure(-i+1,j)= 31.0d0/3.0d0

		
		density(NX+i,j)		= density(NX-i+1,j)
		u_vel(NX+i,j)		= u_vel(NX-i+1,j)
		v_vel(NX+i,j)		= v_vel(NX-i+1,j)
		pressure(NX+i,j)	= pressure(NX-i+1,j)

		enddo
	enddo

	do j = 1,ghostp
		do i = 1, NX

			density(I,-J+1)  = density(i,j)
			u_vel(I,-J+1)    = u_vel(i,j)
			v_vel(I,-J+1)    = v_vel(i,j)
			pressure(I,-J+1) = pressure(i,j)

	    
	      	density(i,NY+j)		=  density(i,NY-j+1)
	      	u_vel(i,NY+j)		=  u_vel(i,NY-j+1)
	      	v_vel(i,NY+j)		=  v_vel(i,NY-j+1)
	      	pressure(i,NY+j)	=  pressure(i,NY-j+1)


	  enddo
	enddo

	elseif(test_case .eq. 3) then

	do i = 1,ghostp
		do j = 1, NY

		density(-i+1,j)	= density(i,j)
		u_vel(-i+1,j)		= u_vel (i,j)
		v_vel(-i+1,j)		= v_vel (i,j)
		pressure(-i+1,j)	= pressure(i,j)
		
		density(NX+i,j)		= density(NX-1+i,j)
		u_vel(NX+i,j)		= u_vel(NX-1+i,j)
		v_vel(NX+i,j)		= v_vel(NX-1+i,j)
		pressure(NX+i,j)	= pressure(NX-1+i,j)

		enddo
	enddo

	do j = 1,ghostp
		do i = 1, NX

			density(i,-j+1) = density(i,NY-j+1)
			u_vel(i,-j+1) = u_vel(i,NY-j+1)
			v_vel(i,-j+1) = v_vel(i,NY-j+1)
			pressure(i,-j+1) = pressure(i,NY-j+1)

	    
	      	density(i,NY+j)		=  density(i,j)
	      	u_vel(i,NY+j)		=  u_vel(i,j)
	      	v_vel(i,NY+j)		=  v_vel(i,j)
	      	pressure(i,NY+j)	=  pressure(i,j)

	  enddo
	enddo

	elseif(test_case .eq. 4) then

	do i = 1,ghostp

		do j = 1, NY

		
		density(-i+1,j)	= density(i,j)
		u_vel(-i+1,j)		= u_vel (i,j)
		v_vel(-i+1,j)		= v_vel (i,j)
		pressure(-i+1,j)	= pressure(i,j)

		
		density(NX+i,j)		= 1.0d0*(gamma-1.0d0+(gamma+1.0d0)*1.3d0)/(gamma+1.0d0+(gamma-1.0d0)*1.3d0)
		u_vel(NX+i,j)		= (gamma**0.5d0)+dsqrt(2.0d0)*((1.0d0-1.3d0)/dsqrt((gamma-1.0d0)+1.30d0*(gamma+1.0d0)))
		v_vel(NX+i,j)		= 0.0d0
		pressure(NX+i,j)	= 1.3d0 

		enddo
	enddo

	! Set top and bottom boundary condition

	do j = 1,ghostp

		do i = 1, NX

			density(I,-J+1)  = density(i,j)
			u_vel(I,-J+1)    = u_vel(i,j)
			v_vel(I,-J+1)    = v_vel(i,j)
			pressure(I,-J+1) = pressure(i,j)

	    
	      	density(i,NY+j)		=  density(i,NY-j+1)
	      	u_vel(i,NY+j)		=  u_vel(i,NY-j+1)
	      	v_vel(i,NY+j)		=  v_vel(i,NY-j+1)
	      	pressure(i,NY+j)	=  pressure(i,NY-j+1)

	  enddo
	enddo
     

	else

	write(*,*) 'Oi, Wrong'


	endif
     
      end

 	!***********************************************************************
	!*****                       Time step, TVD- Runge Kutta           *****
	!***********************************************************************

	subroutine rungekutta(density, u_vel, v_vel, pressure, cons, residual, gamma, rk_step, dt,NX,NY,ghostp,n_eqn)
	implicit none

	integer 				:: i, j, k, rk_step, rk

	integer 				:: NX, NY, ghostp,n_eqn

	double precision 		:: dt, uv, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

	double precision		:: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn,0:2) ! 0:2 is for the three Runge-Kutta time steps

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)


	if(rk_step .EQ.0) then

	do k=1, n_eqn

	    do j = 1, NY
	
		  do i = 1, NX
	
		    cons(i,j,k,1) = cons(i,j,k,0) + dt*residual(i,j,k)
	
		  enddo
	
		enddo
	enddo

	elseif(rk_step .EQ.1) then

	do k=1, n_eqn
	    
	    do j = 1, NY
		
		  do i = 1, NX
		
		    cons(i,j,k,2) = (3.0d0/4.0d0)*cons(i,j,k,0) + (1.0/4.0d0)*(cons(i,j,k,1) + dt*residual(i,j,k))
		
		  enddo
		
		enddo
	enddo

	else

	do k=1, n_eqn
	    
	    do j = 1, NY
		
		  do i = 1, NX
		
		    cons(i,j,k,0) = (1.0d0/3.0d0)*cons(i,j,k,0) + (2.0d0/3.0d0)*(cons(i,j,k,2) + dt*residual(i,j,k))
		
		  enddo
		
		enddo
	
	enddo

	endif
	
	rk = MOD(rk_step +1, 3)

			do i = 1, NX
				do j = 1, NY

			    density(i,j)		= cons(i,j,1,rk)
		       	u_vel(i,j)			= cons(i,j,2,rk)/density(i,j)
		        v_vel(i,j)			= cons(i,j,3,rk)/density(i,j)
		        uv 					= u_vel(i,j)**2 + v_vel(i,j)**2
			    pressure(i,j)		= (gamma-1.0)*(cons(i,j,4,rk)-0.5*cons(i,j,1,rk)*uv)

				enddo
			enddo

	end subroutine rungekutta



 	!***********************************************************************
	!*****                      Flux in X-direction                    *****
	!***********************************************************************


 	subroutine FX(density, u_vel, v_vel, pressure, residual, gamma,NX,NY,ghostp)
 	implicit none


 	integer     		:: i, j, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4
 	double precision    ::	 dx, dy, gamma



	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: primitive(-ghostp:NX+ghostp,4) 
	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)


	double precision		:: consl(4),consr(4),fright(4), fleft(4)
	double precision		:: priml(-ghostp:NX+ghostp,4),primr(-ghostp:NX+ghostp,4)
	double precision		:: flux_half(-ghostp:NX+ghostp,4)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average,sound_left,sound_right
	double precision		:: enthalpyleft, enthalpyright
	double precision		:: mx,my
	double precision		:: sqrt_rho_R,sqrt_rho_L,rho,divisor
	double precision 		:: SL,SR,SP,EL,ER

    

	common /grid/ dx, dy


	mx=1.0d0;my=0.0d0
	

	do j = 1, NY

		do i = -ghostp, NX + ghostp

			    primitive(i,1) = density(i,j)
			    primitive(i,2) = u_vel(i,j)
			    primitive(i,3) = v_vel(i,j)
			    primitive(i,4) = pressure(i,j)


		enddo
	  
	  call HOCUS(NX, primitive, priml,primr,mx,my,ghostp,n_eqn)


	  do i =0,NX



	  		sound_left		    = (gamma*priml(i,4)/priml(i,1))**0.5d0
	  		sound_right		  	= (gamma*primr(i,4)/primr(i,1))**0.5d0



	  		! ****************************************************************************

			consl(1)			= priml(i,1)
			consl(2)			= priml(i,1)*priml(i,2)
			consl(3)			= priml(i,1)*priml(i,3)
			consl(4)			= priml(i,4)/(gamma-1.0d0) + 0.5d0*priml(i,1)*(priml(i,2)**2.0d0+priml(i,3)**2.0d0)

			consr(1)			= primr(i,1)
			consr(2)			= primr(i,1)*primr(i,2)
			consr(3)			= primr(i,1)*primr(i,3)
			consr(4)			= primr(i,4)/(gamma-1.0d0) + 0.5d0*primr(i,1)*(primr(i,2)**2.0d0+primr(i,3)**2.0d0)

			enthalpyleft	   = (priml(i,4) + consl(4))/priml(i,1)
			enthalpyright    = (primr(i,4) + consr(4))/primr(i,1) 

	  		fleft(1)	   = consl(2)
			fleft(2)	   = consl(2) * priml(i,2) + priml(i,4)
			fleft(3)	   = consl(2) * priml(i,3) 
			fleft(4)	   = priml(i,2) * (consl(4) + priml(i,4))

			fright(1)	   = consr(2)
			fright(2)	   = consr(2) * primr(i,2)+ primr(i,4)
			fright(3)	   = consr(2) * primr(i,3) 
			fright(4)	   = primr(i,2) * (consr(4) + primr(i,4))

		sqrt_rho_L = dsqrt(priml(i,1))
		sqrt_rho_R = dsqrt(primr(i,1))

		rho 	   = dsqrt(primr(i,1)/priml(i,1))*priml(i,1)
		divisor	   = 1.0d0/(sqrt_rho_R+sqrt_rho_L)

		u_average 		= (  (priml(i,2)*sqrt_rho_L)    + (   primr(i,2)*sqrt_rho_R))*divisor
		v_average 		= (  (priml(i,3)*sqrt_rho_L)    + (   primr(i,3)*sqrt_rho_R))*divisor
		enthalpy_average= (  (enthalpyleft*sqrt_rho_L) + (enthalpyright*sqrt_rho_R))*divisor
		uv_average 		= 	0.5d0 * (u_average**2.0d0 + v_average**2.0d0)
		sound_average 	= 	dsqrt((gamma - 1.0d0) * (enthalpy_average - uv_average))


	    SL 			= MIN(priml(i,2) -sound_left , u_average-sound_average)
		SR 			= MAX(primr(i,2) +sound_right, u_average+sound_average)



		SP = (primr(i,4) - priml(i,4) + (priml(i,1)*priml(i,2)*(SL- priml(i,2))& 
			     -primr(i,1)*primr(i,2)*(SR- primr(i,2))) )/&
			     (priml(i,1)*(SL-priml(i,2)) - primr(i,1)*(SR-  primr(i,2)) )


		if(SL.gt.0.0d0) then

			do k= 1,n_eqn
				flux_half(i,k) =fleft(k)
			enddo
			
		else if (SL.le. 0.0d0 .and. 0.0d0 .lt. SP ) then 

		EL = (consl(4)/priml(i,1)) + (SP-priml(i,2))*(SP + (priml(i,4)/(priml(i,1)*(SL-priml(i,2)))) )

		flux_half(i,1) =fleft(1) + SL* (( (priml(i,1))*(SL-priml(i,2))/(SL-SP)) * 1.0d0       -consl(1))
		flux_half(i,2) =fleft(2) + SL* (( (priml(i,1))*(SL-priml(i,2))/(SL-SP)) * SP          -consl(2))
		flux_half(i,3) =fleft(3) + SL* (( (priml(i,1))*(SL-priml(i,2))/(SL-SP)) * priml(i,3)  -consl(3))
		flux_half(i,4) =fleft(4) + SL* (( (priml(i,1))*(SL-priml(i,2))/(SL-SP)) * EL          -consl(4))
    		
    	else if (SP.le.0.0d0 .and. 0.0d0 .le. SR) then
      	
      	ER = (consr(4)/primr(i,1)) + (SP-primr(i,2))*(SP + (primr(i,4)/(primr(i,1)*(SR-primr(i,2)))) )

		flux_half(i,1) =fright(1) + SR* (((primr(i,1))*(SR-primr(i,2))/(SR-SP)) * 1.0d0        -consr(1))
		flux_half(i,2) =fright(2) + SR* (((primr(i,1))*(SR-primr(i,2))/(SR-SP)) * SP           -consr(2))
		flux_half(i,3) =fright(3) + SR* (((primr(i,1))*(SR-primr(i,2))/(SR-SP)) * primr(i,3)   -consr(3))
		flux_half(i,4) =fright(4) + SR* (((primr(i,1))*(SR-primr(i,2))/(SR-SP)) * ER           -consr(4))

		elseif(SR .lt. 0.0d0) then

			do k= 1,n_eqn
				flux_half(i,k) =fright(k)
			enddo

		end if


		enddo


	  do M = 1, n_eqn
	    do i = 1, NX
		  residual(i,j,M) = -(flux_half(i,M) - flux_half(i-1,M))/dx
		enddo
	  enddo
	
	enddo
 	

 	end subroutine FX

 	!***********************************************************************
	!*****                      Flux in Y-direction                    *****
	!***********************************************************************

	subroutine GY(density, u_vel, v_vel, pressure, residual, gamma,NX,NY,ghostp)
	implicit none

 	integer     		::  i, j, M, NX, NY, ghostp,k
 	integer, parameter	::  n_eqn =4
 	double precision    ::	dx, dy, gamma

	double precision		:: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp), pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp), v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
	double precision		:: primitive(-ghostp:NY+ghostp,4)

	double precision		:: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,4)
	double precision		:: consl(4),consr(4),fright(4), fleft(4)
	double precision		:: priml(-ghostp:NY+ghostp,4),primr(-ghostp:NY+ghostp,4)
	double precision		:: flux_half(-ghostp:NY+ghostp,4)
	double precision   		:: enthalpy_average, uv_average, sound_average, u_average, v_average, sound_right,sound_left	

	double precision		:: enthalpyleft, enthalpyright
	double precision		:: mx,my
	double precision		:: sqrt_rho_R,sqrt_rho_L,rho,divisor
	double precision 		:: SL,SR,SP,EL,ER



	common /grid/ dx, dy


	mx=0.0d0;my=1.0d0
	

	do j = 1, NX

	 		do i = -ghostp, NY + ghostp

			   	primitive(i,1) = density(j,i)
			    primitive(i,2) = u_vel(j,i)
			    primitive(i,3) = v_vel(j,i)
			    primitive(i,4) = pressure(j,i)

			enddo	  

! Call interpolation

	  call HOCUS(NY, primitive, priml,primr,mx,my,ghostp,n_eqn)

! Calculations for the Riemann solver

	  do i =0, NY


	  		sound_left		    = (gamma*priml(i,4)/priml(i,1))**0.5d0
	  		sound_right		  	= (gamma*primr(i,4)/primr(i,1))**0.5d0



! ***********v******************************************************************
			consl(1)			= priml(i,1)
			consl(2)			= priml(i,1)*priml(i,2)
			consl(3)			= priml(i,1)*priml(i,3)
			consl(4)			= priml(i,4)/(gamma-1.0) + 0.5d0*priml(i,1)*(priml(i,2)**2.0d0+priml(i,3)**2.0d0)

			consr(1)			= primr(i,1)
			consr(2)			= primr(i,1)*primr(i,2)
			consr(3)			= primr(i,1)*primr(i,3)
			consr(4)			= primr(i,4)/(gamma-1.0) + 0.5d0*primr(i,1)*(primr(i,2)**2.0d0+primr(i,3)**2.0d0)


			enthalpyleft	   = (priml(i,4) + consl(4))/priml(i,1)
			enthalpyright    = (primr(i,4) + consr(4))/primr(i,1) 

	  		fleft(1)	   = consl(3)
			fleft(2)	   = consl(2) * priml(i,3)
			fleft(3)	   = consl(3) * priml(i,3) + priml(i,4)
			fleft(4)	   = consl(3)*enthalpyleft

			fright(1)	   = consr(3)
			fright(2)	   = consr(2) * primr(i,3)
			fright(3)	   = consr(3) * primr(i,3) + primr(i,4)
			fright(4)	   = consr(3)*enthalpyright



	   	sqrt_rho_L = dsqrt(priml(i,1))
		sqrt_rho_R = dsqrt(primr(i,1))

		rho 	   = dsqrt(primr(i,1)/priml(i,1))*priml(i,1)
		divisor	   = 1.0d0/(sqrt_rho_R+sqrt_rho_L)

		u_average 		= (  (priml(i,2)*sqrt_rho_L)    + (   primr(i,2)*sqrt_rho_R))*divisor
		v_average 		= (  (priml(i,3)*sqrt_rho_L)    + (   primr(i,3)*sqrt_rho_R))*divisor
		enthalpy_average= (  (enthalpyleft*sqrt_rho_L) + (enthalpyright*sqrt_rho_R))*divisor
		uv_average 		= 	0.5d0 * (u_average**2.0d0 + v_average**2.0d0)
		sound_average 	= 	dsqrt((gamma - 1.0d0) * (enthalpy_average - uv_average))		

	   SL 			= MIN(priml(i,3) -sound_left , v_average-sound_average)
	   SR 			= MAX(primr(i,3) +sound_right, v_average+sound_average)

		SP = (primr(i,4) - priml(i,4) + (priml(i,1)*priml(i,3)*(SL- priml(i,3))& 
			     -primr(i,1)*primr(i,3)*(SR- primr(i,3))) )/&
			     (priml(i,1)*(SL-priml(i,3)) - primr(i,1)*(SR-  primr(i,3)) )


			if(SL.gt.0.0d0) then

			do k= 1,n_eqn
				flux_half(i,k) =fleft(k)
			enddo


		else if (SL.le. 0.0d0 .and. 0.0d0 .lt. SP ) then 

		EL = (consl(4)/priml(i,1)) + (SP-priml(i,3))*(SP + (priml(i,4)/(priml(i,1)*(SL-priml(i,3)))) )

		flux_half(i,1) =fleft(1) + SL* (( (priml(i,1))*(SL-priml(i,3))/(SL-SP)) * 1.0d0       -consl(1))
		flux_half(i,2) =fleft(2) + SL* (( (priml(i,1))*(SL-priml(i,3))/(SL-SP)) * priml(i,2)  -consl(2))
		flux_half(i,3) =fleft(3) + SL* (( (priml(i,1))*(SL-priml(i,3))/(SL-SP)) * SP 	      -consl(3))
		flux_half(i,4) =fleft(4) + SL* (( (priml(i,1))*(SL-priml(i,3))/(SL-SP)) * EL          -consl(4))
    		
    	else if (SP.le.0.0d0 .and. 0.0d0 .le. SR) then
      	
      	ER = (consr(4)/primr(i,1)) + (SP-primr(i,3))*(SP + (primr(i,4)/(primr(i,1)*(SR-primr(i,3)))) )

		flux_half(i,1) =fright(1) + SR* (((primr(i,1))*(SR-primr(i,3))/(SR-SP)) * 1.0d0        -consr(1))
		flux_half(i,2) =fright(2) + SR* (((primr(i,1))*(SR-primr(i,3))/(SR-SP)) * primr(i,2)   -consr(2))
		flux_half(i,3) =fright(3) + SR* (((primr(i,1))*(SR-primr(i,3))/(SR-SP)) * SP 		   -consr(3))
		flux_half(i,4) =fright(4) + SR* (((primr(i,1))*(SR-primr(i,3))/(SR-SP)) * ER           -consr(4))
      
      		elseif(SR .lt. 0.0d0) then

			do k= 1,n_eqn
				flux_half(i,k) =fright(k)
			enddo
		
			end if

		enddo


	  do M = 1, n_eqn
	    do i = 1, NY
	      residual(j,i,M) = residual(j,i,M) - (flux_half(i,M) - flux_half(i-1,M))/dy
		enddo
	  enddo
	enddo


	end subroutine GY
	

	!***********************************************************************
	!*****                      HOCUS interpolation	                   *****
	!***********************************************************************
	subroutine HOCUS (NS, un, ulnew,urnew,mx,my,ghostp,n_eqn)

	implicit none
	integer				:: ix,NS, ghostp, n_eqn,i,k,total
	double precision	:: un(-ghostp:NS+ghostp,n_eqn)
	double precision	:: ulnew(-ghostp:NS+ghostp,n_eqn),urnew(-ghostp:NS+ghostp,n_eqn)
	double precision 	:: mx,my, lx, ly
	double precision 	:: v(-2:3,n_eqn),charstencil(-2:2)
    double precision 	:: lefteigen(4,4), righteigen(4,4)
    double precision, parameter	:: epsilon = 1.0d-40, gamma = 1.4d0
    double precision	:: a(0:NS), b(0:NS), c1(0:NS), r(0:NS),temp(0:NS)


   	double precision	:: ulnew_weno(-ghostp:NS+ghostp,n_eqn),urnew_weno(-ghostp:NS+ghostp,n_eqn)
	double precision	:: ulnew_compact(-1:NS,n_eqn),urnew_compact(-1:NS,n_eqn),ul(n_eqn),ur(n_eqn)
  	double precision	:: tbv_weno(-ghostp:NS+ghostp,n_eqn),tbv_compact(-ghostp:NS+ghostp,n_eqn)
  	double precision    :: sqrt_rho, divisor,rho,u_vel,v_vel,p,c,sqrt_rho_L,sqrt_rho_R

  	double precision, PARAMETER ::  B2 = 4.0d0/3.0d0, MP5=7.0d0, EPSM=1.0d-40
	double precision :: VOR,VMP,DJM1,DJ,DJP1,DM4JPH,DM4JMH,VUL,VAV,VMD,VLC,VMIN,VMAX, MINMOD2,MINMOD4,l2norm

	lx = -my ; ly = mx;



	do ix =-1,NS+1


		sqrt_rho_L = sqrt(un(ix+0,1))
		sqrt_rho_R = sqrt(un(ix+1,1))

		rho 	   = sqrt(un(ix+1,1)/un(ix,1))*un(ix,1)

		divisor	   = 1.0d0/(sqrt_rho_R+sqrt_rho_L)

		u_vel 	= ((un(ix,2)*sqrt_rho_L) + (un(ix+1,2)*sqrt_rho_R))*divisor
		v_vel 	= ((un(ix,3)*sqrt_rho_L) + (un(ix+1,3)*sqrt_rho_R))*divisor
		p 		= ((un(ix,4)*sqrt_rho_L) + (un(ix+1,4)*sqrt_rho_R))*divisor
		c       = sqrt(gamma*p/rho)


        righteigen(1,1)=	1.0d0;	 	righteigen(1,2)= 0.0d0;	righteigen(1,3)=	1.0d0;		righteigen(1,4)=	1.0d0;
        righteigen(2,1)=-mx*c/rho;	 	righteigen(2,2)=lx/rho;	righteigen(2,3)=	0.0d0;		righteigen(2,4)= mx*c/rho;
        righteigen(3,1)=-my*c/rho; 		righteigen(3,2)=ly/rho;	righteigen(3,3)= 	0.0d0;		righteigen(3,4)= my*c/rho;
        righteigen(4,1)= c**2.0d0;		righteigen(4,2)= 0.0d0; righteigen(4,3)= 	0.0d0;		righteigen(4,4)= c**2.0d0;

        
        lefteigen(1,1)=0.0d0;	lefteigen(1,2)=-mx*rho*0.5d0/c;	lefteigen(1,3)=-my*rho*0.5d0/c 	; lefteigen(1,4)= 0.5d0/c**2.0d0;
        lefteigen(2,1)=0.0d0;	lefteigen(2,2)=		lx*rho    ;	lefteigen(2,3)=ly*rho 			; lefteigen(2,4)=		   0.0d0;
        lefteigen(3,1)=1.0d0;	lefteigen(3,2)=		0.0d0 	  ;	lefteigen(3,3)=0.0d0  			; lefteigen(3,4)=-1.0d0/c**2.0d0;
        lefteigen(4,1)=0.0d0;	lefteigen(4,2)=mx*rho*0.5d0/c ;	lefteigen(4,3)=my*rho*0.5d0/c   ; lefteigen(4,4)= 0.5d0/c**2.0d0; 
 


        do i=-2,3
		
			v(i,:) =matmul(lefteigen,un(i+ix,:))

		enddo

		do i=1,n_eqn
			    charstencil=v(-2:2,i)

				VOR = (1.0d0/60.0d0)*(2.0d0*charstencil(-2)-13.0d0*charstencil(-1)+47.0d0*charstencil(0)+27.0d0*charstencil(1)-3.0d0*charstencil(2))
				
				VMP = charstencil(0) + MINMOD2(charstencil(1)-charstencil(0),MP5*(charstencil(0)-charstencil(-1)))
				
				if ((VOR-charstencil(0))*(VOR-VMP)<EPSM) then
				   ul(i)=VOR
				else

				   DJM1 = charstencil(-2)-2.0d0*charstencil(-1)+charstencil(0)
				   DJ   = charstencil(-1)-2.0d0*charstencil(0) +charstencil(1)
				   DJP1 = charstencil(0) -2.0d0*charstencil(1) +charstencil(2)

				   DM4JPH = MINMOD4(4.0d0*DJ-DJP1,4.0d0*DJP1-DJ,DJ,DJP1)
				   DM4JMH = MINMOD4(4.0d0*DJ-DJM1,4.0d0*DJM1-DJ,DJ,DJM1)

				   VUL = charstencil(0) + MP5*(charstencil(0)-charstencil(-1))
				   VAV = 0.5d0*(charstencil(0)+charstencil(1))
				   VMD = VAV - 0.5d0*DM4JPH
				   VLC = charstencil(0) + 0.5d0*(charstencil(0)-charstencil(-1)) + B2*DM4JMH

				   VMIN = MAX(MIN(charstencil(0),charstencil(1),VMD),MIN(charstencil(0),VUL,VLC))
				   VMAX = MIN(MAX(charstencil(0),charstencil(1),VMD),MAX(charstencil(0),VUL,VLC))

				   ul(i) = VOR + MINMOD2(VMIN-VOR,VMAX-VOR)
				endif

				charstencil=v(-1:3,i)

				VOR = (1.0d0/60.0d0)*(-3.0d0*charstencil(-2)+27.0d0*charstencil(-1)+47.0d0*charstencil(0)-13.0d0*charstencil(1)+2.0d0*charstencil(2))
				VMP = charstencil(0) + MINMOD2(charstencil(-1)-charstencil(0),MP5*(charstencil(0)-charstencil(1)))
				
				if ((VOR-charstencil(0))*(VOR-VMP)<EPSM) then
				   ur(i)=VOR
				else
				   
				    DJP1 = charstencil(+2)-2.0d0*charstencil(+1)+charstencil(0)
				   DJ   = charstencil(+1)-2.0d0*charstencil(0) +charstencil(-1)
				   DJM1 = charstencil(0) -2.0d0*charstencil(-1) +charstencil(-2)

				   DM4JPH = MINMOD4(4.0d0*DJ-DJP1,4.0d0*DJP1-DJ,DJ,DJP1)
				   DM4JMH = MINMOD4(4.0d0*DJ-DJM1,4.0d0*DJM1-DJ,DJ,DJM1)

				   VUL = charstencil(0) + MP5*(charstencil(0)-charstencil(+1))
				   VAV = 0.5d0*(charstencil(0)+charstencil(-1))
				   VMD = VAV - 0.5d0*DM4JMH

				   VLC = charstencil(0) + 0.5d0*(charstencil(0)-charstencil(+1)) + B2*DM4JPH
				   
				   VMIN = MAX(MIN(charstencil(0),charstencil(-1),VMD),MIN(charstencil(0),VUL,VLC))
				   VMAX = MIN(MAX(charstencil(0),charstencil(-1),VMD),MAX(charstencil(0),VUL,VLC))

			 	   ur(i) = VOR + MINMOD2(VMIN-VOR,VMAX-VOR)


				endif
        enddo

        urnew_weno(ix,:)=matmul(righteigen,ur)
        ulnew_weno(ix,:)=matmul(righteigen,ul)
	enddo
!*****************************!*****************************!*****************************!*****************************!*****************************

 	  ulnew_compact(-1,:) = ulnew_weno(-1,:)
      urnew_compact(-1,:) = urnew_weno(-1,:)

      ulnew_compact(NS+1,:) = ulnew_weno(NS+1,:)
      urnew_compact(NS+1,:) = urnew_weno(NS+1,:)

 	do k = 1,n_eqn
	    	
	    do i = 0, NS
			      a(i) = 1.0d0/2.0d0
			      b(i) = 1.0d0
			      c1(i) = 1.0d0/6.0d0
	   		 end do
	    
	    i = 0  ; a(i) = 0.0d+00
	    i = NS  ; c1(i) = 0.0d+00
	    total = NS+1
	    
	    r(0)		= (1.0d0/18.0d0)*un(-1,k)+(19.0d0/18.0d0)*un(0,k)+(10.0d0/18.0d0)*un(1,k)-(1.0d0/2.0d0)*ulnew_weno(-1,k)
	    r(NS)		= (1.0d0/18.0d0)*un(NS-1,k)+(19.0d0/18.0d0)*un(NS,k)+(10.0d0/18.0d0)*un(NS+1,k)-(1.0d0/6.0d0)*ulnew_weno(NS+1,k)
	    
	    do i=1,NS-1
	    r(i)		=(1.0d0/18.0d0)*un(i-1,k)+(19.0d0/18.0d0)*un(i,k)+(10.0d0/18.0d0)*un(i+1,k)
	    enddo
	    
	    call tridiag(a,b,c1,r,temp,total)
	        
	    ulnew_compact(0:NS,k)=temp(0:NS) 

		    do i = 0, NS
			      a(i) = 1.0d0/6.0d0
			      b(i) = 1.0d0
			      c1(i) = 1.0d0/2.0d0
		    end do

	    i = 0  ; a(i) = 0.0d+00
	    i = NS  ; c1(i) = 0.0d+00

		r(0)		= (10.0d0/18.0d0)*un(0,k)+(19.0d0/18.0d0)*un(1,k)+(1.0d0/18.0d0)*un(2,k)-(1.0d0/6.0d0)*urnew_weno(-1,k)
	    r(NS)		= (10.0d0/18.0d0)*un(NS,k)+(19.0d0/18.0d0)*un(NS+1,k)+(1.0d0/18.0d0)*un(NS+2,k)-(1.0d0/2.0d0)*urnew_weno(NS+1,k)
	   
	    do i=1,NS-1
	    r(i)		= (10.0d0/18.0d0)*un(i,k)+(19.0d0/18.0d0)*un(i+1,k)+(1.0d0/18.0d0)*un(i+2,k)
	    enddo
	    
	    call tridiag(a,b,c1,r,temp,total)
	        
	    urnew_compact(0:NS,k)=temp(0:NS)

	    ulnew(0:NS,k)   = (ulnew_compact(0:NS,k)+urnew_compact(0:NS,k))/2.0d0
	    urnew(0:NS,k)   = (urnew_compact(0:NS,k)+ulnew_compact(0:NS,k))/2.0d0

	enddo


		do k=1,n_eqn
	    	do i=1,NS

	    	tbv_compact(i,k)  = abs(ulnew_compact(i,k)-urnew_compact(i,k))+abs(ulnew_compact(i-1,k)-urnew_compact(i-1,k))

	    	tbv_weno(i,k)  	  = abs(ulnew_weno(i,k)-urnew_weno(i,k))+abs(ulnew_weno(i-1,k)-urnew_weno(i-1,k))
	    	enddo
	    enddo



		    do k=1,n_eqn
		    	do i=1,NS
		    		if(tbv_weno(i,k) .lt. tbv_compact(i,k)) then

		    		ulnew(i-2,k)   = ulnew_weno(i-2,k)
					urnew(i-2,k)   = urnew_weno(i-2,k)

					ulnew(i-1,k)   = ulnew_weno(i-1,k)
					urnew(i-1,k)   = urnew_weno(i-1,k)

					ulnew(i,k)   = ulnew_weno(i,k)
					urnew(i,k)   = urnew_weno(i,k)

					ulnew(i+1,k)   = ulnew_weno(i+1,k)
					urnew(i+1,k)   = urnew_weno(i+1,k)

					endif
		    	enddo
		    enddo


	
	end subroutine HOCUS

  
double precision function minmod2(X,Y)
double precision, INTENT(IN) :: X,Y
minmod2 = 0.5d0*(sign(1.0d0,X)+sign(1.0d0,Y))*min(ABS(X),ABS(Y))
end function minmod2


double precision function minmod4(W,X,Y,Z)
double precision, INTENT(IN) :: W,X,Y,Z

minmod4 = 0.125d0*(sign(1.0d0,W)+sign(1.0d0,X))* &
          ABS( (sign(1.0d0,W)+sign(1.0d0,Y))*(sign(1.0d0,W)+SIGN(1.0d0,Z)) )* &
          min(ABS(W),ABS(X),ABS(Y),ABS(Z))

end function minmod4




  !***********************************************************************
  !*****   Tridiagonal solver for compact scheme                 *****
  !***********************************************************************

  subroutine tridiag(l, d, cons, r, x, n)

    
  ! local variables are not implicit by default
  !
      implicit none

  ! input/output arguments
  !
      integer                   , intent(in)  :: n
      double precision, dimension(n), intent(in)  :: l, d, cons, r
      double precision, dimension(n), intent(out) :: x

  ! local variables
  !
      integer                        :: i, j
      double precision               :: t
      double precision, dimension(n) :: g
  !
  !-------------------------------------------------------------------------------
  !
  ! decomposition and forward substitution
  ! 
      t    = d(1)
      x(1) = r(1) / t      
      do i = 2, n
        j  = i - 1
        g(i) = cons(j) / t
        t    = d(i) - l(i) * g(i)

        if (t == 0.0d0) then
          write (*,*)"algebra::tridiag", "solution failed!"
          stop
          ! This is for gfortran --- a weird bug that doesn't happen with ifort
          ! t=1.0d0-40
          ! write(*,*) t,i
        end if
        x(i) = (r(i) - l(i) * x(j)) / t
      end do

  ! backsubstitution
  !
      do i = n - 1, 1, -1
        j    = i + 1
        x(i) = x(i) - g(j) * x(j)
      end do

  !-------------------------------------------------------------------------------
  !
  end subroutine tridiag





