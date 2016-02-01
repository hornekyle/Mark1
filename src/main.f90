! Draw a bounding box in pymol
!
! http://www.pymolwiki.org/index.php/DrawBoundingBox

program main_prg
	use kinds_mod
	use units_mod
	use settings_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_thermo
		!! I/O unit for thermo report output
	
	call setupSim()
	call runSim()
	call endSim()
	
contains

	subroutine setupSim
		open(file='mark1.xyz',newunit=iou_xyz)
		open(file='mark1.thermo',newunit=iou_thermo)
		
		call initialize()
		
		enableLennardJones = .true.
		call setThermostat(.false.,T0,10.0_wp*dt)
		call setBarostat(.false.,P0,5.0E10_wp*dt)
		call buildSystem(convert(5.181_wp,'A','m'),[7,7,7],T0) !5.26_wp
		
		call doBox()
		call writeStepXYZ(iou_xyz)
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
	end subroutine setupSim

	subroutine runSim
		integer::k
		
		do k=0,N_steps
			
			if(mod(k,skip_thermo)==0) call thermoReport(k)
			if(mod(k,skip_dump  )==0) call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			
			call velocityVerlet(dt)
			call doBox()
			
! 			if(k==Ns/2) call setThermostat(.false.)
			
		end do
	end subroutine runSim

	subroutine endSim
		close(iou_xyz)
		close(iou_thermo)
	end subroutine endSim

	subroutine thermoReport(k)
		integer,intent(in)::k
		integer,save::c = 0
		character(128)::buf
		
		if(mod(c,20)==0) then
			write(buf,'(1A5,6A15)') 'k [#]','t [ps]','T [K]','TE [eV]','KE [eV]','PE [eV]','P [bar]'
			write(stdout,'(2A,1I4)') colorize(trim(buf),[5,5,0]),' Nc = ',nint(averageNeighbors())
			if(c==0) write(iou_thermo,'(1A)') '#'//trim(buf)
		end if
		write(stdout,'(1I5,8G15.7)') k,convert(t,'s','ps'),temperature(),convert(E(),'J','eV'), &
			& convert(KE(),'J','eV'),convert(PE(),'J','eV'),convert(pressure(),'Pa','bar'),Pepsilon,convert(mean(box),'m','A')
		
		write(iou_thermo,'(1I5,6G25.15)') k,convert(t,'s','ps'),temperature(),convert(E(),'J','eV'), &
			& convert(KE(),'J','eV'),convert(PE(),'J','eV'),convert(pressure(),'Pa','bar')
		
		c = c+1
	end subroutine thermoReport

end program main_prg 

