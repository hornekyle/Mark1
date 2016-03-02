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
	use lmpIntegrator_mod
	implicit none
	
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_thermo
		!! I/O unit for thermo report output
	integer::iou_lammps=1
		!! I/O unit to read lammps dump file
	
	call setupSim()
	call runSim()
	call endSim()
	
contains

	subroutine setupSim
		open(file='mark1.xyz',newunit=iou_xyz)
		open(file='mark1.thermo',newunit=iou_thermo)
						
		call initialize_parameters()
					
		enableLennardJones = .true.
		call setThermostat(.false.,T0,10.0_wp*dt)
		call setBarostat(.false.,P0, 5.0E10_wp*dt)
		call buildSystem(convert(lattice_const,'A','m'),[2,2,2],T0)
		
		call doBox()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
	end subroutine setupSim
		
	subroutine runSim
		integer::k
		do k=0,N_steps
			call mullerPlathe()
			if(mod(k,skip_mullerPlathe)==0) call mullerplatheReport(k)
			!if(mod(k,skip_thermo)==0) call thermoReport(k)
			if(mod(k,skip_dump)==0) call writeStepXYZ(iou_xyz)
			!if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			
			call velocityVerlet(dt)
			call doBox()
 		end do
	end subroutine runSim

	subroutine endSim
		close(iou_xyz)
		close(iou_thermo)
	end subroutine endSim
	
	subroutine thermoReport(k)
		integer,intent(in)::k
		integer,save::c = 0
		real(wp), dimension(3)::o
		integer::i,j
		
		o = heatflux()
		
		if (mod(c,50)==0) then
			write(*,*)
			write(stdout,'(1X,1A17,1A5, 3F6.2)')"\x1B[96mSystem size:","\x1B[97m", [(convert(box(i), 'm', 'A'), i=1,3)]
			write(stdout,'(1X,1A21, 1A5, 1I5)') "\x1B[96mNumber of atoms:","\x1B[97m", size(atoms)
			write(stdout,'(1X,1A29, 1A5, 1I3)') "\x1B[96mAv. number of neighbors:", "\x1B[97m", nint(averageNeighbors())
			write(*,*)			
			write(*, '(1X, 1A5, 1A9, 6A13, 1A22)') 'Step', 'Temp', 'KE()', 'PE()', 'TotEng', 'Jx', 'Jy', 'Jz', 'Fnorm'
		end if
		
		write(stdout,'(1X,1I3,1F11.3,3F13.6,3ES17.6, 1F13.6)') k, temperature(), &
			& convert(KE(),'J','eV'), &
			& convert(PE(),'J','eV'), &
			& convert(E(),'J','eV'),  &
			& (convert(o(i),'W/m2','eV/ps/A2')/product(box),i=1,3), &
			& convert(fnorm(),'N','eV/A')
		c = c+1
	end subroutine thermoReport
	
	subroutine mullerplatheReport(k)
		integer,intent(in)::k
		integer::i, j, c
		write(*,*)
		if (mod(c,32)==0) then
			write(stdout,'(1X, 1A7, 1I3)') "Step N:", k
		end if
		
		write(*,'(1X, 1A4, 3A5, 3A10, 3A10, 1A11)') 'k', 'r(x)', 'r(y)', 'r(z)', 'v(x)', 'v(y)', 'v(z)', 'KE(i)', 'KE()', &
			& 'TT(i)','Temp()'
		do i=1, size(atoms)
			write(stdout,'(1X, 1I4, 3F5.1, 3F10.4, 3F10.4, 1F10.4)') atoms(i)%atom_id, &
				& [(convert(atoms(i)%r(j),'m','A'),j=1,3)], &
				& [(convert(atoms(i)%v(j),'m/s', 'A/ps'), j=1,3)], &
				& convert(KEi(i), 'J', 'eV'), &
				& convert(KE(), 'J', 'eV'), &
				& atoms(i)%tt, temperature()
			c = c+1
		end do
			
	end subroutine mullerplatheReport

end program main_prg 
