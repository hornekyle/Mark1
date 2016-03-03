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
		call buildSystem(lattice_const,latM,T0)
		
		call doBox()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
	end subroutine setupSim
		
	subroutine runSim
		integer::k
		do k=0,N_steps
			call sub1()
			if(mod(k,skip_mullerPlathe)==0) call mullerplatheReport(k)
			if(mod(k,skip_dump)==0) call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			!if(mod(k,skip_thermo)==0) call thermoReport(k)
		
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
		integer::i, j, p
		integer, save::c=0
		real(wp),dimension(3)::r
		real(wp)::centre, radius
		character(len=100) writeFormat1, writeFormat2
						
		radius = lj%cutoff+lj%skin
		centre = latM(1)*lattice_const/2.0_wp
		writeFormat1 = '(1X, /, 2A5, 2X, 3A6, 3A10, 3A10, 2A10)'
		writeFormat2 = '(1X, 1I5, 1I5, 2X, 3F6.1, 3F10.4, 3F10.4, 2F10.4)'
		p = 0
				
		write(stdout,'(1X, 1A9, 3F7.2)') "Box size:", [(convert(box(j), 'm', 'A'), j=1,3)]
		write(stdout,'(1X, 1A12, 1F8.4)') "Cutoff+skin:", convert((lj%cutoff+lj%skin),'m','A')
		write(stdout,'(1X, 1A16, 1I4)') "Number of atoms:", size(atoms)
		write(stdout,'(1X, 1A17, /)') "Units used: metal"
		
		if (mod(c,10000)==0) write(stdout,'(1X, 1A7, 1I3)') "Step N:", k
		
		write(*,writeFormat1)'#','id','r(x)', 'r(y)', 'r(z)', 'v(x)', 'v(y)','v(z)','KE(i)', 'KE()', 'TT(i)','Temp()','Distance'
						
		do i=1, size(atoms)
			if (abs(fn2(centre, atoms(i))) <= radius .and. (atoms(i)%r(3)>centre-radius .and. atoms(i)%r(3)<centre+radius)) then
				p = p+1
				!! add to hot list
				!! go with cold condition and add to cold list
				if (mod(p, 30)==0) write(*,writeFormat1)'#', 'id', 'r(x)', 'r(y)', 'r(z)', 'v(x)', 'v(y)', 'v(z)', 'KE(i)', &
					& 'KE()', 'TT(i)','Temp()', 'Distance'
				write(stdout, writeFormat2) p, atoms(i)%atom_id, &
					& [(convert(atoms(i)%r(j),'m','A'),j=1,3)], &
					& [(convert(atoms(i)%v(j),'m/s', 'A/ps'), j=1,3)], &
					& convert(KEi(i), 'J', 'eV'), &
					& convert(KE(), 'J', 'eV'), &
					& atoms(i)%tt, temperature(), &
					& convert(abs(fn2(centre, atoms(i))), 'm', 'A')
			end if
		end do
		write(*,*)
		call fn1(k)
		c = c+1	
	end subroutine mullerplatheReport

end program main_prg 
