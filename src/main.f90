! Draw a bounding box in pymol
!
! http://www.pymolwiki.org/index.php/DrawBoundingBox


! write: step temp1 temp2 temp2
! write: step rx1 ry1 rz1 KE1 rx2 ry2 rz2 KE2
! call swap atoms(k)
! call calculate conduction
! write: step conduction


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
		integer,dimension(:), allocatable::l
		do k=0, N_steps
			
			call rnemd(k)
			write(*,*)
			write(*,'(1X,1F10.7)') regions(k+1)%temps
			
			!call sub1()
			!call mullerplatheReport(k,1) 
				!! 1 - print all, 2 - print cold, 3 - print hot (mid) region

			!l = regionList(24.0_wp*1E-10_wp,28*1E-10_wp)
			!write(*,*) listTemp(l)
			!write(*,*) selectHot(l)
			!write(*,*) selectCold(l)
			
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
	
	subroutine mullerplatheReport(k, region)
		integer,intent(in)::k, region
		integer::i, j, p
		integer, save::c=0
		real(wp),dimension(3)::r
		real(wp)::centre, radius, rstep
		character(len=100) writeFormat1, writeFormat2
						
		radius = lj%cutoff+lj%skin
		centre = latM(3)*lattice_const/2.0_wp
		rstep = lattice_const/2.0_wp
		writeFormat1 = '(1X, /, 2A5, 3X, 3A6, 6A10)'
		writeFormat2 = '(1X, 2I5, 2X, 3F6.1, 6F10.4)'
		p = 0
				
		write(stdout,'(1X, 1A9, 3F7.2)') "Box size:", &
			& [(convert(box(j), 'm', 'A'), j=1,3)]
		write(stdout,'(1X, 1A16, 1I4)') "Number of atoms:", size(atoms)
		write(stdout,'(1X, 1A17, /)') "Units used: metal"
		
		if (mod(c,10000)==0) write(stdout,'(1X, 1A7, 1I3)') "Step N:", k
				
		if (region == 1) then
			!! print all
			write(*,*) "=====Printing whole region along Z axis====="
			write(*,writeFormat1)'#','id','r(x)', 'r(y)', 'r(z)', 'v(x)', & 
				& 'v(y)','v(z)','KE(i)', 'KE()', 'Temp(i)'
			do i=1, size(atoms)
				if 	(atoms(i)%r(3) > 24.0*1E-10_wp) then
					p = p+1
					call printRegion(i,p,writeFormat1,writeFormat2)
				end if
			end do
		
		else if (region == 2) then
			!! print COLD region
			write(*,*) "=====Printing cold region along Z axis====="
			write(*,writeFormat1)'#','id','r(x)', 'r(y)', 'r(z)', 'v(x)', & 
				& 'v(y)','v(z)','KE(i)', 'KE()', 'Temp(i)'
			do i=1,size(atoms)
				if 	(atoms(i)%r(3) < lattice_const) then
					p = p+1
					call printRegion(i,p,writeFormat1,writeFormat2)
				end if
			end do
		
		else if (region == 3) then
			!! print HOT region
			write(*,*) "=====Printing hot region along Z axis====="
			write(*,writeFormat1)'#','id','r(x)', 'r(y)', 'r(z)', 'v(x)', & 
				& 'v(y)','v(z)','KE(i)', 'KE()', 'Temp(i)'
			do i=1, size(atoms)
				if	((atoms(i)%r(3) >= centre-rstep) .and. &
					& atoms(i)%r(3) < centre+rstep) then
					p = p+1
					call printRegion(i,p,writeFormat1,writeFormat2)
				end if
			end do
		end if
		c = c+1	
	end subroutine mullerplatheReport
		
	subroutine printRegion(i,p,writeFormat1,writeFormat2)
		integer, intent(in)::i,p
		character(len=100), intent(in):: writeFormat1, writeFormat2
		real(wp)::centre, radius
		integer::j
		radius = lj%cutoff+lj%skin

		
		if (mod(p, 30)==0) write(*,writeFormat1)'#', 'id', 'r(x)', 'r(y)', &
			& 'r(z)', 'v(x)', 'v(y)', 'v(z)', 'KE(i)', 'KE()', 'Temp(i)'
		write(stdout, writeFormat2) p, atoms(i)%atom_id, &
			& [(convert(atoms(i)%r(j),'m','A'),j=1,3)], &
			& [(convert(atoms(i)%v(j),'m/s', 'A/ps'), j=1,3)], &
			& convert(KEi(i), 'J', 'eV'), &
			& convert(KE(), 'J', 'eV'), &
			& atoms(i)%tt
		
	end subroutine printRegion

end program main_prg 
