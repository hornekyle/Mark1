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
			write(*,*) regions(1)%ttt(1,1)
				
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


end program main_prg 
