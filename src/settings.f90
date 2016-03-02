module settings_mod
	use kinds_mod
	use units_mod
	implicit none
	
	!=======================!
	!= Universal Constants =!
	!=======================!
	
	real(wp),parameter::kB = 1.3806504E-23_wp!1.3806488E-23_wp
		!! Boltzmann constant in SI units
	
	!=========!
	!= Types =!
	!=========!
	
	type::potential_t
		real(wp)::cutoff
			!! Cutoff radius
		real(wp)::skin
			!! Skin beyond cutoff for neighbor lists
		real(wp),dimension(:),allocatable::coeffs
			!! Potential coefficients
	end type
	
	type::nosehoover_t
		real(wp)::set = 10.0_wp
			!! Set value
		real(wp)::tau = 100.0_wp
			!! Characteristic time
	end type
	
	!=============!
	!= Variables =!
	!=============!
	
	type(potential_t)::lj ! coeffs=[E0,S0]
		!! Lennard-Jones potential data
	type(nosehoover_t)::thermostat
		!! Thermostat settings
	type(nosehoover_t)::barostat
		!! Barostat settings
	
	integer::N_steps
		!! Number of total simulation steps
	integer::skip_thermo
		!! Number of steps between thermo reports
	integer::skip_dump
		!! Number of steps between xyz file dumps
	integer::skip_neighbor
		!! Number of steps between neighbor list rebuilds
	integer::skip_mullerPlathe
		!! Number of steps between muller-plathe report
	real(wp)::lattice_const
	integer, dimension(3)::latM
	
	real(wp)::T0
		!! Target temperature [K]
	real(wp)::P0
		!! Target pressure [Pa]
	real(wp)::dt
		!! Timestep [seconds]
	
contains

	subroutine initialize_parameters()
		real(wp)::E0,S0
		
		!= Potential =!
		
		E0 = kB*convert(125.7_wp,'K','K')
		S0 = convert(3.345_wp,'A','m')
		lj%cutoff = 2.5_wp*S0
		lj%skin = 0.5_wp*S0
		
		lj%coeffs = [E0,S0]
		
		!= Simulation =!
		N_steps       = 0
		skip_mullerPlathe = 1
		skip_thermo   = 1
		skip_dump     = 1
		skip_neighbor = 20
		lattice_const = 5.40_wp
		latM = [5,5,5]
		
		T0 = convert(45.0_wp,'K','K')
		P0 = convert(1.0_wp,'bar','Pa')
		dt = convert(10.0_wp,'fs','s')
		
		!= Thermostat =!
		thermostat%tau = convert(100.0_wp*dt,'s','ps')
		
		!= Barostat =!
		barostat%tau = convert(1000.0_wp*dt,'s','ps')
	end subroutine initialize_parameters

	subroutine writeLammpsVars(fn)
		character(*),intent(in)::fn
		integer::iou
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A,1I10)')   'variable N_steps       equal ',N_steps
		write(iou,'(1A,1I10)')   'variable skip_thermo   equal ',skip_thermo
		write(iou,'(1A,1I10)')   'variable skip_dump     equal ',skip_dump
		write(iou,'(1A,1I10)')   'variable skip_neighbor equal ',skip_neighbor
		
		write(iou,'(1A,1EN25.5)') 'variable dt            equal ',convert(dt,'s','ps')
		write(iou,'(1A,1EN25.5)') 'variable tau_T         equal ',convert(thermostat%tau,'s','ps')
		write(iou,'(1A,1EN25.5)') 'variable tau_P         equal ',convert(barostat%tau,'s','ps')
		write(iou,'(1A,1EN25.5)') 'variable T0            equal ',convert(T0,'K','K')
		write(iou,'(1A,1EN25.5)') 'variable P0            equal ',convert(P0,'Pa','bar')
		
		write(iou,'(1A,1EN25.5)') 'variable cutoff        equal ',convert(lj%cutoff,'m','A')
		write(iou,'(1A,1EN25.5)') 'variable skin          equal ',convert(lj%skin,'m','A')
		
		close(iou)
	end subroutine writeLammpsVars

end module settings_mod
