module integrate_mod
	use kinds_mod
	use system_mod
	implicit none
	private
	
	logical::doThermostat = .false.
	logical::doBarostat = .false.
	
	public::setThermostat
	public::setBarostat
	public::velocityVerlet
	public::leapFrog
	public::doBox
	public::rnemd
	
contains

	subroutine setThermostat(state,T,tau)
		!! Turns on/off damping parameter "eta" in the Integrator
		logical,intent(in)::state
		real(wp),intent(in),optional::T,tau
		
		if(state .and. present(T) .and. present(tau) ) then
			thermostat%set = T
			thermostat%tau = tau
			Teta = 0.0_wp
		else if(state) then
			write(stdout,'(1A)') colorize('Error: Called thermostate(true) without T and tau.',[5,0,0])
			stop 1
		else
			Teta = 0.0_wp
		end if
		
		doThermostat = state
	end subroutine setThermostat

	subroutine setBarostat(state,P,tau)
		!! Turns on/off damping parameter "eta" in the Integrator
		logical,intent(in)::state
		real(wp),intent(in),optional::P,tau
		
		if(state .and. present(P) .and. present(tau) ) then
			barostat%set = P
			barostat%tau = tau
			Pepsilon = 0.0_wp
		else if(state) then
			write(stdout,'(1A)') colorize('Error: Called barostate(true) without P and tau.',[5,0,0])
			stop 1
		else
			Pepsilon = 0.0_wp
		end if
		
		doBarostat = state
	end subroutine setBarostat

	subroutine velocityVerlet(dt)
		!! Velocity-Verlet integration
		real(wp),intent(in)::dt
		
		real(wp),dimension(3)::d
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r =  atoms(k)%r+d+Pepsilon*atoms(k)%r
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-(Teta+Pepsilon)*atoms(k)%v
			atoms(k)%f = -delV(k)
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		if(doThermostat) Teta = Teta+DetaDt()*dt
		if(doBarostat) Pepsilon = Pepsilon+DepsilonDt()*dt
		box = box+Pepsilon*box
		t  = t+dt
		ts = ts+1
	end subroutine velocityVerlet

	subroutine leapFrog(dt)
		!! Leap frog integration
		real(wp),intent(in)::dt
		
		real(wp),dimension(3)::d,ao
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r = atoms(k)%r+d+Pepsilon*atoms(k)%r
		end do
		do k=1,size(atoms)
			ao = atoms(k)%a
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-(Teta+Pepsilon)*atoms(k)%v
			atoms(k)%v = atoms(k)%v+0.5_wp*(ao+atoms(k)%a)*dt
		end do
		
		if(doThermostat) Teta = Teta+DetaDt()*dt
		if(doBarostat) Pepsilon = Pepsilon+DepsilonDt()*dt
		box = box+Pepsilon*box
		t  = t+dt
		ts = ts+1
	end subroutine leapFrog

	subroutine doBox
		!! Returns moving atoms into the simulation box
		integer::k
		
		forall(k=1:3)
			where(atoms(:)%r(k)>box(k)) atoms(:)%r(k) = atoms(:)%r(k)-box(k)
			where(atoms(:)%r(k)<0.0_wp) atoms(:)%r(k) = atoms(:)%r(k)+box(k)
		end forall
	end subroutine doBox

	pure function DetaDt() result(o)
		!! Calculates damping parameter("eta") change over time
		real(wp)::o
		
		o = (1.0_wp/thermostat%tau**2)*(temperature()/thermostat%set-1.0_wp)
	end function DetaDt

	function DepsilonDt() result(o)
		!! Calculates damping parameter("eta") change over time.
		real(wp)::o
		
		o = (1.0_wp/barostat%tau**2)*(pressure()/barostat%set-1.0_wp)
	end function DepsilonDt
	
	subroutine rnemd (k)
		integer, intent(in)::k
		integer:: i, j
		integer,dimension(:), allocatable::l
		real(wp):: abc
		
		abc = real(latM(3)*lattice_const/N_slabs, wp)
		
		!! for each step do 1, size(atoms), write atoms into regions    +
		!! do 1, size(regions) calculate av. temp                       +
		!! for region 1 and mid get hot and cold and swap
		!! run and watch
		
		do j=1, N_slabs
			l = regionList(j*abc - abc, j*abc)
			regions(k+1)%temps(j) = listTemp(l)
		end do
		
	end subroutine rnemd
	
	subroutine swapV(k)
		integer,  intent(in)::k
		integer:: i,j, h, c
		
		!h = selectHot(
		
		!end do
		
	end subroutine swapV

end module integrate_mod
