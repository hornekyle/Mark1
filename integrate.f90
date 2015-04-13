module integrate_mod
	use kinds_mod
	use system_mod
	implicit none
	private
	
	logical::doThermostat = .false.
	real(wp)::eta = 0.0_wp
	real(wp)::Tset = 10.0_wp
	real(wp)::tauT = 100.0_wp
	
	public::setThermostat
	public::velocityVerlet
	public::leapFrog
	public::doBox
	
contains

	subroutine setThermostat(state,T,tau)
		logical,intent(in)::state
		real(wp),intent(in),optional::T,tau
		
		if(state .and. present(T) .and. present(tau) ) then
			Tset = T
			tauT = tau
			eta = 0.0_wp
		else if(state) then
			write(*,*) 'Error: Called thermostate(true) without T and tau!'
			stop 1
		else
			eta = 0.0_wp
		end if
		
		doThermostat = state
	end subroutine setThermostat

	subroutine velocityVerlet(dt)
		real(wp),intent(in)::dt
		integer::k
		
		forall(k=1:size(atoms))
			atoms(k)%r =  atoms(k)%r+atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-eta*atoms(k)%v
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end forall
		
		if(doThermostat) eta = eta+DetaDt()*dt
		t  = t+dt
		ts = ts+1
	end subroutine velocityVerlet

	subroutine leapFrog(dt)
		real(wp),intent(in)::dt
		real(wp),dimension(2)::ao
		integer::k
		
		do k=1,size(atoms)
			atoms(k)%r = atoms(k)%r+atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
		end do
		do k=1,size(atoms)
			ao = atoms(k)%a
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-eta*atoms(k)%v
			atoms(k)%v = atoms(k)%v+0.5_wp*(ao+atoms(k)%a)*dt
		end do
		
		if(doThermostat) eta = eta+DetaDt()*dt
		t  = t+dt
		ts = ts+1
	end subroutine leapFrog

	subroutine doBox
		integer::k
		
		forall(k=1:2)
			where(atoms(:)%r(k)>box(k)) atoms(:)%r(k) = atoms(:)%r(k)-box(k)
			where(atoms(:)%r(k)<0.0_wp) atoms(:)%r(k) = atoms(:)%r(k)+box(k)
		end forall
	end subroutine doBox

	pure function DetaDt() result(o)
		real(wp)::o
		
		o = (temperature()/Tset-1.0_wp)/tauT**2
	end function DetaDt

end module integrate_mod