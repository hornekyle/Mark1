module units_mod
	use kinds_mod
	use utilities_mod
	implicit none
	private
	
	integer,parameter::Nu = 25
	
	character(16),dimension(Nu),parameter::names = [ character(16):: &
		& 'm','mm','um','nm','pm','fm','A', &
		& 's','ms','us','ns','ps','fs', &
		& 'J','eV', &
		& 'kg','u', &
		& 'K', 'C', &
		& 'Pa', 'bar', &
		& 'm/s', 'A/ps', &
		& 'W/m2', 'eV/ps/A2'&
		& ]
	
	logical::isSetup = .false.
	
	real(wp),dimension(Nu,Nu)::cf,co
	
	public::convert
	
contains

	function convert(v,iu,ou) result(o)
		real(wp),intent(in)::v
			!! Value to convert
		character(*),intent(in)::iu
			!! Units of value
		character(*),intent(in)::ou
			!! Units for return value
		real(wp)::o
		
		integer::ki,ko
		
		if(.not.isSetup) call setup()
		
		ki = getIndex(iu)
		ko = getIndex(ou)
		
		if(ki<0 .or. ko<0) then
			write(stdout,'(1A)') colorize('Units Error: Cannot find specified units.',[5,0,0])
			if(ki<0) write(stdout,'(1A)') colorize('   Unit: ',[5,5,0])//trim(adjustl(iu))
			if(ko<0) write(stdout,'(1A)') colorize('   Unit: ',[5,5,0])//trim(adjustl(ou))
			stop 1
		end if
		
		o = v*cf(ki,ko)+co(ki,ko)
	end function convert

	subroutine setup()
		integer::i,j,k
		
		!= Conversion Factors =!
		
		cf = -1.0_wp
		
		forall(k=1:Nu) cf(k,k) = 1.0_wp
		cf( getIndex('m' ) , getIndex('A'  ) ) = 1.0E10_wp
		cf( getIndex('m' ) , getIndex('fm' ) ) = 1.0E15_wp
		cf( getIndex('m' ) , getIndex('pm' ) ) = 1.0E12_wp
		cf( getIndex('m' ) , getIndex('nm' ) ) = 1.0E9_wp
		cf( getIndex('m' ) , getIndex('um' ) ) = 1.0E6_wp
		cf( getIndex('m' ) , getIndex('mm' ) ) = 1.0E3_wp
		
		cf( getIndex('s' ) , getIndex('fs' ) ) = 1.0E15_wp
		cf( getIndex('s' ) , getIndex('ps' ) ) = 1.0E12_wp
		cf( getIndex('s' ) , getIndex('ns' ) ) = 1.0E9_wp
		cf( getIndex('s' ) , getIndex('us' ) ) = 1.0E6_wp
		cf( getIndex('s' ) , getIndex('ms' ) ) = 1.0E3_wp
		
		cf( getIndex('J' ) , getIndex('eV' ) ) = 6.24150636309E18_wp
		cf( getIndex('kg') , getIndex('u'  ) ) = 6.02213665168E26_wp
		
		cf( getIndex('K' ) , getIndex('C'  ) ) = 1.0_wp
		cf( getIndex('Pa') , getIndex('bar') ) = 1.0E-5_wp
		
		cf( getIndex('m/s' ) , getIndex('A/ps' ) ) = 1.0E-2_wp
		
		cf( getIndex('W/m2') , getIndex('eV/ps/A2') ) = 6.24150636309E-14
		
		! Assert inverse conversions
		do j=1,Nu
			do i=1,Nu
				if(cf(i,j)<0.0_wp) cf(i,j) = 1.0_wp/cf(j,i)
			end do
		end do
		
		!= Conversion offsets =!
		
		co = 0.0_wp
		
		co( getIndex('K' ) , getIndex('C' ) ) = -273.15_wp
		! Assert reverse offsets
		do j=1,Nu
			do i=1,Nu
				if(abs(co(i,j))<2.0_wp**5*epsilon(1.0_wp)) co(i,j) = -co(j,i)
			end do
		end do
		
		isSetup = .true.
	end subroutine setup

	function getIndex(un) result(k)
		character(*),intent(in)::un
		integer::k
		integer::i
		
		k = -1
		do i=1,Nu
			if(names(i)==un) k = i
		end do
	end function getIndex

end module units_mod