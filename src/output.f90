module output_mod
	use kinds_mod
	use units_mod
	use system_mod
	implicit none
	
contains

	subroutine writeStepXYZ(iou)
		integer,intent(in)::iou
		
		character(32)::buf
		integer::i,k
		
		write(buf,*) size(atoms)
		write(iou,'(1A)') trim(adjustl(buf))
		write(iou,'(1A,1I10)') 'Atoms. Timestep: ',ts
		do k=1,size(atoms)
			write(iou,'(1I4, 3F5.1,6F13.9)') atoms(k)%atom_id, &
				& [(convert(atoms(k)%r(i),'m','A'),i=1,3)], &
				& [(convert(atoms(k)%v(i),'m/s', 'A/ps'), i=1,3)], &
				& [(convert(atoms(k)%f(i), 'N', 'eV/A'), i=1,3)]			
		end do
	end subroutine writeStepXYZ
	

end module output_mod
