program nanotube
	implicit none

	real*4,parameter	:: pi=3.141592654
	integer			:: n,m
	real*4			:: b
	integer			:: ncells
	integer			:: hcd
	integer			:: l,dr
	real*4			:: rm,radius
	integer			:: nc
	real*4,allocatable	:: x(:),y(:),z(:),u(:)	
	integer			:: p1,p2
	real*4			:: t
	real*4			:: phi
	real*4			:: alpha
	real*4			:: h
	integer			:: i
	real*4			:: length
	character*6		:: output

!Input from user
10	write(0,*) '===> Chiral indexes (two integers):'
	read(*,*,err=10) n,m

	if(n<m)then
		write(0,*) '   ERROR: n must be greater than or equal to m!'
		
		goto 10
	elseif(n<1)then
		write(0,*) '   ERROR: n must be an integer greater than 1!'
		
		goto 10
	elseif(m<0)then
		write(0,*) '   ERROR: m must be an integer greater than or equal to zero!'
		
		goto 10
	end if

20	write(0,*) '===> Bond length (in Angstroms):'
	read(*,*,err=20) b

	if(b<=0.)then
		write(0,*) '   ERROR: Bond length must be greater than zero!'
		
		goto 20
	end if

!Get the highest common divisor of both n and m
	if(m==0)then
		hcd=n
	else
		do hcd=m,1,-1
			if(mod(n,hcd)==0)then
				if(mod(m,hcd)==0)then
					exit
				end if
			end if
		end do
	end if 

!Get l
	l=m*m+n*m+n*n

!Get |R|
	rm=b*sqrt(3.*real(l))

!Get the tube radius
	radius=rm/(2.*pi)

!Get dr
	if(mod((n-m),(3*hcd))==0)then
		dr=3*hcd
	else
		dr=hcd
	end if

!Get the tube length in the Z direction
	length=(3.*b*sqrt(real(l)))/real(dr)

	write(0,*) '   Unit cell length: ',length,' Angstroms.'

30	write(0,*) '===> Number of unit cells:'
	read(*,*,err=30) ncells
	
	if(ncells<1)then
		write(0,*) '   ERROR: Number of cells must be greater than or equal to 1!'
	end if
	
	write(0,*) '   Tube length: ',ncells*length,' Angstroms.'
	write(0,*) '   Tube radius: ',radius,' Angstroms.'

40	write(0,*) '===> Output format (xyz or lammps):'
	read(*,*,err=40) output
	
	if(output/='xyz'.and.output/='lammps')then
		write(0,*) '   ERROR: Format not allowed!'
		
		goto 40
	end if

!Get the number of two-atom unit cells in the translational cell
	nc=2*l/dr

!Allocate the size for the used arrays
	allocate(x(2*nc*ncells))
	allocate(y(2*nc*ncells))
	allocate(z(2*nc*ncells))
	allocate(u(2*nc*ncells))

!Get phi
	phi=(pi*real(m+n))/real(l)

!Get t
	t=b*real(n-m)/(2.*sqrt(real(l)))

!Get p1 and p2
	do p1=0,n
		if(mod((hcd+p1*m),n)==0)then
			p2=(hcd+p1*m)/n

			exit
		end if
	end do

!Get alpha
	alpha=pi*(m*(2.*p2+p1)+n*(2.*p1+p2))/real(l)

!Get h
	h=(3.*hcd*b)/(2.*sqrt(real(l)))

!Set the first atom position in the two-atom unit cell
	x(1)=radius
	y(1)=0.
	z(1)=0.
	u(1)=0.

!Set the second atom position in the two-atom unit cell
	x(2)=radius*cos(phi)
	y(2)=radius*sin(phi)
	z(2)=t
	u(2)=u(1)+phi

!Write the file header and the positions of the two-atom unit cell
	if(output=='xyz')then
		write(6,*) 2*nc*ncells
		write(6,*) 'NANOTUBE: d=',2.*radius,'; l=',ncells*length,'; b=',b
        	write(6,*) 'C',x(1),y(1),z(1)
        	write(6,*) 'C',x(2),y(2),z(2)
	elseif(output=='lammps')then
		write(6,*) '# Ready to be used with the read_data command in a LAMMPS script file'
		write(6,*)
		write(6,*) 2*nc*ncells,' atoms'
		write(6,*)
        	write(6,*) '1 atom types'
		write(6,*)
		write(6,*) -radius-5.0,radius+5.0,' xlo xhi'
		write(6,*) -radius-5.0,radius+5.0,' ylo yhi'
		write(6,*) 0.0,ncells*length,' zlo zhi'
		write(6,*)
		write(6,*) 'Atoms'
		write(6,*)
		write(6,*) 1,'1',x(1),y(1),z(1)
        	write(6,*) 2,'1',x(2),y(2),z(2)
	
	end if

!Complete the tube helical motif
	do i=3,2*hcd
		x(i)=radius*cos(u(i-2)+(2.*pi)/real(hcd))
		y(i)=radius*sin(u(i-2)+(2.*pi)/real(hcd))
		z(i)=z(i-2)
		u(i)=u(i-2)+(2.*pi)/real(hcd)

		if(output=='xyz')then
			write(6,*) 'C',x(i),y(i),z(i)
		elseif(output=='lammps')then
			write(6,*) i,'1',x(i),y(i),z(i)
		end if
	end do

!Complete the translational unit cell
	do i=2*hcd+1,2*nc*ncells
                x(i)=radius*cos(u(i-(2*hcd))+alpha)
                y(i)=radius*sin(u(i-(2*hcd))+alpha)
                z(i)=z(i-(2*hcd))+h
                u(i)=u(i-(2*hcd))+alpha

		if(output=='xyz')then
			write(6,*) 'C',x(i),y(i),z(i)
		elseif(output=='lammps')then
			write(6,*) i,'1',x(i),y(i),z(i)
		end if
	end do

	stop

100	stop 'Chiral indexes must be integers!'
200	stop 'Bond length must be a real number!'
300	stop 'Number of cells must be an integer!'
end program nanotube
