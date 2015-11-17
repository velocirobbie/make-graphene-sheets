!Robert Sinclair 23.09.2015
!Build a lammps posiitional file for a graphene sheet
program makegraphene
implicit none

integer::				i,j,k,l,a
integer::				order
integer::				n_hex_wide,n_bond,n_angle,n_dihedral,n_impropers
integer::				n_H,n_C,n_total
double precision::			bond_CC,cos_CC,sin_CC
double precision::			bond_CH,p_HH,c_HH
double precision::			n,m,pi

character(9)::				model
character(40)::				filename
character(1)::				hydrogen

double precision,dimension(:,:),allocatable::xyz
integer,dimension(:,:),allocatable::bonds,angles,dihedrals

pi=3.14159265359
bond_CC = 1.42
bond_CH = 1.09
p_HH = 2*cos(pi/6)*bond_CC
c_HH = bond_CC+2*sin(pi/6)*bond_CH

cos_CC = cos(pi/6)*bond_CC
sin_CC = sin(pi/6)*bond_CC

order = 0
n_hex_wide = 0
n_C = 0
n_H = 0
n_bond = 0
n_angle = 0
n_dihedral = 0
n_impropers = 0

write(*,*)'===BUILD A GRAPHENE SHEET==='
write(*,*)'Shape? (hexagon, square or rectangle)'
read(*,*)model
if (model.eq.'hexagon'.or.model.eq.'square') then
	write(*,*)'How many A accross?'
	read(*,*)n
	else
	write(*,*)'How many A long and wide?'
	read(*,*)n,m
end if
write(*,*)'Hydrogen on edges? (y/n)'
read(*,*)hydrogen


if (model.eq.'hexagon') then
!===========Hexagon=============

if (n.lt.bond_CC*sqrt(3.0)) then
	write(*,*)'ERROR1: dimensions too small'
	stop
end if

!find largest hex that fits into width specified
!number of hexagons on width must be odd, order 
n_hex_wide=n/(bond_CC*sqrt(3.0))
if (mod(order,2).eq.0) n_hex_wide=n_hex_wide-1
order=n_hex_wide/2 +1

!count carbons atoms in hexagon
do i=1,n_hex_wide
	if (mod(i,2).eq.0) then
		else
		n_C=n_C+i*6
	end if
end do
write(*,*)'Number of carbons in sheet = ',n_C

if (hydrogen.eq.'y') n_H=order*6
n_total=n_H+n_C

allocate(xyz(n_total+1,3))
xyz=0

!calculate coordinates for the first third of the hexagon shape
a=1
do i=1,order
	do j=1,i
		xyz(a,1)=cos_CC*(i-1)
		xyz(a,2)=-(sin_CC+bond_CC)*(i-1)+(j-1)*2*(bond_CC+sin_CC)
		xyz(a+1,1)=xyz(a,1)
		xyz(a+1,2)=xyz(a,2)+bond_CC
		a=a+2
	end do
end do

!calculate the coordinates for the central part of the hexagon
do i=1,n_hex_wide
	do j=1,order
		if (mod(i,2).eq.0) then
			xyz(a,1)=cos_CC*(i-1+order)
			xyz(a,2)=-(sin_CC+bond_CC)*(order-1)+(j-1)*2*(bond_CC+sin_CC)
			xyz(a+1,1)=xyz(a,1)
			xyz(a+1,2)=xyz(a,2)+bond_CC
		else
			xyz(a,1)=cos_CC*(i-1+order)
			xyz(a,2)=-((sin_CC+bond_CC)*(order-1)+sin_CC)+(j-1)*2*(bond_CC+sin_CC)
			xyz(a+1,1)=xyz(a,1)
			xyz(a+1,2)=xyz(a,2)+bond_CC+2*sin_CC
		end if
		a=a+2
	end do
end do

!Calculate the coordinates for the final third of the atoms
do i=order,1,-1
        do j=1,i
                xyz(a,1)=cos_CC*(order+n_hex_wide+(order-i))
                xyz(a,2)=-(sin_CC+bond_CC)*(i-1)+(j-1)*2*(bond_CC+sin_CC)
                xyz(a+1,1)=xyz(a,1)
                xyz(a+1,2)=xyz(a,2)+bond_CC
                a=a+2
        end do
end do
write(*,*)'Number of Carbons = ',a-1

!Calculate the coordinates of the hydrogens
if (hydrogen.eq.'y')then
	xyz(a,1)=-cos(pi/6)*bond_CH
	xyz(a,2)=-sin(pi/6)*bond_CH
	do i=1,6
		do j=1,order-1
			xyz(a+1,1)=xyz(a,1)+p_HH*sin((i-0.5)*pi/3)
			xyz(a+1,2)=xyz(a,2)-p_HH*cos((i-0.5)*pi/3)
			a=a+1
		end do
		xyz(a+1,2)=xyz(a,2)-c_HH*cos((i)*pi/3)
		a=a+1
	end do
	write(*,*)'Number of hydrogens = ',n_H,a-1-n_C
end if

!populate a list of bonds
do i=1,order
	n_bond=n_bond+i*6+(i-1)*12
end do

do i=1,n_C
	do j=1,n_total
		



filename='hexagon.xyz'
open(10,file=filename)
write(10,*)n_total
write(10,*)

do i=1,n_C
	write(10,*)'C  ',xyz(i,1),xyz(i,2),xyz(i,3)
end do

do i=1,n_H
	write(10,*)'H  ',xyz(n_C+i,1),xyz(n_C+i,2),xyz(n_C+i,3)
end do

close(10)
!=============Hexagon=============
end if

stop
end program



