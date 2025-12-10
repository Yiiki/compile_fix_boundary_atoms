program fix_boundary_atoms
implicit none

type config
  integer*4 :: natom
  real*8 :: AL(3,3)
  integer*4,allocatable :: tatom(:),matom(:,:)
  real*8,allocatable :: xatom(:,:)
end type config
type(config) :: ac

character*256 :: fcfg,fout,fchk
logical :: iqst

call getarg(1,fcfg)

write(fout,141) "fix.", trim(fcfg)
write(fchk,141) "chk.", trim(fout)
141 format(A4,A)
write(6,*) "out: ", trim(fout)
write(6,*) "chk: ", trim(fchk)

inquire(file=trim(fcfg),exist=iqst)
if(.not.iqst) then
  write(6,*) "file ",trim(fcfg), " not exist"
  stop
endif

call read_config(fcfg,ac)

call fix_atoms()

call write_config(fout,ac)

call chk_atoms()

call write_config(fchk,ac)

contains

subroutine chk_atoms()
implicit none
integer*4 :: iatom, matom
do iatom=1,ac%natom
matom=ac%matom(1,iatom)
select case (matom)
  case (1)  
       ac%tatom(iatom)=1
  case (0)  
       ac%tatom(iatom)=79
  case DEFAULT
       write(6,*) "wrong matom", iatom, ac%matom(:,iatom)
       stop
end select
end do
end subroutine chk_atoms

subroutine fix_atoms()
implicit none
integer*4 :: iatom,i
real*8,parameter :: atol=2.d-1 ! unit: Angstrom
real*8 :: rad(3)
logical :: boundary
do i=1,3
  rad(i)=ac%AL(i,i)*0.5d0-atol
end do
do iatom=1,ac%natom
  ac%matom(1:3,iatom)=1
  do i=1,3
  boundary=abs(dist_min(ac%xatom(i,iatom)-0.5d0))*ac%AL(i,i)>rad(i)
  if(boundary) then
    ! write(6,*) "fix atom ", iatom
    ac%matom(1:3,iatom)=0
    exit
  end if
  end do
end do
end subroutine fix_atoms

function dist_func(cfg,iatom,jatom) result(y)
implicit none
type(config) :: cfg
integer*4 :: iatom,jatom
real*8 :: y(3)
y=cfg%xatom(1:3,iatom)-cfg%xatom(1:3,jatom)
y=matmul(cfg%AL,y)
end function dist_func

function dist_min(x) result(y)
implicit none
real*8,intent(in) :: x
real*8 :: y
y=mod(mod(x,1.d0)+1.5d0,1.d0)-0.5d0
end

function dist_min3(x) result(y)
implicit none
real*8,intent(in) :: x(3)
real*8 :: y(3)
y=mod(mod(x,1.d0)+1.5d0,1.d0)-0.5d0
end


subroutine write_config(file,cfg)
implicit none
integer*4 :: new_unit,iatom
character(len=*) :: file
type(config) :: cfg
open(newunit=new_unit,file=trim(file))
rewind(new_unit)
write(new_unit,*) cfg%natom
write(new_unit,*) "LATTICE"
write(new_unit,'(1x,3F15.8)') cfg%AL(1:3,1)
write(new_unit,'(1x,3F15.8)') cfg%AL(1:3,2)
write(new_unit,'(1x,3F15.8)') cfg%AL(1:3,3)
write(new_unit,'(A)') "POSITION"
do iatom=1,cfg%natom
  write(new_unit,131) cfg%tatom(iatom),cfg%xatom(1:3,iatom),cfg%matom(1:3,iatom) 
end do
close(new_unit)
131 format(1x,I3,3F15.8,3(1x,I1))

end subroutine write_config

subroutine read_config(file,cfg)
implicit none
integer*4 :: new_unit,iatom
character(len=*) :: file
type(config) :: cfg
open(newunit=new_unit,file=trim(file))
rewind(new_unit)
read(new_unit,*) cfg%natom
read(new_unit,*) ! "LATTICE"
read(new_unit,*) cfg%AL(:,1)
read(new_unit,*) cfg%AL(:,2)
read(new_unit,*) cfg%AL(:,3)
read(new_unit,*) ! "POSITION"
allocate(cfg%tatom(cfg%natom))
allocate(cfg%matom(3,cfg%natom))
allocate(cfg%xatom(3,cfg%natom))
do iatom=1,cfg%natom
read(new_unit,*) cfg%tatom(iatom),cfg%xatom(1:3,iatom),cfg%matom(1:3,iatom)
end do
close(new_unit)
end subroutine read_config


end program fix_boundary_atoms
