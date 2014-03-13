module Gas

  implicit none
  private

  public gasType
  type gasType
    private

    integer, public :: n
    real(8), allocatable, public :: pos(:, :), mom(:, :)
! The nest few guys are private!
    real(8), allocatable :: for(:, :)
    real(8) :: BoxSize, LargeRad, SmallRad
    integer, allocatable :: PairList(:,:)
    integer :: PairNum

  contains
    private

    procedure, public :: create, destroy, make_PairList
    procedure, public :: move_part
    procedure, public :: get_kin, get_pot, get_BoxSize
    procedure, public :: set_temp

    procedure :: calc_force !, update_pairs,

  end type

contains

  subroutine create(this, LinCell, density, temperature)
    class(gasType) :: this
    integer, intent(in) :: LinCell
    real(8) :: density, temperature, volume

    this%n = 4 * LinCell**3
    allocate (this%pos(3, this%n))
    allocate (this%mom(3, this%n))
    allocate (this%for(3, this%n))

    call this%set_temp(temperature)
    Volume = this%n/density
    this%BoxSize = Volume**(1.0_8/3)
    call init_pos(this%pos, LinCell, this%BoxSize)
    call init_mom(this%mom, this%n)
    ALLOCATE(this%PairList(this%n*(this%n-1)/2,2))
    this%LargeRad = 5.0_8
    this%SmallRad = 3.7_8
  end subroutine
  
  subroutine destroy(this)
    class(gasType) :: this
    
    deallocate(this%pos, this%mom, this%for)
  end subroutine
  
  subroutine init_pos (pos, LinCell, BoxSize)
    integer :: i, j, k, PartCount
    integer, intent(in) :: LinCell
    real(8) :: CellSize
    real(8), intent(out) :: Pos(3,4*LinCell**3)
    real(8), intent(in) :: BoxSize

    CellSize = BoxSize/LinCell
    PartCount = 0
    do i=0, LinCell - 1
      do j=0, LinCell - 1
        do k = 0, LinCell - 1
          PartCount = PartCount + 1
          Pos(:,PartCount) = [i, j, k]
          PartCount = PartCount + 1
          Pos(:,PartCount) = [DBLE(i),(j+0.5_8),  (k+0.5_8)]
          PartCount = PartCount + 1
          Pos(:,PartCount) = [(i+0.5_8), DBLE(j), (k+0.5_8)]
          PartCount = PartCount + 1
          Pos(:,PartCount) = [(i+0.5_8), (j+0.5_8), DBLE(k)]
        end do
      end do
    end do
    Pos = Pos*CellSize
  end subroutine init_pos

  subroutine init_mom(mom, n)
    integer, intent(in) :: n
    real(8), dimension (N) :: r1, r2
    real(8) :: AvMom(3), PI
    real(8), intent(out) :: mom(3,n) 
    integer :: i

    PI = 4._8*atan(1._8)
    call random_number (r1)
    call random_number (r2)
    Mom(1,:) = sqrt(-log(r1))*sin(2*PI*r2)
    Mom(2,:) = sqrt(-log(r1))*cos(2*PI*r2)
    call random_number (r1)
    call random_number (r2)
    Mom(3,:) = sqrt(-log(r1))*sin(r2)
    AvMom = SUM(Mom, 2)/N
    do i = 1, N
      Mom(:,i) = Mom(:,i) - AvMom
    end do
  end subroutine init_mom


  subroutine make_pairList(this)
    class(gasType) :: this
    integer :: i, j
    real(8) :: RelPos(3), RadiusSq
    
    this%PairNum = 0
    do i = 1, this%n-1
      do j = i+1, this%n
        RelPos = this%Pos(:,j) - this%Pos(:,i)
        RelPos = RelPos - nint(RelPos/this%BoxSize)*this%BoxSize
        RadiusSq = sum(RelPos*RelPos)
        if (RadiusSq<this%LargeRad**2) then
          this%PairNum = this%PairNum + 1
          this%PairList(this%PairNum, :) = [i,j]
        end if
      end do
    end do
    end subroutine make_pairList
   
        
    
  
  subroutine move_part(this, h)
    class(gasType) :: this
    real(8) :: h
    
    this%mom = this%mom + 0.5_8*this%for*h
    this%pos = this%pos + this%mom*h
    call this%calc_force
    this%mom = this%mom + 0.5_8*this%for*h
  end subroutine 
  
  subroutine  calc_force(this)
    class(gasType) :: this
    integer :: i, j, cnt
    real(8) :: relpos(3), r2, forceConst
!    real(8), parameter :: cutoff = 3.5_8**2
    
    this%for = 0_8
!    return
    do cnt = 1, this%PairNum
      i = this%PairList(cnt, 1)
      j = this%PairList(cnt, 2)
      relpos = this%pos(:,i) - this%pos(:,j)
      relpos = relpos - nint(relpos/this%BoxSize)*this%BoxSize
      r2 = sum(relpos*relpos)
      if (r2<this%SmallRad**2) then
        forceConst = 24*(2/r2**7-1/r2**4)
        this%for(:,i) = this%for(:,i) + forceConst*relpos
        this%for(:,j) = this%for(:,j) - forceConst*relpos
      end if
    end do
  end subroutine calc_force

  subroutine get_kin(this, Kinetic)
    class(gasType) :: this
    real(8) :: Kinetic
    
    Kinetic = 0.5_8*sum(this%mom*this%mom)
  end subroutine get_kin
  
  subroutine get_pot(this, Potential)
    class(gasType) :: this
    real(8), intent(out) :: Potential
    integer :: i, j
    real(8) :: relpos(3), r2
    real(8), parameter :: cutoff = 3.5_8**2

    Potential = 0._8
    do i = 1, this%n-1
      do j = i+1, this%n
        relpos = this%pos(:,i) - this%pos(:,j)
        relpos = relpos - nint(relpos/this%BoxSize)*this%BoxSize
        r2 = sum(relpos*relpos)
        if (r2<cutoff) then
          Potential = Potential + 4*(1/r2**6-1/r2**3)
        end if
      end do
    end do
  end subroutine get_pot  
  
  subroutine get_boxsize(this, BoxSize)
    class(gasType) :: this
    real(8), intent(out) :: BoxSize
    
    BoxSize = this%BoxSize
  end subroutine get_boxsize
    
  subroutine set_temp(this, temperature)
    class(gasType) :: this
    real(8) :: temperature, Lambda
    
    Lambda = sqrt(3*this%n*Temperature/SUM(this%Mom*this%Mom))
    this%Mom = Lambda*this%Mom
    print *, 0.5*sum(this%mom*this%mom)
  end subroutine set_temp


end module Gas

