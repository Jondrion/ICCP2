module polymer
  use plplot

  implicit none

  real(8), parameter :: PI = 4 * atan(1d0)
  real(8), allocatable, private :: AvWeight(:,:)

  private

  public polymerType
  type polymerType
    private
   
    integer, private :: Length
    integer, private :: Population, PopulationLim
    integer, private :: NumberAngles
    real(8), private :: temperature, epsilon, sigma
    real(8), allocatable, private :: Position(:,:)

  contains
    private
     
    procedure, public :: init, create, reset, destroy
    procedure, private :: calc_Angles, calc_Energy, calc_Weights, calc_Gyradius, calc_EToE
    procedure, private :: choose_Angle, cumsum
    procedure, public :: get_Position, get_Length
    procedure, public :: plot
    procedure, private :: plot_polymer
     
  end type

contains
  
  ! -- initialize the polymer: allocate position array of certain length and create the first two beads
  subroutine init(this, Length, temperature, epsilon, sigma)
  
    class(polymerType) :: this
    integer, intent(in) :: Length
    real(8), intent(in) :: temperature, epsilon, sigma
    integer :: i
    
    this%Length=Length
    allocate (this%Position(2,this%Length))

    this%temperature=temperature

    this%epsilon=epsilon
    this%sigma=sigma

    ! -- create inital two beads
    this%Position(1,1)=0
    this%Position(2,1)=0
    this%Position(1,2)=0
    this%Position(2,2)=1

  
    ! -- allocate averageweight and initialize
    allocate (AvWeight(2,this%Length))
    do i=1,this%Length
      AvWeight(1,i)=1._8
      AvWeight(2,i)=1._8
    end do

    this%NumberAngles=6

    this%PopulationLim=1

    this%Population=1

    ! -- creating the actual polymers, these might be cloned or killed, creating a population of maximum size this%PopulationLim
    open (100, file="polymerdata.txt", ACTION="write", STATUS="unknown", Position="append")
      
    
  end subroutine


  ! -- creating the polymer by adding enough beads until length is reached, this routine is recursive
  recursive subroutine create(this, PolWeight, Number)

    class(polymerType) :: this
    real(8) :: Angle(this%NumberAngles)
    real(8) :: Energy(this%NumberAngles)
    integer, intent(in) :: Number
    real(8), intent(inout) :: PolWeight
    real(8) :: NewWeight
    real(8) :: Weights(this%NumberAngles)
    real(8) :: Norm_Weights(this%NumberAngles)
    real(8) :: SumWeights
    real(8) :: New_Angle
    real(8) :: R
    real(8) :: EndtoEnd
    real(8) :: Gyradius
    real(8) :: Uplim, Lowlim

    ! -- Calculate weights w_j^L and their product W^L
    call this%calc_Angles(Angle,this%NumberAngles)
    call this%calc_Energy(Angle,Number,Energy)
    call this%calc_Weights(Energy,Weights)

    SumWeights = SUM(Weights)

    Norm_Weights = Weights/SumWeights
    
    ! -- choose angle for new bead

    call this%choose_Angle(Angle, Norm_Weights, New_Angle)

    ! -- Add bead number L
    this%Position(1,Number)=this%Position(1,Number-1)+SIN(New_Angle)
    this%Position(2,Number)=this%Position(2,Number-1)+COS(New_Angle)

    ! -- Update PolWeight with correction for PREM
    PolWeight=PolWeight*SumWeights/ (0.75_8 * this%NumberAngles)

    ! -- Update AvWeight
    AvWeight(1,Number)= AvWeight(1,Number) + PolWeight
    AvWeight(2,Number)=AvWeight(2,Number)+1

    ! -- write End to end and gyradius to a file
    call this%calc_EToE(EndtoEnd, Number)
    call this%calc_Gyradius(Gyradius, Number)
    write (100, "(I10,E18.5E4,F20.4,F20.4)" ) Number, PolWeight, EndtoEnd, Gyradius
    
    ! -- Determining Upper and lower limits for cloning and killing
    Lowlim=(1.2_8) * AvWeight(1,Number)/AvWeight(1,3)
    Uplim=(2.2_8) * AvWeight(1,Number)/AvWeight(1,3)

    Lowlim=0
    !Uplim=0
    ! -- recursive part
    if ( Number < this%Length ) then
      ! -- kill if necessar
      if ( PolWeight < Lowlim ) then
        call RANDOM_NUMBER(R)
        if ( R < 0.5 ) then
         
          this%Population = this%Population
          NewWeight = 2 * PolWeight
          call this%create(NewWeight, Number+1)
        else
        

        end if
      ! -- clone if necessary
      else if ( PolWeight > Uplim .and. this%Population < this%PopulationLim ) then
!         print *, "N : ", Number, "Population: ", this%population
!         print *, "AvWeight1: ", AvWeight
        this%Population = this%Population+1
        NewWeight = 0.5 * PolWeight
        call this%create(NewWeight, Number+1)
!         print *, "N : ", Number, "Population: ", this%population
!         print *, "AvWeight2: ", AvWeight
        NewWeight = 0.5 * PolWeight
        call this%create(NewWeight, Number+1)
      else
        
        call this%create(PolWeight, Number+1)
      end if
    else
      ! -- End of polymer, nothing happens now.
      
    end if

  end subroutine


  subroutine destroy(this)
    class(polymerType) :: this

    close (100)
    
    deallocate(this%Position)
    deallocate(AvWeight)
  end subroutine

  subroutine reset(this)
    class(polymerType) :: this

    
    this%Population=1
    deallocate(this%Position)
    allocate (this%Position(2,this%Length))
    this%Position(1,1)=0
    this%Position(2,1)=0
    this%Position(1,2)=0
    this%Position(2,2)=1
    
  end subroutine

  subroutine calc_Angles(this, Angle, N)

    class(polymerType) :: this
    integer, intent(in) :: N
    integer :: i
    real(8), intent(out) :: Angle(N)
    real(8) :: Interval
    real(8) :: Offset

    Interval = 2 * PI / N
    call RANDOM_NUMBER(Offset)
    Offset = Offset * Interval
    Angle(:) = [(Offset + ((i-1)*Interval), i=1,N)]

  end subroutine


  subroutine calc_Energy(this, Angle, N, E)
    
    class(polymerType) :: this
    real(8), intent(in) :: Angle(this%NumberAngles)
    integer, intent(in) :: N
    real(8) :: Ri(2)
    real(8) :: Rsqi
    real(8) :: rm2,rm6,rm12
    real(8) :: new_pos(2,this%NumberAngles)
    real(8), intent(out) :: E(this%NumberAngles)

    integer :: i,j

    E=0
    
    ! -- calculate position given the angle
    new_pos(1,:) = this%Position(1,N-1) + SIN(Angle(:))
    new_pos(2,:) = this%Position(2,N-1) + COS(Angle(:))

    ! -- loop over all angles
    do j = 1, this%NumberAngles, 1
      ! -- loop over all existing elements
      do i = 1, N-1, 1
        Ri = new_pos(:,j) - this%Position(:,i)
        Rsqi = dot_product(Ri,Ri)
        
        rm2 = 1.d0/Rsqi
        rm6 = rm2**3
        rm12 = rm6**2
        ! Calculate energy with sigma=0.8 and epsilon=0.25
        E(j) = E(j) + 4._8 * this%epsilon * ( this%sigma**12 * rm12 - this%sigma**6 * rm6 )
      end do
    end do

  end subroutine


  subroutine calc_Weights(this, Energy, Weights)
    
    class(polymerType) :: this
    real(8), intent(out) :: Weights(this%NumberAngles)
    real(8), intent(in) :: Energy(this%NumberAngles)

    Weights(:) = EXP(-Energy(:)/this%temperature)

  end subroutine


  subroutine choose_Angle(this, Angles, Norm_Weights, output_Angle)

    class(polymerType) :: this
    real(8), intent(inout) :: Norm_Weights(this%NumberAngles)
    real(8), intent(in) :: Angles(this%NumberAngles)
    real(8), intent(out) :: output_Angle
    real(8) :: Roulette
    integer :: angle_number

    call RANDOM_NUMBER(Roulette)
    call this%cumsum(Norm_Weights)
    
    angle_number=1
    do while (Norm_Weights(angle_number)<Roulette)
      angle_number=angle_number+1
    end do
      
    output_Angle = Angles(angle_number)

  end subroutine


  subroutine cumsum(this, Array)
    class(polymerType) :: this
    real(8), intent(inout) :: Array(this%NumberAngles)
    integer :: i
    
    do i=1, this%NumberAngles-1
      Array(i+1)=Array(i+1)+Array(i)
    end do

  end subroutine


  subroutine calc_Gyradius(this, Gyradius, Number)

    class(polymerType) :: this
    
    real(8), intent(out) :: Gyradius
    integer, intent(in) :: Number
    real(8) :: Rmean(2)
    integer :: i

    Rmean = 0._8
    Gyradius = 0._8

    do i = 1, Number
      Rmean = Rmean + this%Position(:,i)
    end do
    Rmean = Rmean/Number
  
    do i = 1, Number
      Gyradius = Gyradius + dot_product((this%Position(:,i)-Rmean),(this%Position(:,i)-Rmean))
    end do
    Gyradius = Gyradius/Number

  end subroutine calc_Gyradius


  subroutine calc_EToE(this, EToE, Number)

    class(polymerType) :: this
    integer, intent(in) :: Number    
    real(8), intent(out) :: EToE

    EToE = (dot_product((this%Position(:,Number)-this%Position(:,1)),(this%Position(:,Number)-this%Position(:,1))))
    
  end subroutine calc_EToE


  subroutine get_Position(this, Position)

    class(polymerType) :: this
    real(8), allocatable, intent(out) :: Position(:,:)

    Position=this%Position

  end subroutine

  subroutine get_Length(this, Length)

    class(polymerType) :: this
    real(8), intent(out) :: Length

    Length=this%Length

  end subroutine


  subroutine plot(this)

    class(polymerType) :: this

    call plsdev("xcairo")

    ! gnuplot color scheme
    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)      ! red
    call plscol0(2, 0, 255, 0)      ! green
    call plscol0(3, 0, 0, 255)      ! blue
    call plscol0(4, 255, 0, 255)    ! magenta
    call plscol0(5, 0, 255, 255)    ! cyan
    call plscol0(6, 255, 255, 0)    ! yellow
    call plscol0(7, 0, 0, 0)        ! black
    call plscol0(8, 255, 76, 0)     ! orange
    call plscol0(9, 128, 128, 128)  ! gray

    call plinit()

    call this%plot_polymer

    call plend()

  end subroutine
  

  subroutine plot_polymer(this)

    class(polymerType) :: this
    real(8) :: min, max,border
    Character(50)::  weightstring

    min = minval(this%Position)
    max = maxval(this%Position)
  
    border=(max-min)/20
    min = min - border
    max = max + border

    Write( weightstring, '(A7,ES10.3)' )  "Weight: "

    call plcol0(7)
    call plenv(min, max, min, max, 0, 0)
    call pllab("x", "y", weightstring)
  
    call plcol0(1)
    call plline(this%Position(1,:),this%Position(2,:))
    
    call plcol0(2)
    call plpoin(this%Position(1,:),this%Position(2,:),4)

  end subroutine


end module polymer