module plot
use plplot

contains

  subroutine init_plot(BoxSize)
    real(8) :: BoxSize
    call plsdev("xcairo")
    call plscolbg(185,225,255)  
    call plinit()
    call plscol0 (3, 100, 0, 0);
    call pladv(0)
    call plvpor(0d0, 1d0, 0d0, 1d0)
    call plwind(-1d0, 1d0, -2d0 / 3, 4d0 / 3)
    call plw3d(1d0, 1d0, 1d0, 0._8, BoxSize, 0._8, BoxSize, 0._8, BoxSize, 45d0, -45d0)
  end subroutine init_plot


  subroutine plot_points(xyz, BoxSize)

    real(8), intent(in) :: xyz(:, :), BoxSize

    call plclear()
    call plcol0(1)
    call plot_box(0._8, BoxSize, 0._8, BoxSize, 0._8, BoxSize)
    call plcol0(3)
    call plpoin3(xyz(1, 2:), xyz(2, 2:), xyz(3, 2:), 4)
    call plcol0(3)
    call plpoin3(xyz(1, 1:1), xyz(2, 1:1), xyz(3, 1:1), 4)
    call plflush()
  end subroutine

  subroutine plot_box(xmin, xmax, ymin, ymax, zmin, zmax)
    real(8), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(8) :: x(5), y(5), z(5)

    x = (/xmin, xmax, xmax, xmin, xmin/)
    y = (/ymin, ymin, ymax, ymax, ymin/)
    z = (/zmin, zmin, zmin, zmin, zmin/)
    call plline3(x, y, z)

    x = (/xmin, xmax, xmax, xmin, xmin/)
    y = (/ymin, ymin, ymax, ymax, ymin/)
    z = (/zmax, zmax, zmax, zmax, zmax/)
    call plline3(x, y, z)

    x = (/xmin, xmax, xmax, xmin, xmin/)
    y = (/ymin, ymin, ymin, ymin, ymin/)
    z = (/zmin, zmin, zmax, zmax, zmin/)
    call plline3(x, y, z)

    x = (/xmin, xmax, xmax, xmin, xmin/)
    y = (/ymax, ymax, ymax, ymax, ymax/)
    z = (/zmin, zmin, zmax, zmax, zmin/)
    call plline3(x, y, z)
  end subroutine
  
  subroutine end_plot
    call plspause(.false.)
    call plend()
  end subroutine end_plot

   
end module plot
