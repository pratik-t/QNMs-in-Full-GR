module IO
   implicit none
   integer, parameter:: DBPR=kind(0.d0) !double precision
   real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR
   real(DBPR) :: NaN = 0.0_DBPR

contains

!-------------------------------------------------------------------------------------------------------------------------------

subroutine read_eos(ID, filename, a_data, b_data, c_data)
   use iso_fortran_env, only: real64 => real64
   implicit none
 
   character(len=*), intent(in) :: filename
   integer, intent(in) :: ID
   character(len=1000) :: line
   integer :: n_cols, ios, i
   real(real64), allocatable, intent(out) :: a_data(:), b_data(:)
   real(real64), allocatable, intent(out), optional :: c_data(:)
   real(real64) :: a_point, b_point, c_point
 
   allocate(a_data(0), b_data(0))
   if (present(c_data)) allocate(c_data(0))
 
   !--- open and read header line (CSV) ---
   open(unit=ID+7, file=filename, status='old', action='read')
   read(ID+7,'(A)', iostat=ios) line
   if (ios /= 0) then
     print *, 'Error reading header from ', trim(filename)
     stop
   end if
 
   !--- determine number of columns by counting commas ---
   n_cols = count([(line(i:i), i=1,len_trim(line))] == ',') + 1
   close(ID+7)
 
   if (n_cols < 2) then
     print *, 'Input file has fewer than 2 columns:', trim(filename)
     stop
   end if
 
   !--- reopen and skip header, then read data list-directed (commas okay) ---
   open(unit=ID+7, file=filename, status='old', action='read')
   read(ID+7, '(A)', iostat=ios) line     ! skip header row
 
   do
     if (n_cols == 3) then
       read(ID+7, *, iostat=ios) a_point, b_point, c_point
     else
       read(ID+7, *, iostat=ios) a_point, b_point
       c_point = -10.0_DBPR
     end if
     if (ios /= 0) exit
 
     a_data = [a_data, a_point]
     b_data = [b_data, b_point]
     if (present(c_data)) c_data = [c_data, c_point]
   end do
 
   close(ID+7)

 end subroutine read_eos

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine read_TOV_data(MATRIX, r_data, m_data, p_data, e_data, de_dp_data, lamb_data, phi_data, dphi_data)

      real(DBPR), allocatable, dimension(:,:), intent(in) :: MATRIX
      real(DBPR), allocatable, dimension(:), intent(out) :: r_data, m_data, p_data, e_data, de_dp_data, &
         lamb_data, phi_data, dphi_data

      r_data     = MATRIX(1,:)
      m_data     = MATRIX(2,:)
      p_data     = MATRIX(3,:)
      e_data     = MATRIX(4,:)
      de_dp_data = MATRIX(5,:)
      lamb_data  = MATRIX(6,:)
      phi_data   = MATRIX(7,:)
      dphi_data  = MATRIX(8,:)

   end subroutine read_TOV_data

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine make_filename(dir_in, dir_out, fn, out_text, input_fn, output_fn)

      character(len=*), intent(in) :: dir_in, dir_out
      character(len=*), intent(in), optional :: fn, out_text
      character(len=*), intent(out), optional ::  input_fn, output_fn
      character(len=128) :: in_path, out_path
      
      in_path = dir_in(1:len_trim(dir_in))
      
      if (present(fn)) then
         input_fn = './DATA/'//in_path(1:len_trim(in_path))//'/'//fn(1:len_trim(fn))
      end if

      ! Specifying the out_path to it where the data will be stored.
      out_path = dir_out(1:len_trim(dir_out))//'/'

      if(present(fn)) then
         output_fn = './DATA/'//out_path(1:len_trim(out_path))//fn(1:len_trim(fn)-4)//&
            '_'//out_text(1:len_trim(out_text))//'.csv'
      end if

   end subroutine make_filename

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine interpolate(xdata, ydata, xval, yval)

      real(DBPR), intent(in) :: xval
      real(DBPR), intent(in) :: xdata(:), ydata(:)
      real(DBPR), intent(out) :: yval
      integer :: data_size, start, mid, finish, range, low, high
      real(DBPR) :: slope

      data_size = size(xdata)
      start= 1
      finish= data_size
      range= finish-start
      mid= (start+finish)/2

      ! Do a binary search of data to find the bounds for given xval

      if (xval>xdata(finish) .and. xval-xdata(finish)>1E-10) then
         print*, 'Value out of range', xval, '>', xdata(finish)

      else
         do while (xdata(mid) /= xval .and. range>0)
            if (xval> xdata(mid)) then
               start= mid+1
            else
               finish= mid-1
            end if
            range= finish-start
            mid= (start+finish)/2
         end do

         if (xdata(mid-1)<xval .and. xval<xdata(mid)) then
            low= mid-1
         else
            low= mid
         end if

         high= low+1

         ! print*, 'Value', xval, 'is between: ', xdata(low), ', ', xdata(high)

         if (xdata(low)==xval) then
            yval= ydata(low)

         else
            slope = (ydata(high)-ydata(low))/(xdata(high)-xdata(low))
            yval  = ydata(low)+ slope*(xval- xdata(low))
         end if

      end if

   end subroutine interpolate

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine continuity_check(x_data, z_data, disc_start, disc_end)

      real(DBPR), allocatable, intent(inout) :: x_data(:), z_data(:)
      integer, intent(out) :: disc_start, disc_end
      integer :: i

      do i= 1, size(x_data)-1
         if (z_data(i)==0) then
            disc_start  = i
            disc_end    = i+1
            x_data(i+1) = x_data(i)+1E-10_DBPR
            exit
         else
            disc_start = -10
            disc_end   = -10
         end if
      end do

   end subroutine continuity_check

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine matinv4(A, B)
      !! Performs a direct calculation of the inverse of a 4×4 matrix.
      real(DBPR), intent(in)  :: A(4,4)   !! Matrix
      real(DBPR), intent(out) :: B(4,4)   !! Inverse matrix
      real(DBPR)              :: detinv

      ! Calculate the inverse determinant of the matrix
   detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

      ! Calculate the inverse of the matrix
   B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
   B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
   B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
   B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
   B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
   B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
   B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
   B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
   B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
   B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
   B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
   B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))

   end subroutine matinv4

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine cramers(Y1, Y2, Y3, Y4, Y5, solution)

      real(DBPR), intent(in) :: Y1(4), Y2(4), Y3(4), Y4(4), Y5(4)
      real(DBPR), intent(out) :: solution(4)
      real(DBPR) :: detY, detA1, detA2, detA3, detA4, MATY(4,4), MATA1(4,4), MATA2(4,4), MATA3(4,4), MATA4(4,4), &
         a1, a2, a3, a4

      MATY  = reshape([Y1, Y2, Y3, Y4], shape=[4,4])
      MATA1 = reshape([Y5, Y2, Y3, Y4], shape=[4,4])
      MATA2 = reshape([Y1, Y5, Y3, Y4], shape=[4,4])
      MATA3 = reshape([Y1, Y2, Y5, Y4], shape=[4,4])
      MATA4 = reshape([Y1, Y2, Y3, Y5], shape=[4,4])

      detY  = det(MATY)
      detA1 = det(MATA1)
      detA2 = det(MATA2)
      detA3 = det(MATA3)
      detA4 = det(MATA4)

      a1 = detA1/detY
      a2 = detA2/detY
      a3 = detA3/detY
      a4 = detA4/detY

      solution = [a1, a2, a3, a4]

   contains

      function det(A) result(detA)
         real(DBPR), intent(in) :: A(4,4)
         real(DBPR) :: detA

   detA = (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

      end function det

   end subroutine cramers

!-------------------------------------------------------------------------------------------------------------------------------

   function polyfit(vx, vy, d) result(poly_val)

      implicit none
      integer, intent(in) :: d
      double complex, dimension(d+1) :: poly_val
      double complex, intent(in) :: vx(:), vy(:)
      double complex, allocatable :: X(:,:), XT(:,:), XTX(:,:)
      integer :: i, j, n, lda, lwork, info
      integer, allocatable :: ipiv(:)
      double complex, allocatable :: work(:)

      n = d+1
      lda = n
      lwork = n

      allocate(ipiv(n), work(lwork), XT(n, size(vx)), X(size(vx), n), XTX(n, n))

      ! prepare the matrix
      do i = 0, d
         do j = 1, size(vx)
            X(j, i+1) = vx(j)**i
         end do
      end do

      XT = X
      call ZGETRF(n, n, XT, lda, ipiv, info)
      call ZGETRI(n, XT, lda, ipiv, work, lwork, info)
      poly_val = matmul(XT, vy)

   end function

!-------------------------------------------------------------------------------------------------------------------------------

end module IO
