module FULLGR
   use IO
   use COWLING
   use odepack_mod
   implicit none
contains

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine frequency(type, MATRIX, r_disc, mode, damp)

    real(DBPR), intent(in) :: r_disc
    real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
    real(DBPR), intent(out) :: mode, damp
    real(DBPR) :: initial, step, l, Mass, Radius
    real(DBPR) :: g_guess, f_guess, p1_guess, p2_guess, in_guess
    real(DBPR), allocatable :: guesses(:)
    integer, allocatable :: nodes(:)
    logical :: type(3)

      initial = 1.0_DBPR
      step    = 0.1_DBPR
      call find_guess(MATRIX, r_disc, initial, step, 4, guesses, nodes)
      guesses = guesses*0.0014766_DBPR

      if (r_disc<0 .or. (r_disc>0 .and. nodes(1)==0)) then
         g_guess  = 0.0_DBPR
         f_guess  = guesses(1)
         p1_guess = guesses(2)
         p2_guess = guesses(3)
      else if (r_disc>0 .and. nodes(2)==0) then
         g_guess  = guesses(1)
         f_guess  = guesses(2)
         p1_guess = guesses(3)
         p2_guess = guesses(4)
      else if (r_disc>0 .and. nodes(1)==1 .and. nodes(2)==2) then
         f_guess  = NaN/NaN
         p1_guess = guesses(1)
      end if

      if (type(1)) then

         if (isnan(f_guess)) then
            mode = NaN/NaN
            damp = NaN/NaN
         else
            l       = 2.0_DBPR
            Mass    = MATRIX(2,size(MATRIX,2))*1000/1.4766_DBPR
            Radius  = MATRIX(1,size(MATRIX,2))*1000
            in_guess = sqrt(((2*l*(l-1.0_DBPR))/(2*l+1.0_DBPR))* Mass/((Radius/1.4766_DBPR)**3))
            call find_freq_GR(10, MATRIX, r_disc, in_guess, 0.0_DBPR, 1E-3_DBPR, 1E-5_DBPR, mode, damp)

            if ((mode>f_guess) .or. (r_disc>0 .and. mode<g_guess)) then
               call find_freq_GR(1, MATRIX, r_disc, f_guess, g_guess, -1E-3_DBPR, 1E-4_DBPR, mode)
               in_guess = mode
               call find_freq_GR(2, MATRIX, r_disc, in_guess, 0.0_DBPR, 1E-6_DBPR, 1E-4_DBPR, mode, damp)
               if (isnan(mode)) mode = in_guess
            end if
         end if

      end if

      if (type(2)) then

         call find_freq_GR(1, MATRIX, r_disc, p1_guess, f_guess, -1E-3_DBPR, 1E-4_DBPR, mode)
         in_guess = mode
         if (isnan(mode)) then
            damp = NaN/NaN
         else
            call find_freq_GR(2, MATRIX, r_disc, in_guess, 0.0_DBPR, 1E-6_DBPR, 1E-4_DBPR, mode, damp)
            if (isnan(mode)) mode = in_guess
         end if
      end if

      if (type(3)) then
         if (r_disc<0) then
            mode = NaN/NaN
            damp = NaN/NaN
         else
            call find_freq_GR(1, MATRIX, r_disc, g_guess, 0.0_DBPR, -1E-3_DBPR, 1E-4_DBPR, mode)
            in_guess = mode
            call find_freq_GR(2, MATRIX, r_disc, in_guess, 0.0_DBPR, 1E-6_DBPR, 1E-3_DBPR, mode, damp)
            if (isnan(mode)) mode = in_guess
         end if

      end if

      mode = (mode/0.0014766_DBPR)*47.71345159236942258889E-3_DBPR
      damp = damp*1476.6_DBPR/(2.99792458E8_DBPR)

   end subroutine frequency

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine find_freq_GR(type, MATRIX, r_disc, om_guess, low_bound, h_freq, tol_in, om_eigen, damp_eigen)
      integer, intent(in) :: type
      real(DBPR), intent(in) :: om_guess, low_bound, r_disc, h_freq, tol_in
      real(DBPR), intent(in), allocatable :: MATRIX(:,:)
      real(DBPR), intent(out) :: om_eigen
      real(DBPR), intent(out), optional :: damp_eigen
      real(DBPR) :: om, om_prev, h, rtol, atol, tol, error, a_incs(3), slopes(3)
      real(DBPR), allocatable :: om_vals(:), inc_vals(:)
      integer :: iter, flag
      logical :: condition
      double complex :: i, a, b, c, b_prev, poly_vals(3), Ainc_vals(3), root1, root2, om_cmplx_prev, om_cmplx

      atol = 1E-15_DBPR
      iter = 0
      om   = om_guess

      if (type==1) then

         ! open(unit=311, file= 'omcheck.dat', status= 'replace', action= 'write')

         flag   = 0
         slopes = 0.0_DBPR
         a_incs = 0.0_DBPR
         h      = h_freq
         rtol   = tol_in
         tol    = 1E-5_DBPR

111      do

            if (om<=low_bound) then
               if (flag==0) then
                  flag   = flag+1
                  om     = om_guess
                  rtol   = 1E-8_DBPR
                  slopes = 0.0_DBPR
                  a_incs = 0.0_DBPR
                  iter   = 0
                  goto 111
               else if (flag==2) then
                  flag   = 10
                  om     = om_guess
                  h      = h_freq/10.0_DBPR
                  rtol   = 1E-10_DBPR
                  slopes = 0.0_DBPR
                  a_incs = 0.0_DBPR
                  iter   = 0
                  goto 111
               else if(flag==10) then
                  om_eigen = NaN/NaN
                  exit
               end if
            end if

            iter = iter+1

            call outoing_omega(rtol, atol, MATRIX, r_disc, om, a)

            a_incs = [a_incs(2:3), log10(abs(a))]
            slopes = [slopes(2:3), (a_incs(3)-a_incs(2))/abs(h)]

            if (flag==0 .or. flag==2) then
               condition = (a_incs(2) < a_incs(1) .and. a_incs(2) < a_incs(3))
            else
               condition = (slopes(2) < slopes(1) .and. slopes(2) < slopes(3))
            end if

            ! write(311,*) om, a_incs(3)

            ! print*, iter, flag, om, a_incs

            if (iter>2 .and. condition) then
               if (h==h_freq) then
                  om     = om-5*h
                  h      = h/10.0_DBPR
                  slopes = 0.0_DBPR
                  a_incs = 0.0_DBPR
                  rtol   = 1E-8_DBPR
                  iter   = 0
                  flag   = 2
                  goto 111
               end if
               call local_minima(1E-10_DBPR, atol, tol, MATRIX, r_disc, om, om-2*h, om_eigen)
               exit
            else
               om = om+h
            end if

         end do

         ! close(311)

      else if (type==2 .or. type==10) then

         allocate(om_vals(0), inc_vals(0))
         i     = (0.0_DBPR, 1.0_DBPR)
         h     = h_freq !keeping this of the order of the imaginary part works, less/more doesn't work.
         om    = om_guess
         rtol  = 1E-10_DBPR
         tol   = tol_in
         error = 10
         b     = 10

         do while (error>tol)

            iter          = iter + 1
            om_prev       = om
            om_cmplx_prev = om_cmplx
            b_prev = b

            call outoing_omega(rtol, atol, MATRIX, r_disc, om-h, a)
            call outoing_omega(rtol, atol, MATRIX, r_disc, om, b)
            call outoing_omega(rtol, atol, MATRIX, r_disc, om+h, c)

            error = abs(b-b_prev)

            om_vals   = [om-h, om, om+h]
            Ainc_vals = [a,b,c]
            poly_vals = polyfit(om_vals+ 0*i, Ainc_vals, 2)

            if (isnan(realpart(poly_vals(1)))) then
               ! print*, 'PROBLEM IN CURVE FITTING'
               om   = NaN/NaN
               exit
            end if

            ! ! a = poly_vals(3); b = poly_vals(2); c = poly_vals(1); root = (-b +- sqrt(b^2- 4ac))/2a
            root1 = (-poly_vals(2)+sqrt((poly_vals(2)**2)- 4*poly_vals(3)*poly_vals(1)))/(2*poly_vals(3))
            root2 = (-poly_vals(2)-sqrt((poly_vals(2)**2)- 4*poly_vals(3)*poly_vals(1)))/(2*poly_vals(3))

            ! print*, root1, root2, abs(b), error

            if (abs(om_vals(1)-abs(realpart(root1))) < abs(om_vals(1)-abs(realpart(root2)))) then
               om       = realpart(root1)
               om_cmplx = root1
            else
               om       = realpart(root2)
               om_cmplx = root2
            end if

            if ((abs(om-om_guess)>1E-3_DBPR) .and. (type.ne.10)) then
               om_cmplx = NaN/NaN
               om       = NaN/NaN
               exit
            end if

            if (om<0 .or. iter>10) then
               om_cmplx = NaN/NaN
               om       = NaN/NaN
               exit
            end if

         end do

         if(.not. isnan(om) .and. (abs(b)<1E-2_DBPR .or. type==10)) then
            om_eigen   = om
            damp_eigen = 1.0_DBPR/abs((imagpart(om_cmplx)))
         else
            om_eigen   = NaN/NaN
            damp_eigen = NaN/NaN
         end if

      end if

   end subroutine find_freq_GR

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine local_minima(rtol, atol, root_tol, MATRIX, r_disc, a, b, root)

      ! Brent's method to find local minima within an interval [a,b] (Guaranteed convergence)
      ! using function cache to avoid recomputation of function values

      real(DBPR), intent(in) :: rtol, atol, root_tol, r_disc
      real(DBPR), intent(in), allocatable :: MATRIX(:,:)
      real(DBPR), intent(out) :: root
      real(DBPR), allocatable :: cache(:)
      real(DBPR) :: invphi, ratio, a, b, x, v, w, e, d, mid, fx, fv, fw, q, p, dx, u, fu
      integer :: i, calls

      calls = 0
      allocate(cache(0))

      invphi = (1.0_DBPR+ sqrt(5.0_DBPR))/2.0_DBPR - 1.0_DBPR
      ratio = (3.0_DBPR- sqrt(5.0_DBPR))/2.0_DBPR

      x = b+invphi*(a-b)
      v = x
      w = x
      e = 0.0_DBPR
      d = 0.0_DBPR

      do i=1,100

         mid = 0.5*(a+b)

         if (abs(b-a) < root_tol) exit

         fx = func(x, cache)
         fv = func(v, cache)
         fw = func(w, cache)

         q = 2.0_DBPR*((x-v)*(fx-fw)- (x-w)*(fx-fv))
         p = -sign(1.0_DBPR,q)*((x-v)*(x-v)*(fx-fw) - (fx-fv)*(x-w)*(x-w))
         q = abs(q)

         if (q .ne. 0.0_DBPR) then
            dx = p/q
         else
            dx = 0.0_DBPR
         end if

         if ((q.ne.0) .and. (a<x+dx) .and. (x+dx<b) .and. (abs(dx)<0.5*abs(e))) then
            e = d
            d = dx
         else if (x<mid) then
            e = b-x
            d = ratio * e
         else
            e = a-x
            d = ratio * e
         end if

         u = x+d
         fu = func(u, cache)

         ! print*, i, a, b, u

         if (fu<=fx) then
            if (u<x) then
               b = x
            else
               a = x
            end if
            v = w
            w = x
            x = u
         else
            if (u<x) then
               a = u
            else
               b = u
            end if

            if (fu<=fw .or. w==x) then
               v = w
               w = u
            else if (fu<=fv .or. v==x .or. v==w) then
               v = u
            end if
         end if
      end do

      root = mid
      ! print*, 'No. of function calls = ', calls, 'Root = ', root

      !---------------------------------------------------------------------------------------------------------------------------

   contains

      function func(om_in, cache_in) result(a_inc)
         real(DBPR), intent(in) :: om_in
         real(DBPR), allocatable, intent(inout) :: cache_in(:)
         real(DBPR) :: a_inc
         real(DBPR), allocatable :: cache_mat(:,:)
         integer :: indx(1), n_rows
         double complex :: a_cmplx

         if (size(cache_in) > 0) then
            n_rows = size(cache_in) / 2
            allocate(cache_mat(n_rows, 2))
            cache_mat = reshape(cache_in, [n_rows, 2], order=[2,1])
            indx = findloc(cache_mat(:,1), om_in)
         else
            indx(1) = 0
         end if

         if (indx(1) > 0) then
            a_inc = cache_mat(indx(1),2)
         else
            calls = calls+1
            call outoing_omega(rtol, atol, MATRIX, r_disc, om_in, a_cmplx)
            a_inc = log10(abs(a_cmplx))
            cache_in = [cache_in, [om_in, a_inc]]
         end if

      end function func

   end subroutine local_minima

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine outoing_omega(rtol, atol, MATRIX, r_disc, omega, B_out)

      real(DBPR), intent(in) :: rtol, atol, r_disc, omega
      real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
      double complex, intent(out) :: B_out
      real(DBPR) :: rc, taylor_end, W_init, K_init, sol(4)
      real(DBPR), allocatable :: r(:), m(:), p(:), e(:)
      real(DBPR), allocatable :: Y1(:,:), Y2(:,:), Y3(:,:), Y4(:,:), Y5(:,:), Y_front(:,:), Y_back(:,:), final_Y(:,:)
      real(DBPR) :: step, Radius, Mass, H1_surf, K_surf, r_far
      real(DBPR), allocatable :: r_far_dat(:), Z_dat(:), dZ_dat(:)
      ! real(DBPR), allocatable :: r_grid(:), H0(:), V(:), H1(:), K(:), W(:), X(:)
      ! integer :: i

      allocate(r(0), m(0), p(0), e(0))

      r = MATRIX(1,:)
      m = MATRIX(2,:)
      p = MATRIX(3,:)
      e = MATRIX(4,:)

      W_init = 1.0_DBPR
      K_init = e(1)+p(1)

      rc     = r(size(r))/(0.0014766_DBPR*2.0_DBPR)
      if (rc<r_disc/(0.0014766_DBPR)) rc = r_disc/(0.0014766_DBPR)+ &
         0.5_DBPR*((r(size(r))/(0.0014766_DBPR))-(r_disc/(0.0014766_DBPR)))

      step = 5E-3_DBPR

      if (r_disc>0) then
         taylor_end = 0.5_DBPR
      else
         taylor_end = 0.5_DBPR
      end if

      call fluid_solve(step, rtol, atol, 'Y1', r_disc, omega, MATRIX, Y1, rc, taylor_end, &
         H1_init=0.004_DBPR, K_init=0.01_DBPR, W_init=1.0_DBPR)
      call fluid_solve(step, rtol, atol, 'Y2', r_disc, omega, MATRIX, Y2, rc, taylor_end, &
         H1_init=0.001_DBPR, K_init=0.1_DBPR, W_init=0.01_DBPR)
      call fluid_solve(step, rtol, atol, 'Y3', r_disc, omega, MATRIX, Y3, rc, taylor_end, &
         H1_init=0.06_DBPR, K_init=0.04_DBPR, W_init=0.1_DBPR )

      call fluid_solve(step, rtol, atol, 'Y4', r_disc, omega, MATRIX, Y4, rc, taylor_end, K_init=e(1)+p(1))
      call fluid_solve(step, rtol, atol, 'Y5', r_disc, omega, MATRIX, Y5, rc, taylor_end, K_init=-(e(1)+p(1)))

      call cramers(Y1(size(Y1,1), 2:5), Y2(size(Y2,1), 2:5), Y3(size(Y3,1), 2:5), -Y4(size(Y4,1), 2:5), Y5(size(Y5,1), 2:5), sol)

      allocate(Y_back(size(Y5,1), size(Y5, 2)), Y_front(size(Y1,1), size(Y1, 2)))

      Y_front = sol(1)*Y1+ sol(2)*Y2+ sol(3)*Y3
      Y_back  = sol(4)*Y4+ Y5
      Y_front(:,1) = Y1(:,1)
      Y_back(:,1)  = Y4(:,1)

      allocate(final_Y(size(Y1,1)+size(Y4,1)-1, 5))

      Y_front = Y_front(size(Y_front,1):1:-1, :)

      ! call all_fluid_functions(Y_back, Y_front, MATRIX, omega, r_grid, H0, V, H1, K, W, X)

      ! open(unit=343, file= 'fluidcheck.dat', status= 'replace', action= 'write')
      ! do i= 1, size(H1)
      !     write(343, *) r_grid(i), H0(i), V(i), H1(i), K(i), W(i), X(i)
      ! end do
      ! close(343)

      H1_surf = Y_front(size(Y_front, 1), 2)
      K_surf  = Y_front(size(Y_front, 1), 3)
      Radius  = r(size(r))/(0.0014766_DBPR)
      Mass    = m(size(m))/(0.0014766_DBPR)

      step  = 5E-2_DBPR
      r_far = 50.0_DBPR/omega

      call zerilli_solve(step, rtol, atol, r_far, Mass, Radius, omega, H1_surf, K_surf, r_far_dat, Z_dat, dZ_dat)

      ! open(unit=343, file= 'farcheck.dat', status= 'replace', action= 'write')
      ! do i= 1, size(Z_dat)
      !     write(343, *) r_far_dat(i), Z_dat(i), dZ_dat(i)
      ! end do
      ! close(343)

      call zerilli_far_lindblom(omega, Mass, r_far_dat(size(r_far_dat)), Z_dat(size(Z_dat)), dZ_dat(size(dZ_dat)), B_out)

   end subroutine outoing_omega

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine all_fluid_functions(Y_back, Y_front, MATRIX, omega, r_grid, H0, V, H1, K, W, X)

      real(DBPR), intent(in) :: omega
      real(DBPR), allocatable, intent(in) ::  Y_back(:,:), Y_front(:,:), MATRIX(:,:)
      real(DBPR), allocatable, intent(out) :: r_grid(:), H0(:), V(:), H1(:), K(:), W(:), X(:)
      real(DBPR), allocatable :: r_dat(:), m_dat(:), p_dat(:), e_dat(:), de_dp_dat(:), lamb_dat(:), &
         phi_dat(:), dphi_dat(:)
      real(DBPR) :: r, m, p, e, dedp, lamb, phi, dphi, l
      integer :: i

      l = 2.0_DBPR

      call read_TOV_data(MATRIX, r_dat, m_dat, p_dat, e_dat, de_dp_dat, lamb_dat, phi_dat, dphi_dat)
      r_dat    = r_dat/(0.0014766_DBPR)
      m_dat    = m_dat/(0.0014766_DBPR)
      p_dat    = p_dat*((0.0014766_DBPR)**2)
      e_dat    = e_dat*((0.0014766_DBPR)**2)
      dphi_dat = dphi_dat*(0.0014766_DBPR)

      r_grid  = [Y_back(1:size(Y_back,1)-1:1,1),Y_front(:,1)]
      H1      = [Y_back(1:size(Y_back,1)-1:1,2),Y_front(:,2)]
      K       = [Y_back(1:size(Y_back,1)-1:1,3),Y_front(:,3)]
      W       = [Y_back(1:size(Y_back,1)-1:1,4),Y_front(:,4)]
      X       = [Y_back(1:size(Y_back,1)-1:1,5),Y_front(:,5)]

      allocate(H0(size(r_grid)), V(size(r_grid)))

      do i=1, size(r_grid)

         r = r_grid(i)
         call interpolate(r_dat, m_dat,     r, m)
         call interpolate(r_dat, p_dat,     r, p)
         call interpolate(r_dat, e_dat,     r, e)
         call interpolate(r_dat, de_dp_dat, r, dedp)
         call interpolate(r_dat, lamb_dat,  r, lamb)
         call interpolate(r_dat, phi_dat,   r, phi)
         call interpolate(r_dat, dphi_dat,  r, dphi)

         H0(i) = (8*PI*(r**3)*exp(-phi)*X(i)- (0.5_DBPR*l*(l+1)*(m+4*PI*(r**3)*p)- &
            (omega**2)*(r**3)*exp(-2*phi-2*lamb))*H1(i) + &
            (0.5_DBPR*(l+2)*(l-1)*r- (omega**2)*(r**3)*exp(-2*phi)- (1/r)*exp(2*lamb)*(m+4*PI*(r**3)*p)* &
            (3*m- r+ 4*PI*(r**3)*p))*K(i))/ (3*m+ 0.5_DBPR*(l+2)*(l-1)*r+ 4*PI*(r**3)*p)

         V(i)  = (X(i) + (1/r)*(-(p+e)*dphi)*exp(phi-lamb)*W(i) - 0.5_DBPR*(p+e)*exp(phi)*H0(i))/ ((omega**2)*(p+e)*exp(-phi))

      end do

      V(1:5) = V(6)

   end subroutine all_fluid_functions

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine fluid_solve(in_step, rtol, atol, type, r_disc_in, omega, MATRIX, Y, rc, taylor_end, H1_init, K_init, W_init)

      type(lsoda_class) :: lsoda
      character(len=*), intent(in) :: type
      real(DBPR), intent(in) :: in_step, rtol, atol, r_disc_in, omega, rc, taylor_end
      real(DBPR), intent(in), allocatable :: MATRIX (:,:)
      real(DBPR), intent(out), allocatable :: Y(:,:)
      real(DBPR), intent(in), optional :: H1_init, W_init, K_init
      integer :: i, j, checkpoint, flag, r_disc_flag
      real(DBPR) :: r_disc, step, rout_prev, rout, m, p, e, dedp, lamb, phi, dphi, d_lamb_phi_r
      real(DBPR), allocatable :: r_dat(:), m_dat(:), p_dat(:), e_dat(:), de_dp_dat(:), lamb_dat(:), &
         phi_dat(:), dphi_dat(:), temp(:), d_lamb_phi_r_dat(:), temp_dat(:)
      integer :: istate, itask, neq
      real(DBPR) :: r, x(4)
      real(DBPR) :: l, W_0, K_0, H1_0, X_0, x_in_0(4), x_in_2(4)

      !---------------------------------------------------------------------------------------------------------------------------

      call read_TOV_data(MATRIX, r_dat, m_dat, p_dat, e_dat, de_dp_dat, lamb_dat, phi_dat, dphi_dat)

      r_dat    = r_dat/(0.0014766_DBPR)
      r_disc   = r_disc_in/(0.0014766_DBPR)
      m_dat    = m_dat/(0.0014766_DBPR)
      p_dat    = p_dat*((0.0014766_DBPR)**2)
      e_dat    = e_dat*((0.0014766_DBPR)**2)
      dphi_dat = dphi_dat*(0.0014766_DBPR)

      r_disc_flag = 0

      allocate(temp(size(e_dat)), d_lamb_phi_r_dat(size(e_dat)))
      temp = exp(-lamb_dat)*dphi_dat/(r_dat**2)
      do i= 2, size(e_dat)
         d_lamb_phi_r_dat(i) = (temp(i)-temp(i-1))/(r_dat(i)-r_dat(i-1))
         if (isnan(d_lamb_phi_r_dat(i))) d_lamb_phi_r_dat(i)= d_lamb_phi_r_dat(i-1)
      end do
      d_lamb_phi_r_dat(1) = d_lamb_phi_r_dat(2)

      !---------------------------------------------------------------------------------------------------------------------------

      neq  = 4
      call lsoda%initialize(func_fluid_lsoda, neq, istate=istate, mxstep= 50000000)

      itask  = 1
      istate = 1

      do i= 1, size(r_dat)
         if (r_dat(i).ge.taylor_end) then
            checkpoint = i
            exit
         end if
      end do

      allocate(temp_dat(0), Y(0,0))

      step = in_step/(0.0014766_DBPR)
      l = 2.0_DBPR
      j = 0

      !---------------------------------------------------------------------------------------------------------------------------

      if (type == 'Y4' .or. type == 'Y5') then

         step = step
         W_0  = 1.0_DBPR
         K_0  = K_init*((0.0014766_DBPR)**2)
         H1_0 = (1/(l*(l+1)))*(2*l*K_0+ 16*PI*(p_dat(1)+e_dat(1))*W_0)
         X_0  = (e_dat(1)+p_dat(1))*exp(phi_dat(1))*(((4*PI/3)*(e_dat(1)+3*p_dat(1))- &
            (omega**2)*exp(-2*phi_dat(1))/l)*W_0+ 0.5_DBPR*K_0)

         x_in_0 = [H1_0, K_0, W_0, X_0]

         call taylor_2(x_in_0, p_dat(1), e_dat(1), phi_dat(1), de_dp_dat(1), omega, x_in_2)

         do i= 1, checkpoint
            x = x_in_0 + 0.5_DBPR*(r_dat(i)**2)*x_in_2
            temp_dat = [temp_dat, [r_dat(i), x]]
            r = r_dat(i)
            rout = r_dat(i+1)
            j = j+1
         end do

      else if (type == 'Y1' .or. type == 'Y2' .or. type == 'Y3') then

         step = -step
         x    = [H1_init, K_init, W_init, 0.0_DBPR]
         r    = r_dat(size(r_dat))
         rout = r_dat(size(r_dat)-1)
         temp_dat = [temp_dat, [r, x]]
         j = j+1

      end if

      flag = 0

890   itask = 1
      istate = 1

      do

         if (type == 'Y1' .or. type == 'Y2' .or. type == 'Y3') then

            if (rout<rc .and. flag==0) then
               flag = 1
               rout = rc
            else if (rout<rc .and. flag==1) then
               exit
            end if

         else if (type == 'Y4' .or. type == 'Y5') then

            if (rout>rc .and. flag==0) then
               flag = 1
               rout = rc
            else if (rout>=rc .and. flag==1) then
               exit
            end if

         end if

         call lsoda%integrate(x, r, rout, rtol, [atol], itask, istate)
         temp_dat = [temp_dat, [rout, x]]

         if ( istate<0 ) then
            write (6,11) istate
11          format (///' Error halt in fluid solve... ISTATE =',i3)
            stop 1
         else
            rout = rout + step
         end if
         j = j+1

         if (r_disc>0) then
            if (type == 'Y4' .or. type == 'Y5') then
               if (rout>r_disc .and. r<r_disc .and. r_disc_flag == 0) then
                  rout        = r_disc - 1E-5_DBPR/(0.0014766_DBPR)
                  if (rout>r) rout_prev = rout
                  r_disc_flag = r_disc_flag + 1
               end if
               if (r_disc_flag == 1 .and. (rout_prev .ne. rout)) then
                  r           = r_disc + 1E-5_DBPR/(0.0014766_DBPR)
                  rout        = r + step
                  x           = x
                  temp_dat    = [temp_dat, [r, x]]
                  j           = j+1
                  r_disc_flag = 0
                  goto 890
               end if
            end if
         end if

      end do

      Y = reshape(temp_dat, [j,5], order= [2,1])

      !---------------------------------------------------------------------------------------------------------------------------

   contains

      subroutine func_fluid_lsoda(self_, Neq_, r_, x_, Ydash_, ierr_)
         implicit none
         class(lsoda_class), intent(inout) :: self_
         integer, intent(in) :: Neq_
         real(DBPR), intent(in) :: r_
         real(DBPR), intent(in) :: x_(Neq_)
         real(DBPR), intent(out) :: Ydash_(Neq_)
         integer, intent(out) :: ierr_
         real(DBPR) :: HH0, HH1, VV, WW, KK, XX
         ierr_ = 0
         self_%iwork = self_%iwork

         HH1 = x_(1)
         KK  = x_(2)
         WW  = x_(3)
         XX  = x_(4)

         call interpolate(r_dat, m_dat,            r_, m)
         call interpolate(r_dat, p_dat,            r_, p)
         call interpolate(r_dat, e_dat,            r_, e)
         call interpolate(r_dat, de_dp_dat,        r_, dedp)
         call interpolate(r_dat, lamb_dat,         r_, lamb)
         call interpolate(r_dat, phi_dat,          r_, phi)
         call interpolate(r_dat, dphi_dat,         r_, dphi)
         call interpolate(r_dat, d_lamb_phi_r_dat, r_, d_lamb_phi_r)

         HH0 = (8*PI*(r**3)*exp(-phi)*XX- (0.5_DBPR*l*(l+1)*(m+4*PI*(r**3)*p)- (omega**2)*(r**3)*exp(-2*phi-2*lamb))*HH1 + &
            (0.5_DBPR*(l+2)*(l-1)*r- (omega**2)*(r**3)*exp(-2*phi)- (1/r)*exp(2*lamb)*(m+4*PI*(r**3)*p)* &
            (3*m- r+ 4*PI*(r**3)*p))*KK)/ (3*m+ 0.5_DBPR*(l+2)*(l-1)*r+ 4*PI*(r**3)*p)

         VV  = (XX + (1/r)*(-(p+e)*dphi)*exp(phi-lamb)*WW - 0.5_DBPR*(p+e)*exp(phi)*HH0)/ ((omega**2)*(p+e)*exp(-phi))

         Ydash_ = fluid_func(x_, omega, r_, m, p, e, dedp, lamb, phi, dphi, d_lamb_phi_r)

      end subroutine func_fluid_lsoda

   end subroutine fluid_solve

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine taylor_2(x_in_0, p0, e0, phi0, de_dp_0, omega, x_in_2)

      real(DBPR), intent(in) ::  x_in_0(:), p0, e0, phi0, omega, de_dp_0
      real(DBPR), intent(out) :: x_in_2(4)
      real(DBPR) :: l, Q0_X, Q0_K, Q0_H1, Q1_X, Q1_K, Q1_W, p2, e2, phi2, phi4, p4 !, H0_2, V_2
      real(DBPR) :: T_MAT(4,4), U_MAT(4,4), T_INV_MAT(4,4), Y0(4,1), Y2(4,1)

      l = 2.0_DBPR

      ! Y   = [H1; K; W; X]
      !x_in = [H1, K, W, X]

      Y0(1,:) = [x_in_0(1)]
      Y0(2,:) = [x_in_0(2)]
      Y0(3,:) = [x_in_0(3)]
      Y0(4,:) = [x_in_0(4)]

      ! 1/(γp)= de/dp* 1/(p+e)

      p2   = -(4*PI/3)*(e0+ p0)*(e0+ 3*p0)
      e2   = de_dp_0*p2
      phi2 = (4*PI/3)*(e0+ 3*p0)
      p4   = -(2*PI/5)*(e0+ p0)*(e2+ 5*p2)- (2*PI/3)*(e2+ p2)*(e0+ 3*p0)- (32*(PI**2)/9)*e0*(e0+ p0)*(e0+ 3*p0)
      phi4 = (2*PI/5)*(e2+ 5*p2)+ (32*(PI**2)/9)*e0*(e0+ 3*p0)

      Q0_X  = (4/((l+2)*(l-1)))*(8*PI*exp(-phi0))
      Q0_K  = -(4/((l+2)*(l-1)))*((omega**2)*exp(-2*phi0)+ 8*PI*e0/3)
      Q0_H1 = -(4/((l+2)*(l-1)))*(-(omega**2)*exp(-2*phi0)+ (2*PI/3)*l*(l+1)*(e0+ 3*p0))

      Q1_X = (2/(l*(l+1)))*(de_dp_0*exp(-phi0)/(p0+e0))
      Q1_K = (2/(l*(l+1)))*(3./2.)
      Q1_W = (2/(l*(l+1)))*(4*PI*(l+1)*e0/3)

      T_MAT(1,:) = [ (l+3)/2, -1.0_DBPR, -8*PI*(p0+e0)*(l+3)/(l*(l+1)), 0.0_DBPR ]

      U_MAT(1,:) = [ (4*PI*(2*l+3)*e0/3) - 4*PI*p0 + 0.5_DBPR*Q0_H1, -8*PI*(p0+e0)*Q1_K + 0.5_DBPR*Q0_K, &
         -8*PI*(p0+e0)*Q1_W+ 8*PI*(p2+e2)/l, -8*PI*(p0+e0)*Q1_X + 0.5_DBPR*Q0_X ]

      T_MAT(2,:) = [ -l*(l+1)/4, (l+2)/2, 4*PI*(p0+e0), 0.0_DBPR ]

      U_MAT(2,:) = [ 0.5_DBPR*Q0_H1, (4*PI/3)*(e0+ 3*p0)+ 0.5_DBPR*Q0_K, -4*PI*(p2+e2+ 8*PI*e0*(e0+p0)/3), 0.5_DBPR*Q0_X ]

      T_MAT(3,:) = [ -l*(l+1)*exp(phi0)*(p0+e0)/8, 0.0_DBPR, -0.5_DBPR*exp(phi0)*(p0+e0)*(-(omega**2)*exp(-2*phi0)- &
         4*PI*(p0+e0)+ (l+2)*phi2), (l+2)/2 ]

      U_MAT(3,:) = [ exp(phi0)*(p0+e0)* ((omega**2)*exp(-2*phi0)/2+ Q0_H1/4), &
         exp(phi0)*(p0+e0)* (phi2+ Q0_K/4- l*(l+1)*phi2*Q1_K/2), &
         exp(phi0)*(p0+e0)* ((l+1)*phi4- 2*PI*(p2+e2)- 16*(PI**2)*e0*(e0+p0)/3+ (phi4- 4*PI*e0*phi2/3)+ &
         (omega**2)*exp(-2*phi0)*(phi2-4*PI*e0/3)- l*(l+1)*phi2*Q1_W/2), &
         exp(phi0)*(p0+e0)* (Q0_X/4 -l*(l+1)*phi2*Q1_X/2)+ 0.5_DBPR*l*(phi2+ (p2+e2)/(p0+e0)) ]

      T_MAT(4,:) = [ 0.0_DBPR, -(p0+e0)/4, 0.5_DBPR*(p2+ (l+3)*(omega**2)*exp(-2*phi0)*(p0+e0)/(l*(l+1))), &
         0.5_DBPR*exp(-phi0) ]

      U_MAT(4,:) = [ (p0+e0)*Q0_H1/4, (p2+e2)/4+ (p0+e0)*Q0_K/4+ 0.5_DBPR*(omega**2)*exp(-2*phi0)*(p0+e0)*Q1_K, &
         -p4+ 4*PI*e0*p2/3- (omega**2)*exp(-2*phi0)*(p2+e2- 2*(p0+e0)*phi2)/(2*l)+ &
         0.5_DBPR*(omega**2)*exp(-2*phi0)*(p0+e0)*Q1_W, 0.5_DBPR*exp(-phi0)*phi2+ (p0+e0)*Q0_X/4+ &
         0.5_DBPR*(omega**2)*exp(-2*phi0)*(p0+e0)*Q1_X ]

      call matinv4(T_MAT, T_INV_MAT)

      Y2 = matmul(matmul(T_INV_MAT, U_MAT), Y0)

      ! H0_2 = Q0_X+Q0_K+Q0_H1 + Y2(2,1)
      ! V_2  = Q1_X+Q1_K+Q1_W  - (l+3)*Y2(3,1)/(l*(l+1))

      ! x_in_2 = [H1_2, K_2, W_2, X_2]

      x_in_2 = [Y2(1,1), Y2(2,1), Y2(3,1), Y2(4,1)]

   end subroutine taylor_2

!-------------------------------------------------------------------------------------------------------------------------------

   function fluid_func(x_in, omega, r, m, p, e, de_dp, lamb, phi, dphi, d_lamb_phi_r) result(d_H1_K_W_X)

      real(DBPR), intent(in) :: omega, r, m, p, e, de_dp, lamb, phi, dphi, d_lamb_phi_r, x_in(:)
      real(DBPR) :: d_H1_K_W_X(4), H0, V, H1, K, W, X, l, dH1, dK, dW, dX

      H1 = x_in(1)
      K  = x_in(2)
      W  = x_in(3)
      X  = x_in(4)
      l = 2.0_DBPR

      H0 = (8*PI*(r**3)*exp(-phi)*X- (0.5_DBPR*l*(l+1)*(m+4*PI*(r**3)*p)- (omega**2)*(r**3)*exp(-2*phi-2*lamb))*H1 + &
         (0.5_DBPR*(l+2)*(l-1)*r- (omega**2)*(r**3)*exp(-2*phi)- (1/r)*exp(2*lamb)*(m+4*PI*(r**3)*p)* &
         (3*m- r+ 4*PI*(r**3)*p))*K)/ (3*m+ 0.5_DBPR*(l+2)*(l-1)*r+ 4*PI*(r**3)*p)

      V  = (X + (1/r)*(-(p+e)*dphi)*exp(phi-lamb)*W - 0.5_DBPR*(p+e)*exp(phi)*H0)/ ((omega**2)*(p+e)*exp(-phi))

      dH1 = (-1/r)*(l+1+ 2*m*exp(2*lamb)/r+ 4*PI*(r**2)*(p-e)*exp(2*lamb))*H1+ (1/r)*exp(2*lamb)*(H0+ K- 16*PI*(p+e)*V)

      dK  = H0/r+ l*(l+1)*H1/(2*r)- ((l+1)/r- dphi)*K- 8*PI*(p+e)*exp(lamb)*W/r

      dW  = -(l+1)*W/r+ r*exp(lamb)*(de_dp*exp(-phi)*X/(p+e) -l*(l+1)*V/(r**2)+ H0/2+ K)

      dX  = -l*X/r + (p+e)*exp(phi)*( 0.5_DBPR*((1/r)-dphi)*H0+ 0.5_DBPR*((omega**2)*r*exp(-2*phi)+ l*(l+1)/(2*r))*H1+ &
         0.5_DBPR*(3_DBPR*dphi-(1/r))*K- l*(l+1)*dphi*V/(r**2)- (1/r)*(4*PI*(p+e)*exp(lamb)+ &
         (omega**2)*exp(lamb-2*phi)- (r**2)*d_lamb_phi_r)*W )

      d_H1_K_W_X = [dH1, dK, dW, dX]

   end function fluid_func

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine zerilli_solve(in_step, rtol, atol, rfar, Mass, Radius, omega, H1_surf, K_surf, r_out_dat, Z_dat, dZ_dat)

      type(lsoda_class) :: lsoda
      real(DBPR), intent(in) :: in_step, rtol, atol, rfar, Mass, Radius, omega, H1_surf, K_surf
      real(DBPR), intent(out), allocatable :: r_out_dat(:), Z_dat(:), dZ_dat(:)
      logical :: halt
      integer :: flag
      real(DBPR), allocatable :: b_, g_, h_, k_
      integer :: istate,itask,neq
      real(DBPR) :: l, n, step, r, rout, Z_surf, dZ_surf, x(2)

      !---------------------------------------------------------------------------------------------------------------------------

      neq  = size(x)
      call lsoda%initialize(func_zerilli_lsoda, neq, istate=istate, mxstep= 500000000)

      itask  = 1
      istate = 1

      !---------------------------------------------------------------------------------------------------------------------------

      allocate(r_out_dat(0), Z_dat(0), dZ_dat(0))

      halt = .false.
      flag = 0

      l = 2.0_DBPR
      n = 0.5_DBPR*(l-1)*(l+2)

      step = (in_step/(0.0014766_DBPR))
      r    = Radius

      b_ = Mass/Radius
      g_ = (n*(n+1)+ 3*n*b_+ 6*(b_**2))/(n+ 3*b_)
      h_ = (n- 3*n*b_- 3*(b_**2))/((1.0_DBPR-2*b_)*(n+ 3*b_))
      k_ = (1.0_DBPR)/(1.0_DBPR- 2*b_)

      Z_surf  = Radius*(k_*K_surf- H1_surf)/(k_*g_- h_)
      dZ_surf = (-h_*K_surf+ g_*H1_surf)/(k_*g_- h_)

      r_out_dat = [r_out_dat, r]
      Z_dat     = [Z_dat, Z_surf]
      dZ_dat    = [dZ_dat, dZ_surf]

      x    = [Z_surf, dZ_surf]
      rout = r + step

      do while (halt .eqv. .false.)

         call lsoda%integrate(x, r, rout, rtol, [atol], itask, istate)
         r_out_dat = [r_out_dat, rout]
         Z_dat     = [Z_dat, x(1)]
         dZ_dat    = [dZ_dat, x(2)]

         if (rout>rfar .and. flag==0) then
            flag = 1
            rout = rfar
         else if (rout>rfar .and. flag==1) then
            halt = .true.
         end if

         if ( istate<0 ) then
            write (6,11) istate
11          format (///' Error halt in Zerilli Solve... ISTATE =',i3)
            stop 1
         else
            rout = rout + step
         end if

      end do

      !----------------------------------------------------------------------------------------------------------------------------

   contains

      subroutine func_zerilli_lsoda(self_, Neq_, r_, x_, Ydash_, ierr_)
         implicit none
         class(lsoda_class), intent(inout) :: self_
         integer, intent(in) :: Neq_
         real(DBPR), intent(in) :: r_
         real(DBPR), intent(in) :: x_(Neq_)
         real(DBPR), intent(out) :: Ydash_(Neq_)
         integer, intent(out) :: ierr_
         ierr_ = 0
         self_%iwork = self_%iwork

         Ydash_ = zerilli_func(x_, omega, r_, Mass)

      end subroutine func_zerilli_lsoda

   end subroutine zerilli_solve

!-------------------------------------------------------------------------------------------------------------------------------

   function zerilli_func(x_in, omega, r, m) result(d_Z_dZ)

      real(DBPR), intent(in) :: x_in(:), omega, r, m
      real(DBPR) :: d_Z_dZ(2)
      real(DBPR) :: l, n, V, Z_in, dZ_in, dZ_dr, d2Z_dr2

      Z_in  = x_in(1)
      dZ_in = x_in(2)

      l = 2.0_DBPR
      n  = 0.5*(l-1)*(l+2)

      V= ((2*(1- 2*m/r))/(r**3* (n*r+3*m)**2))*(n**2*(n+1)*r**3+ 3*n**2*m*r**2+9*n*m**2*r+9*m**3)

      dZ_dr   = dZ_in/(1.0_DBPR- 2*m/r)
      d2Z_dr2 = (V-(omega**2))*Z_in/(1.0_DBPR- 2*m/r)

      d_Z_dZ = [dZ_dr, d2Z_dr2]

   end function zerilli_func

!-------------------------------------------------------------------------------------------------------------------------------

   subroutine zerilli_far_lindblom(omega, M, r_far, Z_far, dZ_far, A_inc)

      real(DBPR), intent(in) :: omega, M, Z_far, dZ_far
      double complex, intent(out) :: A_inc
      real(DBPR) :: n, l, r_tort, r_far
      double complex :: i, Z_out, Z_inc, dZ_out, dZ_inc, a1, a2, a0, A_out

      integer ::j
      real(DBPR):: r(1000), rt(1000)
      double complex :: ZZ(1000), dZZ(1000), ZZi(1000), dZZi(1000), ZZo(1000), dZZo(1000)

      i = (0.0_DBPR, 1.0_DBPR)

      l = 2.0_DBPR
      n = (l-1)*(l+2)/2.0_DBPR

      a0 = (0.01,0.0001)
      a1 = -i*(n+1)*a0/omega
      a2 = (-n*(n+1) + i*M*omega*(1.5 + 3/n))*a0/(2*(omega**2))

      r_tort = r_far + 2*M*log(r_far/(2*M)- 1)

      Z_inc = exp(-i*omega*r_tort)*(a0 + a1/r_far + a2/(r_far**2))
      dZ_inc = -i*omega*exp(-i*omega*r_tort)*(a0 + a1/r_far + (a2-i*a1*(1-2*M/r_far)/omega)/(r_far**2))

      ! Z_out = exp(i*omega*r_tort)*(conjg(a0) + conjg(a1)/r_far + conjg(a2)/(r_far**2))
      ! dZ_out = i*omega*exp(i*omega*r_tort)*(conjg(a0) + conjg(a1)/r_far + (conjg(a2)+i*conjg(a1)*(1-2*M/r_far)/omega)/(r_far**2))

      Z_out = conjg(Z_inc)
      dZ_out = conjg(dZ_inc)

      A_out = (Z_far*dZ_inc- dZ_far*Z_inc)/(Z_out*dZ_inc- Z_inc*dZ_out)
      A_inc = (-Z_far*dZ_out+ dZ_far*Z_out)/(Z_out*dZ_inc- Z_inc*dZ_out)

      do j=1, size(r)
         r(j)= r_far+j
      end do

      rt = r + 2*M*log(r/(2*M)- 1)

      ZZi = exp(-i*omega*rt)*(a0 + a1/r + a2/(r**2))
      dZZi = -i*omega*exp(-i*omega*rt)*(a0 + a1/r + (a2-i*a1*(1-2*M/r)/omega)/(r**2))

      ZZo = exp(i*omega*rt)*(conjg(a0) + conjg(a1)/r + conjg(a2)/(r**2))
      dZZo = i*omega*exp(i*omega*rt)*(conjg(a0) + conjg(a1)/r + (conjg(a2)+i*conjg(a1)*(1-2*M/r)/omega)/(r**2))

      ZZ  = A_inc*ZZi+ A_out*ZZo
      dZZ = A_inc*dZZi + A_out*dZZo

      ! open(unit=1111, file= 'farcheck.dat', status= 'old', action= 'write', position='append')
      ! do j= 1, size(r)
      !     write(1111, *) r(j), DREAL(ZZ(j)), DREAL(dZZ(j))
      ! end do
      ! close(1111)

      ! A_inc = abs(A_inc)

   end subroutine zerilli_far_lindblom

!-------------------------------------------------------------------------------------------------------------------------------

end module FULLGR
