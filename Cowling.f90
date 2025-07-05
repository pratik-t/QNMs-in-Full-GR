module COWLING
    use IO
    use TOV
    implicit none
contains

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine cowling_frequency(MATRIX, in_guess, r_disc, fmode, pmode, gmode)

    real(DBPR), intent(in) :: r_disc
    real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
    real(DBPR), intent(out) :: fmode, pmode, gmode
    real(DBPR) :: i, step, f_guess, p_guess, g_guess, in_guess(:), flag, temp
    integer :: f_node, p_node, g_node, tempi
    real(DBPR), allocatable :: r(:), m(:), p(:), e(:), de_dp(:), lamb(:), phi(:), dphi(:), guesses(:)
    integer, allocatable :: nodes(:)

    call read_TOV_data(MATRIX, r, m, p, e, de_dp, lamb, phi, dphi)
    
    flag = 0

    ! This do loop finds the first (exit if flag_sign== n), n points where the function cuts the x axis
    ! These points will be the initial guesses for f, p and g modes. 
    ! We will then do newton raphson for these initial guesses. 
            
    222 if (in_guess(1)<0.or.in_guess(2)<0) then
        i    = 1
        step = 0.1
        call find_guess(MATRIX, r_disc, i, step, 5, guesses, nodes)

        f_guess = guesses(1)
        p_guess = guesses(2)
        g_guess = guesses(3)
        
        call find_freq(f_guess, r, p, e, de_dp, lamb, phi, dphi, fmode, f_node, r_disc)
        call find_freq(p_guess, r, p, e, de_dp, lamb, phi, dphi, pmode, p_node, r_disc)
        
        ! print*, '$ f-node= ', f_node, fmode, 'p1-node= ',p_node, pmode
            
    else

        f_guess = in_guess(1)
        p_guess = in_guess(2)
        g_guess = in_guess(3)
        
        call find_freq(f_guess, r, p, e, de_dp, lamb, phi, dphi, fmode, f_node, r_disc)
        call find_freq(p_guess, r, p, e, de_dp, lamb, phi, dphi, pmode, p_node, r_disc)
        
        ! print*, 'f-node= ', f_node, fmode, 'p1-node= ',p_node, pmode
        
    end if
    
    gmode = -10.0_DBPR
    
    if (r_disc>0.) then
        if (g_guess<0) then
            in_guess = (/-10.0_DBPR,-10.0_DBPR,-10.0_DBPR/)
            goto 222
        end if
        call find_freq(g_guess, r, p, e, de_dp, lamb, phi, dphi, gmode, g_node, r_disc)
    end if

    if (gmode>fmode .and. g_node==1.0_DBPR) then
        ! gmode has the highest frequency and is actually the pmode. fmode has the lowest frequency. So fix the labels. 
        temp  = gmode
        gmode = fmode
        fmode = pmode
        pmode = temp

        tempi  = g_node
        g_node = f_node
        f_node = p_node
        p_node = tempi
    end if

    if ((f_node.ne.0.0_DBPR).or.(p_node.ne.1.0_DBPR).or.(gmode>0.and.g_node.ne.1)) then
        
        if (flag == 1) then
            if (f_node.ne.0.0_DBPR) fmode = NaN/NaN
            if (p_node.ne.1.0_DBPR) pmode = NaN/NaN
            if (g_node.ne.1.0_DBPR) gmode = NaN/NaN
        else
            in_guess = (/-10.0_DBPR,-10.0_DBPR,-10.0_DBPR/)
            flag = flag+1
            goto 222
        end if

    end if

    ! print*, 'f-node= ', f_node, fmode, 'p1-node= ',p_node, pmode, 'g-node= ', g_node, gmode 

    fmode = fmode*47.71345159236942258889E-3_DBPR
    pmode = pmode*47.71345159236942258889E-3_DBPR
    if (gmode>0) gmode = gmode*47.71345159236942258889E-3_DBPR

end subroutine cowling_frequency

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine find_guess(MATRIX, r_disc, i, step, stop, guesses, nodes)
    
    integer, intent(in) :: stop
    real(DBPR), intent(in) :: r_disc
    real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
    real(DBPR), allocatable, intent(out) :: guesses(:)
    integer, allocatable, intent(out) :: nodes(:)
    integer :: flag_sign, node, node_prev
    real(DBPR) :: i, step, func_val, sign_prev, sign_
    real(DBPR), allocatable :: r(:), m(:), p(:), e(:), de_dp(:), lamb(:), phi(:), dphi(:)
    
    if (allocated(nodes)) deallocate(nodes)
    allocate(guesses(0), nodes(0))

    call read_TOV_data(MATRIX, r, m, p, e, de_dp, lamb, phi, dphi)
    
    flag_sign = 0
    call mode_rk4(i, r, p, e, de_dp, lamb, phi, dphi, r_disc, func_val, node)
    sign_ = sign(real(1, DBPR), func_val)
    
    666 do

        sign_prev = sign_
        node_prev = node
        call mode_rk4(i, r, p, e, de_dp, lamb, phi, dphi, r_disc, func_val, node)
        
        if (isnan(func_val)) then
            print*, 'ERROR IN COWLING FIND_GUESS: NaN Values'
            stop
        end if

        sign_ = sign(real(1, DBPR), func_val)
        
        if (sign_/= sign_prev) then 
            guesses = [guesses, i-step]
            nodes   = [nodes, node_prev]
            flag_sign = flag_sign +1
            goto 666
        end if

        i=i+step

        if (flag_sign== stop) then 
            exit
        end if

    end do

end subroutine find_guess

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine find_freq(guess, r, p, e, de_dp, lamb, phi, dphi, root, nodes, r_disc, not_Newton)

    integer, intent(in), optional :: not_Newton
    real(DBPR), intent(in) :: r_disc
    real(DBPR), intent(in) :: r(:), p(:), e(:), de_dp(:), lamb(:), phi(:), dphi(:)
    real(DBPR), intent(out) ::  root
    integer, intent(out) :: nodes
    real(DBPR) :: tol, error, h, original_guess, guess, new_guess, f, fh, fdash

    111 original_guess = guess

    h     = 1e-10_DBPR     ! This is the order of error on derivative for Newton-Raphson method
    tol   = 1e-10_DBPR     ! Tolerance for closeness of values in Newton-Raphson
    error = 10             ! Temporary initial value for error
    
    do while (error>tol)
    
        call mode_rk4(guess, r, p, e, de_dp, lamb, phi, dphi, r_disc, f, nodes)
        call mode_rk4(guess-h, r, p, e, de_dp, lamb, phi, dphi, r_disc, fh, nodes)

        fdash     = (f- fh)/h
        new_guess = guess- (f/fdash)           
        error     = abs(new_guess-guess)
        guess = new_guess
        
    end do

    root = guess

    if (present(not_Newton)) then
        call mode_rk4(guess, r, p, e, de_dp, lamb, phi, dphi, r_disc, f, nodes)
        h     = 0.01
        error = 10
        do while (nodes>1 .and. error>tol)
            guess= guess+h
            call mode_rk4(guess, r, p, e, de_dp, lamb, phi, dphi, r_disc, f, nodes)
            error = abs(f)
            if (nodes<1) then
                guess = guess-h
                h     = h/10.0_DBPR
            end if
        end do
    end if

    if (root<0.0_DBPR .or. isnan(root)) then
        guess = original_guess+0.01_DBPR
        goto 111
    end if

end subroutine find_freq

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine mode_rk4(omega, r, p, e, de_dp, lamb, phi, dphi, r_disc, fomega, nodes) 

    real(DBPR), intent(in) :: omega, r_disc, r(:), p(:), e(:), de_dp(:), lamb(:), phi(:), dphi(:)
    real(DBPR), intent(out) :: fomega 
    integer, intent(out) :: nodes
    integer :: i
    real(DBPR) :: h, l, W, V, Wm, Vm, Wp, Vp, c1, c2, x(2), k1(2), k2(2), k3(2), k4(2), sign_W_p, sign_V_p, sign_W, sign_V

    l = 2.0_DBPR
    W = r(1)**(l+1)
    V = -(r(1)**l)/l
    x = [W, V]
    
    nodes= 0
    sign_W_p= sign(real(1, DBPR), W)
    sign_V_p= sign(real(1, DBPR), V)
    
    do i= 1, size(r)-1, 2
        
        Wm = x(1)
        Vm = x(2)
        
        h  = r(i+2)-r(i)
        
        k1 = h* mode_func(r(i), x, omega, de_dp(i), lamb(i), phi(i), dphi(i))
        k2 = h* mode_func(r(i+1), x+0.5_DBPR*k1, omega, de_dp(i+1), lamb(i+1), phi(i+1), dphi(i+1))
        k3 = h* mode_func(r(i+1), x+0.5_DBPR*k2, omega, de_dp(i+1), lamb(i+1), phi(i+1), dphi(i+1))
        k4 = h* mode_func(r(i+2), x+k3, omega, de_dp(i+2), lamb(i+2), phi(i+2), dphi(i+2))
        
        x  = x+ (1.0_DBPR/6.0_DBPR)*(k1+ 2*k2+ 2*k3+ k4)
        
        if (abs(r(i)-r_disc)<1E-10) then
            
            Wp   = x(1)
            Vp   = x(2)
            c1   = exp(lamb(i))*exp(-2*phi(i))*(omega**2)*(r_disc**2)                          
            c2   = ((e(i) + p(i))/(e(i+2) + p(i+2)))*((c1*Vm) + (dphi(i)*Wm))
            x(1) = Wm
            x(2) = (c1**(-1.))*(c2 - dphi(i)*Wp)

        end if
        
        sign_W= sign(real(1, DBPR), x(1))
        sign_V= sign(real(1, DBPR), x(2))
        if ((sign_W.ne.sign_W_p) .and. (sign_V.ne.sign_V_p)) then
            nodes= nodes+1
            sign_W_p = sign_W
            sign_V_p = sign_V
        end if
        
    end do
    
    W = x(1)
    V = x(2)
    
    fomega = (omega**2)*exp(lamb(i))*exp(-2*phi(i))*V+ (1/(r(i)**2))*dphi(i)*W
    
end subroutine mode_rk4

!-----------------------------------------------------------------------------------------------------------------------------------

function mode_func(r, x, omega, de_dp, lamb, phi, dphi) result(dW_dV)

    real(DBPR), intent(in) :: r, omega, de_dp, lamb, phi, dphi
    real(DBPR), intent(in) :: x(:) 
    real(DBPR) :: dW_dv(2), l, W, V, dW, dV

    W = x(1)
    V = x(2)
    l = 2.0_DBPR
    
    dW = de_dp* ((omega**2)* (r**2)* exp(lamb)* exp(-2*phi)* V+ (dphi* W))- (l*(l+1)* exp(lamb)* V)  
    dV = 2* dphi* V- (exp(lamb)* W)/ (r**2)

    dW_dV = [dW, dV]
    
end function mode_func

!-----------------------------------------------------------------------------------------------------------------------------------

end module COWLING