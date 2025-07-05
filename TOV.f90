module TOV
    use IO
    use odepack_mod
    implicit none
contains

!TOV ODEs-----------------------------------------------------------------------------------------------------------------------

function TOV_func(e, r, x) result(dm_dp_dphi)  ! meff is the effective nucleon mass
    real(DBPR), intent(in) :: e, r
    real(DBPR), intent(in) :: x(3)
    real(DBPR) :: dm_dp_dphi(3)
    real(DBPR) :: m, p, phi, dm, dp, dphi
    
    m   = x(1)
    p   = x(2)
    phi = x(3)

    dm   = 4*PI*e*r**2
    dp   = -(e + p)*(m + 4*PI*p*r**3)/(r*(r - 2*m))
    dphi = (m + 4*PI*p*r**3)/(r*(r - 2*m))

    dm_dp_dphi= [dm, dp, dphi]
    
end function TOV_func

!Tidal ODEs---------------------------------------------------------------------------------------------------------------------

function tidal_func(r, y, e1, p1, de_dp_1, lamb, dphi) result(dy) 

    real(DBPR), intent(in) :: y, r, e1, p1, de_dp_1, lamb, dphi
    real(DBPR) :: Q, dy, p, e
    
    p     = p1
    e     = e1

    Q  = 4*PI*exp(2*lamb)*(5*e+ 9*p+ ((p1+e1)*de_dp_1))- (6*exp(2*lamb))/(r**2)- 4*(dphi)**2
    dy = (-1./r)*(y**2+ y*exp(2*lamb)*(1.+ 4*PI*r**2*(p- e))+ Q*r**2)
    
end function tidal_func

!Solving TOV--------------------------------------------------------------------------------------------------------------------

subroutine SolveTOV(p_data, e_data, h_in, ec, r_lim, r_fin, m_fin, e_disc, r_disc, MATRIX, tide) 
    
    real(DBPR), intent(in) :: h_in, ec, r_lim, e_disc
    real(DBPR), intent(in) :: p_data(:), e_data(:)
    real(DBPR), intent(out) :: m_fin, r_fin, r_disc
    real(DBPR), allocatable, intent(out) :: MATRIX(:,:)
    real(DBPR), intent(out), optional :: tide
    integer :: i, flag
    real(DBPR) :: h, r0, e, p, p_limit, mc, pc, phic, phi_fin, phi_far, x0(3), x0_prev(3), b, yR, k2
    real(DBPR), allocatable :: rdisc_dat(:), mdisc_dat(:), pdisc_dat(:), phidisc_dat(:), edisc_dat(:), &
                            r_array(:), m_array(:), p_array(:), e_array(:), de_dp(:), phi_array(:), &
                            dphi_array(:), lambda_array(:)

    
    call interpolate(e_data, p_data, ec, pc)
    
    h    = h_in
    mc   = 0.0_DBPR
    phic = 0.0_DBPR
    
    p_limit = p_data(1)         ! The limit replacing p<0

    flag = 0
    i    = 1
    r0   = 1E-15_DBPR           ! Starting point of radius. This increases by h in each step.
    x0   = [mc, pc, phic]       ! [Initial mass, Initial pressure, initial phi]
    p    = x0(2)                ! p= initial pressure. This will change in the loop.
    e    = ec                   ! e= initial energy density. This will change in the loop.
    
    allocate(r_array(0), m_array(0), p_array(0), phi_array(0), e_array(0))

    r_array   = [r0]
    m_array   = [mc]
    p_array   = [pc]
    phi_array = [phic]
    e_array   = [ec]
    
    if (e_disc>0.0_DBPR .and. ec>e_disc) then
        
        call init_discont(h_in, p_data, e_data, ec, e_disc, r_disc, rdisc_dat, &
                            mdisc_dat, pdisc_dat, phidisc_dat, edisc_dat)
    
        r_array   = [r_array, rdisc_dat]
        m_array   = [m_array, mdisc_dat]
        p_array   = [p_array, pdisc_dat]
        phi_array = [phi_array, phidisc_dat]
        e_array   = [e_array, edisc_dat]
        
        i  = size(r_array)              
        r0 = r_array(i)                
        x0 = [m_array(i), p_array(i), phi_array(i)]       
        p  = x0(2)               
        e  = e_array(i)

    else
        r_disc = -100.0_DBPR
        
    end if
    
    !Main RK4 Loop--------------------------------------------------------------------------------------------------------------

    do while(p>0)
        
        x0_prev = x0
        m_fin   = x0(1)
        r_fin   = r0
        phi_fin = x0(3)
        
        call rk4_(h, e, r0, x0)

        if (r_fin*1000 > r_lim) then
            ! Exit if radius is greater than r_lim since these stars won't be considered anyways.
            ! Flag for this condition is r_fin= r_lim+ 100
            r_fin= r_lim+1000.0_DBPR
            exit
        end if
        
        if (p<p_limit) then
            ! p < the first crust pressure point defines the surface
            flag= 1
            exit   
        end if
        
        r0 = r0+h
        p  = x0(2)
        i  = i+1
        call interpolate(p_data, e_data, p, e)
        
        if (p.ge.p_limit) then
            r_array   = [r_array, r0]
            m_array   = [m_array, x0(1)]
            p_array   = [p_array, p]
            phi_array = [phi_array, x0(3)]
            e_array   = [e_array, e]
        end if
        
    end do
    
    !---------------------------------------------------------------------------------------------------------------------------
    !If the do while loop exits by its while condition, the last element will be p<0. 
    !p will be negative so this last element must be removed.
    
    if (flag==0) then

        i         = size(r_array)
        r_array   = r_array(1:i-1)
        m_array   = m_array(1:i-1)
        p_array   = p_array(1:i-1)
        phi_array = phi_array(1:i-1)
        e_array   = e_array(1:i-1)

    end if

    !dρ/dp array----------------------------------------------------------------------------------------------------------------
    
    allocate(de_dp(size(r_array)))
    do i= 1, size(r_array)-1
        de_dp(i)= (e_array(i+1)-e_array(i))/(p_array(i+1)-(p_array(i)))
    end do
    de_dp(size(r_array))= de_dp(size(r_array)-1)

    !φ array--------------------------------------------------------------------------------------------------------------------
    
    phi_far   = (0.5_DBPR)*log(1.0_DBPR- (2.0_DBPR*m_fin/r_fin))
    phi_array = phi_array+ (phi_far-phi_fin)

    !dφ/dr, Λ array-------------------------------------------------------------------------------------------------------------

    dphi_array= (m_array + 4*PI*p_array*r_array**3)/(r_array*(r_array - 2*m_array))
    
    lambda_array= log( (1.0_DBPR- 2.0_DBPR*m_array/r_array)**(-0.5_DBPR) )

    !---------------------------------------------------------------------------------------------------------------------------
    !If the do while loop exits by its while condition then after the last element is removed, the number of elements
    !might be even. I need the number of elements to be odd so I repeat the last element

    if (mod(size(r_array),2)==0) then

        i            = size(r_array)
        r_array      = [r_array      , r_array(i)]
        m_array      = [m_array      , m_array(i)]
        p_array      = [p_array      , p_array(i)]
        e_array      = [e_array      , e_array(i)]
        de_dp        = [de_dp        , de_dp(i)]
        lambda_array = [lambda_array , lambda_array(i)]
        phi_array    = [phi_array    , phi_array(i)]
        dphi_array   = [dphi_array   , dphi_array(i)]

    end if

    !---------------------------------------------------------------------------------------------------------------------------

    if (present(tide)) then

        call rk4_tide(r_array, e_array, p_array, de_dp, lambda_array, dphi_array, yR)

        b  = m_fin/r_fin   
        
        k2 = ( (8.0_DBPR/5.0_DBPR)*(b**5)*((1-2.0_DBPR*b)**2)*(2.0_DBPR-yR+2.0_DBPR*b*(yR-1)) )/&
                ( 2.0_DBPR*b*(6-3.0_DBPR*yR+3.0_DBPR*b*(5.0_DBPR*yR-8)) + &
                4.0_DBPR*(b**3)*(13-11.0_DBPR*yR+b*(3.0_DBPR*yR-2)+2.0_DBPR*b**2*(1+yR)) + &
                3.0_DBPR*((1-2.0_DBPR*b)**2)*(2-yR+2.0_DBPR*b*(yR-1))*log(1-2.0_DBPR*b) )
        
        tide = (2.0_DBPR/3.0_DBPR)*k2*(r_fin/m_fin)**5

    end if

    !Returning------------------------------------------------------------------------------------------------------------------
    
    allocate(MATRIX(8,size(r_array)))
    
    ! do i=1, size(r_array)
    ! print*, r_array(i), p_array(i)/1.323423_DBPR, e_array(i)/1.323423_DBPR, &
    ! m_array(i)*1000/1.4766_DBPR
    ! end do

    MATRIX(1,:) = r_array
    MATRIX(2,:) = m_array
    MATRIX(3,:) = p_array
    MATRIX(4,:) = e_array
    MATRIX(5,:) = de_dp
    MATRIX(6,:) = lambda_array
    MATRIX(7,:) = phi_array
    MATRIX(8,:) = dphi_array

    deallocate(r_array, m_array, phi_array, e_array, de_dp, dphi_array, lambda_array)
    
    m_fin  = m_fin*1000/1.4766_DBPR    ! Output mass in units of solar mass
    r_fin  = r_fin*1000                ! Output radius in kilometers
    
end subroutine SolveTOV

!Solving TOV with LSODAR--------------------------------------------------------------------------------------------------------

subroutine SolveTOV_lsoda(p_data, e_data, h, ec, r_lim, r_fin, m_fin, e_disc, r_disc, MATRIX, tide)
    
    type(lsoda_class) :: lsoda, lsoda_disc
    real(DBPR), intent(in) :: h, ec, r_lim, e_disc
    real(DBPR), intent(in) :: p_data(:), e_data(:)
    real(DBPR), intent(out) :: m_fin, r_fin, r_disc
    real(DBPR), allocatable, intent(out) :: MATRIX(:,:)
    real(DBPR), intent(out), optional :: tide
    integer :: i, r_flag, p_flag, d_flag, istate, itask, neq, istate_disc, itask_disc
    real(DBPR) :: p_disc, r, rout, e, p, p_limit, mc, pc, phic, phi_fin, phi_far, x(3), b, yR, k2
    real(DBPR), allocatable :: r_array(:), m_array(:), p_array(:), e_array(:), de_dp(:), phi_array(:), &
                            dphi_array(:), lambda_array(:)
    real(DBPR) :: rtol, atol
    
    !---------------------------------------------------------------------------------------------------------------------------
    
    neq  = size(x)
    call lsoda%initialize(TOV_lsoda, neq, g=surface, ng=1, istate=istate, mxstep= 500000000)

    rtol   = 1E-10_DBPR
    atol   = 1E-10_DBPR   
    itask  = 1
    istate = 1
    
    !---------------------------------------------------------------------------------------------------------------------------

    call interpolate(e_data, p_data, ec, pc)
    
    mc   = 0.0_DBPR
    phic = 0.0_DBPR
    
    p_limit = p_data(1)         ! The limit replacing p<0

    r_flag = 0
    p_flag = 0
    d_flag = 0

    i    = 1
    r    = 1E-15_DBPR           ! Starting point of radius. This increases by h in each step.
    x    = [mc, pc, phic]       ! [Initial mass, Initial pressure, initial phi]
    p    = x(2)                 ! p= initial pressure. This will change in the loop.
    e    = ec                   ! e= initial energy density. This will change in the loop.
    rout = r + h
    
    allocate(r_array(0), m_array(0), p_array(0), phi_array(0), e_array(0))

    r_array   = [r]
    m_array   = [mc]
    p_array   = [pc]
    phi_array = [phic]
    e_array   = [ec]
    
    if (e_disc>0.0_DBPR .and. ec>e_disc) then

        call interpolate(e_data, p_data, e_disc, p_disc)

        call lsoda_disc%initialize(TOV_lsoda, neq, g= discont, ng=1, istate=istate_disc, mxstep= 5000000)
        rtol        = rtol
        atol        = atol
        itask_disc  = 1
        istate_disc = 1
        
        do 
            call lsoda_disc%integrate(x, r, rout, rtol, [atol], itask_disc, istate_disc)
            call interpolate(p_data, e_data, x(2), e)
            
            if (abs(e-e_disc)<1E-6_DBPR .and. x(2)>0) d_flag = 1
            
            r_array   = [r_array, r]
            m_array   = [m_array, x(1)]
            p_array   = [p_array, x(2)]
            phi_array = [phi_array, x(3)]
            e_array   = [e_array, e]

            if ( istate<0 ) then
                write (6,11) istate
                11 format (///' Error halt in TOV discontinuity... ISTATE = ', i3)
                stop 1
            else
                rout = rout + h
            end if

            if (d_flag==1) exit

        end do
        
        call back_interp(h, r_array, m_array, p_array, e_array, phi_array)
        
        r_array   = [1E-15_DBPR, r_array] 
        m_array   = [mc, m_array] 
        p_array   = [pc, p_array] 
        e_array   = [ec, e_array] 
        phi_array = [phic, phi_array]
        
        i = size(r_array)
        r_disc = r_array(i)
        r      = r_disc
        x      = [m_array(i), p_array(i), phi_array(i)]
        
        do while(p_disc-x(2).le.0.01)
            call rk4_(h, e_array(i), r, x)
            call interpolate(p_data, e_data, x(2), e)
            r_array   = [r_array, r+h] 
            m_array   = [m_array, x(1)] 
            p_array   = [p_array, x(2)] 
            e_array   = [e_array, e] 
            phi_array = [phi_array, x(3)]
            r         = r+h
            rout      = r+h
        end do
        
    else

        r_disc = -100.0_DBPR

    end if
    
    do while (p.ge.p_limit)
        
        call lsoda%integrate(x, r, rout, rtol, [atol], itask, istate)

        r_fin   = r
        m_fin   = x(1)
        phi_fin = x(3)
        p       = x(2)
        
        if (abs(x(2)-p_limit)<atol .and. p>0) p_flag = 1
        
        if (p<0) then

            p      = p_limit
            p_flag = 1
        end if
        
        call interpolate(p_data, e_data, p, e)
        
        if (r_fin*1000 > r_lim) then
            ! Exit if radius is greater than r_lim since these stars won't be considered anyways.
            ! Flag for this condition is r_fin= r_lim+ 1000
            r_fin = r_lim+1000.0_DBPR
            r_flag  = 1
            exit
        end if
        
        r_array   = [r_array, r]
        m_array   = [m_array, x(1)]
        p_array   = [p_array, p]
        phi_array = [phi_array, x(3)]
        e_array   = [e_array, e]
        
        if (p_flag ==1) exit

        if ( istate<0 ) then
            write (6,12) istate
            12 format (///' Error halt in TOV... ISTATE = ', i3)
            stop 1
        else
            rout = rout + h
        end if
        
    end do
    
    if ((istate < 0 .or. istate /= 3 .or. lsoda%jroot(1) /= -1) .and. r_flag==0 .and. p_flag==0) then
      print*, 'ISTATE = ', istate
      error stop 'Root finding Failed'
    endif

    !dρ/dp array----------------------------------------------------------------------------------------------------------------
    
    allocate(de_dp(size(r_array)))
    do i= 1, size(r_array)-1
        de_dp(i)= (e_array(i+1)-e_array(i))/(p_array(i+1)-(p_array(i)))
    end do
    de_dp(size(r_array))= de_dp(size(r_array)-1)

    !φ array--------------------------------------------------------------------------------------------------------------------
    
    phi_far   = (0.5_DBPR)*log(1.0_DBPR- (2.0_DBPR*m_fin/r_fin))
    phi_array = phi_array+ (phi_far-phi_fin)

    !dφ/dr, Λ array-------------------------------------------------------------------------------------------------------------

    dphi_array= (m_array + 4*PI*p_array*r_array**3)/(r_array*(r_array - 2*m_array))
    
    lambda_array= log( (1.0_DBPR- 2.0_DBPR*m_array/r_array)**(-0.5_DBPR) )

    !---------------------------------------------------------------------------------------------------------------------------
    ! If the do while loop exits by its while condition then after the last element is removed, the number of elements
    ! might be even. I need the number of elements to be odd so I repeat the last element

    if (mod(size(r_array), 2)==0) then

        i            = size(r_array)
        r_array      = [r_array      , r_array(i)]
        m_array      = [m_array      , m_array(i)]
        p_array      = [p_array      , p_array(i)]
        e_array      = [e_array      , e_array(i)]
        de_dp        = [de_dp        , de_dp(i)]
        lambda_array = [lambda_array , lambda_array(i)]
        phi_array    = [phi_array    , phi_array(i)]
        dphi_array   = [dphi_array   , dphi_array(i)]

    end if

    !Tidal Deformability--------------------------------------------------------------------------------------------------------

    if (present(tide)) then

        call rk4_tide(r_array, e_array, p_array, de_dp, lambda_array, dphi_array, yR)

        b  = m_fin/r_fin   
        
        k2 = ( (8.0_DBPR/5.0_DBPR)*(b**5)*((1-2.0_DBPR*b)**2)*(2.0_DBPR-yR+2.0_DBPR*b*(yR-1)) )/&
                ( 2.0_DBPR*b*(6-3.0_DBPR*yR+3.0_DBPR*b*(5.0_DBPR*yR-8)) + &
                4.0_DBPR*(b**3)*(13-11.0_DBPR*yR+b*(3.0_DBPR*yR-2)+2.0_DBPR*b**2*(1+yR)) + &
                3.0_DBPR*((1-2.0_DBPR*b)**2)*(2-yR+2.0_DBPR*b*(yR-1))*log(1-2.0_DBPR*b) )
        
        tide = (2.0_DBPR/3.0_DBPR)*k2*(r_fin/m_fin)**5

    end if

    !Returning------------------------------------------------------------------------------------------------------------------
    
    allocate(MATRIX(8,size(r_array)))
    
    MATRIX(1,:) = r_array
    MATRIX(2,:) = m_array
    MATRIX(3,:) = p_array
    MATRIX(4,:) = e_array
    MATRIX(5,:) = de_dp
    MATRIX(6,:) = lambda_array
    MATRIX(7,:) = phi_array
    MATRIX(8,:) = dphi_array

    deallocate(r_array, m_array, phi_array, e_array, de_dp, dphi_array, lambda_array)
    
    m_fin  = m_fin*1000/1.4766_DBPR    ! Output mass in units of solar mass
    r_fin  = r_fin*1000                ! Output radius in kilometers

!-------------------------------------------------------------------------------------------------------------------------------

contains

    subroutine TOV_lsoda(self_, Neq_, r_, x_, Ydash_, ierr_)
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: Neq_
        real(DBPR), intent(in) :: r_
        real(DBPR), intent(in) :: x_(Neq_)
        real(DBPR), intent(out) :: Ydash_(Neq_)
        integer, intent(out) :: ierr_
        real(DBPR) :: e_, p_
        ierr_ = 0
        self_%iwork = self_%iwork
        
        if(x_(2)> p_data(size(p_data))) then
            p_ = p_data(size(p_data))
        else
            p_ = x_(2)
        end if

        call interpolate(p_data, e_data, p_, e_)
        
        Ydash_ = TOV_func(e_, r_, x_)
        
    end subroutine TOV_lsoda
    
    !---------------------------------------------------------------------------------------------------------------------------

    subroutine surface(self_, neq_, r_, x_, ng, gout, ierr_) !The root function, ng is number of roots
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: neq_, ng
        real(DBPR), intent(in) :: r_, x_(neq_)
        integer, intent(out) :: ierr_
        real(DBPR), intent(out) :: gout(ng)
        real(DBPR) :: temp_
        ierr_ = 0
        temp_ = r_
        self_%iwork = self_%iwork

        gout(1) = x_(2)-p_limit
        
    end subroutine surface

    !---------------------------------------------------------------------------------------------------------------------------

    subroutine discont(self_, neq_, r_, x_, ng, gout, ierr_) !The root function, ng is number of roots
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: neq_, ng
        real(DBPR), intent(in) :: r_, x_(neq_)
        integer, intent(out) :: ierr_
        real(DBPR), intent(out) :: gout(ng)
        real(DBPR) :: temp_, e_
        ierr_ = 0
        temp_ = r_
        self_%iwork = self_%iwork

        call interpolate(p_data, e_data, x_(2), e_)
        gout(1) = e_-e_disc
        
    end subroutine discont

    !---------------------------------------------------------------------------------------------------------------------------

    subroutine back_interp(h_in_, r_array_, m_array_, p_array_, e_array_, phi_array_)
        implicit none
        real(DBPR), intent(in) :: h_in_
        real(DBPR), allocatable, intent(inout) :: r_array_(:), m_array_(:), p_array_(:), e_array_(:), phi_array_(:)
        integer :: i_
        real(DBPR) :: h_, r_, m_, p_, e_, phi_
        real(DBPR), allocatable ::  terp_r(:), terp_m(:), terp_p(:), terp_e(:), terp_phi(:)
        
        allocate(terp_r(0), terp_m(0), terp_p(0), terp_e(0), terp_phi(0))

        terp_r   = r_array_
        terp_m   = m_array_
        terp_p   = p_array_
        terp_e   = e_array_
        terp_phi = phi_array_
        i_       = size(terp_r)

        if (allocated(r_array_)) deallocate(r_array_, m_array_, p_array_, e_array_, phi_array_)
        allocate(r_array_(0), m_array_(0), p_array_(0), e_array_(0), phi_array_(0))

        h_ = -h_in_
        r_ = terp_r(i_)
    
        r_array_   = [r_]
        m_array_   = [terp_m(i_)]
        p_array_   = [terp_p(i_)]
        phi_array_ = [terp_phi(i_)]
        e_array_   = [terp_e(i_)]

        do while(r_>h_in_)
            
            r_ = r_+h_
            call interpolate(terp_r, terp_m,   r_, m_)
            call interpolate(terp_r, terp_p,   r_, p_)
            call interpolate(terp_r, terp_phi, r_, phi_)
            call interpolate(terp_r, terp_e,   r_, e_)
            
            if (isnan(m_*p_*phi_*e_)) then
            else    
                r_array_   = [r_array_, r_]
                m_array_   = [m_array_, m_]
                p_array_   = [p_array_, p_]
                phi_array_ = [phi_array_, phi_]
                e_array_   = [e_array_, e_]
            end if

        end do

        if (mod(size(r_array_),2)==0) then
            i_ = size(r_array_)
        else
            i_ = size(r_array_)-1
        end if
        
        r_array_   = r_array_(i_:1:-1)
        m_array_   = m_array_(i_:1:-1)
        p_array_   = p_array_(i_:1:-1)
        phi_array_ = phi_array_(i_:1:-1)
        e_array_   = e_array_(i_:1:-1)

    end subroutine back_interp

end subroutine SolveTOV_lsoda

!Discontinuity Correction-------------------------------------------------------------------------------------------------------

subroutine init_discont(h_in, p_data, e_data, ec, e_disc, r_disc, r_array, m_array, p_array, phi_array, e_array)

    real(DBPR), intent(in) :: h_in, e_disc, ec, p_data(:), e_data(:)
    real(DBPR), intent(out) :: r_disc
    real(DBPR), allocatable, intent(out) :: r_array(:), m_array(:), p_array(:), phi_array(:), e_array(:)
    integer :: i
    real(DBPR) :: h, r, pc, p, p_disc, mb_fin, e, e_prev, m, phi, x_prev(4), x(4), x_disc(4)
    real(DBPR), allocatable :: terp_r_data(:), terp_e_data(:), terp_p_data(:), terp_m_data(:), terp_phi_data(:)
    
    !---------------------------------------------------------------------------------------------------------------------------
    !This part finds radius at which discontinuity occurs
    
    call interpolate(e_data, p_data, ec, pc)
    
    h= 1E-5_DBPR
    r= 1E-15_DBPR
    x= [real(0, DBPR), pc, real(0, DBPR), real(0, DBPR)]
    p= x(2)
    e=ec
    777 do while(e>e_disc)

        e_prev = e
        call interpolate(p_data, e_data, p, e)
        
        if (e<e_disc) then
            e = e_prev
            r = r-h
            x = x_prev
            p = x(2)
            h = h/10.
            if (abs(e-e_disc)<1E-10_DBPR) then
                exit
            else 
                goto 777
            end if
        end if
        
        x_prev = x
        call rk4_(h, e, r, x)
        r      = r+h
        p      = x(2)
        mb_fin = x(4)

    end do
    
    r_disc  = r
    p_disc  = p
    x_disc  = x
    
    !---------------------------------------------------------------------------------------------------------------------------
    !This part does a forward RK4 from r=1E-15 to radius of discontinuity to generate
    !the data as a function of r. This data will be used for in'terp'olation later

    h = h_in
    r = 1E-15_DBPR
    x = [real(0, DBPR), pc, real(0, DBPR), real(0, DBPR)]
    p = x(2)
    e = ec

    terp_r_data   = [r]
    terp_m_data   = [x(1)]
    terp_p_data   = [p]
    terp_phi_data = [x(3)]
    terp_e_data   = [e]

    do while(r<r_disc)
        
        call rk4_(h, e, r, x)
        r = r+h
        p = x(2)
        
        call interpolate(p_data, e_data, p, e)

        terp_r_data   = [terp_r_data, r]
        terp_m_data   = [terp_m_data, x(1)]
        terp_p_data   = [terp_p_data, p]
        terp_phi_data = [terp_phi_data, x(3)]
        terp_e_data   = [terp_e_data, e]

    end do
    
    i                = size(terp_r_data)
    terp_r_data(i)   = r_disc
    terp_m_data(i)   = x_disc(1)
    terp_p_data(i)   = p_disc
    terp_phi_data(i) = x_disc(3)
    terp_e_data(i)   = e_disc 

    !---------------------------------------------------------------------------------------------------------------------------
    !This part does a backwards interpolation from discontinuity to a point close to next point 
    !after r=1E-15, using the previously generated data to interpolate. It then
    !returns all the initial values, which will then be used in SolveTOV
    !Doing this ensures we have a point, where radius point is the discontinuity point.
    
    h = -h_in
    r = r_disc
    x = x_disc
    p = x_disc(2)
    e = e_disc
    
    r_array   = [r]
    m_array   = [x(1)]
    p_array   = [p]
    phi_array = [x(3)]
    e_array   = [e]

    do while(r>h_in)
        
        r = r+h
        call interpolate(terp_r_data, terp_m_data, r, m)
        call interpolate(terp_r_data, terp_p_data, r, p)
        call interpolate(terp_r_data, terp_phi_data, r, phi)
        call interpolate(terp_r_data, terp_e_data, r, e)
        
        if (isnan(m*p*phi*e)) then
        else    
            r_array   = [r_array, r]
            m_array   = [m_array, m]
            p_array   = [p_array, p]
            phi_array = [phi_array, phi]
            e_array   = [e_array, e]
        end if

    end do
    
    !---------------------------------------------------------------------------------------------------------------------------
    !This part ensures that the r_disc point, starting from r= 1E-15, will be an odd point, and f-mode rk4
    !will be able to catch the point. It then inverts all the arrays and passes it on to SolveTOV
    
    i = size(r_array)
    if (mod(size(r_array),2)==0) then
        i = size(r_array)
    else
        i = size(r_array)-1
    end if
    
    r_array   = r_array(i:1:-1)
    m_array   = m_array(i:1:-1)
    p_array   = p_array(i:1:-1)
    phi_array = phi_array(i:1:-1)
    e_array   = e_array(i:1:-1)
    
end subroutine init_discont

!RK4 TOV------------------------------------------------------------------------------------------------------------------------

subroutine rk4_(h_, e_, r_, x_)

    real(DBPR), intent(in) :: h_, e_, r_
    real(DBPR), dimension(3), intent(inout) :: x_
    real(DBPR), dimension(3):: k1_, k2_, k3_, k4_
    
    k1_ = h_* TOV_func (e_, r_, x_)
    k2_ = h_* TOV_func (e_, r_+0.5_DBPR*h_, x_+0.5_DBPR*k1_)
    k3_ = h_* TOV_func (e_, r_+0.5_DBPR*h_, x_+0.5_DBPR*k2_)
    k4_ = h_* TOV_func (e_, r_+h_, x_+k3_)

    x_  = x_+ (1.0_DBPR/6.0_DBPR)*(k1_+2*k2_+2*k3_+k4_)
    
end subroutine rk4_

!RK4 TIDE-----------------------------------------------------------------------------------------------------------------------
    
subroutine rk4_tide(r_, e_, p_, de_dp_, lamb_, dphi_, y_)

    real(DBPR), intent(in) :: r_(:), e_(:), p_(:), de_dp_(:), lamb_(:), dphi_(:)
    real(DBPR), intent(out) :: y_
    integer :: i
    real(DBPR) :: h_, k1_, k2_, k3_, k4_

    y_ = 2.0_DBPR

    do i= 1, size(r_)-1, 2

        h_ = r_(i+2)-r_(i)

        k1_ = h_* tidal_func (r_(i), y_, e_(i), p_(i), de_dp_(i), lamb_(i), dphi_(i))
        k2_ = h_* tidal_func (r_(i+1), y_+0.5_DBPR*k1_, e_(i+1), p_(i+1), de_dp_(i+1), lamb_(i+1), dphi_(i+1))
        k3_ = h_* tidal_func (r_(i+1), y_+0.5_DBPR*k2_, e_(i+1), p_(i+1), de_dp_(i+1), lamb_(i+1), dphi_(i+1))
        k4_ = h_* tidal_func (r_(i+2), y_+k3_, e_(i+2), p_(i+2), de_dp_(i+2), lamb_(i+2), dphi_(i+2))

        y_  = y_+ (1.0_DBPR/6.0_DBPR)*(k1_+2.0_DBPR*k2_+2.0_DBPR*k3_+k4_)

    end do
    
end subroutine rk4_tide

!-------------------------------------------------------------------------------------------------------------------------------

end module TOV