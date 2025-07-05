program main
   use IO
   use TOV
   use COWLING
   use FULLGR
   use omp_lib
   implicit none

!------------------------------------------------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !                                                                                       !
   ! EOS file should have the first column as pressure and second one as energy density    !
   ! third column should be number density                                                 !
   !                                                                                       !
   ! 1 mev/fm^3 = 1.7827 × 10^12 g/cm^3                                                    !
   ! 1 mev/fm^3 = 1.6022 × 10^33 dyne/cm^2                                                 !
   ! 1 g/cm^3= 7.4237 × 10^(-19) /km^2                                                     !
   ! 1 dyne/cm^3= 8.2601 × 10^(-40) /km^2                                                  !
   !                                                                                       !
   ! From these we have the following relations: (/Mm^2= 10^6 /km^2)                       !
   ! 1 mev/fm^3= 1.323423 /Mm^2 (energy density in Megameters= 10^3 km)                    !
   ! 1 mev/fm^3= 1.323423 /Mm^2 (pressure in Megameters)                                   !
   ! Mass of sun in km= 1.4766 km                                                          !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!

   !----------------------------------------------!
   !                                              !
   ! ω [1/Mm]= 299.792458 ω [1/sec]               !
   ! ν= ω/2π                                      !
   ! ν [kHz]= 299.792458/2000π ω [1/Mm]           !
   ! ν [kHz]= 47.71345159236942258889E-3 ω [1/Mm] !
   !                                              !
   !----------------------------------------------!

   ! The ODE integrator can only give accurate results within the relative tolerance. That is, say rel tol is 1E-6. Then
   ! if I find the incoming GW amplitude for frequency w and w+h where h<1E-6, I will get inaccurate results
   ! since h is lower than the integration tolerance. So I can only safely find results for freqeuncies separated by more
   ! than the relative tolerance.

!------------------------------------------------------------------------------------------------------------------------------!

!$ integer :: nThreads
   character(len=1) :: creturn
   character(len= 128) :: directory, filename, input_filename, output_filename, fmt
   character(len=250) :: header
   integer :: ID, t_start, t_stop, clock_max, clock_rate, file_counter, &
      first_file, total_files, i, j, start_, end_, disc_start, disc_end
   real(DBPR), allocatable :: p_data(:), e_data(:), cs_data(:), MATRIX(:,:), OUT_MAT(:,:)
   real(DBPR) :: h, Mass, Radius, r_lim, hcalc, ecalc(100), pcalc(100), ec, pc, e_disc, p_disc, r_disc, tide, masses(2)
   real(DBPR) :: in_guess_cowl(3), f_cowl, p_cowl, g_cowl, f_GR, f_damp, p1_GR, p1_damp, g_GR, g_damp
   character(len=128), allocatable :: filenames(:)
   character(len=128) :: input_line
   logical :: IS_PE, IS_MEV, IS_STOP, IS_TIDE, IS_COWL, IS_GRG, IS_GRF, IS_GRP1

!------------------------------------------------------------------------------------------------------------------------------!
!$ nThreads = omp_get_max_threads()
!$ call omp_set_num_threads(nThreads)

   directory = 'EOS'
   
   allocate(filenames(0))
   open(1, file='EOS_inputs.txt', status='old', action='read')
   do
      read(1, '(A)', iostat=i) input_line
      if (i /= 0) exit
      filenames = [filenames, input_line(1:len_trim(input_line))]
   end do
   close(1)
   first_file  = 1
   total_files = size(filenames)
   
   file_counter = 0

   ! IS_PE:  EOS is Pressure-Energy density (TRUE) or Energy Density-Pressure
   ! IS_MEV: EOS is in Mev/fm^3 (TRUE) or CGS units (FALSE)
   ! IS_STOP: Will stop calculating beyond maximum mass (TRUE)

   open(10, file='Params.txt', status='OLD', action='READ')
      read(10, *) IS_PE
      read(10, *) IS_MEV
      read(10, *) IS_STOP
      read(10, *) IS_TIDE
      read(10, *) IS_COWL
      read(10, *) IS_GRG
      read(10, *) IS_GRF
      read(10, *) IS_GRP1
   close(10)

!------------------------------------------------------------------------------------------------------------------------------!

   call system_clock (t_start, clock_rate, clock_max)

   write(*,*)

   !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(directory, filenames, file_counter, first_file, total_files, &
   !$OMP& IS_PE, IS_MEV, IS_STOP, IS_TIDE, IS_COWL, IS_GRG, IS_GRF, IS_GRP1) SCHEDULE(DYNAMIC)

   do ID = first_file, total_files

!------------------------------------------------------------------------------------------------------------------------------!

      !$OMP CRITICAL
      file_counter = file_counter + 1
      !$OMP END CRITICAL

      creturn = achar(13)  ! Carriage return

      filename = filenames(ID)(1:len_trim(filenames(ID)))

      call make_filename(directory, directory(1:len_trim(directory))//'_DATA', filename, 'data', input_filename, output_filename)

!------------------------------------------------------------------------------------------------------------------------------!

      e_disc = -10.0_DBPR
      p_disc = -10.0_DBPR

      if (IS_PE) then
         call read_eos(ID, input_filename, p_data, e_data, cs_data)
      else
         call read_eos(ID, input_filename, e_data, p_data, cs_data)
      end if
      
      if (.not. IS_MEV) then ! EOS is in CGS units. So convert to MeV/fm^3
         e_data = e_data/(1.78266E12_DBPR)
         p_data = p_data/(1.60218E33_DBPR)
      end if

      if (cs_data(1)>0) then
         call continuity_check(p_data, cs_data, disc_start, disc_end)
         if (disc_start>0 .and. disc_end>0) then
            e_disc = e_data(disc_end)
            p_disc = p_data(disc_end)
         else
            e_disc = -10.0_DBPR
            p_disc = -10.0_DBPR
         end if
      end if

      if (e_disc>0 .and. p_disc>0) then

         ! write(*, '(/a, a, f10.5, a, f10.5, a, f10.5, a,/a)') '', 'EOS has density discontinuity at: Energy Density= ', &
         !    e_data(disc_start), ' <-> ', e_data(disc_end),' [MeV/fm^3] and Pressure= ', p_disc, ' [MeV/fm^3]'

         e_disc = e_disc*1.323423
         p_disc = p_disc*1.323423

      else

         r_disc = -10.0_DBPR

      end if

      ! Converting from Mev/fm^3 to /Mm^2
      p_data    = p_data*1.323423_DBPR
      e_data    = e_data*1.323423_DBPR

      r_lim = 20         ! The desired limiting radius in [Km] beyond which MR won't be calculated
      h     = 5E-6_DBPR  !step-size for TOV Grid

      ! if (ID==0) write(*,*)

!------------------------------------------------------------------------------------------------------------------------------!

      start_ = 0
      do i = 1, size(e_data)
         if (e_data(i)/1.323423>150.) then
            start_= i
            exit
         end if
      end do

      end_ = 0
      do i = 1, size(e_data)
         if (e_data(i)/1.323423>3000.) then
            end_= i
            exit
         end if
      end do

      if(end_==0) end_ = size(e_data)

      if (e_disc==-10.0_DBPR) then
         hcalc = (log10(e_data(end_))-log10(e_data(start_)))/99.0_DBPR
         do i=1,100
            ecalc(i) = log10(e_data(start_)) + hcalc*(i-1)
         end do
         ecalc = 10.0_DBPR**ecalc
         do i=1,100
            call interpolate(e_data, p_data, ecalc(i), pcalc(i))
         end do
      else
         hcalc = (log10(e_data(disc_start))-log10(e_data(start_)))/49.0_DBPR
         do i=1,50
            ecalc(i) = log10(e_data(start_)) + hcalc*(i-1)
         end do
         ecalc(51) = log10(0.5_DBPR*(e_data(disc_start)+e_data(disc_end)))
         hcalc = (log10(e_data(end_))-log10(e_data(disc_end)))/48.0_DBPR
         do i=52,100
            ecalc(i) = log10(e_data(disc_end)) + hcalc*(i-52)
         end do
         ecalc = 10.0_DBPR**ecalc
         do i=1,100
            call interpolate(e_data, p_data, ecalc(i), pcalc(i))
         end do
      end if

      start_ = 1
      end_   = 100
      
!------------------------------------------------------------------------------------------------------------------------------!

      j      = 4
      header = 'ε_c [MeV/fm^3],p_c [MeV/fm^3],R [KM],M [M_sun]'

      if(IS_TIDE) then
         j      = j+1
         header = trim(header)//',Λ'
      end if
      if(IS_COWL) then
         j      = j+3
         header = trim(header)//',g-cowl [kHz],f-cowl [kHz],p1-cowl [kHz]'
      end if
      if(IS_GRG) then
         j      = j+2
         header = trim(header)//',g-GR [kHz],τ_g [s]'
      end if
      if(IS_GRF) then
         j      = j+2
         header = trim(header)//',f-GR [kHz],τ_f [s]'
      end if
      if(IS_GRP1) then
         j      = j+2
         header = trim(header)//',p1-GR [kHz],τ_p1 [s]'
      end if
      
      allocate(OUT_MAT(end_-start_, j))
      
      masses        = 0.0_DBPR
      in_guess_cowl = [-10, -10, -10]

      do i= start_, end_

         write(*, fmt= 1, advance= 'no') creturn, 'Files : ', file_counter, ' of ', (total_files-first_file+1)

         ec = ecalc(i)
         pc = pcalc(i)

         if (IS_TIDE) then
            call SolveTOV_lsoda(p_data, e_data, h, ec, r_lim, Radius, Mass, e_disc, r_disc, MATRIX, tide= tide)
         else
            call SolveTOV_lsoda(p_data, e_data, h, ec, r_lim, Radius, Mass, e_disc, r_disc, MATRIX)
         end if
         
         if (Radius> r_lim) then
            OUT_MAT(i-start_+1, :)= NaN/Nan

!------------------------------------------------------------------------------------------------------------------------------!

         else if (Radius>0.0_DBPR .and. Radius<=r_lim .and. Mass>0.0_DBPR) then

            if (IS_STOP) then
               if (masses(2)<masses(1) .and. Mass>1.80_DBPR) then
                  end_ = i-1
                  exit
               end if
            end if

            masses = [masses(2), Mass]

            if (IS_COWL) then

               call cowling_frequency(MATRIX, in_guess_cowl, r_disc, f_cowl, p_cowl, g_cowl)

               if (isnan(f_cowl) .or. isnan(p_cowl) .or. isnan(g_cowl)) then
                  in_guess_cowl = [-10, -10, -10]/47.71345159236942258889E-3_DBPR
               else
                  in_guess_cowl = [f_cowl, p_cowl, g_cowl]/47.71345159236942258889E-3_DBPR
               end if

               if (g_cowl<0.0_DBPR) then
                  g_cowl= NaN/NaN
               end if

            end if

!------------------------------------------------------------------------------------------------------------------------------!

            if (Mass>0.5) then

               if (IS_GRF)  call frequency([.TRUE.,.FALSE.,.FALSE.], MATRIX, r_disc, f_GR, f_damp)
               if (IS_GRP1) call frequency([.FALSE.,.TRUE.,.FALSE.], MATRIX, r_disc, p1_GR, p1_damp)
               if (IS_GRG)  call frequency([.FALSE.,.FALSE.,.TRUE.], MATRIX, r_disc, g_GR, g_damp)
            else
               f_GR    = NaN/NaN
               f_damp  = NaN/NaN
               p1_GR   = NaN/NaN
               p1_damp = NaN/NaN
            end if

            ! print*, ec/1.323423, Radius, Mass, f_GR, g_GR

            write(*,fmt='(a,f20.10)', advance='no') ' :: Running: ', ec/1.323423_DBPR

            OUT_MAT(i-start_+1, 1:4)= [ec/1.323423_DBPR, pc/1.323423_DBPR, Radius, Mass]
            j = 4
            if(IS_TIDE) then
               j = j+1
               OUT_MAT(i-start_+1, j)= tide
            end if
            if(IS_COWL) then
               j = j+3
               OUT_MAT(i-start_+1, j-2:j)= [g_cowl, f_cowl, p_cowl]
            end if
            if(IS_GRG) then
               j = j+2
               OUT_MAT(i-start_+1, j-1:j)= [g_GR, g_damp]
            end if
            if(IS_GRF) then
               j = j+2
               OUT_MAT(i-start_+1, j-1:j)= [f_GR, f_damp]
            end if
            if(IS_GRP1) then
               j = j+2
               OUT_MAT(i-start_+1, j-1:j)= [p1_GR, p1_damp]
            end if

         end if

      end do

!------------------------------------------------------------------------------------------------------------------------------!

1     format (a, a, i9, a, i9)

      open(unit=ID+7, file= output_filename, status= 'replace', action= 'write')

      write(ID+7, '(a)') trim(header)

      do i= 2, end_-start_
         if (IS_TIDE) then
            write(fmt, '(i2)') size(OUT_MAT, 2)-6
            fmt = '(4(f20.10,","), E20.10,",", '//trim(fmt)//'(f20.10,","), f20.10)'
            if (.not.isnan(OUT_MAT(i,1))) write(ID+7,fmt) OUT_MAT(i, :)
         else
            write(fmt, '(i2)') size(OUT_MAT, 2)-1
            fmt = '('//trim(fmt)//'(f20.10,","), f20.10)'
            if (.not.isnan(OUT_MAT(i,1))) write(ID+7,fmt) OUT_MAT(i, :)
         end if
      end do

      close(ID+7)

      deallocate(OUT_MAT)
!------------------------------------------------------------------------------------------------------------------------------!

   end do

   !$OMP END PARALLEL DO

   call system_clock (t_stop, clock_rate, clock_max)
   write(*,*)
   write(*,*)
   write(*,*) 'Elapsed time is: ', real(t_stop-t_start)/real(clock_rate), 'seconds'

!------------------------------------------------------------------------------------------------------------------------------!

end program main
