subroutine da_setup_be_regional(xb, be, grid)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   !
   ! Updates:
   !
   !       Implementation of multi-variate BE
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type (xb_type), intent(in)    :: xb                    ! First guess structure.
   type (be_type), intent(inout) :: be                    ! Back. errors structure.
   type (domain),  intent(in)    :: grid

   logical                do_normalize1      ! Test do_normalize.
   integer             :: i, j, k, m         ! Loop counters.
   integer, allocatable:: bin(:,:,:)         ! Bin assigned to each 3D point
   integer, allocatable:: bin2d(:,:)         ! Bin assigned to each 2D point
   integer             :: bin_type           ! Type of bin to average over.
   integer             :: num_bins           ! Number of bins (3D fields).
   integer             :: num_bins2d         ! Number of bins (3D fields).
   real*8  :: junk                           ! For reading truncated sd, wsd
   real*8  :: lat_min, lat_max, binwidth_lat ! Used if bin_type = 2 (degrees)..
   real*8  :: hgt_min, hgt_max, binwidth_hgt ! Used if bin_type = 2 (m). .

   integer :: num_cv_3d_basic  ! number of 3D mandatory control variables
   integer :: num_cv_3d_extra  ! number of 3D optional control variables
   integer :: num_cv_2d        ! number of 2D control variables
   integer :: num_cv_qn        ! number of number concentration control variables   

   real*8, allocatable           :: be1_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be2_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be3_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be4_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be5_eval_loc(:,:)     ! Temp arrays.

   real*8, allocatable           :: be6_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be7_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be8_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be9_eval_loc(:,:)     ! Temp arrays.
   real*8, allocatable           :: be10_eval_loc(:,:)    ! Temp arrays.
   real*8, allocatable           :: be11_eval_loc(:,:)    ! Temp arrays.
   real*8, allocatable           :: be12_eval_loc(:,:)    ! Temp arrays.  !zhengh20221014

   real*8, allocatable           :: be1_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be2_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be3_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be4_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be5_eval_glo(:)       ! Global Eigenvalues.

   real*8, allocatable           :: be6_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be7_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be8_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be9_eval_glo(:)       ! Global Eigenvalues.
   real*8, allocatable           :: be10_eval_glo(:)      ! Global Eigenvalues.
   real*8, allocatable           :: be11_eval_glo(:)      ! Global Eigenvalues.
   real*8, allocatable           :: be12_eval_glo(:)      ! Global Eigenvalues. !zheng20221014

   real*8, allocatable           :: alpha_val(:)          ! Global Eigenvalues.

   real*8, allocatable           :: be1_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be2_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be3_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be4_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be5_evec_loc(:,:,:)   ! Local Eigenvectors.

   real*8, allocatable           :: be6_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be7_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be8_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be9_evec_loc(:,:,:)   ! Local Eigenvectors.
   real*8, allocatable           :: be10_evec_loc(:,:,:)  ! Local Eigenvectors.
   real*8, allocatable           :: be11_evec_loc(:,:,:)  ! Local Eigenvectors.
   real*8, allocatable           :: be12_evec_loc(:,:,:)  ! Local Eigenvectors. !zhengh20221014

   real*8, allocatable           :: be1_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be2_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be3_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be4_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be5_evec_glo(:,:)     ! Global Eigenvectors.

   real*8, allocatable           :: be6_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be7_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be8_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be9_evec_glo(:,:)     ! Global Eigenvectors.
   real*8, allocatable           :: be10_evec_glo(:,:)    ! Global Eigenvectors.
   real*8, allocatable           :: be11_evec_glo(:,:)    ! Global Eigenvectors.
   real*8, allocatable           :: be12_evec_glo(:,:)    ! Global Eigenvectors. !zhengh20221014

   real*8, allocatable           :: alpha_evec(:,:)       ! Global Eigenvectors.

   real*8, allocatable           :: be1_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be2_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be3_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be4_rf_lengthscale(:) ! RF lengthscale.
   real*8, allocatable           :: be5_rf_lengthscale(:)

   real*8, allocatable           :: be6_rf_lengthscale(:)
   real*8, allocatable           :: be7_rf_lengthscale(:)
   real*8, allocatable           :: be8_rf_lengthscale(:)
   real*8, allocatable           :: be9_rf_lengthscale(:)
   real*8, allocatable           :: be10_rf_lengthscale(:)
   real*8, allocatable           :: be11_rf_lengthscale(:)
   real*8, allocatable           :: be12_rf_lengthscale(:) !zhengh20221014

   real*8, allocatable           :: alpha_rf_lengthscale(:)
   real*8, allocatable           :: alpha_rf_scale_factor(:)

   real*8, allocatable           :: evec_loc(:,:,:)        ! Latitudinally varying eigenvectors.
   real*8, allocatable           :: eval_loc(:,:)          ! Latitudinally varying eigenvalues.

   character*5                 :: namv(5)=(/"psi  ","chi_u","t_u  ","rh   ","ps_u "/)
   character*10                :: variable
   character*80                :: var80c
   character(len=15)           :: write_fmt
   integer                     :: ni, nj, nk, nk_2d, b
   integer                     :: ix, jy, kz, kzs(5), mz, mzs(5)
   real, allocatable   :: regcoeff_psi_chi(:)        ! psi/chi    regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! psi/t      regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! psi/ps     regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! psi/rh     regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! chi_u/t    regression coefficient
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! chi_u/ps   regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! chi_u/rh   regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! t_u/rh     regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! ps_u/rh    regression coefficient

   !Set other regcoeffs, Zheng, 20200513
   real, allocatable   :: regcoeff_psi_qcloud(:,:,:)             ! psi/qcloud         regression cooefficient.
   real, allocatable   :: regcoeff_psi_qice(:,:,:)               ! psi/qice           regression cooefficient.
   real, allocatable   :: regcoeff_psi_qrain(:,:,:)              ! psi/qrain          regression cooefficient.
   real, allocatable   :: regcoeff_psi_qsnow(:,:,:)              ! psi/qsnow          regression cooefficient.
   real, allocatable   :: regcoeff_psi_qgraup(:,:,:)             ! psi/qgraup         regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qcloud(:,:,:)           ! chi_u/qcloud       regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qice(:,:,:)             ! chi_u/qice         regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qrain(:,:,:)            ! chi_u/qrain        regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qsnow(:,:,:)            ! chi_u/qsnow        regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qgraup(:,:,:)           ! chi_u/qgraup       regression cooefficient.
   real, allocatable   :: regcoeff_t_u_ps(:,:)                   ! t_u/ps             regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qcloud(:,:,:)             ! t_u/qcloud         regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qice(:,:,:)               ! t_u/qice           regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qrain(:,:,:)              ! t_u/qrain          regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qsnow(:,:,:)              ! t_u/qsnow          regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qgraup(:,:,:)             ! t_u/qgraup         regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qcloud(:,:)              ! ps_u/qcloud        regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qice(:,:)                ! ps_u/qice          regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qrain(:,:)               ! ps_u/qrain         regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qsnow(:,:)               ! ps_u/qsnow         regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qgraup(:,:)              ! ps_u/qgraup        regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qcloud(:,:,:)            ! rh_u/qcloud        regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qice(:,:,:)              ! rh_u/qice          regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qrain(:,:,:)             ! rh_u/qrain         regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qsnow(:,:,:)             ! rh_u/qsnow         regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qgraup(:,:,:)            ! rh_u/qgraup        regression cooefficient.
   real, allocatable   :: regcoeff_qcloud_u_qice(:,:,:)          ! qcloud_u/qice      regression cooefficient.
   real, allocatable   :: regcoeff_qcloud_u_qrain(:,:,:)         ! qcloud_u/qrain     regression cooefficient.
   real, allocatable   :: regcoeff_qcloud_u_qsnow(:,:,:)         ! qcloud_u/qsnow     regression cooefficient.
   real, allocatable   :: regcoeff_qcloud_u_qgraup(:,:,:)        ! qcloud_u/qgraup    regression cooefficient.
   real, allocatable   :: regcoeff_qice_u_qrain(:,:,:)           ! qice_u/qrain       regression cooefficient.
   real, allocatable   :: regcoeff_qice_u_qsnow(:,:,:)           ! qice_u/qsnow       regression cooefficient.
   real, allocatable   :: regcoeff_qice_u_qgraup(:,:,:)          ! qice_u/qgraup      regression cooefficient.
   real, allocatable   :: regcoeff_qrain_u_qsnow(:,:,:)          ! qrain_u/qsnow      regression cooefficient.
   real, allocatable   :: regcoeff_qrain_u_qgraup(:,:,:)         ! qrain_u/qgraup     regression cooefficient.
   real, allocatable   :: regcoeff_qsnow_u_qgraup(:,:,:)         ! qsnow_u/qgraup     regression cooefficient.
   real, allocatable   :: regcoeff_u_v(:,:,:)                    ! u/v                regression cooefficient.
   real, allocatable   :: regcoeff_u_t(:,:,:)                    ! u/t                regression cooefficient.
   real, allocatable   :: regcoeff_u_ps(:,:)                     ! u/ps               regression cooefficient.
   real, allocatable   :: regcoeff_u_rh(:,:,:)                   ! u/rh               regression cooefficient.
   real, allocatable   :: regcoeff_u_qcloud(:,:,:)               ! u/qcloud           regression cooefficient.
   real, allocatable   :: regcoeff_u_qice(:,:,:)                 ! u/qice             regression cooefficient.
   real, allocatable   :: regcoeff_u_qrain(:,:,:)                ! u/qrain            regression cooefficient.
   real, allocatable   :: regcoeff_u_qsnow(:,:,:)                ! u/qsnow            regression cooefficient.
   real, allocatable   :: regcoeff_u_qgraup(:,:,:)               ! u/qgraup           regression cooefficient.
   real, allocatable   :: regcoeff_v_u_t(:,:,:)                  ! v_u/t              regression cooefficient.
   real, allocatable   :: regcoeff_v_u_ps(:,:)                   ! v_u/ps             regression cooefficient.
   real, allocatable   :: regcoeff_v_u_rh(:,:,:)                 ! v_u/rh             regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qcloud(:,:,:)             ! v_u/qcloud         regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qice(:,:,:)               ! v_u/qice           regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qrain(:,:,:)              ! v_u/qrain          regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qsnow(:,:,:)              ! v_u/qsnow          regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qgraup(:,:,:)             ! v_u/qgraup         regression cooefficient.


   !Set other regcoeffs, zhengh, 20221014  add qnrain
   real, allocatable   :: regcoeff_psi_qnrain(:,:,:)             ! psi/qnrain         regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_qnrain(:,:,:)           ! chi_u/qnrain       regression cooefficient.
   real, allocatable   :: regcoeff_t_u_qnrain(:,:,:)             ! t_u/qnrain         regression cooefficient.
   real, allocatable   :: regcoeff_ps_u_qnrain(:,:)              ! ps_u/qnrain        regression cooefficient.
   real, allocatable   :: regcoeff_rh_u_qnrain(:,:,:)            ! rh_u/qnrain        regression cooefficient.
   real, allocatable   :: regcoeff_qcloud_u_qnrain(:,:,:)        ! qcloud_u/qnrain    regression cooefficient.
   real, allocatable   :: regcoeff_qice_u_qnrain(:,:,:)          ! qice_u/qnrain      regression cooefficient.
   real, allocatable   :: regcoeff_qrain_u_qnrain(:,:,:)         ! qrain_u/qnrain     regression cooefficient.
   real, allocatable   :: regcoeff_qsnow_u_qnrain(:,:,:)         ! qsnow_u/qnrain     regression cooefficient.
   real, allocatable   :: regcoeff_qgraup_u_qnrain(:,:,:)        ! qgraup_u/qnrain     regression cooefficient.
   real, allocatable   :: regcoeff_u_qnrain(:,:,:)               ! u/qnrain           regression cooefficient.
   real, allocatable   :: regcoeff_v_u_qnrain(:,:,:)             ! v_u/qnrain         regression cooefficient.



   !real                :: qrain_th_low, qrain_th_high

   integer                     :: be_unit, ier, be_rf_unit, be_print_unit, it, idummy

!-----------for interpolating CV5--------------------------------------------------------------
   REAL, ALLOCATABLE           :: reg_psi_ps0(:,:), reg_psi_chi0(:,:), reg_psi_t0(:,:,:), &
                                  reg_psi_ps   (:), reg_psi_chi   (:), reg_psi_t   (:,:), &
                                  reg_psi_ps_be(:), reg_psi_chi_be(:), reg_psi_t_be(:,:)
   REAL, DIMENSION(:,:), ALLOCATABLE:: covm1_be, covm2_be, covm3_be, covm4_be,&
                                       covm1   , covm2   , covm3   , covm4
   REAL, DIMENSION(:),   ALLOCATABLE:: rfls1_be, rfls2_be, rfls3_be, rfls4_be
   REAL, DIMENSION(:),   ALLOCATABLE:: eta_be
! For write out the Interpolated CV5 BE (bin_type=5) file:
   integer              :: bin_type_out, num_bins_out, num_bins2d_out
   integer, allocatable :: bin_out(:,:,:), bin2d_out(:,:)
   integer, parameter  :: be_out_unit   = 135

!===========for interpolating CV5==============================================================

   if (trace_use) call da_trace_entry("da_setup_be_regional -- Zhengh 20221014")
!   if( (cv_options == 5) .or. (cv_options == 6) .or. (cv_options == 7) ) then
   if( (cv_options == 5) .or. (cv_options == 6) .or. (cv_options == 7) .or. (cv_options == 9) ) then		!Zheng,20200513
      write (unit=message(1),fmt='(A,I3)') &
         'Set up background errors for regional application for cv_options = ',cv_options
      write (unit=message(2),fmt='(a)') ' '
      call da_message(message(1:2))
      if ( (cv_options == 7 .or. cv_options == 6) .and. interpolate_stats ) then
         write (unit=message(1),fmt='(A,I3)') 'interpolate_stats is not implemented for cv_options=',cv_options
         call da_error(__FILE__,__LINE__,message(1:1))
      end if
      if ( ( cloud_cv_options >=2 ) .and. interpolate_stats ) then
         write (unit=message(1),fmt='(A,I3)') 'interpolate_stats is not implemented for cloud_cv_options=',cloud_cv_options
         call da_error(__FILE__,__LINE__,message(1:1))
      end if
      if ( ( qnr_option >=1 ) .and. interpolate_stats ) then                                                                     !zhengh20224014
         write (unit=message(1),fmt='(A,I3)') 'interpolate_stats is not implemented for qnr_options=',qnr_option
         call da_error(__FILE__,__LINE__,message(1:1))
      end if
   else 
!      write (unit=message(1),fmt='(a)') 'cv_options should be 5, 6, or 7'
      write (unit=message(1),fmt='(a)') 'cv_options should be 5, 6, 7, or 9'			!Zheng,20200521
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   ix = xb % mix
   jy = xb % mjy
   kz = xb % mkz

   be_rf_unit    = unit_end + 1
   be_print_unit = unit_end + 2

   call da_get_unit(be_unit)
   open(unit=be_unit,file="be.dat", status="old",form="unformatted")

   rewind (be_unit)
   read (be_unit, iostat=ier) ni, nj, nk
   print *, 'ni, nj, nk = ', ni, nj, nk
   if (ier /= 0) then
      write (unit=message(1),fmt='(a,i3)') 'Error in reading be.dat, unit= ',be_unit
      call da_error(__FILE__,__LINE__,message(1:1))
   end if
   read (be_unit) bin_type
   print *, 'bin_type = ', bin_type

!-----------for interpolating CV5--------------------------------------------------------------
   if ( .not. interpolate_stats ) then 
      if (ni /= ix .or. nj /= jy .or. nk /= kz) then
         write(unit=message(1),fmt='(a)')  &
            'Dimensions of the assimilation domain do not match those in the BE file.'
         write(unit=message(2),fmt='(3x,a,3i4)') "in be.dat: ni, nj, nk = ",ni, nj, nk
         write(unit=message(3),fmt='(3x,a,3i4)') "in fg:     ix, jy, kz = ",ix, jy, kz
         write(message(4),'(a)') 'be.dat need to be generated using data consistent with the assimilation domain.'
         call da_error(__FILE__,__LINE__,message(1:4))
      end if
   else  ! when interpolate_stats = true
      ! Must use the domain averaged Reg. Coeff.:
      lat_stats_option = .FALSE.
      ! Must use the global Eigenvector/eigenvalue
      vert_evalue = 1
      write (unit=message(1),fmt='(A,I3,A)') 'interpolate_stats is .true., BE for cv_options = ',cv_options,' will be interpolated.'
      call da_message(message(1:1))
      if (rootproc) write(be_out_unit) ix, jy, kz
      write(6,'(a/("k=",i3,2x,f10.4))') 'xb%sigmah=',(k,xb%sigmah(k),k=1,kz)
      allocate (eta_be(1:nk))
      read(be_unit,iostat=ier) eta_be
      !  If not available from BE file, get it from namelist
      if (ier /= 0) then
         write(message(1),'("ier=",i5,2x,a)') ier, " NO Eta values existed in BE file."
         write(message(2),'(a)') "Will use BE Eta values from the namelist be_eta."
         call da_message(message(1:2))
         rewind(be_unit)
         read(be_unit) idummy, idummy, idummy
         read(be_unit) idummy
         eta_be = be_eta(1:nk)
      endif
      write(6,'(a/("k=",i3,2x,f10.4))') 'eta_be=',(k,eta_be(k),k=1,nk)
      !  If eta values not defined correctly, stop
      if (eta_be(1) == eta_be(nk)) then
         write(message(1),'(a)') 'Wrong eta values found.'
         write(message(2),'(a)') 'The eta levels used in generating be.dat should be specified in namelist &wrfvar10 be_eta.'
         call da_error(__FILE__,__LINE__,message(1:2))
      endif
   end if

!===========for interpolating CV5==============================================================

   allocate (bin(1:ni,1:nj,1:nk))
   allocate (bin2d(1:ni,1:nj))

  !if(cloud_cv_options.eq.2)then
  ! read (be_unit)num_bins, num_bins2d
  ! read (be_unit)lat_min, lat_max, binwidth_lat
  ! read (be_unit)hgt_min, hgt_max, binwidth_hgt
  ! read (be_unit)qrain_th_low, qrain_th_high
  ! read (be_unit)bin(1:ni,1:nj,1:nk)
  ! read (be_unit)bin2d(1:ni,1:nj)
  !else
   read (be_unit)lat_min, lat_max, binwidth_lat
   read (be_unit)hgt_min, hgt_max, binwidth_hgt
   read (be_unit)num_bins, num_bins2d
   read (be_unit)bin(1:ni,1:nj,1:nk)
   read (be_unit)bin2d(1:ni,1:nj)
  !end if

   print *, lat_min, lat_max, binwidth_lat
   print *, hgt_min, hgt_max, binwidth_hgt
   print *, 'num_bins, num_bins2d = ', num_bins, num_bins2d
   print *, 'bin = ',   bin(1:1,1:1,1:1)
   print *, 'bin2d = ', bin2d(1:1,1:1)

   num_cv_3d_basic = 4 ! (Psi, Chi_u, T_u, Pseudo RH) or (U,V,T,Pseudo RH)
   num_cv_3d_extra = 0 ! default: no cloud and w
   num_cv_2d       = 1 ! Ps_u
   num_cv_qn       = 0 !                                        zhengh 20221017
   if ( cloud_cv_options == 2 ) then									!Zheng,20200514    !!!!!!!!!!!!!important!!!!!!!!!!!!!!!!!
!   if ( cloud_cv_options >= 2 ) then
      num_cv_3d_extra = 5 ! 5 hydrometeors, qc,qr,qi,qs,qg
      if (qnr_option >= 1)  then
         num_cv_qn = num_cv_qn + 1  ! qnr
      end if
   end if
   if ( use_cv_w ) then
      num_cv_3d_extra = num_cv_3d_extra + 1 ! clouds+w
   end if

!-----------for interpolating CV5--------------------------------------------------------------
! Write out the interpolated CV5 BE:
    if (interpolate_stats) then
 
! No matter which bin_type BE inputed, after Interpolated, the BE always
! became bin_type = 1 BE, i.e. the global BE.
     bin_type_out   = 5
     num_bins_out   = kz
     num_bins2d_out = 1
     if (rootproc) then 
     write(be_out_unit)bin_type_out
     write(be_out_unit)lat_min, lat_max, binwidth_lat
     write(be_out_unit)hgt_min, hgt_max, binwidth_hgt
     write(be_out_unit)num_bins_out, num_bins2d_out
     endif
      allocate (bin_out(1:ix,1:jy,1:kz))
      allocate (bin2d_out(1:ix,1:jy))
        do k = 1, kz
           bin_out (:,:,k) = k
        enddo
        bin2d_out = 1
      if (rootproc) then
      write(be_out_unit)bin_out
      write(be_out_unit)bin2d_out
      endif
      deallocate (bin_out)
      deallocate (bin2d_out)
    endif
!===========for interpolating CV5==============================================================
!
!   if ( cv_options /= 7 ) then   ! No regression coefficients for cv_options == 7

   ! 1.1 Read in regression coefficients
   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then		!Zheng,20200522
      allocate  (regcoeff_psi_chi(1:num_bins))
      allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
      allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
      allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
      if ( cv_options == 6 ) then
         allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
         allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
         allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
         ! allocate  (regcoeff_t_u_ps(1:nk,1:num_bins2d))
         allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))
      end if
      if ( cv_options == 9 .and. use_regcoeff ) then				!Zheng,20200513
         allocate (regcoeff_psi_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_ps(1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_rh(1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qcloud(1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qice(1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qrain(1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qsnow(1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qgraup(1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qcloud_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qcloud_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qcloud_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qcloud_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qice_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qice_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qice_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qrain_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qrain_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qsnow_u_qgraup(1:nk,1:nk,1:num_bins2d))
         if (qnr_option >= 1)  then                       !zhengh 20221014 add qnrain (reg to others) 10
            allocate (regcoeff_psi_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_chi_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_t_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_ps_u_qnrain(1:nk,1:num_bins2d))
            allocate (regcoeff_rh_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_qcloud_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_qice_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_qrain_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_qsnow_u_qnrain(1:nk,1:nk,1:num_bins2d))
            allocate (regcoeff_qgraup_u_qnrain(1:nk,1:nk,1:num_bins2d))

         end if
      end if
      
   else	if ( cv_options == 7 .and. use_regcoeff ) then						! cv_options == 7
      allocate (regcoeff_u_v(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_t(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_ps(1:nk,1:num_bins2d))
      allocate (regcoeff_u_rh(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_qcloud(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_qice(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_t(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_ps(1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_rh(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_qcloud(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_qice(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_v_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_ps(1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_qcloud(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_qice(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_t_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_rh(1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_qcloud(1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_qice(1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_qrain(1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_qsnow(1:nk,1:num_bins2d))
      allocate (regcoeff_ps_u_qgraup(1:nk,1:num_bins2d))
      allocate (regcoeff_rh_u_qcloud(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_rh_u_qice(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_rh_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_rh_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_rh_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qcloud_u_qice(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qcloud_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qcloud_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qcloud_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qice_u_qrain(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qice_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qice_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qrain_u_qsnow(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qrain_u_qgraup(1:nk,1:nk,1:num_bins2d))
      allocate (regcoeff_qsnow_u_qgraup(1:nk,1:nk,1:num_bins2d))
      if ( qnr_option >= 1) then   !zhengh 20221014 add qnrain (reg to others) 10

         allocate (regcoeff_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_v_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_ps_u_qnrain(1:nk,1:num_bins2d))
         allocate (regcoeff_rh_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qcloud_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qice_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qrain_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qsnow_u_qnrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_qgraup_u_qnrain(1:nk,1:nk,1:num_bins2d))

      end if

   end if

   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then		!Zheng,20200522
      regcoeff_psi_chi = 0.
      regcoeff_psi_t = 0.
      regcoeff_psi_ps = 0.
      regcoeff_psi_rh = 0.
      if ( cv_options == 6 ) then
         regcoeff_chi_u_t = 0.
         regcoeff_chi_u_ps = 0.
         regcoeff_chi_u_rh = 0.
         ! regcoeff_t_u_ps = 0.
         regcoeff_t_u_rh = 0.
         regcoeff_ps_u_rh = 0.
      end if
      if ( cv_options == 9 .and. use_regcoeff ) then				!Zheng,20200522
         regcoeff_psi_qcloud = 0.
         regcoeff_psi_qice = 0.
         regcoeff_psi_qrain = 0.
         regcoeff_psi_qsnow = 0.
         regcoeff_psi_qgraup = 0.
         regcoeff_chi_u_t = 0.
         regcoeff_chi_u_ps = 0.
         regcoeff_chi_u_rh = 0.
         regcoeff_chi_u_qcloud = 0.
         regcoeff_chi_u_qice = 0.
         regcoeff_chi_u_qrain = 0.
         regcoeff_chi_u_qsnow = 0.
         regcoeff_chi_u_qgraup = 0.
         regcoeff_t_u_ps = 0.
         regcoeff_t_u_rh = 0.
         regcoeff_t_u_qcloud = 0.
         regcoeff_t_u_qice = 0.
         regcoeff_t_u_qrain = 0.
         regcoeff_t_u_qsnow = 0.
         regcoeff_t_u_qgraup = 0.
         regcoeff_ps_u_rh = 0.
         regcoeff_ps_u_qcloud = 0.
         regcoeff_ps_u_qice = 0.
         regcoeff_ps_u_qrain = 0.
         regcoeff_ps_u_qsnow = 0.
         regcoeff_ps_u_qgraup = 0.
         regcoeff_rh_u_qcloud = 0.
         regcoeff_rh_u_qice = 0.
         regcoeff_rh_u_qrain = 0.
         regcoeff_rh_u_qsnow = 0.
         regcoeff_rh_u_qgraup = 0.
         regcoeff_qcloud_u_qice = 0.
         regcoeff_qcloud_u_qrain = 0.
         regcoeff_qcloud_u_qsnow = 0.
         regcoeff_qcloud_u_qgraup = 0.
         regcoeff_qice_u_qrain = 0.
         regcoeff_qice_u_qsnow = 0.
         regcoeff_qice_u_qgraup = 0.
         regcoeff_qrain_u_qsnow = 0.
         regcoeff_qrain_u_qgraup = 0.
         regcoeff_qsnow_u_qgraup = 0.
         if (qnr_option >= 1) then   !zhengh 20221014 add qnrain (reg to others) 10
            regcoeff_psi_qnrain = 0.
            regcoeff_chi_u_qnrain = 0.
            regcoeff_t_u_qnrain = 0.
            regcoeff_ps_u_qnrain = 0.
            regcoeff_rh_u_qnrain = 0.
            regcoeff_qcloud_u_qnrain = 0.
            regcoeff_qice_u_qnrain = 0.
            regcoeff_qrain_u_qnrain = 0.
            regcoeff_qsnow_u_qnrain = 0.
            regcoeff_qgraup_u_qnrain = 0.
         end if
      end if
   else	if ( cv_options == 7 .and. use_regcoeff ) then						! cv_options == 7
      regcoeff_u_v = 0.
      regcoeff_u_t = 0.
      regcoeff_u_ps = 0.
      regcoeff_u_rh = 0.
      regcoeff_u_qcloud = 0.
      regcoeff_u_qice = 0.
      regcoeff_u_qrain = 0.
      regcoeff_u_qsnow = 0.
      regcoeff_u_qgraup = 0.
      regcoeff_v_u_t = 0.
      regcoeff_v_u_ps = 0.
      regcoeff_v_u_rh = 0.
      regcoeff_v_u_qcloud = 0.
      regcoeff_v_u_qice = 0.
      regcoeff_v_u_qrain = 0.
      regcoeff_v_u_qsnow = 0.
      regcoeff_v_u_qgraup = 0.
      regcoeff_t_u_ps = 0.
      regcoeff_t_u_rh = 0.
      regcoeff_t_u_qcloud = 0.
      regcoeff_t_u_qice = 0.
      regcoeff_t_u_qrain = 0.
      regcoeff_t_u_qsnow = 0.
      regcoeff_t_u_qgraup = 0.
      regcoeff_ps_u_rh = 0.
      regcoeff_ps_u_qcloud = 0.
      regcoeff_ps_u_qice = 0.
      regcoeff_ps_u_qrain = 0.
      regcoeff_ps_u_qsnow = 0.
      regcoeff_ps_u_qgraup = 0.
      regcoeff_rh_u_qcloud = 0.
      regcoeff_rh_u_qice = 0.
      regcoeff_rh_u_qrain = 0.
      regcoeff_rh_u_qsnow = 0.
      regcoeff_rh_u_qgraup = 0.
      regcoeff_qcloud_u_qice = 0.
      regcoeff_qcloud_u_qrain = 0.
      regcoeff_qcloud_u_qsnow = 0.
      regcoeff_qcloud_u_qgraup = 0.
      regcoeff_qice_u_qrain = 0.
      regcoeff_qice_u_qsnow = 0.
      regcoeff_qice_u_qgraup = 0.
      regcoeff_qrain_u_qsnow = 0.
      regcoeff_qrain_u_qgraup = 0.
      regcoeff_qsnow_u_qgraup = 0.
      if ( qnr_option >= 1) then   !zhengh 20221014 add qnrain (reg to others) 10
         regcoeff_u_qnrain = 0.
         regcoeff_v_u_qnrain = 0.
         regcoeff_t_u_qnrain = 0.
         regcoeff_ps_u_qnrain = 0.
         regcoeff_rh_u_qnrain = 0.
         regcoeff_qcloud_u_qnrain = 0.
         regcoeff_qice_u_qnrain = 0.
         regcoeff_qrain_u_qnrain = 0.
         regcoeff_qsnow_u_qnrain = 0.
         regcoeff_qgraup_u_qnrain = 0.
      end if
   end if

   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then			!Zheng,20200522
      if ( cv_options == 5 ) then
         read (be_unit) regcoeff_psi_chi
         read (be_unit) regcoeff_psi_ps 
         read (be_unit) regcoeff_psi_t
!      else
      else if ( cv_options == 6 ) then							!Zheng,20200513
         do i = 1 , 9
            read (be_unit) var80c  
            select case( trim(adjustl(var80c)) )
            case ('regcoeff_psi_chi')
               read (be_unit) regcoeff_psi_chi
            case ('regcoeff_psi_t')
               read (be_unit) regcoeff_psi_t
            case ('regcoeff_psi_ps')
               read (be_unit) regcoeff_psi_ps
            case ('regcoeff_psi_rh')
               read (be_unit) regcoeff_psi_rh
            case ('regcoeff_chi_u_t')
               read (be_unit) regcoeff_chi_u_t
            case ('regcoeff_chi_u_ps')
               read (be_unit) regcoeff_chi_u_ps
            case ('regcoeff_chi_u_rh')
               read (be_unit) regcoeff_chi_u_rh
            case ('regcoeff_t_u_rh')
               read (be_unit) regcoeff_t_u_rh
            case ('regcoeff_ps_u_rh')
               read (be_unit) regcoeff_ps_u_rh
            case default;
               message(1)=' Read problem in regression coefficients in BE file '
               write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
               call da_error(__FILE__,__LINE__,message(1:2))
            end select
         end do
      else if ( cv_options == 9 .and. (use_regcoeff .and. (qnr_option == 0))) then				!Zheng,20200513
         do while ( 1 == 1 )
            read (be_unit) var80c  
            select case( trim(adjustl(var80c)) )
            case ('regcoeff_psi_chi')
	       write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_psi_chi
               if ( remain_regcoeff2(1) .eq. 0 ) then
                  regcoeff_psi_chi = 0.
               end if
            case ('regcoeff_psi_t')
	       write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_psi_t
               if ( remain_regcoeff3(1) .eq. 0 ) then
                  regcoeff_psi_t = 0.
               end if
            case ('regcoeff_psi_ps')
	       write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_psi_ps
               if ( remain_regcoeff4(1) .eq. 0 ) then
                  regcoeff_psi_ps = 0.
               end if
            case ('regcoeff_psi_rh')
	       write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_psi_rh
               if ( remain_regcoeff5(1) .eq. 0 ) then
                  regcoeff_psi_rh = 0.
               end if
            case ('regcoeff_psi_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_psi_qcloud
                     if ( remain_regcoeff6(1) .eq. 0 ) then
                        regcoeff_psi_qcloud = 0.
                     end if
            case ('regcoeff_psi_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_psi_qice
                     if ( remain_regcoeff7(1) .eq. 0 ) then
                        regcoeff_psi_qice = 0.
                     end if
            case ('regcoeff_psi_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_psi_qrain
                     if ( remain_regcoeff8(1) .eq. 0 ) then
                        regcoeff_psi_qrain = 0.
                     end if
            case ('regcoeff_psi_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_psi_qsnow
                     if ( remain_regcoeff9(1) .eq. 0 ) then
                        regcoeff_psi_qsnow = 0.
                     end if
            case ('regcoeff_psi_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_psi_qgraup
                     if ( remain_regcoeff10(1) .eq. 0 ) then
                        regcoeff_psi_qgraup = 0.
                     end if
                  case ('regcoeff_chi_u_t','regcoeff_chi_t')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_t
                     if ( remain_regcoeff3(2) .eq. 0 ) then
                        regcoeff_chi_u_t = 0.
                     end if
                  case ('regcoeff_chi_u_ps','regcoeff_chi_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_ps
                     if ( remain_regcoeff4(2) .eq. 0 ) then
                        regcoeff_chi_u_ps = 0.
                     end if
                  case ('regcoeff_chi_u_rh','regcoeff_chi_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_rh
                     if ( remain_regcoeff5(2) .eq. 0 ) then
                        regcoeff_chi_u_rh = 0.
                     end if
            case ('regcoeff_chi_u_qcloud','regcoeff_chi_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_qcloud
                     if ( remain_regcoeff6(2) .eq. 0 ) then
                        regcoeff_chi_u_qcloud = 0.
                     end if
            case ('regcoeff_chi_u_qice','regcoeff_chi_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_qice
                     if ( remain_regcoeff7(2) .eq. 0 ) then
                        regcoeff_chi_u_qice = 0.
                     end if
            case ('regcoeff_chi_u_qrain','regcoeff_chi_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_qrain
                     if ( remain_regcoeff8(2) .eq. 0 ) then
                        regcoeff_chi_u_qrain = 0.
                     end if
            case ('regcoeff_chi_u_qsnow','regcoeff_chi_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_qsnow
                     if ( remain_regcoeff9(2) .eq. 0 ) then
                        regcoeff_chi_u_qsnow = 0.
                     end if
            case ('regcoeff_chi_u_qgraup','regcoeff_chi_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_chi_u_qgraup
                     if ( remain_regcoeff10(2) .eq. 0 ) then
                        regcoeff_chi_u_qgraup = 0.
                     end if
                  case ('regcoeff_t_u_ps','regcoeff_t_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_ps
                     if ( remain_regcoeff4(3) .eq. 0 ) then
                        regcoeff_t_u_ps = 0.
                     end if
                  case ('regcoeff_t_u_rh','regcoeff_t_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_rh
                     if ( remain_regcoeff5(3) .eq. 0 ) then
                        regcoeff_t_u_rh = 0.
                     end if
            case ('regcoeff_t_u_qcloud','regcoeff_t_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qcloud
                     if ( remain_regcoeff6(3) .eq. 0 ) then
                        regcoeff_t_u_qcloud = 0.
                     end if
            case ('regcoeff_t_u_qice','regcoeff_t_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qice
                     if ( remain_regcoeff7(3) .eq. 0 ) then
                        regcoeff_t_u_qice = 0.
                     end if
            case ('regcoeff_t_u_qrain','regcoeff_t_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qrain
                     if ( remain_regcoeff8(3) .eq. 0 ) then
                        regcoeff_t_u_qrain = 0.
                     end if
            case ('regcoeff_t_u_qsnow','regcoeff_t_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qsnow
                     if ( remain_regcoeff9(3) .eq. 0 ) then
                        regcoeff_t_u_qsnow = 0.
                     end if
            case ('regcoeff_t_u_qgraup','regcoeff_t_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qgraup
                     if ( remain_regcoeff10(3) .eq. 0 ) then
                        regcoeff_t_u_qgraup = 0.
                     end if
                  case ('regcoeff_ps_u_rh','regcoeff_ps_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_rh
                     if ( remain_regcoeff5(4) .eq. 0 ) then
                        regcoeff_ps_u_rh = 0.
                     end if
            case ('regcoeff_ps_u_qcloud','regcoeff_ps_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qcloud
                     if ( remain_regcoeff6(4) .eq. 0 ) then
                        regcoeff_ps_u_qcloud = 0.
                     end if
            case ('regcoeff_ps_u_qice','regcoeff_ps_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qice
                     if ( remain_regcoeff7(4) .eq. 0 ) then
                        regcoeff_ps_u_qice = 0.
                     end if
            case ('regcoeff_ps_u_qrain','regcoeff_ps_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qrain
                     if ( remain_regcoeff8(4) .eq. 0 ) then
                        regcoeff_ps_u_qrain = 0.
                     end if
            case ('regcoeff_ps_u_qsnow','regcoeff_ps_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qsnow
                     if ( remain_regcoeff9(4) .eq. 0 ) then
                        regcoeff_ps_u_qsnow = 0.
                     end if
            case ('regcoeff_ps_u_qgraup','regcoeff_ps_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qgraup
                     if ( remain_regcoeff10(4) .eq. 0 ) then
                        regcoeff_ps_u_qgraup = 0.
                     end if
            case ('regcoeff_rh_u_qcloud','regcoeff_rh_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qcloud
                     if ( remain_regcoeff6(5) .eq. 0 ) then
                        regcoeff_rh_u_qcloud = 0.
                     end if
            case ('regcoeff_rh_u_qice','regcoeff_rh_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qice
                     if ( remain_regcoeff7(5) .eq. 0 ) then
                        regcoeff_rh_u_qice = 0.
                     end if
            case ('regcoeff_rh_u_qrain','regcoeff_rh_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qrain
                     if ( remain_regcoeff8(5) .eq. 0 ) then
                        regcoeff_rh_u_qrain = 0.
                     end if
            case ('regcoeff_rh_u_qsnow','regcoeff_rh_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qsnow
                     if ( remain_regcoeff9(5) .eq. 0 ) then
                        regcoeff_rh_u_qsnow = 0.
                     end if
            case ('regcoeff_rh_u_qgraup','regcoeff_rh_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qgraup
                     if ( remain_regcoeff10(5) .eq. 0 ) then
                        regcoeff_rh_u_qgraup = 0.
                     end if
            case ('regcoeff_qcloud_u_qice','regcoeff_qcloud_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qice
                     if ( remain_regcoeff7(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qice = 0.
                     end if
            case ('regcoeff_qcloud_u_qrain','regcoeff_qcloud_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qrain
                     if ( remain_regcoeff8(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qrain = 0.
                     end if
            case ('regcoeff_qcloud_u_qsnow','regcoeff_qcloud_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qsnow
                     if ( remain_regcoeff9(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qsnow = 0.
                     end if
            case ('regcoeff_qcloud_u_qgraup','regcoeff_qcloud_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qgraup
                     if ( remain_regcoeff10(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qgraup = 0.
                     end if
            case ('regcoeff_qice_u_qrain','regcoeff_qice_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qrain
                     if ( remain_regcoeff8(7) .eq. 0 ) then
                        regcoeff_qice_u_qrain = 0.
                     end if
            case ('regcoeff_qice_u_qsnow','regcoeff_qice_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qsnow
                     if ( remain_regcoeff9(7) .eq. 0 ) then
                        regcoeff_qice_u_qsnow = 0.
                     end if
            case ('regcoeff_qice_u_qgraup','regcoeff_qice_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qgraup
                     if ( remain_regcoeff10(7) .eq. 0 ) then
                        regcoeff_qice_u_qgraup = 0.
                     end if
            case ('regcoeff_qrain_u_qsnow','regcoeff_qrain_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qsnow
                     if ( remain_regcoeff9(8) .eq. 0 ) then
                        regcoeff_qrain_u_qsnow = 0.
                     end if
            case ('regcoeff_qrain_u_qgraup','regcoeff_qrain_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qgraup
                     if ( remain_regcoeff10(8) .eq. 0 ) then
                        regcoeff_qrain_u_qgraup = 0.
                     end if
            case ('regcoeff_qsnow_u_qgraup','regcoeff_qsnow_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qsnow_u_qgraup
                     if ( remain_regcoeff10(9) .eq. 0 ) then
                        regcoeff_qsnow_u_qgraup = 0.
                     end if
            case ('regcoeff_completed')
               write(*,*) "---Reading regcoeff completed---"
	            exit
            case default;
               message(1)=' Read problem in regression coefficients in BE file '
               write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
               call da_error(__FILE__,__LINE__,message(1:2))
            end select
         end do
      else if ( cv_options == 9 .and. (use_regcoeff .and. (qnr_option >= 1)))  then  !zhengh 20221014 add qnrain (reg to others) 10
         ! regcoeff_psi_qnrain = 0.
         ! regcoeff_chi_u_qnrain = 0.
         ! regcoeff_t_u_qnrain = 0.
         ! regcoeff_ps_u_qnrain = 0.
         ! regcoeff_rh_u_qnrain = 0.
         ! regcoeff_qcloud_u_qnrain = 0.
         ! regcoeff_qice_u_qnrain = 0.
         ! regcoeff_qrain_u_qnrain = 0.
         ! regcoeff_qsnow_u_qnrain = 0.
         ! regcoeff_qgraup_u_qnrain = 0.
         do while ( 1 == 1 )
            read (be_unit) var80c  
            select case( trim(adjustl(var80c)) )
                     case ('regcoeff_psi_chi')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_psi_chi
                  if ( remain_regcoeff2(1) .eq. 0 ) then
                     regcoeff_psi_chi = 0.
                  end if
               case ('regcoeff_psi_t')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_psi_t
                  if ( remain_regcoeff3(1) .eq. 0 ) then
                     regcoeff_psi_t = 0.
                  end if
               case ('regcoeff_psi_ps')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_psi_ps
                  if ( remain_regcoeff4(1) .eq. 0 ) then
                     regcoeff_psi_ps = 0.
                  end if
               case ('regcoeff_psi_rh')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_psi_rh
                  if ( remain_regcoeff5(1) .eq. 0 ) then
                     regcoeff_psi_rh = 0.
                  end if
               case ('regcoeff_psi_qcloud')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_psi_qcloud
                        if ( remain_regcoeff6(1) .eq. 0 ) then
                           regcoeff_psi_qcloud = 0.
                        end if
               case ('regcoeff_psi_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_psi_qice
                        if ( remain_regcoeff7(1) .eq. 0 ) then
                           regcoeff_psi_qice = 0.
                        end if
               case ('regcoeff_psi_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_psi_qrain
                        if ( remain_regcoeff8(1) .eq. 0 ) then
                           regcoeff_psi_qrain = 0.
                        end if
               case ('regcoeff_psi_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_psi_qsnow
                        if ( remain_regcoeff9(1) .eq. 0 ) then
                           regcoeff_psi_qsnow = 0.
                        end if
               case ('regcoeff_psi_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_psi_qgraup
                        if ( remain_regcoeff10(1) .eq. 0 ) then
                           regcoeff_psi_qgraup = 0.
                        end if
                     case ('regcoeff_chi_u_t','regcoeff_chi_t')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_t
                        if ( remain_regcoeff3(2) .eq. 0 ) then
                           regcoeff_chi_u_t = 0.
                        end if
                     case ('regcoeff_chi_u_ps','regcoeff_chi_ps')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_ps
                        if ( remain_regcoeff4(2) .eq. 0 ) then
                           regcoeff_chi_u_ps = 0.
                        end if
                     case ('regcoeff_chi_u_rh','regcoeff_chi_rh')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_rh
                        if ( remain_regcoeff5(2) .eq. 0 ) then
                           regcoeff_chi_u_rh = 0.
                        end if
               case ('regcoeff_chi_u_qcloud','regcoeff_chi_qcloud')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_qcloud
                        if ( remain_regcoeff6(2) .eq. 0 ) then
                           regcoeff_chi_u_qcloud = 0.
                        end if
               case ('regcoeff_chi_u_qice','regcoeff_chi_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_qice
                        if ( remain_regcoeff7(2) .eq. 0 ) then
                           regcoeff_chi_u_qice = 0.
                        end if
               case ('regcoeff_chi_u_qrain','regcoeff_chi_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_qrain
                        if ( remain_regcoeff8(2) .eq. 0 ) then
                           regcoeff_chi_u_qrain = 0.
                        end if
               case ('regcoeff_chi_u_qsnow','regcoeff_chi_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_qsnow
                        if ( remain_regcoeff9(2) .eq. 0 ) then
                           regcoeff_chi_u_qsnow = 0.
                        end if
               case ('regcoeff_chi_u_qgraup','regcoeff_chi_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_chi_u_qgraup
                        if ( remain_regcoeff10(2) .eq. 0 ) then
                           regcoeff_chi_u_qgraup = 0.
                        end if
                     case ('regcoeff_t_u_ps','regcoeff_t_ps')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_ps
                        if ( remain_regcoeff4(3) .eq. 0 ) then
                           regcoeff_t_u_ps = 0.
                        end if
                     case ('regcoeff_t_u_rh','regcoeff_t_rh')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_rh
                        if ( remain_regcoeff5(3) .eq. 0 ) then
                           regcoeff_t_u_rh = 0.
                        end if
               case ('regcoeff_t_u_qcloud','regcoeff_t_qcloud')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_qcloud
                        if ( remain_regcoeff6(3) .eq. 0 ) then
                           regcoeff_t_u_qcloud = 0.
                        end if
               case ('regcoeff_t_u_qice','regcoeff_t_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_qice
                        if ( remain_regcoeff7(3) .eq. 0 ) then
                           regcoeff_t_u_qice = 0.
                        end if
               case ('regcoeff_t_u_qrain','regcoeff_t_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_qrain
                        if ( remain_regcoeff8(3) .eq. 0 ) then
                           regcoeff_t_u_qrain = 0.
                        end if
               case ('regcoeff_t_u_qsnow','regcoeff_t_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_qsnow
                        if ( remain_regcoeff9(3) .eq. 0 ) then
                           regcoeff_t_u_qsnow = 0.
                        end if
               case ('regcoeff_t_u_qgraup','regcoeff_t_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_t_u_qgraup
                        if ( remain_regcoeff10(3) .eq. 0 ) then
                           regcoeff_t_u_qgraup = 0.
                        end if
                     case ('regcoeff_ps_u_rh','regcoeff_ps_rh')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_rh
                        if ( remain_regcoeff5(4) .eq. 0 ) then
                           regcoeff_ps_u_rh = 0.
                        end if
               case ('regcoeff_ps_u_qcloud','regcoeff_ps_qcloud')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_qcloud
                        if ( remain_regcoeff6(4) .eq. 0 ) then
                           regcoeff_ps_u_qcloud = 0.
                        end if
               case ('regcoeff_ps_u_qice','regcoeff_ps_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_qice
                        if ( remain_regcoeff7(4) .eq. 0 ) then
                           regcoeff_ps_u_qice = 0.
                        end if
               case ('regcoeff_ps_u_qrain','regcoeff_ps_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_qrain
                        if ( remain_regcoeff8(4) .eq. 0 ) then
                           regcoeff_ps_u_qrain = 0.
                        end if
               case ('regcoeff_ps_u_qsnow','regcoeff_ps_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_qsnow
                        if ( remain_regcoeff9(4) .eq. 0 ) then
                           regcoeff_ps_u_qsnow = 0.
                        end if
               case ('regcoeff_ps_u_qgraup','regcoeff_ps_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_ps_u_qgraup
                        if ( remain_regcoeff10(4) .eq. 0 ) then
                           regcoeff_ps_u_qgraup = 0.
                        end if
               case ('regcoeff_rh_u_qcloud','regcoeff_rh_qcloud')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_rh_u_qcloud
                        if ( remain_regcoeff6(5) .eq. 0 ) then
                           regcoeff_rh_u_qcloud = 0.
                        end if
               case ('regcoeff_rh_u_qice','regcoeff_rh_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_rh_u_qice
                        if ( remain_regcoeff7(5) .eq. 0 ) then
                           regcoeff_rh_u_qice = 0.
                        end if
               case ('regcoeff_rh_u_qrain','regcoeff_rh_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_rh_u_qrain
                        if ( remain_regcoeff8(5) .eq. 0 ) then
                           regcoeff_rh_u_qrain = 0.
                        end if
               case ('regcoeff_rh_u_qsnow','regcoeff_rh_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_rh_u_qsnow
                        if ( remain_regcoeff9(5) .eq. 0 ) then
                           regcoeff_rh_u_qsnow = 0.
                        end if
               case ('regcoeff_rh_u_qgraup','regcoeff_rh_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_rh_u_qgraup
                        if ( remain_regcoeff10(5) .eq. 0 ) then
                           regcoeff_rh_u_qgraup = 0.
                        end if
               case ('regcoeff_qcloud_u_qice','regcoeff_qcloud_qice')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qcloud_u_qice
                        if ( remain_regcoeff7(6) .eq. 0 ) then
                           regcoeff_qcloud_u_qice = 0.
                        end if
               case ('regcoeff_qcloud_u_qrain','regcoeff_qcloud_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qcloud_u_qrain
                        if ( remain_regcoeff8(6) .eq. 0 ) then
                           regcoeff_qcloud_u_qrain = 0.
                        end if
               case ('regcoeff_qcloud_u_qsnow','regcoeff_qcloud_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qcloud_u_qsnow
                        if ( remain_regcoeff9(6) .eq. 0 ) then
                           regcoeff_qcloud_u_qsnow = 0.
                        end if
               case ('regcoeff_qcloud_u_qgraup','regcoeff_qcloud_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qcloud_u_qgraup
                        if ( remain_regcoeff10(6) .eq. 0 ) then
                           regcoeff_qcloud_u_qgraup = 0.
                        end if
               case ('regcoeff_qice_u_qrain','regcoeff_qice_qrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qice_u_qrain
                        if ( remain_regcoeff8(7) .eq. 0 ) then
                           regcoeff_qice_u_qrain = 0.
                        end if
               case ('regcoeff_qice_u_qsnow','regcoeff_qice_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qice_u_qsnow
                        if ( remain_regcoeff9(7) .eq. 0 ) then
                           regcoeff_qice_u_qsnow = 0.
                        end if
               case ('regcoeff_qice_u_qgraup','regcoeff_qice_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qice_u_qgraup
                        if ( remain_regcoeff10(7) .eq. 0 ) then
                           regcoeff_qice_u_qgraup = 0.
                        end if
               case ('regcoeff_qrain_u_qsnow','regcoeff_qrain_qsnow')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qrain_u_qsnow
                        if ( remain_regcoeff9(8) .eq. 0 ) then
                           regcoeff_qrain_u_qsnow = 0.
                        end if
               case ('regcoeff_qrain_u_qgraup','regcoeff_qrain_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qrain_u_qgraup
                        if ( remain_regcoeff10(8) .eq. 0 ) then
                           regcoeff_qrain_u_qgraup = 0.
                        end if
               case ('regcoeff_qsnow_u_qgraup','regcoeff_qsnow_qgraup')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                        read (be_unit) regcoeff_qsnow_u_qgraup
                        if ( remain_regcoeff10(9) .eq. 0 ) then
                           regcoeff_qsnow_u_qgraup = 0.
                        end if


               case ('regcoeff_psi_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_psi_qnrain
                  if ( remain_regcoeff12(1) .eq. 0 ) then
                     regcoeff_psi_qnrain = 0.
                  end if
               case ('regcoeff_chi_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_chi_u_qnrain
                  if ( remain_regcoeff12(2) .eq. 0 ) then
                     regcoeff_chi_u_qnrain = 0.
                  end if
               case ('regcoeff_t_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_t_u_qnrain
                  if ( remain_regcoeff12(3) .eq. 0 ) then
                     regcoeff_t_u_qnrain = 0.
                  end if
               case ('regcoeff_ps_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_ps_u_qnrain
                  if ( remain_regcoeff12(4) .eq. 0 ) then
                     regcoeff_ps_u_qnrain = 0.
                  end if
               case ('regcoeff_rh_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_rh_u_qnrain
                  if ( remain_regcoeff12(5) .eq. 0 ) then
                     regcoeff_rh_u_qnrain = 0.
                  end if
               case ('regcoeff_qcloud_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_qcloud_u_qnrain
                  if ( remain_regcoeff12(6) .eq. 0 ) then
                     regcoeff_qcloud_u_qnrain = 0.
                  end if
               case ('regcoeff_qice_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_qice_u_qnrain
                  if ( remain_regcoeff12(7) .eq. 0 ) then
                     regcoeff_qice_u_qnrain = 0.
                  end if
               case ('regcoeff_qrain_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_qrain_u_qnrain
                  if ( remain_regcoeff12(8) .eq. 0 ) then
                     regcoeff_qrain_u_qnrain = 0.
                  end if
               case ('regcoeff_qsnow_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_qsnow_u_qnrain
                  if ( remain_regcoeff12(9) .eq. 0 ) then
                     regcoeff_qsnow_u_qnrain = 0.
                  end if
               case ('regcoeff_qgraup_u_qnrain')
                  write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                  read (be_unit) regcoeff_qgraup_u_qnrain
                  if ( remain_regcoeff12(10) .eq. 0 ) then
                     regcoeff_qgraup_u_qnrain = 0.
                  end if

                  case ('regcoeff_completed')
                  write(*,*) "---Reading regcoeff completed---"
                  exit
               case default;
                  message(1)=' Read problem in regression coefficients in BE file '
                  write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
                  call da_error(__FILE__,__LINE__,message(1:2))
            end select
         end do
      
      end if
   else if ( cv_options == 7 .and. (use_regcoeff .and. (qnr_option == 0))) then							!Zheng,20200513
      do while ( 1 == 1 )
         read (be_unit) var80c  
         select case( trim(adjustl(var80c)) )
            case ('regcoeff_u_v')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_v
               if ( remain_regcoeff2(1) .eq. 0 ) then
                  regcoeff_u_v = 0.
               end if
            case ('regcoeff_u_t')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_t
               if ( remain_regcoeff3(1) .eq. 0 ) then
                  regcoeff_u_t = 0.
               end if
            case ('regcoeff_u_ps')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_ps
               if ( remain_regcoeff4(1) .eq. 0 ) then
                  regcoeff_u_ps = 0.
               end if
            case ('regcoeff_u_rh')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_rh
               if ( remain_regcoeff5(1) .eq. 0 ) then
                  regcoeff_u_rh = 0.
               end if
            case ('regcoeff_u_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qcloud
                     if ( remain_regcoeff6(1) .eq. 0 ) then
                        regcoeff_u_qcloud = 0.
                     end if
            case ('regcoeff_u_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qice
                     if ( remain_regcoeff7(1) .eq. 0 ) then
                        regcoeff_u_qice = 0.
                     end if
            case ('regcoeff_u_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qrain
                     if ( remain_regcoeff8(1) .eq. 0 ) then
                        regcoeff_u_qrain = 0.
                     end if
            case ('regcoeff_u_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qsnow
                     if ( remain_regcoeff9(1) .eq. 0 ) then
                        regcoeff_u_qsnow = 0.
                     end if
            case ('regcoeff_u_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qgraup
                     if ( remain_regcoeff10(1) .eq. 0 ) then
                        regcoeff_u_qgraup = 0.
                     end if
                  case ('regcoeff_v_u_t','regcoeff_v_t')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_t
                     if ( remain_regcoeff3(2) .eq. 0 ) then
                        regcoeff_v_u_t = 0.
                     end if
                  case ('regcoeff_v_u_ps','regcoeff_v_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_ps
                     if ( remain_regcoeff4(2) .eq. 0 ) then
                        regcoeff_v_u_ps = 0.
                     end if
                  case ('regcoeff_v_u_rh','regcoeff_v_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_rh
                     if ( remain_regcoeff5(2) .eq. 0 ) then
                        regcoeff_v_u_rh = 0.
                     end if
            case ('regcoeff_v_u_qcloud','regcoeff_v_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qcloud
                     if ( remain_regcoeff6(2) .eq. 0 ) then
                        regcoeff_v_u_qcloud = 0.
                     end if
            case ('regcoeff_v_u_qice','regcoeff_v_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qice
                     if ( remain_regcoeff7(2) .eq. 0 ) then
                        regcoeff_v_u_qice = 0.
                     end if
            case ('regcoeff_v_u_qrain','regcoeff_v_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qrain
                     if ( remain_regcoeff8(2) .eq. 0 ) then
                        regcoeff_v_u_qrain = 0.
                     end if
            case ('regcoeff_v_u_qsnow','regcoeff_v_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qsnow
                     if ( remain_regcoeff9(2) .eq. 0 ) then
                        regcoeff_v_u_qsnow = 0.
                     end if
            case ('regcoeff_v_u_qgraup','regcoeff_v_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qgraup
                     if ( remain_regcoeff10(2) .eq. 0 ) then
                        regcoeff_v_u_qgraup = 0.
                     end if
                  case ('regcoeff_t_u_ps','regcoeff_t_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_ps
                     if ( remain_regcoeff4(3) .eq. 0 ) then
                        regcoeff_t_u_ps = 0.
                     end if
                  case ('regcoeff_t_u_rh','regcoeff_t_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_rh
                     if ( remain_regcoeff5(3) .eq. 0 ) then
                        regcoeff_t_u_rh = 0.
                     end if
            case ('regcoeff_t_u_qcloud','regcoeff_t_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qcloud
                     if ( remain_regcoeff6(3) .eq. 0 ) then
                        regcoeff_t_u_qcloud = 0.
                     end if
            case ('regcoeff_t_u_qice','regcoeff_t_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qice
                     if ( remain_regcoeff7(3) .eq. 0 ) then
                        regcoeff_t_u_qice = 0.
                     end if
            case ('regcoeff_t_u_qrain','regcoeff_t_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qrain
                     if ( remain_regcoeff8(3) .eq. 0 ) then
                        regcoeff_t_u_qrain = 0.
                     end if
            case ('regcoeff_t_u_qsnow','regcoeff_t_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qsnow
                     if ( remain_regcoeff9(3) .eq. 0 ) then
                        regcoeff_t_u_qsnow = 0.
                     end if
            case ('regcoeff_t_u_qgraup','regcoeff_t_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qgraup
                     if ( remain_regcoeff10(3) .eq. 0 ) then
                        regcoeff_t_u_qgraup = 0.
                     end if
                  case ('regcoeff_ps_u_rh','regcoeff_ps_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_rh
                     if ( remain_regcoeff5(4) .eq. 0 ) then
                        regcoeff_ps_u_rh = 0.
                     end if
            case ('regcoeff_ps_u_qcloud','regcoeff_ps_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qcloud
                     if ( remain_regcoeff6(4) .eq. 0 ) then
                        regcoeff_ps_u_qcloud = 0.
                     end if
            case ('regcoeff_ps_u_qice','regcoeff_ps_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qice
                     if ( remain_regcoeff7(4) .eq. 0 ) then
                        regcoeff_ps_u_qice = 0.
                     end if
            case ('regcoeff_ps_u_qrain','regcoeff_ps_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qrain
                     if ( remain_regcoeff8(4) .eq. 0 ) then
                        regcoeff_ps_u_qrain = 0.
                     end if
            case ('regcoeff_ps_u_qsnow','regcoeff_ps_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qsnow
                     if ( remain_regcoeff9(4) .eq. 0 ) then
                        regcoeff_ps_u_qsnow = 0.
                     end if
            case ('regcoeff_ps_u_qgraup','regcoeff_ps_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qgraup
                     if ( remain_regcoeff10(4) .eq. 0 ) then
                        regcoeff_ps_u_qgraup = 0.
                     end if
            case ('regcoeff_rh_u_qcloud','regcoeff_rh_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qcloud
                     if ( remain_regcoeff6(5) .eq. 0 ) then
                        regcoeff_rh_u_qcloud = 0.
                     end if
            case ('regcoeff_rh_u_qice','regcoeff_rh_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qice
                     if ( remain_regcoeff7(5) .eq. 0 ) then
                        regcoeff_rh_u_qice = 0.
                     end if
            case ('regcoeff_rh_u_qrain','regcoeff_rh_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qrain
                     if ( remain_regcoeff8(5) .eq. 0 ) then
                        regcoeff_rh_u_qrain = 0.
                     end if
            case ('regcoeff_rh_u_qsnow','regcoeff_rh_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qsnow
                     if ( remain_regcoeff9(5) .eq. 0 ) then
                        regcoeff_rh_u_qsnow = 0.
                     end if
            case ('regcoeff_rh_u_qgraup','regcoeff_rh_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qgraup
                     if ( remain_regcoeff10(5) .eq. 0 ) then
                        regcoeff_rh_u_qgraup = 0.
                     end if
            case ('regcoeff_qcloud_u_qice','regcoeff_qcloud_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qice
                     if ( remain_regcoeff7(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qice = 0.
                     end if
            case ('regcoeff_qcloud_u_qrain','regcoeff_qcloud_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qrain
                     if ( remain_regcoeff8(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qrain = 0.
                     end if
            case ('regcoeff_qcloud_u_qsnow','regcoeff_qcloud_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qsnow
                     if ( remain_regcoeff9(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qsnow = 0.
                     end if
            case ('regcoeff_qcloud_u_qgraup','regcoeff_qcloud_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qgraup
                     if ( remain_regcoeff10(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qgraup = 0.
                     end if
            case ('regcoeff_qice_u_qrain','regcoeff_qice_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qrain
                     if ( remain_regcoeff8(7) .eq. 0 ) then
                        regcoeff_qice_u_qrain = 0.
                     end if
            case ('regcoeff_qice_u_qsnow','regcoeff_qice_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qsnow
                     if ( remain_regcoeff9(7) .eq. 0 ) then
                        regcoeff_qice_u_qsnow = 0.
                     end if
            case ('regcoeff_qice_u_qgraup','regcoeff_qice_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qgraup
                     if ( remain_regcoeff10(7) .eq. 0 ) then
                        regcoeff_qice_u_qgraup = 0.
                     end if
            case ('regcoeff_qrain_u_qsnow','regcoeff_qrain_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qsnow
                     if ( remain_regcoeff9(8) .eq. 0 ) then
                        regcoeff_qrain_u_qsnow = 0.
                     end if
            case ('regcoeff_qrain_u_qgraup','regcoeff_qrain_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qgraup
                     if ( remain_regcoeff10(8) .eq. 0 ) then
                        regcoeff_qrain_u_qgraup = 0.
                     end if
            case ('regcoeff_qsnow_u_qgraup','regcoeff_qsnow_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qsnow_u_qgraup
                     if ( remain_regcoeff10(9) .eq. 0 ) then
                        regcoeff_qsnow_u_qgraup = 0.
                     end if
            case ('regcoeff_completed')
               write(*,*) "---Reading regcoeff completed---"
         exit
            case default;
               message(1)=' Read problem in regression coefficients in BE file '
               write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
               call da_error(__FILE__,__LINE__,message(1:2))
            end select
      end do

   else if ( cv_options == 7 .and. (use_regcoeff .and. (qnr_option >= 1))) then  !zhengh 20221014 add qnrain (reg to others) 10
            ! regcoeff_u_qnrain = 0.
            ! regcoeff_v_u_qnrain = 0.
            ! regcoeff_t_u_qnrain = 0.
            ! regcoeff_ps_u_qnrain = 0.
            ! regcoeff_rh_u_qnrain = 0.
            ! regcoeff_qcloud_u_qnrain = 0.
            ! regcoeff_qice_u_qnrain = 0.
            ! regcoeff_qrain_u_qnrain = 0.
            ! regcoeff_qsnow_u_qnrain = 0.
            ! regcoeff_qgraup_u_qnrain = 0.
      do while ( 1 == 1 )
         read (be_unit) var80c  
         select case( trim(adjustl(var80c)) )
            case ('regcoeff_u_v')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_v
               if ( remain_regcoeff2(1) .eq. 0 ) then
                  regcoeff_u_v = 0.
               end if
            case ('regcoeff_u_t')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_t
               if ( remain_regcoeff3(1) .eq. 0 ) then
                  regcoeff_u_t = 0.
               end if
            case ('regcoeff_u_ps')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_ps
               if ( remain_regcoeff4(1) .eq. 0 ) then
                  regcoeff_u_ps = 0.
               end if
            case ('regcoeff_u_rh')
         write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
               read (be_unit) regcoeff_u_rh
               if ( remain_regcoeff5(1) .eq. 0 ) then
                  regcoeff_u_rh = 0.
               end if
            case ('regcoeff_u_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qcloud
                     if ( remain_regcoeff6(1) .eq. 0 ) then
                        regcoeff_u_qcloud = 0.
                     end if
            case ('regcoeff_u_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qice
                     if ( remain_regcoeff7(1) .eq. 0 ) then
                        regcoeff_u_qice = 0.
                     end if
            case ('regcoeff_u_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qrain
                     if ( remain_regcoeff8(1) .eq. 0 ) then
                        regcoeff_u_qrain = 0.
                     end if
            case ('regcoeff_u_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qsnow
                     if ( remain_regcoeff9(1) .eq. 0 ) then
                        regcoeff_u_qsnow = 0.
                     end if
            case ('regcoeff_u_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_u_qgraup
                     if ( remain_regcoeff10(1) .eq. 0 ) then
                        regcoeff_u_qgraup = 0.
                     end if
                  case ('regcoeff_v_u_t','regcoeff_v_t')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_t
                     if ( remain_regcoeff3(2) .eq. 0 ) then
                        regcoeff_v_u_t = 0.
                     end if
                  case ('regcoeff_v_u_ps','regcoeff_v_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_ps
                     if ( remain_regcoeff4(2) .eq. 0 ) then
                        regcoeff_v_u_ps = 0.
                     end if
                  case ('regcoeff_v_u_rh','regcoeff_v_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_rh
                     if ( remain_regcoeff5(2) .eq. 0 ) then
                        regcoeff_v_u_rh = 0.
                     end if
            case ('regcoeff_v_u_qcloud','regcoeff_v_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qcloud
                     if ( remain_regcoeff6(2) .eq. 0 ) then
                        regcoeff_v_u_qcloud = 0.
                     end if
            case ('regcoeff_v_u_qice','regcoeff_v_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qice
                     if ( remain_regcoeff7(2) .eq. 0 ) then
                        regcoeff_v_u_qice = 0.
                     end if
            case ('regcoeff_v_u_qrain','regcoeff_v_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qrain
                     if ( remain_regcoeff8(2) .eq. 0 ) then
                        regcoeff_v_u_qrain = 0.
                     end if
            case ('regcoeff_v_u_qsnow','regcoeff_v_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qsnow
                     if ( remain_regcoeff9(2) .eq. 0 ) then
                        regcoeff_v_u_qsnow = 0.
                     end if
            case ('regcoeff_v_u_qgraup','regcoeff_v_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_v_u_qgraup
                     if ( remain_regcoeff10(2) .eq. 0 ) then
                        regcoeff_v_u_qgraup = 0.
                     end if
                  case ('regcoeff_t_u_ps','regcoeff_t_ps')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_ps
                     if ( remain_regcoeff4(3) .eq. 0 ) then
                        regcoeff_t_u_ps = 0.
                     end if
                  case ('regcoeff_t_u_rh','regcoeff_t_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_rh
                     if ( remain_regcoeff5(3) .eq. 0 ) then
                        regcoeff_t_u_rh = 0.
                     end if
            case ('regcoeff_t_u_qcloud','regcoeff_t_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qcloud
                     if ( remain_regcoeff6(3) .eq. 0 ) then
                        regcoeff_t_u_qcloud = 0.
                     end if
            case ('regcoeff_t_u_qice','regcoeff_t_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qice
                     if ( remain_regcoeff7(3) .eq. 0 ) then
                        regcoeff_t_u_qice = 0.
                     end if
            case ('regcoeff_t_u_qrain','regcoeff_t_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qrain
                     if ( remain_regcoeff8(3) .eq. 0 ) then
                        regcoeff_t_u_qrain = 0.
                     end if
            case ('regcoeff_t_u_qsnow','regcoeff_t_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qsnow
                     if ( remain_regcoeff9(3) .eq. 0 ) then
                        regcoeff_t_u_qsnow = 0.
                     end if
            case ('regcoeff_t_u_qgraup','regcoeff_t_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_t_u_qgraup
                     if ( remain_regcoeff10(3) .eq. 0 ) then
                        regcoeff_t_u_qgraup = 0.
                     end if
                  case ('regcoeff_ps_u_rh','regcoeff_ps_rh')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_rh
                     if ( remain_regcoeff5(4) .eq. 0 ) then
                        regcoeff_ps_u_rh = 0.
                     end if
            case ('regcoeff_ps_u_qcloud','regcoeff_ps_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qcloud
                     if ( remain_regcoeff6(4) .eq. 0 ) then
                        regcoeff_ps_u_qcloud = 0.
                     end if
            case ('regcoeff_ps_u_qice','regcoeff_ps_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qice
                     if ( remain_regcoeff7(4) .eq. 0 ) then
                        regcoeff_ps_u_qice = 0.
                     end if
            case ('regcoeff_ps_u_qrain','regcoeff_ps_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qrain
                     if ( remain_regcoeff8(4) .eq. 0 ) then
                        regcoeff_ps_u_qrain = 0.
                     end if
            case ('regcoeff_ps_u_qsnow','regcoeff_ps_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qsnow
                     if ( remain_regcoeff9(4) .eq. 0 ) then
                        regcoeff_ps_u_qsnow = 0.
                     end if
            case ('regcoeff_ps_u_qgraup','regcoeff_ps_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_ps_u_qgraup
                     if ( remain_regcoeff10(4) .eq. 0 ) then
                        regcoeff_ps_u_qgraup = 0.
                     end if
            case ('regcoeff_rh_u_qcloud','regcoeff_rh_qcloud')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qcloud
                     if ( remain_regcoeff6(5) .eq. 0 ) then
                        regcoeff_rh_u_qcloud = 0.
                     end if
            case ('regcoeff_rh_u_qice','regcoeff_rh_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qice
                     if ( remain_regcoeff7(5) .eq. 0 ) then
                        regcoeff_rh_u_qice = 0.
                     end if
            case ('regcoeff_rh_u_qrain','regcoeff_rh_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qrain
                     if ( remain_regcoeff8(5) .eq. 0 ) then
                        regcoeff_rh_u_qrain = 0.
                     end if
            case ('regcoeff_rh_u_qsnow','regcoeff_rh_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qsnow
                     if ( remain_regcoeff9(5) .eq. 0 ) then
                        regcoeff_rh_u_qsnow = 0.
                     end if
            case ('regcoeff_rh_u_qgraup','regcoeff_rh_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_rh_u_qgraup
                     if ( remain_regcoeff10(5) .eq. 0 ) then
                        regcoeff_rh_u_qgraup = 0.
                     end if
            case ('regcoeff_qcloud_u_qice','regcoeff_qcloud_qice')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qice
                     if ( remain_regcoeff7(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qice = 0.
                     end if
            case ('regcoeff_qcloud_u_qrain','regcoeff_qcloud_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qrain
                     if ( remain_regcoeff8(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qrain = 0.
                     end if
            case ('regcoeff_qcloud_u_qsnow','regcoeff_qcloud_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qsnow
                     if ( remain_regcoeff9(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qsnow = 0.
                     end if
            case ('regcoeff_qcloud_u_qgraup','regcoeff_qcloud_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qcloud_u_qgraup
                     if ( remain_regcoeff10(6) .eq. 0 ) then
                        regcoeff_qcloud_u_qgraup = 0.
                     end if
            case ('regcoeff_qice_u_qrain','regcoeff_qice_qrain')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qrain
                     if ( remain_regcoeff8(7) .eq. 0 ) then
                        regcoeff_qice_u_qrain = 0.
                     end if
            case ('regcoeff_qice_u_qsnow','regcoeff_qice_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qsnow
                     if ( remain_regcoeff9(7) .eq. 0 ) then
                        regcoeff_qice_u_qsnow = 0.
                     end if
            case ('regcoeff_qice_u_qgraup','regcoeff_qice_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qice_u_qgraup
                     if ( remain_regcoeff10(7) .eq. 0 ) then
                        regcoeff_qice_u_qgraup = 0.
                     end if
            case ('regcoeff_qrain_u_qsnow','regcoeff_qrain_qsnow')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qsnow
                     if ( remain_regcoeff9(8) .eq. 0 ) then
                        regcoeff_qrain_u_qsnow = 0.
                     end if
            case ('regcoeff_qrain_u_qgraup','regcoeff_qrain_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qrain_u_qgraup
                     if ( remain_regcoeff10(8) .eq. 0 ) then
                        regcoeff_qrain_u_qgraup = 0.
                     end if
            case ('regcoeff_qsnow_u_qgraup','regcoeff_qsnow_qgraup')
               write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
                     read (be_unit) regcoeff_qsnow_u_qgraup
                     if ( remain_regcoeff10(9) .eq. 0 ) then
                        regcoeff_qsnow_u_qgraup = 0.
                     end if



         case ('regcoeff_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_u_qnrain
            if ( remain_regcoeff12(1) .eq. 0 ) then
               regcoeff_u_qnrain = 0.
            end if
         case ('regcoeff_v_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_v_u_qnrain
            if ( remain_regcoeff12(2) .eq. 0 ) then
               regcoeff_v_u_qnrain = 0.
            end if
         case ('regcoeff_t_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_t_u_qnrain
            if ( remain_regcoeff12(3) .eq. 0 ) then
               regcoeff_t_u_qnrain = 0.
            end if
         case ('regcoeff_ps_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_ps_u_qnrain
            if ( remain_regcoeff12(4) .eq. 0 ) then
               regcoeff_ps_u_qnrain = 0.
            end if
         case ('regcoeff_rh_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_rh_u_qnrain
            if ( remain_regcoeff12(5) .eq. 0 ) then
               regcoeff_rh_u_qnrain = 0.
            end if
         case ('regcoeff_qcloud_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_qcloud_u_qnrain
            if ( remain_regcoeff12(6) .eq. 0 ) then
               regcoeff_qcloud_u_qnrain = 0.
            end if
         case ('regcoeff_qice_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_qice_u_qnrain
            if ( remain_regcoeff12(7) .eq. 0 ) then
               regcoeff_qice_u_qnrain = 0.
            end if
         case ('regcoeff_qrain_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_qrain_u_qnrain
            if ( remain_regcoeff12(8) .eq. 0 ) then
               regcoeff_qrain_u_qnrain = 0.
            end if
         case ('regcoeff_qsnow_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_qsnow_u_qnrain
            if ( remain_regcoeff12(9) .eq. 0 ) then
               regcoeff_qsnow_u_qnrain = 0.
            end if
         case ('regcoeff_qgraup_u_qnrain')
            write(*,*) "---Reading "//trim(adjustl(var80c))//"---"
            read (be_unit) regcoeff_qgraup_u_qnrain
            if ( remain_regcoeff12(10) .eq. 0 ) then
               regcoeff_qgraup_u_qnrain = 0.
            end if

         case ('regcoeff_completed')
               write(*,*) "---Reading regcoeff completed---"
               exit
         case default;
            message(1)=' Read problem in regression coefficients in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read regression coefficients for variable: ',trim(adjustl(var80c))
            call da_error(__FILE__,__LINE__,message(1:2))
         end select
      end do
   end if

      ! 1.2 Fill regression coeff. array for model BE, vertical dimension is kz:

   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then	!Zheng,20200522
      allocate (be%reg_psi_chi  (1:jy,1:kz))
      allocate (be%reg_psi_t    (1:jy,1:kz,1:kz))
      allocate (be%reg_psi_ps   (1:jy,1:kz))
      if ( cv_options == 6 ) then
         allocate (be%reg_psi_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_t  (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_ps (1:jy,1:kz))
         allocate (be%reg_chi_u_rh (1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_ps_u_rh  (1:jy,1:kz))
      end if
      if ( cv_options == 9 .and. use_regcoeff ) then				!Zheng,20200513
         allocate (be%reg_psi_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_psi_qcloud(1:jy,1:kz,1:kz))
         allocate (be%reg_psi_qice(1:jy,1:kz,1:kz))
         allocate (be%reg_psi_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_psi_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_psi_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_t  (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_ps (1:jy,1:kz))
         allocate (be%reg_chi_u_rh (1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_qcloud(1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_qice(1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_chi_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_ps(1:jy,1:kz))
         allocate (be%reg_t_u_rh   (1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qcloud(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qice(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_ps_u_rh  (1:jy,1:kz))
         allocate (be%reg_ps_u_qcloud(1:jy,1:kz))
         allocate (be%reg_ps_u_qice(1:jy,1:kz))
         allocate (be%reg_ps_u_qrain(1:jy,1:kz))
         allocate (be%reg_ps_u_qsnow(1:jy,1:kz))
         allocate (be%reg_ps_u_qgraup(1:jy,1:kz))
         allocate (be%reg_rh_u_qcloud(1:jy,1:kz,1:kz))
         allocate (be%reg_rh_u_qice(1:jy,1:kz,1:kz))
         allocate (be%reg_rh_u_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_rh_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_rh_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_qcloud_u_qice(1:jy,1:kz,1:kz))
         allocate (be%reg_qcloud_u_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qcloud_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_qcloud_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_qice_u_qrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qice_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_qice_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_qrain_u_qsnow(1:jy,1:kz,1:kz))
         allocate (be%reg_qrain_u_qgraup(1:jy,1:kz,1:kz))
         allocate (be%reg_qsnow_u_qgraup(1:jy,1:kz,1:kz))
         if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
            allocate (be%reg_psi_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_chi_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_t_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_ps_u_qnrain(1:jy,1:kz))
            allocate (be%reg_rh_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_qcloud_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_qice_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_qrain_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_qsnow_u_qnrain(1:jy,1:kz,1:kz))
            allocate (be%reg_qgraup_u_qnrain(1:jy,1:kz,1:kz))
         end if
      end if
   else if ( cv_options == 7 .and. use_regcoeff ) then							!Zheng,20200513
      allocate (be%reg_u_v(1:jy,1:kz,1:kz))
      allocate (be%reg_u_t(1:jy,1:kz,1:kz))
      allocate (be%reg_u_ps(1:jy,1:kz))
      allocate (be%reg_u_rh(1:jy,1:kz,1:kz))
      allocate (be%reg_u_qcloud(1:jy,1:kz,1:kz))
      allocate (be%reg_u_qice(1:jy,1:kz,1:kz))
      allocate (be%reg_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_t(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_ps(1:jy,1:kz))
      allocate (be%reg_v_u_rh(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_qcloud(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_qice(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_v_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_ps(1:jy,1:kz))
      allocate (be%reg_t_u_rh(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_qcloud(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_qice(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_t_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_ps_u_rh(1:jy,1:kz))
      allocate (be%reg_ps_u_qcloud(1:jy,1:kz))
      allocate (be%reg_ps_u_qice(1:jy,1:kz))
      allocate (be%reg_ps_u_qrain(1:jy,1:kz))
      allocate (be%reg_ps_u_qsnow(1:jy,1:kz))
      allocate (be%reg_ps_u_qgraup(1:jy,1:kz))
      allocate (be%reg_rh_u_qcloud(1:jy,1:kz,1:kz))
      allocate (be%reg_rh_u_qice(1:jy,1:kz,1:kz))
      allocate (be%reg_rh_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_rh_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_rh_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_qcloud_u_qice(1:jy,1:kz,1:kz))
      allocate (be%reg_qcloud_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_qcloud_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_qcloud_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_qice_u_qrain(1:jy,1:kz,1:kz))
      allocate (be%reg_qice_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_qice_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_qrain_u_qsnow(1:jy,1:kz,1:kz))
      allocate (be%reg_qrain_u_qgraup(1:jy,1:kz,1:kz))
      allocate (be%reg_qsnow_u_qgraup(1:jy,1:kz,1:kz))
      if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
         allocate (be%reg_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_v_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_t_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_ps_u_qnrain(1:jy,1:kz))     ! ps 2D
         allocate (be%reg_rh_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qcloud_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qice_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qrain_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qsnow_u_qnrain(1:jy,1:kz,1:kz))
         allocate (be%reg_qgraup_u_qnrain(1:jy,1:kz,1:kz))
      end if
   end if

      !-----------for interpolating CV5--------------------------------------------------------------
   if ( (cv_options == 5 .or. cv_options == 6) .or. (cv_options == 9 .and. use_regcoeff) ) then				!Zheng,20200513
      if ( interpolate_stats ) then
         allocate (reg_psi_chi0(1:nj,1:nk))
         allocate (reg_psi_ps0 (1:nj,1:nk))
         allocate (reg_psi_t0  (1:nj,1:nk,1:nk))
      end if

      be%reg_psi_chi    = 0.
      be%reg_psi_t      = 0.
      be%reg_psi_ps     = 0.
      if ( cv_options == 6 ) then
         be%reg_psi_rh     = 0.
         be%reg_chi_u_t    = 0.
         be%reg_chi_u_ps   = 0.
         be%reg_chi_u_rh   = 0.
         be%reg_t_u_rh     = 0.
         be%reg_ps_u_rh     = 0.								!Zheng,20200513
      end if
      if ( cv_options == 9 .and. use_regcoeff ) then						!Zheng,20200514
         be%reg_psi_rh     = 0.
         be%reg_psi_qcloud = 0.
         be%reg_psi_qice = 0.
         be%reg_psi_qrain = 0.
         be%reg_psi_qsnow = 0.
         be%reg_psi_qgraup = 0.
         be%reg_chi_u_t    = 0.
         be%reg_chi_u_ps   = 0.
         be%reg_chi_u_rh   = 0.
         be%reg_chi_u_qcloud = 0.
         be%reg_chi_u_qice = 0.
         be%reg_chi_u_qrain = 0.
         be%reg_chi_u_qsnow = 0.
         be%reg_chi_u_qgraup = 0.
         be%reg_t_u_ps = 0.
         be%reg_t_u_rh     = 0.
         be%reg_t_u_qcloud = 0.
         be%reg_t_u_qice = 0.
         be%reg_t_u_qrain = 0.
         be%reg_t_u_qsnow = 0.
         be%reg_t_u_qgraup = 0.
         be%reg_ps_u_rh = 0.
         be%reg_ps_u_qcloud = 0.
         be%reg_ps_u_qice = 0.
         be%reg_ps_u_qrain = 0.
         be%reg_ps_u_qsnow = 0.
         be%reg_ps_u_qgraup = 0.
         be%reg_rh_u_qcloud = 0.
         be%reg_rh_u_qice = 0.
         be%reg_rh_u_qrain = 0.
         be%reg_rh_u_qsnow = 0.
         be%reg_rh_u_qgraup = 0.
         be%reg_qcloud_u_qice = 0.
         be%reg_qcloud_u_qrain = 0.
         be%reg_qcloud_u_qsnow = 0.
         be%reg_qcloud_u_qgraup = 0.
         be%reg_qice_u_qrain = 0.
         be%reg_qice_u_qsnow = 0.
         be%reg_qice_u_qgraup = 0.
         be%reg_qrain_u_qsnow = 0.
         be%reg_qrain_u_qgraup = 0.
         be%reg_qsnow_u_qgraup = 0.
         if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
            be%reg_psi_qnrain = 0.
            be%reg_chi_u_qnrain = 0.
            be%reg_t_u_qnrain = 0.
            be%reg_ps_u_qnrain = 0.
            be%reg_rh_u_qnrain = 0.
            be%reg_qcloud_u_qnrain = 0.
            be%reg_qice_u_qnrain = 0.
            be%reg_qrain_u_qnrain = 0.
            be%reg_qsnow_u_qnrain = 0.
            be%reg_qgraup_u_qnrain = 0.
         end if
      end if
   else if ( cv_options == 7 .and. use_regcoeff ) then						!Zheng,20200514
      be%reg_u_v = 0.
      be%reg_u_t = 0.
      be%reg_u_ps = 0.
      be%reg_u_rh = 0.
      be%reg_u_qcloud = 0.
      be%reg_u_qice = 0.
      be%reg_u_qrain = 0.
      be%reg_u_qsnow = 0.
      be%reg_u_qgraup = 0.
      be%reg_v_u_t = 0.
      be%reg_v_u_ps = 0.
      be%reg_v_u_rh = 0.
      be%reg_v_u_qcloud = 0.
      be%reg_v_u_qice = 0.
      be%reg_v_u_qrain = 0.
      be%reg_v_u_qsnow = 0.
      be%reg_v_u_qgraup = 0.
      be%reg_t_u_ps = 0.
      be%reg_t_u_rh     = 0.
      be%reg_t_u_qcloud = 0.
      be%reg_t_u_qice = 0.
      be%reg_t_u_qrain = 0.
      be%reg_t_u_qsnow = 0.
      be%reg_t_u_qgraup = 0.
      be%reg_ps_u_rh = 0.
      be%reg_ps_u_qcloud = 0.
      be%reg_ps_u_qice = 0.
      be%reg_ps_u_qrain = 0.
      be%reg_ps_u_qsnow = 0.
      be%reg_ps_u_qgraup = 0.
      be%reg_rh_u_qcloud = 0.
      be%reg_rh_u_qice = 0.
      be%reg_rh_u_qrain = 0.
      be%reg_rh_u_qsnow = 0.
      be%reg_rh_u_qgraup = 0.
      be%reg_qcloud_u_qice = 0.
      be%reg_qcloud_u_qrain = 0.
      be%reg_qcloud_u_qsnow = 0.
      be%reg_qcloud_u_qgraup = 0.
      be%reg_qice_u_qrain = 0.
      be%reg_qice_u_qsnow = 0.
      be%reg_qice_u_qgraup = 0.
      be%reg_qrain_u_qsnow = 0.
      be%reg_qrain_u_qgraup = 0.
      be%reg_qsnow_u_qgraup = 0.
      if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
         be%reg_u_qnrain = 0.
         be%reg_v_u_qnrain = 0.
         be%reg_t_u_qnrain = 0.
         be%reg_ps_u_qnrain = 0.
         be%reg_rh_u_qnrain = 0.
         be%reg_qcloud_u_qnrain = 0.
         be%reg_qice_u_qnrain = 0.
         be%reg_qrain_u_qnrain = 0.
         be%reg_qsnow_u_qnrain = 0.
         be%reg_qgraup_u_qnrain = 0.
      end if
   end if

   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then		!Zheng,20200522
      if ( interpolate_stats ) then
         do k=1,nk
            do j=1,nj
               reg_psi_chi0(j,k) = regcoeff_psi_chi(bin(1,j,k))
               b = bin2d(1,j)
               reg_psi_ps0(j,k)  = regcoeff_psi_ps(k,b)
            end do
         end do
         do k=1,nk
            do i=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  reg_psi_t0(j,i,k) = regcoeff_psi_t(i,k,b)
               end do
            end do
         end do
      end if   ! if interpolate_stats

      if ( .not. interpolate_stats ) then  !default
         do k=1,nk
            do j=1,nj
               be%reg_psi_chi(j,k) = psi_chi_factor * regcoeff_psi_chi(bin(1,j,k))
               b = bin2d(1,j)
               be%reg_psi_ps(j,k)  = psi_ps_factor  * regcoeff_psi_ps(k,b)
            end do
         end do
         do k=1,nk
            do i=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  be%reg_psi_t(j,i,k) = psi_t_factor * regcoeff_psi_t(i,k,b)
               end do
            end do
         end do

         if ( cv_options == 6 ) then
            do k=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  be%reg_ps_u_rh(j,k)  = ps_u_rh_factor  * regcoeff_ps_u_rh(k,b)
                  be%reg_chi_u_ps(j,k) = chi_u_ps_factor * regcoeff_chi_u_ps(k,b)
               end do
            end do
            do k=1,nk
               do i=1,nk
                  do j=1,nj
                     b = bin2d(1,j)
                     be%reg_psi_rh(j,i,k)   = psi_rh_factor   * regcoeff_psi_rh(i,k,b)
                     be%reg_chi_u_t(j,i,k)  = chi_u_t_factor  * regcoeff_chi_u_t(i,k,b)
                     be%reg_chi_u_rh(j,i,k) = chi_u_rh_factor * regcoeff_chi_u_rh(i,k,b)
                     be%reg_t_u_rh(j,i,k)   = t_u_rh_factor   * regcoeff_t_u_rh(i,k,b)
                  end do
               end do
            end do
         end if   ! if cv_options 6

         if ( cv_options == 9 .and. use_regcoeff ) then					!Zheng,20200513
            do k=1,nk
               do j=1,nj
                  b = bin2d(1,j)
                  be%reg_chi_u_ps(j,k)        = chi_u_ps_factor          * regcoeff_chi_u_ps(k,b)
                  be%reg_t_u_ps(j,k)          = t_u_ps_factor            * regcoeff_t_u_ps(k,b)
                  be%reg_ps_u_rh(j,k)         = ps_u_rh_factor           * regcoeff_ps_u_rh(k,b)
                  be%reg_ps_u_qcloud(j,k)     = ps_u_qcloud_factor       * regcoeff_ps_u_qcloud(k,b)
                  be%reg_ps_u_qice(j,k)       = ps_u_qice_factor         * regcoeff_ps_u_qice(k,b)
                  be%reg_ps_u_qrain(j,k)      = ps_u_qrain_factor        * regcoeff_ps_u_qrain(k,b)
                  be%reg_ps_u_qsnow(j,k)      = ps_u_qsnow_factor        * regcoeff_ps_u_qsnow(k,b)
                  be%reg_ps_u_qgraup(j,k)     = ps_u_qgraup_factor       * regcoeff_ps_u_qgraup(k,b)
                  if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
                     be%reg_ps_u_qnrain(j,k)  = ps_u_qnrain_factor  * regcoeff_ps_u_qnrain(k,b)  !2D
                  end if
               end do
            end do
            do k=1,nk
               do i=1,nk
                  do j=1,nj
                     b = bin2d(1,j)
                     be%reg_psi_rh(j,i,k)            = psi_rh_factor            *  regcoeff_psi_rh(i,k,b)
                     be%reg_psi_qcloud(j,i,k)        = psi_qcloud_factor        *  regcoeff_psi_qcloud(i,k,b)
                     be%reg_psi_qice(j,i,k)          = psi_qice_factor          *  regcoeff_psi_qice(i,k,b)
                     be%reg_psi_qrain(j,i,k)         = psi_qrain_factor         *  regcoeff_psi_qrain(i,k,b)
                     be%reg_psi_qsnow(j,i,k)         = psi_qsnow_factor         *  regcoeff_psi_qsnow(i,k,b)
                     be%reg_psi_qgraup(j,i,k)        = psi_qgraup_factor        *  regcoeff_psi_qgraup(i,k,b)
                     be%reg_chi_u_t(j,i,k)           = chi_u_t_factor           *  regcoeff_chi_u_t(i,k,b)
                     be%reg_chi_u_rh(j,i,k)          = chi_u_rh_factor          *  regcoeff_chi_u_rh(i,k,b)
                     be%reg_chi_u_qcloud(j,i,k)      = chi_u_qcloud_factor      *  regcoeff_chi_u_qcloud(i,k,b)
                     be%reg_chi_u_qice(j,i,k)        = chi_u_qice_factor        *  regcoeff_chi_u_qice(i,k,b)
                     be%reg_chi_u_qrain(j,i,k)       = chi_u_qrain_factor       *  regcoeff_chi_u_qrain(i,k,b)
                     be%reg_chi_u_qsnow(j,i,k)       = chi_u_qsnow_factor       *  regcoeff_chi_u_qsnow(i,k,b)
                     be%reg_chi_u_qgraup(j,i,k)      = chi_u_qgraup_factor      *  regcoeff_chi_u_qgraup(i,k,b)
                     be%reg_t_u_rh(j,i,k)            = t_u_rh_factor            *  regcoeff_t_u_rh(i,k,b)
                     be%reg_t_u_qcloud(j,i,k)        = t_u_qcloud_factor        *  regcoeff_t_u_qcloud(i,k,b)
                     be%reg_t_u_qice(j,i,k)          = t_u_qice_factor          *  regcoeff_t_u_qice(i,k,b)
                     be%reg_t_u_qrain(j,i,k)         = t_u_qrain_factor         *  regcoeff_t_u_qrain(i,k,b)
                     be%reg_t_u_qsnow(j,i,k)         = t_u_qsnow_factor         *  regcoeff_t_u_qsnow(i,k,b)
                     be%reg_t_u_qgraup(j,i,k)        = t_u_qgraup_factor        *  regcoeff_t_u_qgraup(i,k,b)
                     be%reg_rh_u_qcloud(j,i,k)       = rh_u_qcloud_factor       *  regcoeff_rh_u_qcloud(i,k,b)
                     be%reg_rh_u_qice(j,i,k)         = rh_u_qice_factor         *  regcoeff_rh_u_qice(i,k,b)
                     be%reg_rh_u_qrain(j,i,k)        = rh_u_qrain_factor        *  regcoeff_rh_u_qrain(i,k,b)
                     be%reg_rh_u_qsnow(j,i,k)        = rh_u_qsnow_factor        *  regcoeff_rh_u_qsnow(i,k,b)
                     be%reg_rh_u_qgraup(j,i,k)       = rh_u_qgraup_factor       *  regcoeff_rh_u_qgraup(i,k,b)
                     be%reg_qcloud_u_qice(j,i,k)     = qcloud_u_qice_factor     *  regcoeff_qcloud_u_qice(i,k,b)
                     be%reg_qcloud_u_qrain(j,i,k)    = qcloud_u_qrain_factor    *  regcoeff_qcloud_u_qrain(i,k,b)
                     be%reg_qcloud_u_qsnow(j,i,k)    = qcloud_u_qsnow_factor    *  regcoeff_qcloud_u_qsnow(i,k,b)
                     be%reg_qcloud_u_qgraup(j,i,k)   = qcloud_u_qgraup_factor   *  regcoeff_qcloud_u_qgraup(i,k,b)
                     be%reg_qice_u_qrain(j,i,k)      = qice_u_qrain_factor      *  regcoeff_qice_u_qrain(i,k,b)
                     be%reg_qice_u_qsnow(j,i,k)      = qice_u_qsnow_factor      *  regcoeff_qice_u_qsnow(i,k,b)
                     be%reg_qice_u_qgraup(j,i,k)     = qice_u_qgraup_factor     *  regcoeff_qice_u_qgraup(i,k,b)
                     be%reg_qrain_u_qsnow(j,i,k)     = qrain_u_qsnow_factor     *  regcoeff_qrain_u_qsnow(i,k,b)
                     be%reg_qrain_u_qgraup(j,i,k)    = qrain_u_qgraup_factor    *  regcoeff_qrain_u_qgraup(i,k,b)
                     be%reg_qsnow_u_qgraup(j,i,k)    = qsnow_u_qgraup_factor    *  regcoeff_qsnow_u_qgraup(i,k,b)
                     if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
                        be%reg_psi_qnrain(j,i,k)     = psi_qnrain_factor        *  regcoeff_psi_qnrain(i,k,b)
                        be%reg_chi_u_qnrain(j,i,k)   = chi_u_qnrain_factor      *  regcoeff_chi_u_qnrain(i,k,b)
                        be%reg_t_u_qnrain(j,i,k)     = t_u_qnrain_factor        *  regcoeff_t_u_qnrain(i,k,b)
                        be%reg_rh_u_qnrain(j,i,k)    = rh_u_qnrain_factor       *  regcoeff_rh_u_qnrain(i,k,b)
                        be%reg_qcloud_u_qnrain(j,i,k)= qcloud_u_qnrain_factor   *  regcoeff_qcloud_u_qnrain(i,k,b)
                        be%reg_qice_u_qnrain(j,i,k)  = qice_u_qnrain_factor     *  regcoeff_qice_u_qnrain(i,k,b)
                        be%reg_qrain_u_qnrain(j,i,k) = qrain_u_qnrain_factor    *  regcoeff_qrain_u_qnrain(i,k,b)
                        be%reg_qsnow_u_qnrain(j,i,k) = qsnow_u_qnrain_factor    *  regcoeff_qsnow_u_qnrain(i,k,b)
                        be%reg_qgraup_u_qnrain(j,i,k)= qgraup_u_qnrain_factor   *  regcoeff_qgraup_u_qnrain(i,k,b)
                     end if
                  end do
               end do
            end do
         end if   ! if cv_options 9
      end if   ! if not interpolate_stats

   else if ( cv_options == 7 .and. use_regcoeff ) then					!Zheng,20200514
      do k=1,nk
         do j=1,nj
            b = bin2d(1,j)
            be%reg_u_ps(j,k)            = u_ps_factor              * regcoeff_u_ps(k,b)
            be%reg_v_u_ps(j,k)          = v_u_ps_factor            * regcoeff_v_u_ps(k,b)
            be%reg_t_u_ps(j,k)          = t_u_ps_factor            * regcoeff_t_u_ps(k,b)
            be%reg_ps_u_rh(j,k)         = ps_u_rh_factor           * regcoeff_ps_u_rh(k,b)
            be%reg_ps_u_qcloud(j,k)     = ps_u_qcloud_factor       * regcoeff_ps_u_qcloud(k,b)
            be%reg_ps_u_qice(j,k)       = ps_u_qice_factor         * regcoeff_ps_u_qice(k,b)
            be%reg_ps_u_qrain(j,k)      = ps_u_qrain_factor        * regcoeff_ps_u_qrain(k,b)
            be%reg_ps_u_qsnow(j,k)      = ps_u_qsnow_factor        * regcoeff_ps_u_qsnow(k,b)
            be%reg_ps_u_qgraup(j,k)     = ps_u_qgraup_factor       * regcoeff_ps_u_qgraup(k,b)
            if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
               be%reg_ps_u_qnrain(j,k)  = ps_u_qnrain_factor  * regcoeff_ps_u_qnrain(k,b)  !2D
            end if
         end do
      end do
      do k=1,nk
         do i=1,nk
            do j=1,nj
               b = bin2d(1,j)
               be%reg_u_v(j,i,k)               = u_v_factor               *  regcoeff_u_v(i,k,b)
               be%reg_u_t(j,i,k)               = u_t_factor               *  regcoeff_u_t(i,k,b)
               be%reg_u_rh(j,i,k)              = u_rh_factor              *  regcoeff_u_rh(i,k,b)
               be%reg_u_qcloud(j,i,k)          = u_qcloud_factor          *  regcoeff_u_qcloud(i,k,b)
               be%reg_u_qice(j,i,k)            = u_qice_factor            *  regcoeff_u_qice(i,k,b)
               be%reg_u_qrain(j,i,k)           = u_qrain_factor           *  regcoeff_u_qrain(i,k,b)
               be%reg_u_qsnow(j,i,k)           = u_qsnow_factor           *  regcoeff_u_qsnow(i,k,b)
               be%reg_u_qgraup(j,i,k)          = u_qgraup_factor          *  regcoeff_u_qgraup(i,k,b)
               be%reg_v_u_t(j,i,k)             = v_u_t_factor             *  regcoeff_v_u_t(i,k,b)
               be%reg_v_u_rh(j,i,k)            = v_u_rh_factor            *  regcoeff_v_u_rh(i,k,b)
               be%reg_v_u_qcloud(j,i,k)        = v_u_qcloud_factor        *  regcoeff_v_u_qcloud(i,k,b)
               be%reg_v_u_qice(j,i,k)          = v_u_qice_factor          *  regcoeff_v_u_qice(i,k,b)
               be%reg_v_u_qrain(j,i,k)         = v_u_qrain_factor         *  regcoeff_v_u_qrain(i,k,b)
               be%reg_v_u_qsnow(j,i,k)         = v_u_qsnow_factor         *  regcoeff_v_u_qsnow(i,k,b)
               be%reg_v_u_qgraup(j,i,k)        = v_u_qgraup_factor        *  regcoeff_v_u_qgraup(i,k,b)
               be%reg_t_u_rh(j,i,k)            = t_u_rh_factor            *  regcoeff_t_u_rh(i,k,b)
               be%reg_t_u_qcloud(j,i,k)        = t_u_qcloud_factor        *  regcoeff_t_u_qcloud(i,k,b)
               be%reg_t_u_qice(j,i,k)          = t_u_qice_factor          *  regcoeff_t_u_qice(i,k,b)
               be%reg_t_u_qrain(j,i,k)         = t_u_qrain_factor         *  regcoeff_t_u_qrain(i,k,b)
               be%reg_t_u_qsnow(j,i,k)         = t_u_qsnow_factor         *  regcoeff_t_u_qsnow(i,k,b)
               be%reg_t_u_qgraup(j,i,k)        = t_u_qgraup_factor        *  regcoeff_t_u_qgraup(i,k,b)
               be%reg_rh_u_qcloud(j,i,k)       = rh_u_qcloud_factor       *  regcoeff_rh_u_qcloud(i,k,b)
               be%reg_rh_u_qice(j,i,k)         = rh_u_qice_factor         *  regcoeff_rh_u_qice(i,k,b)
               be%reg_rh_u_qrain(j,i,k)        = rh_u_qrain_factor        *  regcoeff_rh_u_qrain(i,k,b)
               be%reg_rh_u_qsnow(j,i,k)        = rh_u_qsnow_factor        *  regcoeff_rh_u_qsnow(i,k,b)
               be%reg_rh_u_qgraup(j,i,k)       = rh_u_qgraup_factor       *  regcoeff_rh_u_qgraup(i,k,b)
               be%reg_qcloud_u_qice(j,i,k)     = qcloud_u_qice_factor     *  regcoeff_qcloud_u_qice(i,k,b)
               be%reg_qcloud_u_qrain(j,i,k)    = qcloud_u_qrain_factor    *  regcoeff_qcloud_u_qrain(i,k,b)
               be%reg_qcloud_u_qsnow(j,i,k)    = qcloud_u_qsnow_factor    *  regcoeff_qcloud_u_qsnow(i,k,b)
               be%reg_qcloud_u_qgraup(j,i,k)   = qcloud_u_qgraup_factor   *  regcoeff_qcloud_u_qgraup(i,k,b)
               be%reg_qice_u_qrain(j,i,k)      = qice_u_qrain_factor      *  regcoeff_qice_u_qrain(i,k,b)
               be%reg_qice_u_qsnow(j,i,k)      = qice_u_qsnow_factor      *  regcoeff_qice_u_qsnow(i,k,b)
               be%reg_qice_u_qgraup(j,i,k)     = qice_u_qgraup_factor     *  regcoeff_qice_u_qgraup(i,k,b)
               be%reg_qrain_u_qsnow(j,i,k)     = qrain_u_qsnow_factor     *  regcoeff_qrain_u_qsnow(i,k,b)
               be%reg_qrain_u_qgraup(j,i,k)    = qrain_u_qgraup_factor    *  regcoeff_qrain_u_qgraup(i,k,b)
               be%reg_qsnow_u_qgraup(j,i,k)    = qsnow_u_qgraup_factor    *  regcoeff_qsnow_u_qgraup(i,k,b)
               if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
                  be%reg_u_qnrain(j,i,k)       = u_qnrain_factor          *  regcoeff_u_qnrain(i,k,b)
                  be%reg_v_u_qnrain(j,i,k)     = v_u_qnrain_factor        *  regcoeff_v_u_qnrain(i,k,b)
                  be%reg_t_u_qnrain(j,i,k)     = t_u_qnrain_factor        *  regcoeff_t_u_qnrain(i,k,b)
                  be%reg_rh_u_qnrain(j,i,k)    = rh_u_qnrain_factor       *  regcoeff_rh_u_qnrain(i,k,b)
                  be%reg_qcloud_u_qnrain(j,i,k)= qcloud_u_qnrain_factor   *  regcoeff_qcloud_u_qnrain(i,k,b)
                  be%reg_qice_u_qnrain(j,i,k)  = qice_u_qnrain_factor     *  regcoeff_qice_u_qnrain(i,k,b)
                  be%reg_qrain_u_qnrain(j,i,k) = qrain_u_qnrain_factor    *  regcoeff_qrain_u_qnrain(i,k,b)
                  be%reg_qsnow_u_qnrain(j,i,k) = qsnow_u_qnrain_factor    *  regcoeff_qsnow_u_qnrain(i,k,b)
                  be%reg_qgraup_u_qnrain(j,i,k)= qgraup_u_qnrain_factor   *  regcoeff_qgraup_u_qnrain(i,k,b)
               end if
            end do
         end do
      end do
   end if   ! if cv_options 7

   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then			!Zheng,20200514
      deallocate (regcoeff_psi_chi)
      deallocate (regcoeff_psi_t)
      deallocate (regcoeff_psi_ps)
      if ( cv_options == 6 ) then
         deallocate (regcoeff_psi_rh)
         deallocate (regcoeff_chi_u_t)
         deallocate (regcoeff_chi_u_ps)
         deallocate (regcoeff_chi_u_rh)
         deallocate (regcoeff_t_u_rh)
         deallocate (regcoeff_ps_u_rh)
      end if
      if ( cv_options == 9 .and. use_regcoeff ) then								!Zheng,20200514
         deallocate (regcoeff_psi_rh)
         deallocate (regcoeff_psi_qcloud)
         deallocate (regcoeff_psi_qice)
         deallocate (regcoeff_psi_qrain)
         deallocate (regcoeff_psi_qsnow)
         deallocate (regcoeff_psi_qgraup)
         deallocate (regcoeff_chi_u_t)
         deallocate (regcoeff_chi_u_ps)
         deallocate (regcoeff_chi_u_rh)
         deallocate (regcoeff_chi_u_qcloud)
         deallocate (regcoeff_chi_u_qice)
         deallocate (regcoeff_chi_u_qrain)
         deallocate (regcoeff_chi_u_qsnow)
         deallocate (regcoeff_chi_u_qgraup)
         deallocate (regcoeff_t_u_ps)
         deallocate (regcoeff_t_u_rh)
         deallocate (regcoeff_t_u_qcloud)
         deallocate (regcoeff_t_u_qice)
         deallocate (regcoeff_t_u_qrain)
         deallocate (regcoeff_t_u_qsnow)
         deallocate (regcoeff_t_u_qgraup)
         deallocate (regcoeff_ps_u_rh)
         deallocate (regcoeff_ps_u_qcloud)
         deallocate (regcoeff_ps_u_qice)
         deallocate (regcoeff_ps_u_qrain)
         deallocate (regcoeff_ps_u_qsnow)
         deallocate (regcoeff_ps_u_qgraup)
         deallocate (regcoeff_rh_u_qcloud)
         deallocate (regcoeff_rh_u_qice)
         deallocate (regcoeff_rh_u_qrain)
         deallocate (regcoeff_rh_u_qsnow)
         deallocate (regcoeff_rh_u_qgraup)
         deallocate (regcoeff_qcloud_u_qice)
         deallocate (regcoeff_qcloud_u_qrain)
         deallocate (regcoeff_qcloud_u_qsnow)
         deallocate (regcoeff_qcloud_u_qgraup)
         deallocate (regcoeff_qice_u_qrain)
         deallocate (regcoeff_qice_u_qsnow)
         deallocate (regcoeff_qice_u_qgraup)
         deallocate (regcoeff_qrain_u_qsnow)
         deallocate (regcoeff_qrain_u_qgraup)
         deallocate (regcoeff_qsnow_u_qgraup)
         if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
            deallocate (regcoeff_psi_qnrain)
            deallocate (regcoeff_chi_u_qnrain)
            deallocate (regcoeff_t_u_qnrain)
            deallocate (regcoeff_ps_u_qnrain)     ! ps 2D
            deallocate (regcoeff_rh_u_qnrain)
            deallocate (regcoeff_qcloud_u_qnrain)
            deallocate (regcoeff_qice_u_qnrain)
            deallocate (regcoeff_qrain_u_qnrain)
            deallocate (regcoeff_qsnow_u_qnrain)
            deallocate (regcoeff_qgraup_u_qnrain)
         end if
      end if
   else if ( cv_options == 7 .and. use_regcoeff ) then
      deallocate (regcoeff_u_v)
      deallocate (regcoeff_u_t)
      deallocate (regcoeff_u_ps)
      deallocate (regcoeff_u_rh)
      deallocate (regcoeff_u_qcloud)
      deallocate (regcoeff_u_qice)
      deallocate (regcoeff_u_qrain)
      deallocate (regcoeff_u_qsnow)
      deallocate (regcoeff_u_qgraup)
      deallocate (regcoeff_v_u_t)
      deallocate (regcoeff_v_u_ps)
      deallocate (regcoeff_v_u_rh)
      deallocate (regcoeff_v_u_qcloud)
      deallocate (regcoeff_v_u_qice)
      deallocate (regcoeff_v_u_qrain)
      deallocate (regcoeff_v_u_qsnow)
      deallocate (regcoeff_v_u_qgraup)
      deallocate (regcoeff_t_u_ps)
      deallocate (regcoeff_t_u_rh)
      deallocate (regcoeff_t_u_qcloud)
      deallocate (regcoeff_t_u_qice)
      deallocate (regcoeff_t_u_qrain)
      deallocate (regcoeff_t_u_qsnow)
      deallocate (regcoeff_t_u_qgraup)
      deallocate (regcoeff_ps_u_rh)
      deallocate (regcoeff_ps_u_qcloud)
      deallocate (regcoeff_ps_u_qice)
      deallocate (regcoeff_ps_u_qrain)
      deallocate (regcoeff_ps_u_qsnow)
      deallocate (regcoeff_ps_u_qgraup)
      deallocate (regcoeff_rh_u_qcloud)
      deallocate (regcoeff_rh_u_qice)
      deallocate (regcoeff_rh_u_qrain)
      deallocate (regcoeff_rh_u_qsnow)
      deallocate (regcoeff_rh_u_qgraup)
      deallocate (regcoeff_qcloud_u_qice)
      deallocate (regcoeff_qcloud_u_qrain)
      deallocate (regcoeff_qcloud_u_qsnow)
      deallocate (regcoeff_qcloud_u_qgraup)
      deallocate (regcoeff_qice_u_qrain)
      deallocate (regcoeff_qice_u_qsnow)
      deallocate (regcoeff_qice_u_qgraup)
      deallocate (regcoeff_qrain_u_qsnow)
      deallocate (regcoeff_qrain_u_qgraup)
      deallocate (regcoeff_qsnow_u_qgraup)
      if ( qnr_option >= 1 )  then                                       !zhengh 20221014 add qnrain (reg to others) 10
         deallocate (regcoeff_u_qnrain)
         deallocate (regcoeff_v_u_qnrain)
         deallocate (regcoeff_t_u_qnrain)
         deallocate (regcoeff_ps_u_qnrain)     ! ps 2D
         deallocate (regcoeff_rh_u_qnrain)
         deallocate (regcoeff_qcloud_u_qnrain)
         deallocate (regcoeff_qice_u_qnrain)
         deallocate (regcoeff_qrain_u_qnrain)
         deallocate (regcoeff_qsnow_u_qnrain)
         deallocate (regcoeff_qgraup_u_qnrain)
      end if
   end if
      ! 1.3 Domain_averaged regression coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   if ( cv_options /= 7 ) then
   if ( cv_options /= 7 .and. (.not.( cv_options == 9 .and. (.not.use_regcoeff))) ) then	!Zheng,20200522
      if (.not.lat_stats_option) then
         write (unit=message(1), fmt='(a)') ' '
         write (unit=message(2), fmt='(3x, a)') &
            'Using the averaged regression coefficients for unbalanced part'

         !-----------for interpolating CV5--------------------------------------------------------------
         if ( interpolate_stats ) then
            allocate (reg_psi_ps_be (1:nk))
            allocate (reg_psi_chi_be(1:nk))
            allocate (reg_psi_t_be  (1:nk,1:nk))
            do k=1,nk
               reg_psi_ps_be(k)  = sum(reg_psi_ps0 (:,k))/dble(nj)
               reg_psi_chi_be(k) = sum(reg_psi_chi0(:,k))/dble(nj)
            end do
            do m=1,nk
               do k=1,nk
                  reg_psi_t_be(k,m)=sum(reg_psi_t0(:,k,m))/dble(nj)
               end do
            end do
            deallocate (reg_psi_chi0)
            deallocate (reg_psi_ps0 )
            deallocate (reg_psi_t0  )
         end if   ! if interpolate_stats
         !===========for interpolating CV5==============================================================

         if ( (.not. interpolate_stats) .and. (bin_type /= 5) ) then
            do k=1,kz
                be%reg_psi_ps  (:,k) = sum(be%reg_psi_ps (:,k))/dble(jy)
                be%reg_psi_chi (:,k) = sum(be%reg_psi_chi(:,k))/dble(jy)
            end do
            do m=1,kz
               do k=1,kz
                   be%reg_psi_t (:,k,m) = sum(be%reg_psi_t(:,k,m))/dble(jy)
                end do
            end do
            if ( cv_options == 6 ) then
               do k=1,kz
                  be%reg_ps_u_rh (:,k) = sum(be%reg_ps_u_rh (:,k))/dble(jy)
                  be%reg_chi_u_ps(:,k) = sum(be%reg_chi_u_ps(:,k))/dble(jy)
               end do
               do m=1,kz
                  do k=1,kz
                     be%reg_psi_rh  (:,k,m) = sum(be%reg_psi_rh  (:,k,m))/dble(jy)
                     be%reg_chi_u_t (:,k,m) = sum(be%reg_chi_u_t (:,k,m))/dble(jy)
                     be%reg_chi_u_rh(:,k,m) = sum(be%reg_chi_u_rh(:,k,m))/dble(jy)
                     be%reg_t_u_rh  (:,k,m) = sum(be%reg_t_u_rh  (:,k,m))/dble(jy)
                  end do
               end do
            end if   ! if cv_options 6

            if ( cv_options == 9 .and. use_regcoeff ) then							!Zheng,20200514
               do k=1,kz
                  be%reg_chi_u_ps(:,k)        =  sum(be%reg_chi_u_ps(:,k))/dble(jy)
                  be%reg_t_u_ps(:,k)          =  sum(be%reg_t_u_ps(:,k))/dble(jy)
                  be%reg_ps_u_rh(:,k)         =  sum(be%reg_ps_u_rh(:,k))/dble(jy)
                  be%reg_ps_u_qcloud(:,k)     =  sum(be%reg_ps_u_qcloud(:,k))/dble(jy)
                  be%reg_ps_u_qice(:,k)       =  sum(be%reg_ps_u_qice(:,k))/dble(jy)
                  be%reg_ps_u_qrain(:,k)      =  sum(be%reg_ps_u_qrain(:,k))/dble(jy)
                  be%reg_ps_u_qsnow(:,k)      =  sum(be%reg_ps_u_qsnow(:,k))/dble(jy)
                  be%reg_ps_u_qgraup(:,k)     =  sum(be%reg_ps_u_qgraup(:,k))/dble(jy)
                  if( qnr_option >= 1 )  then                                             !zhengh 20221015 add qnrain (reg to others) 10
                     be%reg_ps_u_qnrain(j,k)  = sum(be%reg_ps_u_qnrain(:,k))/dble(jy)  !2D
                  end if                  
               end do
               do m=1,kz
                  do k=1,kz
                     be%reg_psi_rh(:,k,m)            =   sum(be%reg_psi_rh(:,k,m))/dble(jy)
                     be%reg_psi_qcloud(:,k,m)        =   sum(be%reg_psi_qcloud(:,k,m))/dble(jy)
                     be%reg_psi_qice(:,k,m)          =   sum(be%reg_psi_qice(:,k,m))/dble(jy)
                     be%reg_psi_qrain(:,k,m)         =   sum(be%reg_psi_qrain(:,k,m))/dble(jy)
                     be%reg_psi_qsnow(:,k,m)         =   sum(be%reg_psi_qsnow(:,k,m))/dble(jy)
                     be%reg_psi_qgraup(:,k,m)        =   sum(be%reg_psi_qgraup(:,k,m))/dble(jy)
                     be%reg_chi_u_t(:,k,m)           =   sum(be%reg_chi_u_t(:,k,m))/dble(jy)
                     be%reg_chi_u_rh(:,k,m)          =   sum(be%reg_chi_u_rh(:,k,m))/dble(jy)
                     be%reg_chi_u_qcloud(:,k,m)      =   sum(be%reg_chi_u_qcloud(:,k,m))/dble(jy)
                     be%reg_chi_u_qice(:,k,m)        =   sum(be%reg_chi_u_qice(:,k,m))/dble(jy)
                     be%reg_chi_u_qrain(:,k,m)       =   sum(be%reg_chi_u_qrain(:,k,m))/dble(jy)
                     be%reg_chi_u_qsnow(:,k,m)       =   sum(be%reg_chi_u_qsnow(:,k,m))/dble(jy)
                     be%reg_chi_u_qgraup(:,k,m)      =   sum(be%reg_chi_u_qgraup(:,k,m))/dble(jy)
                     be%reg_t_u_rh(:,k,m)            =   sum(be%reg_t_u_rh(:,k,m))/dble(jy)
                     be%reg_t_u_qcloud(:,k,m)        =   sum(be%reg_t_u_qcloud(:,k,m))/dble(jy)
                     be%reg_t_u_qice(:,k,m)          =   sum(be%reg_t_u_qice(:,k,m))/dble(jy)
                     be%reg_t_u_qrain(:,k,m)         =   sum(be%reg_t_u_qrain(:,k,m))/dble(jy)
                     be%reg_t_u_qsnow(:,k,m)         =   sum(be%reg_t_u_qsnow(:,k,m))/dble(jy)
                     be%reg_t_u_qgraup(:,k,m)        =   sum(be%reg_t_u_qgraup(:,k,m))/dble(jy)
                     be%reg_rh_u_qcloud(:,k,m)       =   sum(be%reg_rh_u_qcloud(:,k,m))/dble(jy)
                     be%reg_rh_u_qice(:,k,m)         =   sum(be%reg_rh_u_qice(:,k,m))/dble(jy)
                     be%reg_rh_u_qrain(:,k,m)        =   sum(be%reg_rh_u_qrain(:,k,m))/dble(jy)
                     be%reg_rh_u_qsnow(:,k,m)        =   sum(be%reg_rh_u_qsnow(:,k,m))/dble(jy)
                     be%reg_rh_u_qgraup(:,k,m)       =   sum(be%reg_rh_u_qgraup(:,k,m))/dble(jy)
                     be%reg_qcloud_u_qice(:,k,m)     =   sum(be%reg_qcloud_u_qice(:,k,m))/dble(jy)
                     be%reg_qcloud_u_qrain(:,k,m)    =   sum(be%reg_qcloud_u_qrain(:,k,m))/dble(jy)
                     be%reg_qcloud_u_qsnow(:,k,m)    =   sum(be%reg_qcloud_u_qsnow(:,k,m))/dble(jy)
                     be%reg_qcloud_u_qgraup(:,k,m)   =   sum(be%reg_qcloud_u_qgraup(:,k,m))/dble(jy)
                     be%reg_qice_u_qrain(:,k,m)      =   sum(be%reg_qice_u_qrain(:,k,m))/dble(jy)
                     be%reg_qice_u_qsnow(:,k,m)      =   sum(be%reg_qice_u_qsnow(:,k,m))/dble(jy)
                     be%reg_qice_u_qgraup(:,k,m)     =   sum(be%reg_qice_u_qgraup(:,k,m))/dble(jy)
                     be%reg_qrain_u_qsnow(:,k,m)     =   sum(be%reg_qrain_u_qsnow(:,k,m))/dble(jy)
                     be%reg_qrain_u_qgraup(:,k,m)    =   sum(be%reg_qrain_u_qgraup(:,k,m))/dble(jy)
                     be%reg_qsnow_u_qgraup(:,k,m)    =   sum(be%reg_qsnow_u_qgraup(:,k,m))/dble(jy)

                     if ( qnr_option >= 1 )  then                                       !zhengh 20221015 add qnrain (reg to others) 10
                        be%reg_psi_qnrain(j,i,k)     = sum(be%reg_psi_qnrain(:,k,m))/dble(jy)
                        be%reg_chi_u_qnrain(j,i,k)   = sum(be%reg_chi_u_qnrain(:,k,m))/dble(jy)
                        be%reg_t_u_qnrain(j,i,k)     = sum(be%reg_t_u_qnrain(:,k,m))/dble(jy)
                        be%reg_rh_u_qnrain(j,i,k)    = sum(be%reg_rh_u_qnrain(:,k,m))/dble(jy)
                        be%reg_qcloud_u_qnrain(j,i,k)= sum(be%reg_qcloud_u_qnrain(:,k,m))/dble(jy)
                        be%reg_qice_u_qnrain(j,i,k)  = sum(be%reg_qice_u_qnrain(:,k,m))/dble(jy)
                        be%reg_qrain_u_qnrain(j,i,k) = sum(be%reg_qrain_u_qnrain(:,k,m))/dble(jy)
                        be%reg_qsnow_u_qnrain(j,i,k) = sum(be%reg_qsnow_u_qnrain(:,k,m))/dble(jy)
                        be%reg_qgraup_u_qnrain(j,i,k)= sum(be%reg_qgraup_u_qnrain(:,k,m))/dble(jy)
                     end if
                  end do
               end do
            end if
         end if   ! if not interpolate_stats

      else
         write (unit=message(1), fmt='(a)') ' '
         write (unit=message(2), fmt='(3x, a)') &
            'Using the geographically-dependent regression coefficients for unbalanced part'
      end if

      call da_message(message(1:2))

   else if ( cv_options == 7 .and. use_regcoeff ) then								!Zheng,20200514
      if ( (.not. interpolate_stats) .and. (bin_type /= 5) ) then
         do k=1,kz
            be%reg_u_ps(:,k)            =  sum(be%reg_u_ps(:,k))/dble(jy)
            be%reg_v_u_ps(:,k)          =  sum(be%reg_v_u_ps(:,k))/dble(jy)
            be%reg_t_u_ps(:,k)          =  sum(be%reg_t_u_ps(:,k))/dble(jy)
            be%reg_ps_u_rh(:,k)         =  sum(be%reg_ps_u_rh(:,k))/dble(jy)
            be%reg_ps_u_qcloud(:,k)     =  sum(be%reg_ps_u_qcloud(:,k))/dble(jy)
            be%reg_ps_u_qice(:,k)       =  sum(be%reg_ps_u_qice(:,k))/dble(jy)
            be%reg_ps_u_qrain(:,k)      =  sum(be%reg_ps_u_qrain(:,k))/dble(jy)
            be%reg_ps_u_qsnow(:,k)      =  sum(be%reg_ps_u_qsnow(:,k))/dble(jy)
            be%reg_ps_u_qgraup(:,k)     =  sum(be%reg_ps_u_qgraup(:,k))/dble(jy)
            if( qnr_option >= 1 )  then                                             !zhengh 20221015 add qnrain (reg to others) 10
               be%reg_ps_u_qnrain(j,k)  = sum(be%reg_ps_u_qnrain(:,k))/dble(jy)  !2D
            end if  
         end do
         do m=1,kz
            do k=1,kz
               be%reg_u_v(:,k,m)               =   sum(be%reg_u_v(:,k,m))/dble(jy)
               be%reg_u_t(:,k,m)               =   sum(be%reg_u_t(:,k,m))/dble(jy)
               be%reg_u_rh(:,k,m)              =   sum(be%reg_u_rh(:,k,m))/dble(jy)
               be%reg_u_qcloud(:,k,m)          =   sum(be%reg_u_qcloud(:,k,m))/dble(jy)
               be%reg_u_qice(:,k,m)            =   sum(be%reg_u_qice(:,k,m))/dble(jy)
               be%reg_u_qrain(:,k,m)           =   sum(be%reg_u_qrain(:,k,m))/dble(jy)
               be%reg_u_qsnow(:,k,m)           =   sum(be%reg_u_qsnow(:,k,m))/dble(jy)
               be%reg_u_qgraup(:,k,m)          =   sum(be%reg_u_qgraup(:,k,m))/dble(jy)
               be%reg_v_u_t(:,k,m)             =   sum(be%reg_v_u_t(:,k,m))/dble(jy)
               be%reg_v_u_rh(:,k,m)            =   sum(be%reg_v_u_rh(:,k,m))/dble(jy)
               be%reg_v_u_qcloud(:,k,m)        =   sum(be%reg_v_u_qcloud(:,k,m))/dble(jy)
               be%reg_v_u_qice(:,k,m)          =   sum(be%reg_v_u_qice(:,k,m))/dble(jy)
               be%reg_v_u_qrain(:,k,m)         =   sum(be%reg_v_u_qrain(:,k,m))/dble(jy)
               be%reg_v_u_qsnow(:,k,m)         =   sum(be%reg_v_u_qsnow(:,k,m))/dble(jy)
               be%reg_v_u_qgraup(:,k,m)        =   sum(be%reg_v_u_qgraup(:,k,m))/dble(jy)
               be%reg_t_u_rh(:,k,m)            =   sum(be%reg_t_u_rh(:,k,m))/dble(jy)
               be%reg_t_u_qcloud(:,k,m)        =   sum(be%reg_t_u_qcloud(:,k,m))/dble(jy)
               be%reg_t_u_qice(:,k,m)          =   sum(be%reg_t_u_qice(:,k,m))/dble(jy)
               be%reg_t_u_qrain(:,k,m)         =   sum(be%reg_t_u_qrain(:,k,m))/dble(jy)
               be%reg_t_u_qsnow(:,k,m)         =   sum(be%reg_t_u_qsnow(:,k,m))/dble(jy)
               be%reg_t_u_qgraup(:,k,m)        =   sum(be%reg_t_u_qgraup(:,k,m))/dble(jy)
               be%reg_rh_u_qcloud(:,k,m)       =   sum(be%reg_rh_u_qcloud(:,k,m))/dble(jy)
               be%reg_rh_u_qice(:,k,m)         =   sum(be%reg_rh_u_qice(:,k,m))/dble(jy)
               be%reg_rh_u_qrain(:,k,m)        =   sum(be%reg_rh_u_qrain(:,k,m))/dble(jy)
               be%reg_rh_u_qsnow(:,k,m)        =   sum(be%reg_rh_u_qsnow(:,k,m))/dble(jy)
               be%reg_rh_u_qgraup(:,k,m)       =   sum(be%reg_rh_u_qgraup(:,k,m))/dble(jy)
               be%reg_qcloud_u_qice(:,k,m)     =   sum(be%reg_qcloud_u_qice(:,k,m))/dble(jy)
               be%reg_qcloud_u_qrain(:,k,m)    =   sum(be%reg_qcloud_u_qrain(:,k,m))/dble(jy)
               be%reg_qcloud_u_qsnow(:,k,m)    =   sum(be%reg_qcloud_u_qsnow(:,k,m))/dble(jy)
               be%reg_qcloud_u_qgraup(:,k,m)   =   sum(be%reg_qcloud_u_qgraup(:,k,m))/dble(jy)
               be%reg_qice_u_qrain(:,k,m)      =   sum(be%reg_qice_u_qrain(:,k,m))/dble(jy)
               be%reg_qice_u_qsnow(:,k,m)      =   sum(be%reg_qice_u_qsnow(:,k,m))/dble(jy)
               be%reg_qice_u_qgraup(:,k,m)     =   sum(be%reg_qice_u_qgraup(:,k,m))/dble(jy)
               be%reg_qrain_u_qsnow(:,k,m)     =   sum(be%reg_qrain_u_qsnow(:,k,m))/dble(jy)
               be%reg_qrain_u_qgraup(:,k,m)    =   sum(be%reg_qrain_u_qgraup(:,k,m))/dble(jy)
               be%reg_qsnow_u_qgraup(:,k,m)    =   sum(be%reg_qsnow_u_qgraup(:,k,m))/dble(jy)
               if( qnr_option >= 1 )  then                                       !zhengh 20221015 add qnrain (reg to others) 10
                  be%reg_u_qnrain(j,i,k)       = sum(be%reg_u_qnrain(:,k,m))/dble(jy)
                  be%reg_v_u_qnrain(j,i,k)     = sum(be%reg_v_u_qnrain(:,k,m))/dble(jy)
                  be%reg_t_u_qnrain(j,i,k)     = sum(be%reg_t_u_qnrain(:,k,m))/dble(jy)
                  be%reg_rh_u_qnrain(j,i,k)    = sum(be%reg_rh_u_qnrain(:,k,m))/dble(jy)
                  be%reg_qcloud_u_qnrain(j,i,k)= sum(be%reg_qcloud_u_qnrain(:,k,m))/dble(jy)
                  be%reg_qice_u_qnrain(j,i,k)  = sum(be%reg_qice_u_qnrain(:,k,m))/dble(jy)
                  be%reg_qrain_u_qnrain(j,i,k) = sum(be%reg_qrain_u_qnrain(:,k,m))/dble(jy)
                  be%reg_qsnow_u_qnrain(j,i,k) = sum(be%reg_qsnow_u_qnrain(:,k,m))/dble(jy)
                  be%reg_qgraup_u_qnrain(j,i,k)= sum(be%reg_qgraup_u_qnrain(:,k,m))/dble(jy)
               end if
            end do
         end do
      end if
   end if
!   end if ! cv_options /= 7

   ! 2.0 Load the eigenvector and eigenvalue

   allocate (be1_eval_loc (1:nj,1:nk))
   allocate (be2_eval_loc (1:nj,1:nk))
   allocate (be3_eval_loc (1:nj,1:nk))
   allocate (be4_eval_loc (1:nj,1:nk))
   allocate (be5_eval_loc (1:nj,1:1))

   !Give initial value,Zheng,20200527
   be1_eval_loc = 0.
   be2_eval_loc = 0.
   be3_eval_loc = 0.
   be4_eval_loc = 0.
   be5_eval_loc = 0.

   if ( cloud_cv_options >= 2 ) then
      allocate (be6_eval_loc  (1:nj,1:nk))
      allocate (be7_eval_loc  (1:nj,1:nk))
      allocate (be8_eval_loc  (1:nj,1:nk))
      allocate (be9_eval_loc  (1:nj,1:nk))
      allocate (be10_eval_loc (1:nj,1:nk))

      !Give initial value,Zheng,20200527
      be6_eval_loc = 0.
      be7_eval_loc = 0.
      be8_eval_loc = 0.
      be9_eval_loc = 0.
      be10_eval_loc = 0.

   end if
   if ( qnr_option >= 1 )  then                                                            !zhengh 20221015
      allocate (be12_eval_loc  (1:nj,1:nk))
      !Give initial value,Zhengh,20221015
      be12_eval_loc = 0.
   end if

   if ( use_cv_w ) then
      allocate (be11_eval_loc (1:nj,1:nk))

      !Give initial value,Zheng,20200527
      be11_eval_loc = 0.

   end if

   if (vert_corr == vert_corr_2) then

      allocate (be1_eval_glo(1:nk))
      allocate (be2_eval_glo(1:nk))
      allocate (be3_eval_glo(1:nk))
      allocate (be4_eval_glo(1:nk))
      allocate (be5_eval_glo(1:1))

      !Give initial value,Zheng,20200527
      be1_eval_glo = 0.
      be2_eval_glo = 0.
      be3_eval_glo = 0.
      be4_eval_glo = 0.
      be5_eval_glo = 0.

      if ( cloud_cv_options >= 2 ) then
         allocate (be6_eval_glo(1:nk))
         allocate (be7_eval_glo(1:nk))
         allocate (be8_eval_glo(1:nk))
         allocate (be9_eval_glo(1:nk))
         allocate (be10_eval_glo(1:nk))

         !Give initial value,Zheng,20200527
         be6_eval_glo = 0.
         be7_eval_glo = 0.
         be8_eval_glo = 0.
         be9_eval_glo = 0.
         be10_eval_glo = 0.

      end if

      if ( use_cv_w ) then
         allocate (be11_eval_glo(1:nk))

         !Give initial value,Zheng,20200527
         be11_eval_glo = 0.

      end if
      if ( qnr_option >= 1 )  then                                                      !zhengh 20221015
         allocate (be12_eval_glo(1:nk))
         !Give initial value,Zhengh,20221015
         be12_eval_glo = 0.
      end if

      allocate (be1_evec_loc(1:nj,1:nk,1:nk))
      allocate (be2_evec_loc(1:nj,1:nk,1:nk))
      allocate (be3_evec_loc(1:nj,1:nk,1:nk))
      allocate (be4_evec_loc(1:nj,1:nk,1:nk))
      allocate (be5_evec_loc(1:nj,1: 1,1: 1))

      !Give initial value,Zheng,20200527
      be1_evec_loc = 0.
      be2_evec_loc = 0.
      be3_evec_loc = 0.
      be4_evec_loc = 0.
      be5_evec_loc = 0.

      if ( cloud_cv_options >= 2 ) then
         allocate (be6_evec_loc (1:nj,1:nk,1:nk))
         allocate (be7_evec_loc (1:nj,1:nk,1:nk))
         allocate (be8_evec_loc (1:nj,1:nk,1:nk))
         allocate (be9_evec_loc (1:nj,1:nk,1:nk))
         allocate (be10_evec_loc(1:nj,1:nk,1:nk))

         !Give initial value,Zheng,20200527
         be6_evec_loc = 0.
         be7_evec_loc = 0.
         be8_evec_loc = 0.
         be9_evec_loc = 0.
         be10_evec_loc = 0.

      end if
      if ( use_cv_w ) then
         allocate (be11_evec_loc(1:nj,1:nk,1:nk))

         !Give initial value,Zheng,20200527
         be11_evec_loc = 0.

      end if

      if ( qnr_option >= 1 )  then                                          !zhengh 20221015
         allocate (be12_evec_loc(1:nj,1:nk,1:nk))
         !Give initial value,Zhengh,20221015
         be12_evec_loc = 0.
      end if


      allocate (be1_evec_glo(1:nk,1:nk))
      allocate (be2_evec_glo(1:nk,1:nk))
      allocate (be3_evec_glo(1:nk,1:nk))
      allocate (be4_evec_glo(1:nk,1:nk))
      allocate (be5_evec_glo(1:1,1:1))

      !Give initial value,Zheng,20200527
      be1_evec_glo = 0.
      be2_evec_glo = 0.
      be3_evec_glo = 0.
      be4_evec_glo = 0.
      be5_evec_glo = 0.

      if ( cloud_cv_options >= 2 ) then
         allocate (be6_evec_glo (1:nk,1:nk))
         allocate (be7_evec_glo (1:nk,1:nk))
         allocate (be8_evec_glo (1:nk,1:nk))
         allocate (be9_evec_glo (1:nk,1:nk))
         allocate (be10_evec_glo(1:nk,1:nk))

         !Give initial value,Zheng,20200527
         be6_evec_glo = 0.
         be7_evec_glo = 0.
         be8_evec_glo = 0.
         be9_evec_glo = 0.
         be10_evec_glo = 0.

      end if

      if ( use_cv_w ) then
         allocate (be11_evec_glo(1:nk,1:nk))

         !Give initial value,Zheng,20200527
         be11_evec_glo = 0.

      end if

      if ( qnr_option >= 1 )  then                                          !zhengh 20221015
         allocate (be12_evec_glo(1:nk,1:nk))
         !Give initial value,Zhengh,20221015
         be12_evec_glo = 0.
      end if

   end if

   if ( use_rf ) then
      if ( cloud_cv_options >= 2 ) then
         allocate ( be6_rf_lengthscale(1:kz))
         allocate ( be7_rf_lengthscale(1:kz))
         allocate ( be8_rf_lengthscale(1:kz))
         allocate ( be9_rf_lengthscale(1:kz))
         allocate (be10_rf_lengthscale(1:kz))

         !Give initial value,Zheng,20200527
         be6_rf_lengthscale = 0.
         be7_rf_lengthscale = 0.
         be8_rf_lengthscale = 0.
         be9_rf_lengthscale = 0.
         be10_rf_lengthscale = 0.

      end if
      if ( use_cv_w ) then
         allocate (be11_rf_lengthscale(1:kz))

         !Give initial value,Zheng,20200527
         be11_rf_lengthscale = 0.

      end if
      if ( qnr_option >= 1 )  then                                          !zhengh 20221015
         allocate (be12_rf_lengthscale(1:kz))
         !Give initial value,Zhengh,20221015
         be12_rf_lengthscale = 0.
      end if
      
   end if

   if ( cloud_cv_options == 3 ) then                                              ! ?????????????????????? how to?
      ! hard-coded the v6-v11 BE values here
      be % v6  % name     = "qcloud"
      be % v7  % name     = "qrain"
      be % v8  % name     = "qice"
      be % v9  % name     = "qsnow"
      be % v10 % name     = "qgraup"
!      be6_eval_glo        = 1.0e-6
!      be7_eval_glo        = 1.0e-6
!      be8_eval_glo        = 1.0e-6
!      be9_eval_glo        = 1.0e-6
!      be10_eval_glo       = 1.0e-6
!      be6_eval_glo        = 1.2e-8
!      be7_eval_glo        = 2.0e-8
!      be8_eval_glo        = 2.4e-9
!      be9_eval_glo        = 4.0e-8
!      be10_eval_glo       = 4.0e-8
      be6_eval_glo        = 1.2e-8
      be7_eval_glo        = 3.5e-8
      be8_eval_glo        = 2.4e-9
      be9_eval_glo        = 4.0e-8
      be10_eval_glo       = 4.0e-8
      if ( use_cv_w ) then
         be % v11 % name     = "w"
         be11_eval_glo       = 1.0
      end if
      if ( use_rf ) then
         be6_rf_lengthscale  = 1.0
         be7_rf_lengthscale  = 1.0
         be8_rf_lengthscale  = 1.0
         be9_rf_lengthscale  = 1.0
         be10_rf_lengthscale = 1.0
         if ( use_cv_w ) then
            be11_rf_lengthscale = 1.0
         end if
      end if
   end if

   ! 2.2 Read in the eigenvector and eigenvalue

   print *, '-------- reading eigen vector/value -------'
   if ( cv_options /= 7 ) then						!Zheng,20200628
      do i = 1 , num_cv_3d_basic + num_cv_3d_extra + num_cv_qn + 1			!Zhengh 20221017                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!biao ji!!!!!!!!!!!!!!!!
         read (be_unit) variable
         write(*,*) 'Reading evec/eval/lens for variable:',variable		!test,20200401
         select case( trim(adjustl(variable)) )
         case ('psi')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v1 % name = trim(adjustl(variable))
            read (be_unit)  be1_evec_glo
            read (be_unit)  be1_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be1_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be1_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('chi', 'chi_u')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v2 % name = trim(adjustl(variable))
            read (be_unit)  be2_evec_glo
            read (be_unit)  be2_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be2_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be2_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('t_u', 't')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v3 % name = trim(adjustl(variable))
            read (be_unit)  be3_evec_glo
            read (be_unit)  be3_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be3_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be3_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('rh_u' , 'rh' )
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v4 % name = trim(adjustl(variable))
            read (be_unit)  be4_evec_glo
            read (be_unit)  be4_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be4_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be4_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('ps_u' , 'ps' )
            read (be_unit) nk_2d, num_bins2d
            allocate (evec_loc(1:nk_2d,1:nk_2d,1:num_bins2d))
            allocate (eval_loc(1:nk_2d,     1:num_bins2d))
            be % v5 % name = trim(adjustl(variable))
            read (be_unit)  be5_evec_glo
            read (be_unit)  be5_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be5_evec_loc(j,1:nk_2d,1:nk_2d) = evec_loc(1:nk_2d,1:nk_2d,b)
               be5_eval_loc(j,1:nk_2d  ) = eval_loc(1:nk_2d,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('qcloud_u' , 'qcloud' )
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
               be % v6 % name = trim(adjustl(variable))
               read (be_unit)  be6_evec_glo
               read (be_unit)  be6_eval_glo
               read (be_unit)  evec_loc
               read (be_unit)  eval_loc
               do j=1,nj
                  b = bin2d(1,j)
                  be6_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                  be6_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
               end do
           deallocate (evec_loc)
           deallocate (eval_loc)
   
        case ('qrain_u' , 'qrain' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v7 % name = trim(adjustl(variable))
              read (be_unit)  be7_evec_glo
              read (be_unit)  be7_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be7_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be7_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qice_u' , 'qice' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v8 % name = trim(adjustl(variable))
              read (be_unit)  be8_evec_glo
              read (be_unit)  be8_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be8_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be8_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qsnow_u' , 'qsnow' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v9 % name = trim(adjustl(variable))
              read (be_unit)  be9_evec_glo
              read (be_unit)  be9_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be9_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be9_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qgraup_u' , 'qgraup' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v10 % name = trim(adjustl(variable))
              read (be_unit)  be10_evec_glo
              read (be_unit)  be10_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be10_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be10_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('w')											!Zheng,20200529
           if ( use_cv_w ) then
              read (be_unit) nk, num_bins2d
              allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
              allocate (eval_loc(1:nk,     1:num_bins2d))
                 be % v11 % name = trim(adjustl(variable))
                 read (be_unit)  be11_evec_glo
                 read (be_unit)  be11_eval_glo
                 read (be_unit)  evec_loc
                 read (be_unit)  eval_loc
                 do j=1,nj
                    b = bin2d(1,j)
                    be11_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                    be11_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
                 end do
             deallocate (evec_loc)
             deallocate (eval_loc)
          else
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you have set use_cv_w = true !'
            call da_error(__FILE__,__LINE__,message(1:3))
          end if

        case ('qnrain_u' , 'qnrain' )                           !zhengh 20221015
          if ( qnr_option >= 1 ) then
              read (be_unit) nk, num_bins2d
              allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
              allocate (eval_loc(1:nk,     1:num_bins2d))
                 be % v12 % name = trim(adjustl(variable))
                 read (be_unit)  be12_evec_glo
                 read (be_unit)  be12_eval_glo
                 read (be_unit)  evec_loc
                 read (be_unit)  eval_loc
                 do j=1,nj
                    b = bin2d(1,j)
                    be12_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                    be12_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
                 end do
             deallocate (evec_loc)
             deallocate (eval_loc)
          else
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you have set qnr_option >= 1 !'
            call da_error(__FILE__,__LINE__,message(1:3))
          end if
   
         case default;
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you are using the correct be.dat file for your cv_options setting!'
            call da_error(__FILE__,__LINE__,message(1:3))
         end select
   
      end do !all variables

   else						!CV7,Zheng,20200628

      do i = 1 , num_cv_3d_basic + num_cv_3d_extra + num_cv_qn + 1			!Zhengh 20221017
         read (be_unit) variable
         write(*,*) 'Reading evec/eval/lens for variable:',variable		!test,20200401
         select case( trim(adjustl(variable)) )
         case ('u')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v1 % name = trim(adjustl(variable))
            read (be_unit)  be1_evec_glo
            read (be_unit)  be1_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be1_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be1_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('v', 'v_u')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v2 % name = trim(adjustl(variable))
            read (be_unit)  be2_evec_glo
            read (be_unit)  be2_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be2_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be2_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('t_u', 't')
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v3 % name = trim(adjustl(variable))
            read (be_unit)  be3_evec_glo
            read (be_unit)  be3_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be3_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be3_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('rh_u' , 'rh' )
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
            be % v4 % name = trim(adjustl(variable))
            read (be_unit)  be4_evec_glo
            read (be_unit)  be4_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be4_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
               be4_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('ps_u' , 'ps' )
            read (be_unit) nk_2d, num_bins2d
            allocate (evec_loc(1:nk_2d,1:nk_2d,1:num_bins2d))
            allocate (eval_loc(1:nk_2d,     1:num_bins2d))
            be % v5 % name = trim(adjustl(variable))
            read (be_unit)  be5_evec_glo
            read (be_unit)  be5_eval_glo
            read (be_unit)  evec_loc
            read (be_unit)  eval_loc
            do j=1,nj
               b = bin2d(1,j)
               be5_evec_loc(j,1:nk_2d,1:nk_2d) = evec_loc(1:nk_2d,1:nk_2d,b)
               be5_eval_loc(j,1:nk_2d  ) = eval_loc(1:nk_2d,b)
            end do
   	 deallocate (evec_loc)
   	 deallocate (eval_loc)
   
         case ('qcloud_u' , 'qcloud' )
            read (be_unit) nk, num_bins2d
            allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
            allocate (eval_loc(1:nk,     1:num_bins2d))
               be % v6 % name = trim(adjustl(variable))
               read (be_unit)  be6_evec_glo
               read (be_unit)  be6_eval_glo
               read (be_unit)  evec_loc
               read (be_unit)  eval_loc
               do j=1,nj
                  b = bin2d(1,j)
                  be6_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                  be6_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
               end do
           deallocate (evec_loc)
           deallocate (eval_loc)
   
        case ('qrain_u' , 'qrain' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v7 % name = trim(adjustl(variable))
              read (be_unit)  be7_evec_glo
              read (be_unit)  be7_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be7_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be7_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qice_u' , 'qice' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v8 % name = trim(adjustl(variable))
              read (be_unit)  be8_evec_glo
              read (be_unit)  be8_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be8_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be8_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qsnow_u' , 'qsnow' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v9 % name = trim(adjustl(variable))
              read (be_unit)  be9_evec_glo
              read (be_unit)  be9_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be9_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be9_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('qgraup_u' , 'qgraup' )
           read (be_unit) nk, num_bins2d
           allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
           allocate (eval_loc(1:nk,     1:num_bins2d))
              be % v10 % name = trim(adjustl(variable))
              read (be_unit)  be10_evec_glo
              read (be_unit)  be10_eval_glo
              read (be_unit)  evec_loc
              read (be_unit)  eval_loc
              do j=1,nj
                 b = bin2d(1,j)
                 be10_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                 be10_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
              end do
          deallocate (evec_loc)
          deallocate (eval_loc)
   
        case ('w')											!Zheng,20200529
           if ( use_cv_w ) then
              read (be_unit) nk, num_bins2d
              allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
              allocate (eval_loc(1:nk,     1:num_bins2d))
                 be % v11 % name = trim(adjustl(variable))
                 read (be_unit)  be11_evec_glo
                 read (be_unit)  be11_eval_glo
                 read (be_unit)  evec_loc
                 read (be_unit)  eval_loc
                 do j=1,nj
                    b = bin2d(1,j)
                    be11_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                    be11_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
                 end do
             deallocate (evec_loc)
             deallocate (eval_loc)
          else
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you have set use_cv_w = true !'
            call da_error(__FILE__,__LINE__,message(1:3))
          end if

        case ('qnrain_u' , 'qnrain' )                           !zhengh 20221015
          if ( qnr_option >= 1 ) then
              read (be_unit) nk, num_bins2d
              allocate (evec_loc(1:nk,1:nk,1:num_bins2d))
              allocate (eval_loc(1:nk,     1:num_bins2d))
                 be % v12 % name = trim(adjustl(variable))
                 read (be_unit)  be12_evec_glo
                 read (be_unit)  be12_eval_glo
                 read (be_unit)  evec_loc
                 read (be_unit)  eval_loc
                 do j=1,nj
                    b = bin2d(1,j)
                    be12_evec_loc(j,1:nk,1:nk) = evec_loc(1:nk,1:nk,b)
                    be12_eval_loc(j,1:nk  ) = eval_loc(1:nk,b)
                 end do
             deallocate (evec_loc)
             deallocate (eval_loc)
          else
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you have set qnr_option >= 1 !'
            call da_error(__FILE__,__LINE__,message(1:3))
          end if
   
         case default;
            message(1)=' Read problem in eigen vectors/values in BE file '
            write (unit=message(2),fmt='(A,A)') ' Trying to read Eigenvectors for variable: ',trim(adjustl(variable))
            write (unit=message(3),fmt='(A)') ' Make sure you are using the correct be.dat file for your cv_options setting!'
            call da_error(__FILE__,__LINE__,message(1:3))
         end select
   
      end do !all variables

   end if  !cv_options,Zheng,20200628


   if(use_radarobs .and. use_radar_rf .or. use_rad .and. crtm_cloud) then  
      if ( cloud_cv_options == 1 ) be % v4 % name = 'qt   '
   end if

   message(1:5) = ''
   write (unit=message(1),fmt='(3x,A,i2)') 'cloud_cv_options = ', cloud_cv_options
   write (unit=message(2),fmt='(3x,A,7A)') 'WRFDA dry control variables are: ', &
      trim(be % v1 % name), ', ', trim(be % v2 % name), ', ', &
      trim(be % v3 % name), ' and ', trim(be % v5 % name)
   write (unit=message(3),fmt='(3x,A,A)') &
      'WRFDA Humidity control variable is ', trim(be % v4 % name)
   if ( cloud_cv_options >= 2 ) then
      write (unit=message(4),fmt='(3x,A,10A)') 'WRFDA hydrometeor control variables are: ', &
         trim(be % v6 % name), ', ', trim(be % v7 % name), ', ', &
         trim(be % v8 % name), ', ', trim(be % v9 % name), ', and ', &
         trim(be % v10 % name)
   end if

   if ( use_cv_w ) then
      write (unit=message(5),fmt='(3x,A)') &
         'w control variable is activated'
   end if
   if ( qnr_option >= 1 ) then                                  !zhengh 20221014
      write (unit=message(1),fmt='(3x,A,i2)') 'qnr_option = ', qnr_option
      write (unit=message(6),fmt='(3x,A,A)') 'WRFDA number concentration control variables are: ', &
         trim(be % v12 % name)
   end if
   call da_message(message(1:6))

   if (use_rf) then

! 3.0 Load the scale lengths

      ! 3.1 allocate the array for scale lengths
      allocate (rfls1_be(1:nk))
      allocate (rfls2_be(1:nk))
      allocate (rfls3_be(1:nk))
      allocate (rfls4_be(1:nk))

      allocate (be1_rf_lengthscale(1:kz))
      allocate (be2_rf_lengthscale(1:kz))
      allocate (be3_rf_lengthscale(1:kz))
      allocate (be4_rf_lengthscale(1:kz))
      allocate (be5_rf_lengthscale(1:1))


! ===================================================================
! This is used for do_normalize = .T., I am not sure if it is working 
! properly (???), more testing needed. (YRG, 01/12/2012)
!
      b=0			! pointer to variable*mode
      allocate(nij(0:0,0:1,0:0))

      kzs=(/kz,kz,kz,kz,1/)
      be%v1%mz=1; be%v2%mz=1;be%v3%mz=1;be%v4%mz=1;be%v5%mz=1; 
      mzs=(/be%v1%mz,be%v2%mz,be%v3%mz,be%v4%mz,be%v5%mz/)
! ====================================================================

      ! 3.2 read in the scale lengths

      print *, '----- read lengthscale --------'
      if ( cv_options /= 7 ) then					!Zheng,20200628
         do i = 1 , num_cv_3d_basic + num_cv_3d_extra + num_cv_qn + 1				!Zhengh,20221017
            read (be_unit) variable
            select case( trim(adjustl(variable)) )
            case ('psi')
               read(be_unit) rfls1_be
               if (.not.interpolate_stats) be1_rf_lengthscale = rfls1_be
            case ('chi', 'chi_u')
               read(be_unit) rfls2_be
               if (.not.interpolate_stats) be2_rf_lengthscale = rfls2_be
            case ('t_u', 't')
               read(be_unit) rfls3_be
               if (.not.interpolate_stats) be3_rf_lengthscale = rfls3_be
            case ('rh_u' , 'rh')
               read(be_unit) rfls4_be
               if (.not.interpolate_stats) be4_rf_lengthscale = rfls4_be
            case ('ps_u' , 'ps')
               read(be_unit) be5_rf_lengthscale(1:1) 
            case ('qcloud_u' , 'qcloud')
               read(be_unit) be6_rf_lengthscale
            case ('qrain_u' , 'qrain')
               read(be_unit) be7_rf_lengthscale
            case ('qice_u' , 'qice')
               read(be_unit) be8_rf_lengthscale
            case ('qsnow_u' , 'qsnow')
               read(be_unit) be9_rf_lengthscale
            case ('qgraup_u' , 'qgraup')
               read(be_unit) be10_rf_lengthscale
            case ('w')						!Zheng,20200529
               if ( use_cv_w ) then
                  read(be_unit) be11_rf_lengthscale
               else
                  message(1)='Read problem in lengthscales in be.dat'
                  write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
                  write (unit=message(3),fmt='(A)') ' Make sure you have set use_cv_w = true !'
                  call da_error(__FILE__,__LINE__,message(1:3))
               end if

            case ('qnrain_u' , 'qnrain' )                           !zhengh 20221015
               if ( qnr_option >= 1 ) then
                  read(be_unit) be12_rf_lengthscale
               else
                  message(1)='Read problem in lengthscales in be.dat'
                  write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
                  write (unit=message(3),fmt='(A)') ' Make sure you have set qnr_option >= 1 !'
                  call da_error(__FILE__,__LINE__,message(1:3))
               end if
            
            case default;
               message(1)='Read problem in lengthscales in be.dat'
               write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
               call da_error(__FILE__,__LINE__,message(1:2))
            end select
         end do ! all cvs
      else					!CV7,Zheng,20200628
         do i = 1 , num_cv_3d_basic + num_cv_3d_extra + num_cv_qn + 1				!Zhengh 20221017
            read (be_unit) variable
            select case( trim(adjustl(variable)) )
            case ('u')
               read(be_unit) rfls1_be
               if (.not.interpolate_stats) be1_rf_lengthscale = rfls1_be
            case ('v', 'v_u')
               read(be_unit) rfls2_be
               if (.not.interpolate_stats) be2_rf_lengthscale = rfls2_be
            case ('t_u', 't')
               read(be_unit) rfls3_be
               if (.not.interpolate_stats) be3_rf_lengthscale = rfls3_be
            case ('rh_u' , 'rh')
               read(be_unit) rfls4_be
               if (.not.interpolate_stats) be4_rf_lengthscale = rfls4_be
            case ('ps_u' , 'ps')
               read(be_unit) be5_rf_lengthscale(1:1) 
            case ('qcloud_u' , 'qcloud')
               read(be_unit) be6_rf_lengthscale
            case ('qrain_u' , 'qrain')
               read(be_unit) be7_rf_lengthscale
            case ('qice_u' , 'qice')
               read(be_unit) be8_rf_lengthscale
            case ('qsnow_u' , 'qsnow')
               read(be_unit) be9_rf_lengthscale
            case ('qgraup_u' , 'qgraup')
               read(be_unit) be10_rf_lengthscale
            case ('w')						!Zheng,20200529
               if ( use_cv_w ) then
                  read(be_unit) be11_rf_lengthscale
               else
                  message(1)='Read problem in lengthscales in be.dat'
                  write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
                  write (unit=message(3),fmt='(A)') ' Make sure you have set use_cv_w = true !'
                  call da_error(__FILE__,__LINE__,message(1:3))
               end if

            case ('qnrain_u' , 'qnrain' )                           !zhengh 20221015
               if ( qnr_option >= 1 ) then
                  read(be_unit) be12_rf_lengthscale
               else
                  message(1)='Read problem in lengthscales in be.dat'
                  write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
                  write (unit=message(3),fmt='(A)') ' Make sure you have set qnr_option >= 1 !'
                  call da_error(__FILE__,__LINE__,message(1:3))
               end if
            case default;
               message(1)='Read problem in lengthscales in be.dat'
               write(message(2),'("Trying to read lengthscales for variable ",I0,": ",A)')i,trim(adjustl(variable))
               call da_error(__FILE__,__LINE__,message(1:2))
            end select
         end do ! all cvs
      end if	!cv_options,Zheng,20200628



         if( do_normalize )then
            read(be_unit)do_normalize1
            if( do_normalize.neqv.do_normalize1 ) &
               call da_error(__FILE__,__LINE__,(/namv(i)//": do_normalize.neqv.do_normalize1"/))
            read(be_unit)nij
            if( i==1 )allocate(be%sd(nij(0,1,0),nij(0,0,0),sum(mzs)))
            do k=1,mzs(i)
               read(be_unit)be%sd(:,:,b+k)
            enddo
            do k=mzs(i)+1,kzs(i)! read and discard truncated modes:
               read(be_unit)(junk,j=1,nij(0,1,0)*nij(0,0,0))
            enddo
            write(*,'(A,": |sd[",A,"]|=",es9.3)')__FILE__,namv(i),sqrt(sum(be%sd(:,:,b+1:b+mzs(i))**2))
            b=b+mzs(i)          ! point to next variable.
         endif

      if (print_detail_be) then
! Lengthscale variables (except 5!) are arrays of length "kz", so need to account for this in fmt statement
         write(write_fmt,'(A,I0,A)')'(A,',kz,'F10.5)'
         write(unit=message(1),fmt='(A)') 'BE Length scales:'
         write(unit=message(2),fmt=write_fmt) 'be1= ',be1_rf_lengthscale
         write(unit=message(3),fmt=write_fmt) 'be2= ',be2_rf_lengthscale
         write(unit=message(4),fmt=write_fmt) 'be3= ',be3_rf_lengthscale
         write(unit=message(5),fmt=write_fmt) 'be4= ',be4_rf_lengthscale
         write(unit=message(6),fmt='(A,F10.5)') 'be5= ',be5_rf_lengthscale
         call da_message(message(1:6))

         if ( cloud_cv_options >= 2 ) then
            write(unit=message(7),fmt=write_fmt) 'be6= ',be6_rf_lengthscale
            write(unit=message(8),fmt=write_fmt) 'be7= ',be7_rf_lengthscale
            write(unit=message(9),fmt=write_fmt) 'be8= ',be8_rf_lengthscale
            write(unit=message(10),fmt=write_fmt) 'be9= ',be9_rf_lengthscale
            write(unit=message(11),fmt=write_fmt) 'be10= ',be10_rf_lengthscale
            call da_message(message(7:11))
         end if
         if ( use_cv_w ) then
            write(unit=message(12),fmt=write_fmt) 'be11= ',be11_rf_lengthscale
            call da_message(message(12:12))
         end if
         if ( qnr_option >= 1 ) then                                                         !zhengh
            write(unit=message(13),fmt=write_fmt) 'be12= ',be12_rf_lengthscale
            call da_message(message(13:13))
         end if
      end if ! print_detail_be

   end if ! use_rf
!================== for interpolating CV5 =======================================
   ! 3.3 Interpolate the CV5 BE to the model resolution

   if (interpolate_stats) then
 
   ! 3.3.1 Compute the Vertical covariance matrix arrays for input CV5 BE:
 
   ! 3.3.1.1 Allocate the covariance arrays for input CV5 BE
 
      allocate (covm1_be(1:nk,1:nk))
      allocate (covm2_be(1:nk,1:nk))
      allocate (covm3_be(1:nk,1:nk))
      allocate (covm4_be(1:nk,1:nk))
   ! 3.3.1.2 Compute the vertical covariance matrix from Input CV5 BE
   !       eigenvector/eigenvalue:
      call da_Eigen_to_covmatrix(nk, be1_evec_glo, be1_eval_glo, covm1_be)
      call da_Eigen_to_covmatrix(nk, be2_evec_glo, be2_eval_glo, covm2_be)
      call da_Eigen_to_covmatrix(nk, be3_evec_glo, be3_eval_glo, covm3_be)
      call da_Eigen_to_covmatrix(nk, be4_evec_glo, be4_eval_glo, covm4_be)
   ! Deallocate the input CV5 BE Global eigenvector/eigenvalue:
      DEALLOCATE ( be1_eval_glo )
      DEALLOCATE ( be2_eval_glo )
      DEALLOCATE ( be3_eval_glo )
      DEALLOCATE ( be4_eval_glo )
 
      DEALLOCATE ( be1_evec_glo )
      DEALLOCATE ( be2_evec_glo )
      DEALLOCATE ( be3_evec_glo )
      DEALLOCATE ( be4_evec_glo )

   ! Allocate the model CV5 BE Global eigenvector/eigenvalue:
      ALLOCATE ( be1_eval_glo(1:kz) )
      ALLOCATE ( be2_eval_glo(1:kz) )
      ALLOCATE ( be3_eval_glo(1:kz) )
      ALLOCATE ( be4_eval_glo(1:kz) )
 
      ALLOCATE ( be1_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be2_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be3_evec_glo(1:kz,1:kz) )
      ALLOCATE ( be4_evec_glo(1:kz,1:kz) )
 
   ! 3.3.2 Convert the Vertical resolution from input CV5 BE to Model:

   ! 3.3.2.1 Allocate the Covariance matrix arrays for model
 
      allocate (covm1(1:kz,1:kz))
      allocate (covm2(1:kz,1:kz))
      allocate (covm3(1:kz,1:kz))
      allocate (covm4(1:kz,1:kz))
 
   ! 3.3.2.2 Allocate the regression coeff. arrays for model BE
 
      allocate (reg_psi_ps (1:kz))
      allocate (reg_psi_chi(1:kz))
      allocate (reg_psi_t  (1:kz,1:kz))
   ! 3.3.2 Vertical resolution conversion from input BE(nk) to model(kz) BE:
      call da_chg_be_vres(kz, nk, xb%sigmah, eta_be,&
                        reg_psi_chi_be,reg_psi_ps_be,reg_psi_t_be, &
                        reg_psi_chi   ,reg_psi_ps   ,reg_psi_t   , &
                        covm1_be, covm2_be, covm3_be, covm4_be, &
                        covm1   , covm2   , covm3   , covm4   , &
                        rfls1_be, rfls2_be, rfls3_be, rfls4_be, &
                        be1_rf_lengthscale, be2_rf_lengthscale, &
                        be3_rf_lengthscale, be4_rf_lengthscale)
       
   ! 3.3.3 Deallocate the arrays for input CV5 BE:
      deallocate (reg_psi_ps_be )
      deallocate (reg_psi_chi_be)
      deallocate (reg_psi_t_be  )
 
   ! 3.3.4 Assign the Regression coefficients for model BE:
 
   ! 3.3.4.1 Ps and Chi regression coefficients
      do k=1,kz
        do j=1,jy
          be%reg_psi_ps(j,k)  = reg_psi_ps(k)
          be%reg_psi_chi (j,k)= reg_psi_chi(k)
        enddo
      enddo
!          write(6,*)'be%reg_psi_ps(10,10)=',be%reg_psi_ps(10,10)
!          write(6,*)'be%reg_psi_chi(10,10)=',be%reg_psi_chi(10,10)
   ! 3.3.4.2 Temp. regression coefficients
      do m=1,kz
      do k=1,kz
        do j=1,jy
          be%reg_psi_t(j,k,m) = reg_psi_t(k,m)
        enddo
      enddo
      enddo
!          write(6,*)'be%reg_psi_t(10,10,10)=',be%reg_psi_t(10,10,10)
   ! 3.3.4.3 RWite out the Interpolated BE: Regression coefficients
        if (rootproc) then 
        write(be_out_unit) reg_psi_chi
        write(be_out_unit) reg_psi_ps
        write(be_out_unit) reg_psi_t
        endif
   ! 3.3.4.4 Deallocate the domain-averaged Reg. Coeff. arrays:
 
      deallocate (reg_psi_ps )
      deallocate (reg_psi_chi)
      deallocate (reg_psi_t  )
   ! 3.3.5 Re-compute the Eigenvector/eigenvalue from the model covariance matrix
      call da_gen_eigen(kz, covm1, be1_evec_glo, be1_eval_glo)
      call da_gen_eigen(kz, covm2, be2_evec_glo, be2_eval_glo)
      call da_gen_eigen(kz, covm3, be3_evec_glo, be3_eval_glo)
      call da_gen_eigen(kz, covm4, be4_evec_glo, be4_eval_glo)

      deallocate (covm1)
      deallocate (covm2)
      deallocate (covm3)
      deallocate (covm4)
      deallocate (covm1_be)
      deallocate (covm2_be)
      deallocate (covm3_be)
      deallocate (covm4_be)

   ! 3.3.6 Write out the interpolated CV5 BE
 
   ! 3.3.6.1 Eigenvector/eigenvalue:
      variable = be % v1 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be1_evec_glo
         write(be_out_unit) be1_eval_glo
      ! - Because num_bins2d = 1, be1_evec_loc = be1_evec_glo:
         write(be_out_unit) be1_evec_glo
         write(be_out_unit) be1_eval_glo
      endif
      variable = be % v2 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be2_evec_glo
         write(be_out_unit) be2_eval_glo
      ! - Because num_bins2d = 1, be2_evec_loc = be2_evec_glo:
         write(be_out_unit) be2_evec_glo
         write(be_out_unit) be2_eval_glo
      endif
      variable = be % v3 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be3_evec_glo
         write(be_out_unit) be3_eval_glo
      ! - Because num_bins2d = 1, be3_evec_loc = be3_evec_glo:
         write(be_out_unit) be3_evec_glo
         write(be_out_unit) be3_eval_glo
      endif
      variable = be % v4 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) kz, num_bins2d_out
         write(be_out_unit) be4_evec_glo
         write(be_out_unit) be4_eval_glo
      ! - Because num_bins2d = 1, be4_evec_loc = be4_evec_glo:
         write(be_out_unit) be4_evec_glo
         write(be_out_unit) be4_eval_glo
      endif
      variable = be % v5 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) nk_2d, num_bins2d_out
         write(be_out_unit) be5_evec_glo
         write(be_out_unit) be5_eval_glo
      ! - Because num_bins2d = 1, be5_evec_loc = be5_evec_glo:
         write(be_out_unit) be5_evec_glo
         write(be_out_unit) be5_eval_glo
      endif

   ! 3.3.6.2 Scale-length
      variable = be % v1 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be1_rf_lengthscale
      endif
      variable = be % v2 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be2_rf_lengthscale
      endif
      variable = be % v3 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be3_rf_lengthscale
      endif
      variable = be % v4 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be4_rf_lengthscale
      endif
      variable = be % v5 % name
      if (rootproc) then
         write(be_out_unit) variable
         write(be_out_unit) be5_rf_lengthscale
      endif 
   ! - Close the output unit:
 
      close(be_out_unit)
   ! 3.3.7 Deallocate the arrays for reading the input CV5 BE:
 
      DEALLOCATE ( be1_eval_loc )
      DEALLOCATE ( be2_eval_loc )
      DEALLOCATE ( be3_eval_loc )
      DEALLOCATE ( be4_eval_loc )
      DEALLOCATE ( be5_eval_loc )
      DEALLOCATE ( be1_evec_loc )
      DEALLOCATE ( be2_evec_loc )
      DEALLOCATE ( be3_evec_loc )
      DEALLOCATE ( be4_evec_loc )
      DEALLOCATE ( be5_evec_loc )
 
   ! Allocate the arrays for Model local BE:
      ALLOCATE ( be1_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be2_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be3_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be4_eval_loc (1:jy,1:kz) )
      ALLOCATE ( be5_eval_loc (1:jy,1:1) )
      ALLOCATE ( be1_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be2_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be3_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be4_evec_loc(1:jy,1:kz,1:kz) )
      ALLOCATE ( be5_evec_loc(1:jy,1:1,1:1) )
 
   endif ! if interpolate_stats

   if ( use_rf ) then
      deallocate ( rfls1_be )
      deallocate ( rfls2_be )
      deallocate ( rfls3_be )
      deallocate ( rfls4_be )
   end if

   ! 4.0 Check and get the truncated number of the vertical modes
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (vert_corr == vert_corr_2) then

      ! 4.1 Perform checks on eigenvectors:

      if (test_statistics) then
         call da_check_eof_decomposition(be1_eval_glo(:), be1_evec_glo(:,:), be % v1 % name)
         call da_check_eof_decomposition(be2_eval_glo(:), be2_evec_glo(:,:), be % v2 % name)
         call da_check_eof_decomposition(be3_eval_glo(:), be3_evec_glo(:,:), be % v3 % name)
         call da_check_eof_decomposition(be4_eval_glo(:), be4_evec_glo(:,:), be % v4 % name)

         if ( cloud_cv_options == 2 ) then
            call da_check_eof_decomposition(be6_eval_glo(:),  be6_evec_glo(:,:),  be % v6 % name)
            call da_check_eof_decomposition(be7_eval_glo(:),  be7_evec_glo(:,:),  be % v7 % name)
            call da_check_eof_decomposition(be8_eval_glo(:),  be8_evec_glo(:,:),  be % v8 % name)
            call da_check_eof_decomposition(be9_eval_glo(:),  be9_evec_glo(:,:),  be % v9 % name)
            call da_check_eof_decomposition(be10_eval_glo(:), be10_evec_glo(:,:), be % v10 % name)
            if ( use_cv_w ) then
               call da_check_eof_decomposition(be11_eval_glo(:), be11_evec_glo(:,:), be % v11 % name)
            end if
            if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
                call da_check_eof_decomposition(be12_eval_glo(:), be12_evec_glo(:,:), be % v12 % name)
            end if
         end if
      end if

      ! 4.2 Truncate in vertical:

      call da_get_vertical_truncation(max_vert_var1, be1_eval_glo(:), be % v1)
      call da_get_vertical_truncation(max_vert_var2, be2_eval_glo(:), be % v2)
      call da_get_vertical_truncation(max_vert_var3, be3_eval_glo(:), be % v3)
      call da_get_vertical_truncation(max_vert_var4, be4_eval_glo(:), be % v4)

      if ( cloud_cv_options == 2 ) then
         call da_get_vertical_truncation(max_vert_var6, be6_eval_glo(:), be % v6)
         call da_get_vertical_truncation(max_vert_var7, be7_eval_glo(:), be % v7)
         call da_get_vertical_truncation(max_vert_var8, be8_eval_glo(:), be % v8)
         call da_get_vertical_truncation(max_vert_var9, be9_eval_glo(:), be % v9)
         call da_get_vertical_truncation(max_vert_var10,be10_eval_glo(:),be % v10)
         if ( use_cv_w ) then
            call da_get_vertical_truncation(max_vert_var11,be11_eval_glo(:),be % v11)
         end if
         if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
               call da_get_vertical_truncation(max_vert_var12,be12_eval_glo(:),be % v12)
         end if

      else
         if ( jb_factor > 0.0 ) then
            be % v6  % mz = xb % mkz
            be % v7  % mz = xb % mkz
            be % v8  % mz = xb % mkz
            be % v9  % mz = xb % mkz
            be % v10 % mz = xb % mkz
            if ( use_cv_w ) then
               be % v11 % mz = xb % mkz
            end if
            if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
               be % v12 % mz = xb % mkz
            end if
         else
            be % v6  % mz = 0
            be % v7  % mz = 0
            be % v8  % mz = 0
            be % v9  % mz = 0
            be % v10 % mz = 0
            be % v11 % mz = 0
            be % v12 % mz = 0 !zhengh
         end if
      end if

      if (max_vert_var5 == 0.0) then
         be % v5 % mz = 0
      else
         be % v5 % mz = 1
      end if

      write (unit=stdout,fmt=*) ' '

   else

      ! 4.3 no truncated

      be % v1 % mz = xb % mkz
      be % v2 % mz = xb % mkz
      be % v3 % mz = xb % mkz
      be % v4 % mz = xb % mkz
      be % v5 % mz = xb % mkz

      be % v6 % mz = xb % mkz
      be % v7 % mz = xb % mkz
      be % v8 % mz = xb % mkz
      be % v9 % mz = xb % mkz
      be % v10 % mz = xb % mkz
      be % v11 % mz = xb % mkz

      be % v12 % mz = xb % mkz !zhengh,20221015

   end if

   ! 5.0 Initialise control variable space components of header:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! 5.1 Compute the size of the control variables

   be % mix = ix
   be % mjy = jy

   ! 5.2 Transfer errors to error structure:

   call da_allocate_background_errors(jy, kz, be1_eval_glo, be1_evec_glo, be1_eval_loc, &
                                       be1_evec_loc, be % v1)
   call da_allocate_background_errors(jy, kz, be2_eval_glo, be2_evec_glo, be2_eval_loc, &
                                       be2_evec_loc, be % v2)
   call da_allocate_background_errors(jy, kz, be3_eval_glo, be3_evec_glo, be3_eval_loc, &
                                       be3_evec_loc, be % v3)
   call da_allocate_background_errors(jy, kz, be4_eval_glo, be4_evec_glo, be4_eval_loc, &
                                       be4_evec_loc, be % v4)

   ! 5.2.1 transfer the ps_u variance to be % v5:

   call da_allocate_background_errors(jy,  1, be5_eval_glo, be5_evec_glo, be5_eval_loc, &
                                       be5_evec_loc, be % v5)

   if ( cloud_cv_options >=2 ) then
      call da_allocate_background_errors(jy, kz, be6_eval_glo, be6_evec_glo, be6_eval_loc, &
                                         be6_evec_loc, be % v6)
      call da_allocate_background_errors(jy, kz, be7_eval_glo, be7_evec_glo, be7_eval_loc, &
                                         be7_evec_loc, be % v7)
      call da_allocate_background_errors(jy, kz, be8_eval_glo, be8_evec_glo, be8_eval_loc, &
                                         be8_evec_loc, be % v8)
      call da_allocate_background_errors(jy, kz, be9_eval_glo, be9_evec_glo, be9_eval_loc, &
                                         be9_evec_loc, be % v9)
      call da_allocate_background_errors(jy, kz, be10_eval_glo, be10_evec_glo, be10_eval_loc, &
                                         be10_evec_loc, be % v10)
   end if
   if ( use_cv_w ) then
      call da_allocate_background_errors(jy, kz, be11_eval_glo, be11_evec_glo, be11_eval_loc, &
                                         be11_evec_loc, be % v11)
   end if

   if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
      call da_allocate_background_errors(jy, kz, be12_eval_glo, be12_evec_glo, be12_eval_loc, &
                                         be12_evec_loc, be % v12)
   end if


   if (print_detail_be) then
      write (unit=stdout,fmt='(3x,a,i10)') "b) Finished eigenvector processing!"
   end if

   if( use_rf )then		! use recursive filters:

   ! 5.3 Convert the scale lengths in the real distance (meter)

      be1_rf_lengthscale(1:kz) = be1_rf_lengthscale(1:kz) * xb%ds
      be2_rf_lengthscale(1:kz) = be2_rf_lengthscale(1:kz) * xb%ds
      be3_rf_lengthscale(1:kz) = be3_rf_lengthscale(1:kz) * xb%ds
      be4_rf_lengthscale(1:kz) = be4_rf_lengthscale(1:kz) * xb%ds
      be5_rf_lengthscale(1:1)  = be5_rf_lengthscale(1:1)  * xb%ds

      if ( cloud_cv_options >= 2 ) then
         be6_rf_lengthscale(1:kz) =  be6_rf_lengthscale(1:kz) * xb%ds
         be7_rf_lengthscale(1:kz) =  be7_rf_lengthscale(1:kz) * xb%ds
         be8_rf_lengthscale(1:kz) =  be8_rf_lengthscale(1:kz) * xb%ds
         be9_rf_lengthscale(1:kz) =  be9_rf_lengthscale(1:kz) * xb%ds
         be10_rf_lengthscale(1:kz)= be10_rf_lengthscale(1:kz) * xb%ds
      end if
      if ( use_cv_w ) then
         be11_rf_lengthscale(1:kz)= be11_rf_lengthscale(1:kz) * xb%ds
      end if

      if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
         be12_rf_lengthscale(1:kz)= be12_rf_lengthscale(1:kz) * xb%ds
      end if

   else				! use wavelets:
      ! kzs and mzs are used for wavelets
      kzs=(/kz,kz,kz,kz,1/)
      mzs=(/be%v1%mz,be%v2%mz,be%v3%mz,be%v4%mz,be%v5%mz/)

      read(be_unit)k		! trimmed string length
      if( k>len(variable) )then
         write(message(1),'(i0,">",i0," for use_rf=.false.")')k,len(variable)
         call da_error(__FILE__,__LINE__,message(1:1))
      endif
      read(be_unit)variable(1:k)
      if( variable(1:k) /= 'wavelet' )then
         write(message(1),'(A,"/=''wavelet''")')variable(1:k)
         call da_error(__FILE__,__LINE__,message(1:1))
      endif
      b = 0			! pointer to variable*mode
      do m=1,5			! variables loop:
         read(be_unit)k,nk
         if( k>len(variable) )then
            write(message(1),'(i0,">",i0," for ",A)')k,len(variable),namv(m)
            call da_error(__FILE__,__LINE__,message(1:1))
         elseif( nk /= kzs(m) )then
            write(message(1),'(i0,"=nk/=kzs(",i0,")=",i0," for ",A)')nk,m,kzs(m),namv(m)
            call da_error(__FILE__,__LINE__,message(1:1))
         endif
         read(be_unit)variable(1:k)
         if( variable(1:k) /= trim(namv(m)) )then
            write(message(1),'(A,"/=",A)')variable(1:k),trim(namv(m))
            call da_error(__FILE__,__LINE__,message(1:1))
         endif
         if( m==1 )then
!           Possibly reassign namelist do_normalize value:
            read(be_unit)do_normalize,namw,lf,nb
            allocate(nij(0:nb,0:1,0:2))
            read(be_unit)nij	! wavelet indexes
            write(*,'(A,": ")',advance="no")__FILE__
            do i=0,nb		! wavelet-band loop:
               write(*,'(i2,"{",2(i3,","),i3,";",2(i3,","),i3,"}")',advance="no")&
                  i,transpose(nij(i,:,:))
            enddo
            allocate(be%wsd(nij(0,1,2),nij(0,0,2),sum(mzs)), &
!____________________max of {,i}dwtai wavelet scratch workspace sizes:
                     ws(max(maxval(nij(0,:,2)), &
                            2*floor(.5*(real(maxval(nij(0,:,0)))+lf))+lf-2)))
            if( do_normalize )allocate(be%sd(nij(0,1,0),nij(0,0,0),sum(mzs)))
         endif			! if( m==1 )
         do k=1,mzs(m)
!___________mode-k field & wavelet-coefficient std. devs.:
            read(be_unit)be%wsd(:,:,b+k)
            if( do_normalize )read(be_unit)be%sd(:,:,b+k)
         enddo			! mode-k loop.
         do k=mzs(m)+1,kzs(m)	! read and discard truncated modes:
            read(be_unit)(junk,j=1,nij(0,1,2)*nij(0,0,2))
            if( do_normalize )read(be_unit)(junk,j=1,nij(0,1,0)*nij(0,0,0))
         enddo			! mode-k loop.
         write(*,'(A,": |wsd[",A,"]|=",es9.3)')__FILE__,namv(m),sqrt(sum(be%wsd(:,:,b+1:b+mzs(m))**2))
         if( do_normalize )write(*,'(A,": | sd[",A,"]|=",es9.3)')__FILE__,namv(m),sqrt(sum(be%sd(:,:,b+1:b+mzs(m))**2))
         b = b+mzs(m)		! point to next variable.
      enddo			! variables loop
   endif			! if( use_rf )

   ! 6.0 Perform checks on eigenvectors with be data structure:
   if (jb_factor > 0.0 .and. test_statistics) then
      call da_check_eof_decomposition(be%v1%val_g(:), be%v1%evec_g(:,:),&
                                     be%v1%name)
      call da_check_eof_decomposition(be%v2%val_g(:), be%v2%evec_g(:,:),&
                                     be%v2%name)
      call da_check_eof_decomposition(be%v3%val_g(:), be%v3%evec_g(:,:),&
                                     be%v3%name)
      call da_check_eof_decomposition(be%v4%val_g(:), be%v4%evec_g(:,:),&
                                     be%v4%name)

      if ( cloud_cv_options >= 2 ) then
         call da_check_eof_decomposition(be%v6%val_g(:), be%v6%evec_g(:,:),&
                                         be%v6%name)
         call da_check_eof_decomposition(be%v7%val_g(:), be%v7%evec_g(:,:),&
                                         be%v7%name)
         call da_check_eof_decomposition(be%v8%val_g(:), be%v8%evec_g(:,:),&
                                         be%v8%name)
         call da_check_eof_decomposition(be%v9%val_g(:), be%v9%evec_g(:,:),&
                                         be%v9%name)
         call da_check_eof_decomposition(be%v10%val_g(:), be%v10%evec_g(:,:),&
                                         be%v10%name)
      end if
      if ( use_cv_w ) then
         call da_check_eof_decomposition(be%v11%val_g(:), be%v11%evec_g(:,:),&
                                         be%v11%name)
      end if

      if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
         call da_check_eof_decomposition(be%v12%val_g(:), be%v12%evec_g(:,:),&
                                         be%v12%name)
      end if
   end if

   ! 6.1 Close the be unit

   close(be_unit)
   call da_free_unit(be_unit)

   if( use_rf )then
   ! 6.2 Keep the original be % v1, be % v2,...., and lengthscale in the first loop
   !     for the rescaling in the later loops:

      it = 1

      if (max_ext_its > 1 .and. jb_factor > 0.0) then

        if ( rootproc ) then
          write(unit=message(1),fmt='(A,I4)') '>>> Save the variances and scale-lengths in outer-loop', it
          call da_message(message(1:1))
          write(be_rf_unit)  kz, jy, ix, be % v1 % mz, be % v2 % mz, be% v3 % mz, &
                                     be % v4 % mz, be % v5 % mz, xb % ds
          write(be_rf_unit) be % v1 % val, be % v2 % val, be% v3 % val, &
                                     be % v4 % val, be % v5 % val, &
             be1_rf_lengthscale, be2_rf_lengthscale, be3_rf_lengthscale, &
             be4_rf_lengthscale, be5_rf_lengthscale
          if ( cloud_cv_options >= 2 ) then
             write(be_rf_unit) be%v6%val, be%v7%val, be%v8%val, be%v9%val, be%v10%val
             write(be_rf_unit) be6_rf_lengthscale, be7_rf_lengthscale, be8_rf_lengthscale, &
                               be9_rf_lengthscale, be10_rf_lengthscale
          end if
          if ( use_cv_w ) then
             write(be_rf_unit) be%v11%val
             write(be_rf_unit) be11_rf_lengthscale
          end if

          if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
             write(be_rf_unit) be%v12%val
             write(be_rf_unit) be12_rf_lengthscale
          end if
        end if ! rootproc

          if (print_detail_be ) then
             write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                               it, kz, jy, ix, xb % ds
             write(be_print_unit,'("Original val and rf, and mz:",5i5)') &
                      be % v1 % mz, be % v2 % mz, be% v3 % mz, be % v4 % mz, be % v5 % mz
             write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
             write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)
             write(be_print_unit,'(/"scale-length: kz=",i3)') kz
             do i = 1,kz 
               if (i == 1) then
                 write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                   be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i), &
                   be5_rf_lengthscale(i)
               else
                 write(be_print_unit,'(i3,2x,5e15.5)') i,be1_rf_lengthscale(i), &
                   be2_rf_lengthscale(i), be3_rf_lengthscale(i), be4_rf_lengthscale(i)
               endif
             enddo
    
          endif

      endif

      ! 7.0 Apply empirical and recursive filter rescaling factor:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call da_rescale_background_errors(var_scaling1(1), len_scaling1(1), &
                                         xb % ds, be1_rf_lengthscale, be % v1)
      call da_rescale_background_errors(var_scaling2(1), len_scaling2(1), &
                                         xb % ds, be2_rf_lengthscale, be % v2)
      call da_rescale_background_errors(var_scaling3(1), len_scaling3(1), &
                                         xb % ds, be3_rf_lengthscale, be % v3)
      call da_rescale_background_errors(var_scaling4(1), len_scaling4(1), &
                                         xb % ds, be4_rf_lengthscale, be % v4)
      call da_rescale_background_errors(var_scaling5(1), len_scaling5(1), &
                                         xb % ds, be5_rf_lengthscale, be % v5)
      if ( cloud_cv_options >= 2 ) then
         call da_rescale_background_errors(var_scaling6(1), len_scaling6(1), &
                                           xb % ds, be6_rf_lengthscale, be % v6)
         call da_rescale_background_errors(var_scaling7(1), len_scaling7(1), &
                                           xb % ds, be7_rf_lengthscale, be % v7)
         call da_rescale_background_errors(var_scaling8(1), len_scaling8(1), &
                                           xb % ds, be8_rf_lengthscale, be % v8)
         call da_rescale_background_errors(var_scaling9(1), len_scaling9(1), &
                                           xb % ds, be9_rf_lengthscale, be % v9)
         call da_rescale_background_errors(var_scaling10(1), len_scaling10(1), &
                                           xb % ds, be10_rf_lengthscale, be % v10)
      end if
      if ( use_cv_w ) then
         call da_rescale_background_errors(var_scaling11(1), len_scaling11(1), &
                                           xb % ds, be11_rf_lengthscale, be % v11)
      end if

      if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
         call da_rescale_background_errors(var_scaling12(1), len_scaling12(1), &
                                           xb % ds, be12_rf_lengthscale, be % v12)
      end if

   end if ! use_rf

   if (print_detail_be .and. jb_factor > 0.0) then
       write(be_print_unit,'(/"============================================================")')
       write(be_print_unit,'("For outer loop ",i2)') it
       write(be_print_unit,'("it=",i2,2x,"kz=",i3,2x,"jy=",i4,2x,"ix=",i4,2x,"ds=",e12.5)') &
                                                      it, kz, jy, ix, xb % ds
       write(be_print_unit,'("Namelist options specified for this iteration:")')
       write(be_print_unit,'("var_scaling1(it) = ",e12.5,2x,"len_scaling1(it) = "e12.5)')var_scaling1(it),len_scaling1(it)
       write(be_print_unit,'("var_scaling2(it) = ",e12.5,2x,"len_scaling2(it) = "e12.5)')var_scaling2(it),len_scaling2(it)
       write(be_print_unit,'("var_scaling3(it) = ",e12.5,2x,"len_scaling3(it) = "e12.5)')var_scaling3(it),len_scaling3(it)
       write(be_print_unit,'("var_scaling4(it) = ",e12.5,2x,"len_scaling4(it) = "e12.5)')var_scaling4(it),len_scaling4(it)
       write(be_print_unit,'("var_scaling5(it) = ",e12.5,2x,"len_scaling5(it) = "e12.5)')var_scaling5(it),len_scaling5(it)
       write(be_print_unit,'("Background error statistics for this iteration:")')
       write(be_print_unit,'("mz=",i3,2x,"be%v1%val:"/(10e12.5))') be%v1%mz, be%v1%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v2%val:"/(10e12.5))') be%v2%mz, be%v2%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v3%val:"/(10e12.5))') be%v3%mz, be%v3%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v4%val:"/(10e12.5))') be%v4%mz, be%v4%val(1,:)
       write(be_print_unit,'("mz=",i3,2x,"be%v5%val:"/(10e12.5))') be%v5%mz, be%v5%val(1,:)

       if ( cloud_cv_options >= 2 ) then
          write(be_print_unit,'("mz=",i3,2x,"be%v6%val:"/(10e12.5))') be%v6%mz, be%v6%val(1,:)
          write(be_print_unit,'("mz=",i3,2x,"be%v7%val:"/(10e12.5))') be%v7%mz, be%v7%val(1,:)
          write(be_print_unit,'("mz=",i3,2x,"be%v8%val:"/(10e12.5))') be%v8%mz, be%v8%val(1,:)
          write(be_print_unit,'("mz=",i3,2x,"be%v9%val:"/(10e12.5))') be%v9%mz, be%v9%val(1,:)
          write(be_print_unit,'("mz=",i3,2x,"be%v10%val:"/(10e12.5))') be%v10%mz, be%v10%val(1,:)
       end if
       if ( use_cv_w ) then
          write(be_print_unit,'("mz=",i3,2x,"be%v11%val:"/(10e12.5))') be%v11%mz, be%v11%val(1,:)
       end if

       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       write(be_print_unit,'("be%v1%rf_alpha:"/(10e12.5))') be % v1 % rf_alpha(:)
       write(be_print_unit,'("be%v2%rf_alpha:"/(10e12.5))') be % v2 % rf_alpha(:)
       write(be_print_unit,'("be%v3%rf_alpha:"/(10e12.5))') be % v3 % rf_alpha(:)
       write(be_print_unit,'("be%v4%rf_alpha:"/(10e12.5))') be % v4 % rf_alpha(:)
       write(be_print_unit,'("be%v5%rf_alpha:"/(10e12.5))') be % v5 % rf_alpha(:)

       if ( cloud_cv_options >= 2 ) then
          write(be_print_unit,'("be%v6%rf_alpha:"/(10e12.5))')  be % v6 % rf_alpha(:)
          write(be_print_unit,'("be%v7%rf_alpha:"/(10e12.5))')  be % v7 % rf_alpha(:)
          write(be_print_unit,'("be%v8%rf_alpha:"/(10e12.5))')  be % v8 % rf_alpha(:)
          write(be_print_unit,'("be%v9%rf_alpha:"/(10e12.5))')  be % v9 % rf_alpha(:)
          write(be_print_unit,'("be%v10%rf_alpha:"/(10e12.5))') be % v10 % rf_alpha(:)
       end if
       if ( use_cv_w ) then
          write(be_print_unit,'("be%v11%rf_alpha:"/(10e12.5))') be % v11 % rf_alpha(:)
       end if
       if ( qnr_option >= 1 ) then                                                                  !zhengh,20221015
          write(be_print_unit,'("be%v12%rf_alpha:"/(10e12.5))') be % v12 % rf_alpha(:)
       end if

       write(be_print_unit,'(/"scale-length: kz=",i3)') kz
       if ( cloud_cv_options >= 2 ) then
          do i = 1,kz
             if ( qnr_option >= 1 ) then  ! zhengh 20221015
               if (i == 1 ) then
                  write(be_print_unit,'(i3,2x,10e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                     be3_rf_lengthscale(i), be4_rf_lengthscale(i), be5_rf_lengthscale(i), be6_rf_lengthscale(i), &
                     be7_rf_lengthscale(i), be8_rf_lengthscale(i), be9_rf_lengthscale(i), be10_rf_lengthscale(i), be12_rf_lengthscale(i)
               else
                  write(be_print_unit,'(i3,2x,4e15.5,15x,5e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                     be3_rf_lengthscale(i), be4_rf_lengthscale(i), be6_rf_lengthscale(i), be7_rf_lengthscale(i), &
                     be8_rf_lengthscale(i), be9_rf_lengthscale(i), be10_rf_lengthscale(i), be12_rf_lengthscale(i)
               endif

             else
               if (i == 1 ) then
                  write(be_print_unit,'(i3,2x,10e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                     be3_rf_lengthscale(i), be4_rf_lengthscale(i), be5_rf_lengthscale(i), be6_rf_lengthscale(i), &
                     be7_rf_lengthscale(i), be8_rf_lengthscale(i), be9_rf_lengthscale(i), be10_rf_lengthscale(i)
               else
                  write(be_print_unit,'(i3,2x,4e15.5,15x,5e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                     be3_rf_lengthscale(i), be4_rf_lengthscale(i), be6_rf_lengthscale(i), be7_rf_lengthscale(i), &
                     be8_rf_lengthscale(i), be9_rf_lengthscale(i), be10_rf_lengthscale(i)
               endif
             end if
          enddo
       else
          do i = 1,kz
             if (i == 1) then
                write(be_print_unit,'(i3,2x,5e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                   be3_rf_lengthscale(i), be4_rf_lengthscale(i), be5_rf_lengthscale(i)
             else
                write(be_print_unit,'(i3,2x,4e15.5)') i, be1_rf_lengthscale(i), be2_rf_lengthscale(i), &
                   be3_rf_lengthscale(i), be4_rf_lengthscale(i)
             endif
          enddo
       end if

   endif ! print_detail_be

   ! 8.0 deallocate input model state:
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   deallocate (be1_eval_loc)
   deallocate (be2_eval_loc)
   deallocate (be3_eval_loc)
   deallocate (be4_eval_loc)
   deallocate (be5_eval_loc)

   if ( cloud_cv_options >= 2 ) then
      deallocate (be6_eval_loc)
      deallocate (be7_eval_loc)
      deallocate (be8_eval_loc)
      deallocate (be9_eval_loc)
      deallocate (be10_eval_loc)
   end if
   if ( use_cv_w ) then
      deallocate (be11_eval_loc)
   end if
   if ( qnr_option >= 1 ) then  ! zhengh 20221015
      deallocate (be12_eval_loc)
   end if

   if( use_rf )then
      deallocate (be1_rf_lengthscale)
      deallocate (be2_rf_lengthscale)
      deallocate (be3_rf_lengthscale)
      deallocate (be4_rf_lengthscale)
      deallocate (be5_rf_lengthscale)

      if ( cloud_cv_options >= 2 ) then
         deallocate (be6_rf_lengthscale)
         deallocate (be7_rf_lengthscale)
         deallocate (be8_rf_lengthscale)
         deallocate (be9_rf_lengthscale)
         deallocate (be10_rf_lengthscale)
      end if
      if ( use_cv_w ) then
         deallocate (be11_rf_lengthscale)
      end if
      if ( qnr_option >= 1 ) then  ! zhengh 20221015
         deallocate (be12_rf_lengthscale)
      end if
   endif

   if (vert_corr == vert_corr_2) then
      deallocate (be1_eval_glo)
      deallocate (be2_eval_glo)
      deallocate (be3_eval_glo)
      deallocate (be4_eval_glo)
      deallocate (be5_eval_glo)

      if ( cloud_cv_options >= 2 ) then
         deallocate (be6_eval_glo)
         deallocate (be7_eval_glo)
         deallocate (be8_eval_glo)
         deallocate (be9_eval_glo)
         deallocate (be10_eval_glo)
      end if
      if ( use_cv_w ) then
        deallocate (be11_eval_glo)
      end if

      if ( qnr_option >= 1 ) then  ! zhengh 20221015
         deallocate (be12_eval_glo)
      end if

      deallocate (be1_evec_loc)
      deallocate (be2_evec_loc)
      deallocate (be3_evec_loc)
      deallocate (be4_evec_loc)
      deallocate (be5_evec_loc)

      if ( cloud_cv_options >= 2 ) then
         deallocate (be6_evec_loc)
         deallocate (be7_evec_loc)
         deallocate (be8_evec_loc)
         deallocate (be9_evec_loc)
         deallocate (be10_evec_loc)
      end if
      if ( use_cv_w ) then
         deallocate (be11_evec_loc)
      end if
      if ( qnr_option >= 1 ) then   ! zhengh 20221015
         deallocate (be12_evec_loc)
      end if

      deallocate (be1_evec_glo)
      deallocate (be2_evec_glo)
      deallocate (be3_evec_glo)
      deallocate (be4_evec_glo)
      deallocate (be5_evec_glo)

      if ( cloud_cv_options >= 2 ) then
         deallocate (be6_evec_glo)
         deallocate (be7_evec_glo)
         deallocate (be8_evec_glo)
         deallocate (be9_evec_glo)
         deallocate (be10_evec_glo)
      end if
      if ( use_cv_w ) then
         deallocate (be11_evec_glo)
      end if
      if ( qnr_option >= 1 ) then   ! zhengh 20221015
         deallocate (be12_evec_glo)
      end if

   end if

   deallocate (bin)
   deallocate (bin2d)

   if (be % ne > 0) then
      be % alpha % name = 'alpha'
      allocate (alpha_val(1:nk)) ! Not using jy dimension yet.
      allocate (alpha_evec(1:nk,1:nk)) ! Not using jy dimension yet.

      if ( alpha_vertloc_opt > 0 ) then ! Use vertical localization:
         if ( rootproc ) then
            call da_get_alpha_vertloc(xb, alpha_val, alpha_evec)
         end if
         call wrf_dm_bcast_real(alpha_val,  nk)
         call wrf_dm_bcast_real(alpha_evec, nk*nk)
         be % alpha % mz = nk
         call da_get_vertical_truncation(max_vert_var_alpha, alpha_val, be % alpha)
      else
         be % alpha % mz = 1 ! No vertical localization.
         alpha_val(1) = 1.0
         alpha_val(2:kz) = 0.0
         alpha_evec(:,:) = 1.0
      end if
      mz = be % alpha % mz

!     Alpha eigenvalues and eigenvectors:
      allocate (be % alpha % val(1:jy,1:mz)) ! Not using jy dimension but here for consistency.
      allocate (be % alpha % evec(1:nj,1:nk,1:mz))

      if ( anal_type_hybrid_dual_res ) then
         deallocate(be % alpha % val, be % alpha % evec)
         allocate (be % alpha % val(jds_int:jde_int,1:mz))
         allocate (be % alpha % evec(jds_int:jde_int,1:nk,1:mz))
      endif


      do m = 1, mz
         be % alpha % val(:,m) = sigma_alpha * alpha_val(m)
         do k = 1, nk
            be % alpha % evec(:,k,m) = alpha_evec(k,m)
         end do
      end do

!     Alpha RF lengthscales and variance scaling factors:
      allocate (alpha_rf_lengthscale(1:mz))
      allocate (be % alpha % rf_alpha(1:mz))
      allocate (alpha_rf_scale_factor(1:mz))


      if ( anal_type_hybrid_dual_res ) then
         alpha_rf_lengthscale(1:mz) = 1000.0 * alpha_corr_scale / (grid%intermediate_grid % dx )
      else
         alpha_rf_lengthscale(1:mz) = 1000.0 * alpha_corr_scale / xb % ds ! Convert km to grid spacings.
      endif


      call da_calculate_rf_factors( alpha_rf_lengthscale(:), be % alpha % rf_alpha(:), &
                                    alpha_rf_scale_factor(:) )
      do m = 1, mz
         be % alpha % val(:,m) = alpha_rf_scale_factor(m) * be % alpha % val(:,m)
      end do

      if( .not. use_rf ) then
         allocate(be%alpha%wsd(nij(0,1,2),nij(0,0,2),mz))
         call random_number(be%alpha%wsd)			! need to parallelize
         if( do_normalize )then
            allocate(be%alpha%sd(nij(0,1,0),nij(0,0,0),mz))
            call random_number(be%alpha%sd)			! need to parallelize
         endif
      endif

      deallocate(alpha_val)
      deallocate(alpha_evec)
      deallocate(alpha_rf_lengthscale)
      deallocate(alpha_rf_scale_factor)

   else
      be % alpha % mz = 0
   end if

   be%cv_mz(1) = be%v1%mz
   be%cv_mz(2) = be%v2%mz
   be%cv_mz(3) = be%v3%mz
   be%cv_mz(4) = be%v4%mz
   be%cv_mz(5) = be%v5%mz
   if ( cloud_cv_options >= 2 ) then
      be%cv_mz(6)  = be%v6%mz
      be%cv_mz(7)  = be%v7%mz
      be%cv_mz(8)  = be%v8%mz
      be%cv_mz(9)  = be%v9%mz
      be%cv_mz(10) = be%v10%mz
   end if

! use_cv_w: false + qnr_option>1  ->  num_cv_3d_basic:4  num_cv_3d_extra:5  num_cv_2d:1   num_cv_qn:1
! use_cv_w: true  + qnr_option>1  ->  num_cv_3d_basic:4  num_cv_3d_extra:6  num_cv_2d:1   num_cv_qn:1

   if ( use_cv_w ) then  !false
      be%cv_mz(num_cv_3d_basic+num_cv_3d_extra+num_cv_2d) = be%v11%mz
   end if

   ! be%ncv_mz = 13
   write(unit=stdout,fmt='(A,I3)') 'be%ncv_mz = ', be%ncv_mz  !zhengh  output=13 !!!!  define in da_setup_background_errors.inc
   ! mz: Vertical truncation of errors.
   if ( qnr_option >= 1 ) then   ! zhengh 20221015
      ! be%cv_mz(12) = be%v12%mz
      be%cv_mz(num_cv_3d_basic+num_cv_3d_extra+num_cv_2d+num_cv_qn) = be%v12%mz   ! zhengh 20221017  Avoid conflicts with < use_cv_w >
   end if

   ! 但是他把第十二个和第十三个赋值给了 alpha和 ne
   be%cv_mz(be%ncv_mz-1) = be%alpha%mz    ! 12  hybrid
   be%cv_mz(be%ncv_mz)   = be%ne          
   write(unit=stdout,fmt='(A,I3)') 'be%cv_mz = ', be%cv_mz  !zhengh  output=13 !!!!  define in da_setup_background_errors.inc
   if (trace_use) call da_trace_exit("da_setup_be_regional")

end subroutine da_setup_be_regional


