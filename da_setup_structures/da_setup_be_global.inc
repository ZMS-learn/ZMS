subroutine da_setup_be_global (be)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   ! 
   ! Updates: 
   !
   !       Implementation of multi-variate BE  
   !       Syed RH Rizvi,  MMM/NESL/NCAR,  Date: 02/01/2010
   !---------------------------------------------------------------------------

   implicit none

   type (be_type), intent(out) :: be                 ! Back. errors structure

   integer                     :: nrec
   integer                     :: max_wave_in        ! Dimension of input power
   integer                     :: i, j, k, n ! Loop counters
   ! real, allocatable   :: height(:,:,:)      ! Height field.
   integer, allocatable:: bin(:,:,:)      ! Bin assigned to each 3D point
   integer, allocatable:: bin2d(:,:)      ! Bin assigned to each 2D point
   integer             :: bin_type        ! Type of bin to average over.
   integer             :: num_bins        ! Number of bins (3D fields).
   integer             :: num_bins2d      ! Number of bins (3D fields).
   logical             :: dummy

   real                :: binwidth_lat       ! Used if bin_type = 2 (degrees) 
   real                :: binwidth_hgt       ! Used if bin_type = 2 (m)
   real                :: hgt_min, hgt_max   ! Used if bin_type = 2 (m)
   real                :: lat_min, lat_max   ! Used if bin_type = 2 (degrees)

   character*10         :: variable
   integer              :: ni, nj, nk, b, be_unit
   integer, allocatable :: max_wave(:)
   real*8, allocatable  :: evec_g(:,:)
   real*8, allocatable  :: eval_g(:)  
   real*8, allocatable  :: evec_loc(:,:,:)
   real*8, allocatable  :: eval_loc(:,:)  
   real, allocatable   :: regcoeff_psi_chi(:)        ! psi/chi    regression cooefficient.
   real, allocatable   :: regcoeff_psi_t(:,:,:)      ! psi/t   regression cooefficient.
   real, allocatable   :: regcoeff_psi_ps(:,:)       ! psi/ps     regression cooefficient.
   real, allocatable   :: regcoeff_psi_rh(:,:,:)     ! psi/rh     regression cooefficient.
   real, allocatable   :: regcoeff_chi_u_t(:,:,:)    ! chi_u/t regression coefficient
   real, allocatable   :: regcoeff_chi_u_ps(:,:)     ! chi_u/ps   regression coefficient
   real, allocatable   :: regcoeff_chi_u_rh(:,:,:)   ! chi_u/rh   regression coefficient
   real, allocatable   :: regcoeff_t_u_rh(:,:,:)     ! t_u/rh  regression coefficient
   real, allocatable   :: regcoeff_ps_u_rh(:,:)      ! ps_u/rh    regression coefficient

   !Set other regcoeffs, Zheng, 20200512
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

   real*8, allocatable :: power(:)               ! Temporary power spectrum.
   real, allocatable    :: power2d(:,:)           ! Temporary power spectrum.

   if (trace_use) call da_trace_entry("da_setup_be_global")

   call da_message((/"[3.0] Set up background errors (be) for global WRF-Var"/))

   be % max_wave = ide / 2 - 1

   !---------------------------------------------------------------------
   ! [1] Read in be metadata:
   !---------------------------------------------------------------------

   call da_get_unit(be_unit)
   open (unit=be_unit,file="be.dat",status="old",form="unformatted")

   read (be_unit, end= 99, err = 100) ni, nj, nk 
   read (be_unit, err = 100) bin_type
   read (be_unit, err = 100) lat_min, lat_max, binwidth_lat
   read (be_unit, err = 100) hgt_min, hgt_max, binwidth_hgt
   read (be_unit, err = 100) num_bins, num_bins2d

   allocate (bin(1:ni,1:nj,1:nk))
   allocate (bin2d(1:ni,1:nj))
   read (be_unit, err = 100) bin(1:ni,1:nj,1:nk)
   read (be_unit, err = 100) bin2d(1:ni,1:nj)

   if (ni /= ide .or. nj /= jde .or. nk /= kde) then
      call da_error(__FILE__,__LINE__, &
         (/"Cannot generate BE at this resolution"/))
   end if

   !---------------------------------------------------------------------
   ! [2] Read in be regression coefficients:
   !---------------------------------------------------------------------


   allocate  (regcoeff_psi_chi(1:num_bins))
   allocate  (regcoeff_psi_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_psi_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_t(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_ps(1:nk,1:num_bins2d))
   allocate  (regcoeff_chi_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_t_u_rh(1:nk,1:nk,1:num_bins2d))
   allocate  (regcoeff_ps_u_rh(1:nk,1:num_bins2d))

   if ( use_regcoeff ) then		!Zheng,20200512
      if ( cv_options == 9 ) then
         allocate (regcoeff_psi_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_psi_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_chi_u_qgraup(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_ps(1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qgraup(1:nk,1:nk,1:num_bins2d))
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
      else if ( cv_options == 7 ) then
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
         allocate (regcoeff_t_u_qcloud(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qice(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qrain(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qsnow(1:nk,1:nk,1:num_bins2d))
         allocate (regcoeff_t_u_qgraup(1:nk,1:nk,1:num_bins2d))
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
      end if
   end if

   regcoeff_psi_chi = 0.
   regcoeff_psi_t   = 0.
   regcoeff_psi_ps  = 0.
   regcoeff_psi_rh  = 0.
   regcoeff_chi_u_t = 0.
   regcoeff_chi_u_ps= 0.
   regcoeff_chi_u_rh= 0.
   regcoeff_t_u_rh  = 0.
   regcoeff_ps_u_rh = 0.

   if ( use_regcoeff ) then		!Zheng,20200512
      if ( cv_options == 9 ) then
         regcoeff_psi_qcloud = 0.
         regcoeff_psi_qice = 0.
         regcoeff_psi_qrain = 0.
         regcoeff_psi_qsnow = 0.
         regcoeff_psi_qgraup = 0.
         regcoeff_chi_u_qcloud = 0.
         regcoeff_chi_u_qice = 0.
         regcoeff_chi_u_qrain = 0.
         regcoeff_chi_u_qsnow = 0.
         regcoeff_chi_u_qgraup = 0.
         regcoeff_t_u_ps = 0.
         regcoeff_t_u_qcloud = 0.
         regcoeff_t_u_qice = 0.
         regcoeff_t_u_qrain = 0.
         regcoeff_t_u_qsnow = 0.
         regcoeff_t_u_qgraup = 0.
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
      else if ( cv_options == 7 ) then
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
         regcoeff_t_u_qcloud = 0.
         regcoeff_t_u_qice = 0.
         regcoeff_t_u_qrain = 0.
         regcoeff_t_u_qsnow = 0.
         regcoeff_t_u_qgraup = 0.
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
      end if
   end if

   read (be_unit,end= 99,  err = 100) regcoeff_psi_chi
   read (be_unit,end= 99,  err = 100) regcoeff_psi_ps 
   read (be_unit,end= 99,  err = 100) regcoeff_psi_t

   ! Fill regression coeff. array
   allocate (be%reg_psi_chi(1:jde,1:nk))
   allocate (be%reg_psi_t  (1:jde,1:nk,1:nk))
   allocate (be%reg_psi_ps (1:jde,1:nk))
   allocate (be%reg_psi_rh (1:jde,1:nk,1:nk))
   allocate (be%reg_chi_u_t(1:jde,1:nk,1:nk))
   allocate (be%reg_chi_u_ps(1:jde,1:nk))
   allocate (be%reg_chi_u_rh(1:jde,1:nk,1:nk))
   allocate (be%reg_t_u_rh (1:jde,1:nk,1:nk))
   allocate (be%reg_ps_u_rh(1:jde,1:nk))

   if ( use_regcoeff ) then		!Zheng,20200512
      if ( cv_options == 9 ) then
         allocate (be%reg_psi_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_psi_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_psi_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_psi_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_psi_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_chi_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_chi_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_chi_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_chi_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_chi_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_ps(1:jde,1:nk))
         allocate (be%reg_t_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_ps_u_qcloud(1:jde,1:nk))
         allocate (be%reg_ps_u_qice(1:jde,1:nk))
         allocate (be%reg_ps_u_qrain(1:jde,1:nk))
         allocate (be%reg_ps_u_qsnow(1:jde,1:nk))
         allocate (be%reg_ps_u_qgraup(1:jde,1:nk))
         allocate (be%reg_rh_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qrain_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qrain_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qsnow_u_qgraup(1:jde,1:nk,1:nk))
      else if ( cv_options == 7 ) then
         allocate (be%reg_u_v(1:jde,1:nk,1:nk))
         allocate (be%reg_u_t(1:jde,1:nk,1:nk))
         allocate (be%reg_u_ps(1:jde,1:nk))
         allocate (be%reg_u_rh(1:jde,1:nk,1:nk))
         allocate (be%reg_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_t(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_ps(1:jde,1:nk))
         allocate (be%reg_v_u_rh(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_v_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_ps(1:jde,1:nk))
         allocate (be%reg_t_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_t_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_ps_u_qcloud(1:jde,1:nk))
         allocate (be%reg_ps_u_qice(1:jde,1:nk))
         allocate (be%reg_ps_u_qrain(1:jde,1:nk))
         allocate (be%reg_ps_u_qsnow(1:jde,1:nk))
         allocate (be%reg_ps_u_qgraup(1:jde,1:nk))
         allocate (be%reg_rh_u_qcloud(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_rh_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qice(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qcloud_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qrain(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qice_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qrain_u_qsnow(1:jde,1:nk,1:nk))
         allocate (be%reg_qrain_u_qgraup(1:jde,1:nk,1:nk))
         allocate (be%reg_qsnow_u_qgraup(1:jde,1:nk,1:nk))
      end if
   end if


   be%reg_psi_chi = 0.
   be%reg_psi_t   = 0.
   be%reg_psi_ps  = 0.
   be%reg_psi_rh  = 0.
   be%reg_chi_u_t = 0.
   be%reg_chi_u_ps= 0.
   be%reg_chi_u_rh= 0.
   be%reg_t_u_rh  = 0.
   be%reg_ps_u_rh  = 0.			!Zheng,20200513

   if ( use_regcoeff ) then		!Zheng,20200512
      if ( cv_options == 9 ) then
         be%reg_psi_qcloud = 0.
         be%reg_psi_qice = 0.
         be%reg_psi_qrain = 0.
         be%reg_psi_qsnow = 0.
         be%reg_psi_qgraup = 0.
         be%reg_chi_u_qcloud = 0.
         be%reg_chi_u_qice = 0.
         be%reg_chi_u_qrain = 0.
         be%reg_chi_u_qsnow = 0.
         be%reg_chi_u_qgraup = 0.
         be%reg_t_u_ps = 0.
         be%reg_t_u_qcloud = 0.
         be%reg_t_u_qice = 0.
         be%reg_t_u_qrain = 0.
         be%reg_t_u_qsnow = 0.
         be%reg_t_u_qgraup = 0.
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
      else if ( cv_options == 7 ) then
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
         be%reg_t_u_qcloud = 0.
         be%reg_t_u_qice = 0.
         be%reg_t_u_qrain = 0.
         be%reg_t_u_qsnow = 0.
         be%reg_t_u_qgraup = 0.
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
      end if
   end if


   do k=1,nk
      do j =1, jde
         b = bin(1,j,k)
         be%reg_psi_chi(j,k) = psi_chi_factor * regcoeff_psi_chi(b)
      end do
   end do

   do j=1,jde
      b = bin2d(1,j)
      do k=1,nk
         be%reg_psi_ps(j,k)   = psi_ps_factor   * regcoeff_psi_ps(k,b)
         be%reg_ps_u_rh(j,k)  = ps_u_rh_factor  * regcoeff_ps_u_rh(k,b)
         be%reg_chi_u_ps(j,k) = chi_u_ps_factor * regcoeff_chi_u_ps(k,b)
      end do
   end do

   do j=1,jde
      b = bin2d(1,j)
      do i=1,nk
         do k=1,nk
            be%reg_psi_t(j,i,k)   = psi_t_factor      *  regcoeff_psi_t(i,k,b)
            be%reg_psi_rh(j,i,k)  = psi_rh_factor     *  regcoeff_psi_rh(i,k,b)
            be%reg_chi_u_t(j,i,k) = chi_u_t_factor    *  regcoeff_chi_u_t(i,k,b)
            be%reg_chi_u_rh(j,i,k)= chi_u_rh_factor   *  regcoeff_chi_u_rh(i,k,b)
            be%reg_t_u_rh(j,i,k)  = t_u_rh_factor     *  regcoeff_t_u_rh(i,k,b)
         end do
      end do
   end do

   if ( use_regcoeff ) then					!Zheng,20200513
      if ( cv_options == 9 ) then				!Zheng,20200513
         do j=1,jde
            b = bin2d(1,j)
            do k=1,nk
               be%reg_t_u_ps(j,k)          = t_u_ps_factor            * regcoeff_t_u_ps(k,b)
               be%reg_ps_u_qcloud(j,k)     = ps_u_qcloud_factor       * regcoeff_ps_u_qcloud(k,b)
               be%reg_ps_u_qice(j,k)       = ps_u_qice_factor         * regcoeff_ps_u_qice(k,b)
               be%reg_ps_u_qrain(j,k)      = ps_u_qrain_factor        * regcoeff_ps_u_qrain(k,b)
               be%reg_ps_u_qsnow(j,k)      = ps_u_qsnow_factor        * regcoeff_ps_u_qsnow(k,b)
               be%reg_ps_u_qgraup(j,k)     = ps_u_qgraup_factor       * regcoeff_ps_u_qgraup(k,b)
            end do
         end do
         do j=1,jde
            b = bin2d(1,j)
            do i=1,nk
               do k=1,nk
                  be%reg_psi_qcloud(j,i,k)        = psi_qcloud_factor        *  regcoeff_psi_qcloud(i,k,b)
                  be%reg_psi_qice(j,i,k)          = psi_qice_factor          *  regcoeff_psi_qice(i,k,b)
                  be%reg_psi_qrain(j,i,k)         = psi_qrain_factor         *  regcoeff_psi_qrain(i,k,b)
                  be%reg_psi_qsnow(j,i,k)         = psi_qsnow_factor         *  regcoeff_psi_qsnow(i,k,b)
                  be%reg_psi_qgraup(j,i,k)        = psi_qgraup_factor        *  regcoeff_psi_qgraup(i,k,b)
                  be%reg_chi_u_qcloud(j,i,k)      = chi_u_qcloud_factor      *  regcoeff_chi_u_qcloud(i,k,b)
                  be%reg_chi_u_qice(j,i,k)        = chi_u_qice_factor        *  regcoeff_chi_u_qice(i,k,b)
                  be%reg_chi_u_qrain(j,i,k)       = chi_u_qrain_factor       *  regcoeff_chi_u_qrain(i,k,b)
                  be%reg_chi_u_qsnow(j,i,k)       = chi_u_qsnow_factor       *  regcoeff_chi_u_qsnow(i,k,b)
                  be%reg_chi_u_qgraup(j,i,k)      = chi_u_qgraup_factor      *  regcoeff_chi_u_qgraup(i,k,b)
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
               end do
            end do
         end do
      else if ( cv_options == 7 ) then				!Zheng,20200513
         do j=1,jde
            b = bin2d(1,j)
            do k=1,nk
               be%reg_u_ps(j,k)            = u_ps_factor              * regcoeff_u_ps(k,b)
               be%reg_v_u_ps(j,k)          = v_u_ps_factor            * regcoeff_v_u_ps(k,b)
               be%reg_t_u_ps(j,k)          = t_u_ps_factor            * regcoeff_t_u_ps(k,b)
               be%reg_ps_u_qcloud(j,k)     = ps_u_qcloud_factor       * regcoeff_ps_u_qcloud(k,b)
               be%reg_ps_u_qice(j,k)       = ps_u_qice_factor         * regcoeff_ps_u_qice(k,b)
               be%reg_ps_u_qrain(j,k)      = ps_u_qrain_factor        * regcoeff_ps_u_qrain(k,b)
               be%reg_ps_u_qsnow(j,k)      = ps_u_qsnow_factor        * regcoeff_ps_u_qsnow(k,b)
               be%reg_ps_u_qgraup(j,k)     = ps_u_qgraup_factor       * regcoeff_ps_u_qgraup(k,b)
            end do
         end do
         do j=1,jde
            b = bin2d(1,j)
            do i=1,nk
               do k=1,nk
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
               end do
            end do
         end do
      end if
   end if

   deallocate (regcoeff_psi_chi)
   deallocate (regcoeff_psi_t)
   deallocate (regcoeff_psi_ps)
   deallocate (regcoeff_psi_rh)
   deallocate (regcoeff_chi_u_t)
   deallocate (regcoeff_chi_u_ps)
   deallocate (regcoeff_chi_u_rh)
   deallocate (regcoeff_t_u_rh)
   deallocate (regcoeff_ps_u_rh)

   if ( use_regcoeff ) then		!Zheng,20200512
      if ( cv_options == 9 ) then
         deallocate (regcoeff_psi_qcloud)
         deallocate (regcoeff_psi_qice)
         deallocate (regcoeff_psi_qrain)
         deallocate (regcoeff_psi_qsnow)
         deallocate (regcoeff_psi_qgraup)
         deallocate (regcoeff_chi_u_qcloud)
         deallocate (regcoeff_chi_u_qice)
         deallocate (regcoeff_chi_u_qrain)
         deallocate (regcoeff_chi_u_qsnow)
         deallocate (regcoeff_chi_u_qgraup)
         deallocate (regcoeff_t_u_ps)
         deallocate (regcoeff_t_u_qcloud)
         deallocate (regcoeff_t_u_qice)
         deallocate (regcoeff_t_u_qrain)
         deallocate (regcoeff_t_u_qsnow)
         deallocate (regcoeff_t_u_qgraup)
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
      else if ( cv_options == 7 ) then
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
         deallocate (regcoeff_t_u_qcloud)
         deallocate (regcoeff_t_u_qice)
         deallocate (regcoeff_t_u_qrain)
         deallocate (regcoeff_t_u_qsnow)
         deallocate (regcoeff_t_u_qgraup)
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
      end if
   end if

   !---------------------------------------------------------------------
   ! [3] Read in be vertical eigenmodes:
   !---------------------------------------------------------------------

   do nrec = 1, 4
      read (be_unit,end= 99,  err = 100) variable   
      read (be_unit,end= 99,  err = 100) nk, num_bins2d 

      allocate (evec_g(1:nk,1:nk))
      allocate (eval_g(1:nk))
      allocate (evec_loc(1:nk,1:nk,num_bins2d))
      allocate (eval_loc(1:nk,num_bins2d))

      read (be_unit,end= 99, err = 100) evec_g     
      read (be_unit,end= 99, err = 100) eval_g     
      read (be_unit,end= 99, err = 100) evec_loc     
      read (be_unit,end= 99, err = 100) eval_loc    

      if (nrec == 1) then
         be % v1 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var1, var_scaling1(1), be%v1)

      else if (nrec == 2) then
         be % v2 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var2, var_scaling2(1), be%v2)

      else if (nrec == 3) then
         be % v3 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var3, var_scaling3(1), be%v3)

      else if (nrec == 4) then
         be % v4 % name = variable               
         call da_get_bins_info(nj, nk, bin2d, evec_g, eval_g, &
            evec_loc, eval_loc, max_vert_var4, var_scaling4(1), be%v4)
      end if 

      deallocate (evec_g)     
      deallocate (eval_g)     
      deallocate (evec_loc)     
      deallocate (eval_loc)     

   end do ! loop nrec

   deallocate (bin)
   deallocate (bin2d)

   !---------------------------------------------------------------------
   ! [4] Read in be power spectra:
   !---------------------------------------------------------------------

   allocate(max_wave(1:nk))

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format 
      if (k == 1) then
          allocate (power(0:max_wave_in))                      ! Temporary.
          allocate (power2d(0:max_wave_in,1:nk))               ! Temporary.
      end if
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra(be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v1 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(/3x,3a,i6)') &
      'Horizontal truncation for ', be % v1 % name, ' = ', be % v1 % max_wave
   allocate (be % v1 % power(0:be % v1 % max_wave,1:nk))
   be % v1 % power(0:be % v1 % max_wave,1:nk) = power2d(0:be % v1 % max_wave,1:nk)
   be % v1 % power(0,1:nk) = len_scaling1(1) * be % v1 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format    
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v2 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v2 % name, ' = ', be % v2 % max_wave
   allocate (be % v2 % power(0:be % v2 % max_wave,1:nk))
   be % v2 % power(0:be % v2 % max_wave,1:nk) = power2d(0:be % v2 % max_wave,1:nk)
   be % v2 % power(0,1:nk) = len_scaling2(1) * be % v2 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format    
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v3 % max_wave = maxval(max_wave(1:nk))
   write(unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v3 % name, ' = ', be % v3 % max_wave
   allocate (be % v3 % power(0:be % v3 % max_wave,1:nk))
   be % v3 % power(0:be % v3 % max_wave,1:nk) = power2d(0:be % v3 % max_wave,1:nk)
   be % v3 % power(0,1:nk) = len_scaling3(1) * be % v3 % power(0,1:nk) 

   do k = 1, nk
      read (be_unit) variable
      read (be_unit) max_wave_in, nrec
      read (be_unit) dummy ! use to preserve file format
      read (be_unit) power(0:max_wave_in) 
      power2d(:,k) = power(:) 

      ! Truncate power spectra:
      call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                                power, max_wave(k))
   end do

   be % v4 % max_wave = maxval(max_wave(1:nk))
   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v4 % name, ' = ', be % v4 % max_wave
   allocate (be % v4 % power(0:be % v4 % max_wave,1:nk))
   be % v4 % power(0:be % v4 % max_wave,1:nk) = power2d(0:be % v4 % max_wave,1:nk)
   be % v4 % power(0,1:nk) = len_scaling4(1) * be % v4 % power(0,1:nk) 

   ! ps_u:
   read (be_unit) variable
   be % v5 % name = variable
   be % v5 % mz = 1
   if (max_vert_var5 <=  0.0) be % v5 % mz = 0                         
   read (be_unit) max_wave_in, nrec
   read (be_unit) dummy ! use to preserve file format
   read (be_unit) power(0:max_wave_in) 

   ! Truncate power spectra:
   call da_truncate_spectra (be % max_wave, max_wave_in, power_truncation, &
                             power, be % v5 % max_wave)

   write (unit=stdout,fmt='(3x,3a,i6)') &
      'Horizontal truncation for ', be % v5 % name, ' = ', be % v5 % max_wave
   allocate (be % v5 % power(0:be % v5 % max_wave,1))
   be % v5 % power(0:be % v5 % max_wave,1) = power(0:be % v5 % max_wave)
   be % v5 % power(0,1) = len_scaling5(1) * be%v5%power(0,1) 

   deallocate(power)

   !--------------------------------------------------------------
   ! [5] Perform checks on eigenvectors:
   !--------------------------------------------------------------

   if (test_statistics) then
      call da_check_eof_decomposition(be%v1%val_g(:), be%v1%evec_g(:,:),&
                                     be%v1%name)
      call da_check_eof_decomposition(be%v2%val_g(:), be%v2%evec_g(:,:),&
                                     be%v2%name)
      call da_check_eof_decomposition(be%v3%val_g(:), be%v3%evec_g(:,:),&
                                     be%v3%name)
      call da_check_eof_decomposition(be%v4%val_g(:), be%v4%evec_g(:,:),&
                                     be%v4%name)
   end if
 
   !--------------------------------------------------------------
   ! [6] Set up alpha (flow-dependent) control variable "errors": 
   !--------------------------------------------------------------

   if (be % ne > 0) then
      be % alpha % mz = be % ne
      be % alpha % name = 'alpha'
      be % alpha % max_wave = alpha_truncation
      if (print_detail_be) then
         write(unit=stdout,fmt='(3x,3a,i6)') &
            'Horizontal truncation for ', be % alpha % name, ' = ', &
            be % alpha % max_wave
      end if
      allocate (power(0:be % alpha % max_wave))
      call da_calc_power_spectrum(be % alpha % max_wave, power)

      allocate (be % alpha % power(0:be % alpha % max_wave, be % ne))
      do n = 1, be % ne
         be % alpha % power(0:be % alpha % max_wave,n) = power(0:be % alpha % max_wave)
      end do
   end if

   deallocate (max_wave)

   write (unit=stdout,fmt='(A)') " "

   close (be_unit)
   call da_free_unit(be_unit)

   if (trace_use) call da_trace_exit("da_setup_be_global")

   return

99 write (unit=message(1),fmt='(a, i5)')' Unexpected end on BE-unit = ',be_unit
   call da_error(__FILE__,__LINE__,message(1:1))

100 write (unit=message(1),fmt='(a, i5)')' Read error on BE-unit = ',be_unit
   call da_error(__FILE__,__LINE__,message(1:1))

   ! FIX? daft having these here

   deallocate(power)
   deallocate(power2d)

end subroutine da_setup_be_global


