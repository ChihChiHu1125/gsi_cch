module radiance_mod
!$$$   module documentation block
!                .      .    .                                       .
! module:    radiance_mod
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  This module contains variables and routines related
!            to cloud and aerosol usages for radiance assimilation
!
! program history log:
!   2015-07-20 Yanqiu Zhu
!   2016-10-27 Yanqiu - add ATMS
!
! subroutines included:
!   sub radiance_mode_init           -  guess init
!   radiance_mode_destroy
!   radiance_obstype_init
!   radiance_obstype_search
!   radiance_obstype_destroy
!   radiance_parameter_cloudy_init
!   radiance_parameter_aerosol_init
!   radiance_ex_obserr
!   radiance_ex_biascor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

! !USES:

  use kinds, only: r_kind,i_kind
  use constants, only: zero,half
  use mpimod, only: mype
  implicit none
  save

! set subroutines to public
  public :: radiance_mode_init
  public :: radiance_mode_destroy
  public :: radiance_obstype_init
  public :: radiance_obstype_search
  public :: radiance_obstype_destroy
  public :: radiance_parameter_cloudy_init
  public :: radiance_parameter_aerosol_init
  public :: radiance_ex_obserr
  public :: radiance_ex_obserr_gmi

! CCH:
  public :: radiance_ex_obserr_evolve_gauss
  public :: radiance_ex_obserr_evolve_gauss_interp
  public :: read_nongauss_table 
  public :: interpolate

  public :: radiance_ex_biascor
  public :: radiance_ex_biascor_gmi

  public :: icloud_fwd,icloud_cv,iallsky,cw_cv,ql_cv
  public :: n_actual_clouds,n_clouds_fwd,n_clouds_jac
  public :: cloud_names,cloud_names_jac,cloud_names_fwd
  public :: idx_cw,idx_ql,idx_qi,idx_qr,idx_qs,idx_qg,idx_qh

  public :: iaerosol_fwd,iaerosol_cv,iaerosol
  public :: n_actual_aerosols,n_aerosols_fwd,n_aerosols_jac
  public :: aerosol_names,aerosol_names_fwd,aerosol_names_jac

  public :: total_rad_type
  public :: rad_type_info
  public :: total_aod_type
  public :: aod_type_info

  public :: rad_obs_type

  interface radiance_ex_obserr
     module procedure radiance_ex_obserr_1
     module procedure radiance_ex_obserr_2
!     module procedure radiance_ex_obserr_3 !for GMI
  end interface 
 
  interface radiance_ex_biascor
     module procedure radiance_ex_biascor_1
     module procedure radiance_ex_biascor_2
!     module procedure radiance_ex_biascor_3 ! for GMI
  end interface

  character(len=20),save,allocatable,dimension(:) :: cloud_names
  character(len=20),save,allocatable,dimension(:) :: cloud_names_fwd
  character(len=20),save,allocatable,dimension(:) :: cloud_names_jac
  character(len=20),save,allocatable,dimension(:) :: aerosol_names
  character(len=20),save,allocatable,dimension(:) :: aerosol_names_fwd
  character(len=20),save,allocatable,dimension(:) :: aerosol_names_jac
  logical :: icloud_fwd,icloud_cv,iallsky,cw_cv,ql_cv
  logical :: iaerosol_fwd,iaerosol_cv,iaerosol,iprecip
  integer(i_kind) :: n_actual_clouds,n_clouds_jac,n_clouds_fwd
  integer(i_kind) :: n_actual_aerosols,n_aerosols_fwd,n_aerosols_jac
  integer(i_kind) :: idx_cw,idx_ql,idx_qi,idx_qr,idx_qs,idx_qg,idx_qh

  integer(i_kind) :: total_rad_type
  integer(i_kind) :: total_aod_type

  type rad_obs_type
    character(len=10) :: rtype            ! instrument
    integer(i_kind) :: nchannel       ! total channel number
!   character(len=8) :: cfoption          ! cloud fraction option: gmao_lcf4crtm, emc_lcf4crtm 
    character(len=10) :: ex_obserr        ! indicator for special obs error assignment: ex_obserr1 or ex_obserr2
    logical :: cld_sea_only           ! .true. only perform all-sky over ocean
    logical :: ex_biascor             ! .true. for special bias correction
    logical :: cld_effect             ! .true. additional cloud effect quality control
    logical :: lcloud_fwd,lallsky
    logical :: lprecip 
    integer(i_kind),pointer,dimension(:) :: lcloud4crtm=> NULL()    ! -1 clear-sky; 0 forwad operator only; 1 iallsky
    logical :: laerosol_fwd,laerosol
    integer(i_kind),pointer,dimension(:) :: laerosol4crtm => NULL() ! -1 no aero used; 0 forwad operator only; 1 iaerosol 
    real(r_kind),pointer,dimension(:) :: cclr    => NULL()
    real(r_kind),pointer,dimension(:) :: ccld    => NULL()
    real(r_kind),pointer,dimension(:) :: cldval1 => NULL()
  end type rad_obs_type

  type(rad_obs_type),save,dimension(:),allocatable :: rad_type_info
  type(rad_obs_type),save,dimension(:),allocatable :: aod_type_info

contains

  subroutine radiance_mode_init
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_mode_init
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  This routine sets default values for variables used in
!            the radiance processing routines.
!
! program history log:
!   2015-07-20  zhu     
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use gsi_metguess_mod, only: gsi_metguess_get
    use gsi_chemguess_mod, only: gsi_chemguess_get
    use mpeu_util, only: getindex
    use control_vectors, only: cvars3d
    implicit none

    integer(i_kind) icw_av,iql_av,iqi_av,iqtotal,ier
    integer(i_kind) iqr,iqs,iqg 
    integer(i_kind) indx_p25,indx_dust1,indx_dust2,ip25_av,idust1_av,idust2_av

!   initialize variables
    icloud_fwd=.false.
    icloud_cv=.false.
    iallsky=.false.
    cw_cv=.false.
    ql_cv=.false.

    n_actual_clouds=0
    n_clouds_fwd=0 
    n_clouds_jac=0

    iaerosol_fwd=.false.
    iaerosol_cv=.false.
    iaerosol=.false.

    iprecip=.false.

    n_actual_aerosols=0
    n_aerosols_fwd=0
    n_aerosols_jac=0

!   inquire number of clouds 
    call gsi_metguess_get ( 'clouds::3d', n_actual_clouds, ier )
    if (n_actual_clouds>0) then
       allocate(cloud_names(n_actual_clouds))
       call gsi_metguess_get ('clouds::3d', cloud_names, ier)
       call gsi_metguess_get ('clouds_4crtm_fwd::3d', n_clouds_fwd, ier)
       n_clouds_fwd=max(0,n_clouds_fwd)
       if (n_clouds_fwd>0) then
          icloud_fwd=.true.
          allocate(cloud_names_fwd(max(n_clouds_fwd,1)))
          call gsi_metguess_get ('clouds_4crtm_fwd::3d', cloud_names_fwd, ier)

          call gsi_metguess_get ('clouds_4crtm_jac::3d', n_clouds_jac, ier )
          n_clouds_jac=max(0,n_clouds_jac)
          if (n_clouds_jac>0) then
             allocate(cloud_names_jac(max(n_clouds_jac,1)))
             call gsi_metguess_get ('clouds_4crtm_jac::3d', cloud_names_jac, ier)
          end if
          iqr = getindex(cloud_names_fwd,'qr')
          iqs = getindex(cloud_names_fwd,'qs')
          iqg = getindex(cloud_names_fwd,'qg')
          if (iqr>0 .or. iqs>0 .or. iqg>0) iprecip = .true.
       end if

!      inquire number of clouds to participate in CRTM calculations
       call gsi_metguess_get ( 'i4crtm::ql', idx_ql, ier )
       call gsi_metguess_get ( 'i4crtm::qi', idx_qi, ier )
       call gsi_metguess_get ( 'i4crtm::qr', idx_qr, ier )
       call gsi_metguess_get ( 'i4crtm::qs', idx_qs, ier )
       call gsi_metguess_get ( 'i4crtm::qg', idx_qg, ier )
       call gsi_metguess_get ( 'i4crtm::qh', idx_qh, ier )
!      if (idx_ql>10 .or. idx_qi>10 .or. idx_qr>10 .or. idx_qs>10 &
!         .or. idx_qg>10 .or. idx_qh>10) icloud_fwd=.true.

!      Determine whether or not cloud-condensate is the control variable
!      (ges_cw=ges_ql+ges_qi)
       icw_av=getindex(cvars3d,'cw')
       iql_av=getindex(cvars3d,'ql')
       iqi_av=getindex(cvars3d,'qi')

!      Determine whether or not total moisture (water vapor+total cloud
!      condensate) is the control variable
       iqtotal=getindex(cvars3d,'qt')

       if (icw_av>0) cw_cv=.true.
       if (iql_av>0) ql_cv=.true.
       if (icw_av>0 .or. iql_av>0 .or. iqi_av>0 .or. iqtotal>0) icloud_cv=.true.
       if (icloud_cv .and. icloud_fwd) iallsky=.true.

    end if  ! end of (n_actual_clouds>0)


!   inquire number of aerosols
    call gsi_chemguess_get ( 'aerosols::3d', n_actual_aerosols, ier )
    if (n_actual_aerosols > 0) then
       iaerosol_fwd=.true.
       allocate(aerosol_names(n_actual_aerosols))
       call gsi_chemguess_get ('aerosols::3d',aerosol_names,ier)
       indx_p25   = getindex(aerosol_names,'p25')
       indx_dust1 = getindex(aerosol_names,'dust1')
       indx_dust2 = getindex(aerosol_names,'dust2')

       call gsi_chemguess_get ( 'aerosols_4crtm::3d', n_aerosols_fwd, ier )
       if (n_aerosols_fwd >0) then
          allocate(aerosol_names_fwd(n_aerosols_fwd))
          call gsi_chemguess_get ( 'aerosols_4crtm::3d', aerosol_names_fwd, ier)  
       end if
       call gsi_chemguess_get ( 'aerosols_4crtm_jac::3d', n_aerosols_jac, ier )
       if (n_aerosols_jac >0) then
          allocate(aerosol_names_jac(n_aerosols_jac))
          call gsi_chemguess_get ( 'aerosols_4crtm_jac::3d', aerosol_names_jac, ier)  
       end if
    endif

!   Determine whether aerosols are control variables
    ip25_av=getindex(cvars3d,'p25')
    idust1_av=getindex(cvars3d,'dust1')
    idust2_av=getindex(cvars3d,'dust2')
    if (ip25_av>0 .or. idust1_av>0 .or. idust2_av>0) iaerosol_cv=.true.

    if (iaerosol_cv .and. iaerosol_fwd) iaerosol=.true.

    if (mype==0) then
       write(6,*) 'radiance_mode_init: icloud_fwd=',icloud_fwd,' iallsky=',iallsky, &
                  ' cw_cv=',cw_cv,' ql_cv=',ql_cv,' iaerosol_fwd=',iaerosol_fwd,' iaerosol=',iaerosol, &
                  ' iprecip=',iprecip
       write(6,*) 'radiance_mode_init: n_actual_clouds=',n_actual_clouds
       if (n_actual_clouds>0) write(6,*) 'radiance_mode_init: cloud_names=',cloud_names  
       write(6,*) 'radiance_mode_init: n_clouds_fwd=',n_clouds_fwd
       if (n_clouds_fwd>0) write(6,*) 'radiance_mode_init: cloud_names_fwd=',cloud_names_fwd
       write(6,*) 'radiance_mode_init: n_clouds_jac=',n_clouds_jac
       if (n_clouds_jac>0) write(6,*) 'radiance_mode_init: cloud_names_jac=',cloud_names_jac
       write(6,*) 'radiance_mode_init: n_actual_aerosols=',n_actual_aerosols
       if (n_actual_aerosols>0) write(6,*) 'radiance_mode_init: aerosol_names=',aerosol_names
       write(6,*) 'radiance_mode_init: n_aerosols_fwd=',n_aerosols_fwd
       if (n_aerosols_fwd>0) write(6,*) 'radiance_mode_init: aerosol_names_fwd=',aerosol_names_fwd
       write(6,*) 'radiance_mode_init: n_aerosols_jac=',n_aerosols_jac
       if (n_aerosols_jac>0) write(6,*) 'radiance_mode_init: aerosol_names_jac=',aerosol_names_jac
    end if
    
  end subroutine radiance_mode_init

  subroutine radiance_mode_destroy
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_mode_destroy
!
!   prgrmmr:     yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  This routine deallocate arrays
!
! program history log:
!   2015-07-20  zhu     
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

    implicit none
 
    if(allocated(cloud_names)) deallocate(cloud_names)
    if(allocated(cloud_names_fwd)) deallocate(cloud_names_fwd)
    if(allocated(cloud_names_jac)) deallocate(cloud_names_jac)
  
    if(allocated(aerosol_names)) deallocate(aerosol_names)
    if(allocated(aerosol_names_fwd)) deallocate(aerosol_names_fwd)
    if(allocated(aerosol_names_jac)) deallocate(aerosol_names_jac)

  end subroutine radiance_mode_destroy

  subroutine radiance_obstype_init
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_obstype_init
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  This routine sets default values for variables used in
!            the cloudy/with aerosol radiance processing routines.
!
! program history log:
!   2015-07-20  zhu -- initial code    
!   2018-04-04  zhu -- move rad_type_info(k)%cclr and rad_type_info(k)%ccld to this subroutine
!   2018-04-06  derber -- change rad_type_info(k)%cclr default value from zero to a large number
!   2019-03-21  Wei/Martin - fix capabilities for AOD assimilation
!   2021-03-05  X.Li - add viirs (viirs-m with sstviirs data set)
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use radinfo, only: nusis,jpch_rad,icloud4crtm,iaerosol4crtm
    use obsmod, only: ndat,dtype,dsis
    use gsi_io, only: verbose
    use chemmod, only: laeroana_gocart, laeroana_fv3cmaq
    implicit none

    logical :: first,diffistr,found
    integer(i_kind) :: i,j,k,ii,nn1,nn2,nn
    integer(i_kind),dimension(ndat) :: k2i
    character(10),dimension(ndat) :: rtype,rrtype,drtype
    logical print_verbose

    print_verbose=.false.
    if(verbose)print_verbose=.true.
!   Cross-check 
    do j=1,jpch_rad
       if (icloud4crtm(j)>=0) then
          if (.not. iallsky) icloud4crtm(j)=0
          if (.not. icloud_fwd) icloud4crtm(j)=-1
       end if
       if (iaerosol4crtm(j)>=0) then
          if (.not. iaerosol) iaerosol4crtm(j)=0
          if (.not. iaerosol_fwd) iaerosol4crtm(j)=-1
       end if
    end do

    if (icloud_fwd .and. all(icloud4crtm<0)) then 
       icloud_fwd=.false.
       iallsky=.false.
       n_clouds_fwd=0
       n_clouds_jac=0       
       cloud_names_fwd=' '
       cloud_names_jac=' '
    end if
    if (iaerosol_fwd .and. all(iaerosol4crtm<0)) then
       iaerosol=.false.
       if ( .not. laeroana_gocart .and. .not. laeroana_fv3cmaq) then
          iaerosol_fwd=.false.
          n_aerosols_fwd=0
          n_aerosols_jac=0   
          aerosol_names_fwd=' '
          aerosol_names_jac=' '
       end if
    end if

    if (iallsky .and. all(icloud4crtm<1)) then
       iallsky=.false.
       n_clouds_jac=0
       cloud_names_jac=' '
    end if

    if (iaerosol .and. all(iaerosol4crtm<1)) then
       iaerosol=.false.
       n_aerosols_jac=0
       aerosol_names_jac=' '
    end if

!   determine rads type
    drtype='other'
    do i=1,ndat
       rtype(i)=dtype(i)                   !     rtype  - observation types to process
       if (index(dtype(i),'amsre') /= 0)  rtype(i)='amsre'
       if (index(dtype(i),'ssmis') /= 0)  rtype(i)='ssmis'
       if (index(dtype(i),'sndr') /= 0)   rtype(i)='sndr'
       if (index(dtype(i),'hirs') /= 0)   rtype(i)='hirs'
       if (index(dtype(i),'avhrr') /= 0)  rtype(i)='avhrr'
       if (index(dtype(i),'viirs-m') /= 0)  rtype(i)='viirs'
       if (index(dtype(i),'modis') /= 0)  rtype(i)='modis'
       if (index(dtype(i),'seviri') /= 0) rtype(i)='seviri'

       if(rtype(i) == 'hirs'   .or. rtype(i) == 'sndr'     .or.  rtype(i) == 'seviri' .or. &
          rtype(i) == 'airs'   .or. rtype(i) == 'amsua'    .or.  rtype(i) == 'msu'    .or. & 
          rtype(i) == 'iasi'   .or. rtype(i) == 'amsub'    .or.  rtype(i) == 'mhs'    .or. &
          rtype(i) == 'hsb'    .or. rtype(i) == 'goes_img' .or.  rtype(i) == 'ahi'    .or. &
          rtype(i) == 'avhrr'  .or. rtype(i) == 'amsre'    .or.  rtype(i) == 'ssmis'  .or. & 
          rtype(i) == 'ssmi'   .or. rtype(i) == 'atms'     .or.  rtype(i) == 'cris'   .or. & 
          rtype(i) == 'amsr2'  .or. rtype(i) == 'gmi'      .or.  rtype(i) == 'saphir' .or. &
          rtype(i) == 'cris-fsr' .or. rtype(i) == 'abi'    .or.  rtype(i) == 'viirs' ) then
          drtype(i)='rads'
       end if
    end do
 
!   Determine total rad types
    k=0
    k2i=0
    first=.true.
    rrtype=''
    do i=1,ndat
       if (drtype(i) /= 'rads') cycle
 
       found=.false.
       if (first) then
          k=k+1
          rrtype(k)=rtype(i) 
          k2i(k)=i
          first=.false.
       else
          do j=1,k
             if (trim(rtype(i)) == trim(rrtype(j))) then 
                found=.true.
                exit
             end if
          end do
          if (.not. found) then
             k=k+1
             rrtype(k)=rtype(i)
             k2i(k)=i
          end if
       end if
    end do
    total_rad_type=k
    if (mype==0) write(6,*) 'radiance_obstype_init: total_rad_type=', k,' types are: ', rrtype(1:total_rad_type)

    if (total_rad_type<=0) return

    allocate(rad_type_info(total_rad_type)) 

    do k=1, total_rad_type
       rad_type_info(k)%rtype=rrtype(k)
       rad_type_info(k)%cld_sea_only=.false.
!      rad_type_info(k)%ex_obserr=.false.
       rad_type_info(k)%ex_obserr=' '
       rad_type_info(k)%ex_biascor=.false.
       rad_type_info(k)%cld_effect=.false.
       rad_type_info(k)%lcloud_fwd=.false.
       rad_type_info(k)%lprecip=.false.
       rad_type_info(k)%lallsky=.false.
       rad_type_info(k)%laerosol_fwd=.false.
       rad_type_info(k)%laerosol=.false.

       if (iprecip) rad_type_info(k)%lprecip=.true.

       ii=k2i(k)
       first=.true.
       nn1=0
       nn2=0
       do j=1,jpch_rad
          if (j==jpch_rad) then
             diffistr = .true.
          else
             diffistr = trim(nusis(j))/=trim(nusis(j+1))
          end if
          if (trim(dsis(ii))==trim(nusis(j))) then
!         if (index(trim(nusis(j)),trim(rrtype(k))) /= 0) then
             if (first) then
                nn1=j
                first=.false.
             else
                nn2=j
             end if
             if (diffistr) exit
          end if
       end do
       if (nn1/=0 .and. nn2/=0) then
          rad_type_info(k)%nchannel=nn2-nn1+1
       else
          cycle
       end if

!      determine usages of cloud and aerosol in each type
       allocate(rad_type_info(k)%lcloud4crtm(rad_type_info(k)%nchannel)) 
       allocate(rad_type_info(k)%laerosol4crtm(rad_type_info(k)%nchannel)) 
       nn=0
       do j=nn1,nn2
          nn=nn+1
          rad_type_info(k)%lcloud4crtm(nn)=icloud4crtm(j)
          rad_type_info(k)%laerosol4crtm(nn)=iaerosol4crtm(j)

          if (icloud4crtm(j)<0 .and. iaerosol4crtm(j)<0) cycle
          if (.not. rad_type_info(k)%lallsky) then
             if (icloud4crtm(j)==1) then 
                rad_type_info(k)%lallsky=.true.
                rad_type_info(k)%lcloud_fwd=.true.
             end if
          end if
          if (.not. rad_type_info(k)%lcloud_fwd) then
             if (icloud4crtm(j)==0) rad_type_info(k)%lcloud_fwd=.true.
          end if
          if (.not. rad_type_info(k)%laerosol) then
             if (iaerosol4crtm(j)==1) then 
                rad_type_info(k)%laerosol=.true.
                rad_type_info(k)%laerosol_fwd=.true.
             end if
          end if
          if (.not. rad_type_info(k)%laerosol_fwd) then
             if (iaerosol4crtm(j)==0) rad_type_info(k)%laerosol_fwd=.true.
          end if
       end do

       if (mype==0 .and. print_verbose)  &
                               write(6,*) 'radiance_obstype_init: type=', rad_type_info(k)%rtype, &
                               ' nch=',rad_type_info(k)%nchannel, &
                               ' lcloud_fwd=',rad_type_info(k)%lcloud_fwd, &
                               ' lprecip=',rad_type_info(k)%lprecip, &    
                               ' lallsky=',rad_type_info(k)%lallsky, &
                               ' laerosol_fwd=',rad_type_info(k)%laerosol_fwd, &
                               ' laerosol=',rad_type_info(k)%laerosol

       allocate(rad_type_info(k)%cclr(rad_type_info(k)%nchannel)) 
       allocate(rad_type_info(k)%ccld(rad_type_info(k)%nchannel)) 
       allocate(rad_type_info(k)%cldval1(rad_type_info(k)%nchannel)) 
       rad_type_info(k)%cclr(:)=9999.9_r_kind
       rad_type_info(k)%ccld(:)=zero
       rad_type_info(k)%cldval1(:)=zero

    end do ! end total_rad_type

  end subroutine radiance_obstype_init

  subroutine radiance_obstype_search(obstype,radmod)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_obstype_search   find the rad_type_info(i) that
!                                          matches the input obstype
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-08-20
!
! abstract:
!
! program history log:
!   2015-08-20  zhu
!   2019-03-21  Wei/Martin - added in AOD type support
!
!   input argument list:
!         obstype
!
!   output argument list:
!         radmod
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    implicit none
    character(10) :: obstype
    type(rad_obs_type) :: radmod
    logical match
    integer(i_kind) i

    if (total_rad_type<=0 .and. total_aod_type<=0) return
    
    match=.false.
    do i=1,total_rad_type
       if (trim(rad_type_info(i)%rtype)=='msu') then
          match=trim(obstype)==trim(rad_type_info(i)%rtype)
       else
          match=index(trim(obstype),trim(rad_type_info(i)%rtype)) /= 0
       end if
       if (match) then
!         if (mype==0) write(6,*) 'radiance_obstype_search: obstype=',obstype, &
!                                 ' rtype=',rad_type_info(i)%rtype
          radmod%rtype = rad_type_info(i)%rtype
          radmod%nchannel = rad_type_info(i)%nchannel
          radmod%cld_sea_only = rad_type_info(i)%cld_sea_only
          radmod%cld_effect = rad_type_info(i)%cld_effect
          radmod%ex_obserr = rad_type_info(i)%ex_obserr
          radmod%ex_biascor = rad_type_info(i)%ex_biascor

          radmod%lcloud_fwd = rad_type_info(i)%lcloud_fwd
          radmod%lallsky = rad_type_info(i)%lallsky
          radmod%lcloud4crtm => rad_type_info(i)%lcloud4crtm

          radmod%laerosol_fwd = rad_type_info(i)%laerosol_fwd
          radmod%laerosol = rad_type_info(i)%laerosol
          radmod%laerosol4crtm => rad_type_info(i)%laerosol4crtm

          radmod%cclr => rad_type_info(i)%cclr
          radmod%ccld => rad_type_info(i)%ccld
          radmod%cldval1 => rad_type_info(i)%cldval1
          radmod%lprecip = radmod%lcloud_fwd .and. rad_type_info(i)%lprecip 

          return
       end if
    end do
    do i=1,total_aod_type
       match=index(trim(obstype),trim(aod_type_info(i)%rtype)) /= 0
       if (match) then
          radmod%rtype = aod_type_info(i)%rtype
          radmod%nchannel = aod_type_info(i)%nchannel
          radmod%cld_sea_only = aod_type_info(i)%cld_sea_only
          radmod%cld_effect = aod_type_info(i)%cld_effect
          radmod%ex_obserr = aod_type_info(i)%ex_obserr
          radmod%ex_biascor = aod_type_info(i)%ex_biascor
 
          radmod%lcloud_fwd = aod_type_info(i)%lcloud_fwd
          radmod%lallsky = aod_type_info(i)%lallsky
          radmod%lcloud4crtm => aod_type_info(i)%lcloud4crtm

          radmod%laerosol_fwd = aod_type_info(i)%laerosol_fwd
          radmod%laerosol = aod_type_info(i)%laerosol
          radmod%laerosol4crtm => aod_type_info(i)%laerosol4crtm

          radmod%cclr => aod_type_info(i)%cclr
          radmod%ccld => aod_type_info(i)%ccld
          return
       end if ! match
    end do
    if (mype==0) write(6,*) 'radiance_obstype_search type not found: obstype=',obstype

    if (.not. match) then
       if (mype==0) write(6,*) 'radiance_obstype_search: #WARNING# obstype=',obstype,' not found in rtype'
    end if

  end subroutine radiance_obstype_search


  subroutine radiance_obstype_destroy
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_obstype_destroy
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  
!
! program history log:
!   2015-07-20  zhu
!   2019-03-21 Wei/Martin - added in aod_type support
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    implicit none

    integer(i_kind) :: k
    if (total_rad_type>0) then
       do k=1, total_rad_type
          if(associated(rad_type_info(k)%lcloud4crtm)) deallocate(rad_type_info(k)%lcloud4crtm)
          if(associated(rad_type_info(k)%laerosol4crtm)) deallocate(rad_type_info(k)%laerosol4crtm)
          if(associated(rad_type_info(k)%cclr)) deallocate(rad_type_info(k)%cclr)
          if(associated(rad_type_info(k)%ccld)) deallocate(rad_type_info(k)%ccld)
       end do
    end if
    if (total_aod_type>0) then
       do k=1, total_aod_type
          if(associated(aod_type_info(k)%lcloud4crtm)) deallocate(aod_type_info(k)%lcloud4crtm)
          if(associated(aod_type_info(k)%laerosol4crtm)) deallocate(aod_type_info(k)%laerosol4crtm)
          if(associated(aod_type_info(k)%cclr)) deallocate(aod_type_info(k)%cclr)
          if(associated(aod_type_info(k)%ccld)) deallocate(aod_type_info(k)%ccld)
       end do
    end if

    do k=1, total_rad_type
       if(associated(rad_type_info(k)%lcloud4crtm)) deallocate(rad_type_info(k)%lcloud4crtm)
       if(associated(rad_type_info(k)%laerosol4crtm)) deallocate(rad_type_info(k)%laerosol4crtm)
       if(associated(rad_type_info(k)%cclr)) deallocate(rad_type_info(k)%cclr)
       if(associated(rad_type_info(k)%ccld)) deallocate(rad_type_info(k)%ccld)
       if(associated(rad_type_info(k)%cldval1)) deallocate(rad_type_info(k)%cldval1)
    end do
    if(allocated(rad_type_info)) deallocate(rad_type_info)
    if(allocated(aod_type_info)) deallocate(aod_type_info)
  end subroutine radiance_obstype_destroy


  subroutine radiance_parameter_cloudy_init
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_parameter_cloudy_init
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-07-20
!
! abstract:  This routine sets default values for variables used in
!            the cloudy radiance processing routines.
!
! program history log:
!   2015-07-20  zhu
!   2016-10-27  zhu - add ATMS
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    use mpeu_util, only: gettablesize, gettable
    implicit none

    character(len=*),parameter:: fixfilename='cloudy_radiance_info.txt'
    character(len=*),parameter:: toptablename='radiance_mod_instr_input'
    integer(i_kind) :: lunin
    character(len=20) :: tablename
    character(len=10) :: obsname
    character(len=10) :: ex_obserr
    character(len=8)  :: obsloc   ! global, sea, or, land ...
    logical :: ex_biascor,cld_effect
    logical :: pcexist
    logical :: obs_found

    integer(i_kind) i,ii,istr,ntot,nrows
    character(len=256),allocatable,dimension(:):: utable

    if (.not. icloud_fwd .or. total_rad_type<=0) return

    inquire(file=fixfilename,exist=pcexist)
    if (.not. pcexist) then 
       write(6,*)'radiance_parameter_cloudy_init: cloudy_radiance_info.txt is missing' 
       call stop2(79)
    end if
    lunin=11
    open(lunin,file=fixfilename,form='formatted')

!   Scan file for desired table first and get size of table
    call gettablesize(toptablename,lunin,ntot,nrows)
    if (mype==0) write(6,*) 'radiance_parameter_cloudy_init: ',toptablename, ntot, nrows
    if(nrows==0) then
       return
    endif

!   Get contents of table
    allocate(utable(nrows))
    call gettable(toptablename,lunin,ntot,nrows,utable)

    do ii=1,nrows
       read(utable(ii),*) obsname,obsloc,ex_obserr,ex_biascor,cld_effect
       if (mype==0) write(6,*) obsname,obsloc,ex_obserr,ex_biascor,cld_effect

       obs_found=.false.
       do i=1,total_rad_type
          if (index(trim(rad_type_info(i)%rtype),trim(obsname)) /= 0) then
             obs_found=.true.
             istr=i
             if (trim(obsloc)=='sea') rad_type_info(i)%cld_sea_only=.true.
             rad_type_info(i)%ex_obserr=ex_obserr
             rad_type_info(i)%ex_biascor=ex_biascor
             rad_type_info(i)%cld_effect=cld_effect

             if (.not. rad_type_info(i)%lcloud_fwd) then 
                rad_type_info(i)%cld_sea_only=.false.
                rad_type_info(i)%cld_effect=.false.
                rad_type_info(i)%ex_obserr=' '
                rad_type_info(i)%ex_biascor=.false.
             end if

             if (mype==0) write(6,*) 'cloudy_radiance_info for ', trim(obsname),&
                 ' cld_sea_only=', rad_type_info(i)%cld_sea_only, &
                 ' ex_obserr=', rad_type_info(i)%ex_obserr, &
                 ' ex_biascor=', rad_type_info(i)%ex_biascor

!            allocate space for entries from table, Obtain table contents
             tablename='obs_'//trim(obsname)
             if ( rad_type_info(i)%ex_obserr == 'ex_obserr3' ) then
                call sensor_parameter_table(trim(tablename),lunin,rad_type_info(i)%nchannel,rad_type_info(i)%cclr,rad_type_info(i)%ccld,rad_type_info(i)%cldval1)
             else
                call sensor_parameter_table(trim(tablename),lunin,rad_type_info(i)%nchannel,rad_type_info(i)%cclr,rad_type_info(i)%ccld)
             endif
             exit
          end if
       end do
       if (.not. obs_found) cycle
    enddo ! end of nrows
    deallocate(utable)
    close(lunin)
  end subroutine radiance_parameter_cloudy_init


  subroutine sensor_parameter_table(filename,lunin,nchal,cclr,ccld,cldval1)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    sensor_parameter_table
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-09-10
!
! abstract:  This routine retrieves parameters used for AMSUA all-sky radiance
!
! program history log:
!   2015-09-10  zhu
!   2018-08-10  mkim : added a column of data 'cldval1' in the all-sky sensor parameter table
!   2019-07-108 todling : turn cldval1 into optional to avoid forcing another column in info file when not needed
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
    use mpeu_util, only: gettablesize
    use mpeu_util, only: gettable
    use gsi_io, only: verbose
    implicit none

    character(len=*), intent(in) :: filename
    integer(i_kind) , intent(in) :: lunin
    integer(i_kind) , intent(in) :: nchal
    real(r_kind)    , dimension(nchal), intent(inout) :: cclr,ccld
    real(r_kind)    , dimension(nchal), optional, intent(inout) :: cldval1

    integer(i_kind) ii,ntot,nrows,ich0
    real(r_kind) cclr0,ccld0,cldval1_0
    character(len=256),allocatable,dimension(:):: utable
    logical print_verbose

    print_verbose=.false.
    if(verbose .and. mype == 0)print_verbose=.true.
!   Initialize the arrays
    cclr(:)=zero
    ccld(:)=zero
    if ( present(cldval1) ) then
       cldval1(:)=zero
    endif

!   Scan file for desired table first and get size of table
    call gettablesize(filename,lunin,ntot,nrows)
    if (print_verbose) write(6,*) 'sensor_parameter_table: ',filename, nrows
    if(nrows==0) then
       return
    endif

!   Get contents of table
    allocate(utable(nrows))
    call gettable(filename,lunin,ntot,nrows,utable)

!   Retrieve each token of interest from table
    do ii=1,nrows
       if (present(cldval1)) then
         read(utable(ii),*) ich0,cclr0,ccld0,cldval1_0
         cldval1(ich0)=cldval1_0
       else
         read(utable(ii),*) ich0,cclr0,ccld0
       endif
       cclr(ich0)=cclr0
       ccld(ich0)=ccld0
    enddo
    deallocate(utable)

    if (print_verbose) then
       if (present(cldval1)) then
          write(6,*) 'sensor_parameter_table: ich  cclr  ccld  cldval1'
          do ii=1,nchal
             write(6,*) ii,cclr(ii),ccld(ii),cldval1(ii)
          end do
       else
          write(6,*) 'sensor_parameter_table: ich  cclr  ccld'
          do ii=1,nchal
             write(6,*) ii,cclr(ii),ccld(ii)
          end do
       endif
    end if

  end subroutine sensor_parameter_table

  subroutine radiance_parameter_aerosol_init
    ! History:   2018-10-31 Wei/Martin - fix stub
    use aeroinfo, only: jpch_aero,nusis_aero
    use obsmod, only: ndat,dtype,dsis
    use gsi_io, only: verbose
    implicit none
    logical :: first,diffistr,found
    integer(i_kind) :: i,j,k,ii,nn1,nn2
    integer(i_kind),dimension(ndat) :: k2i
    character(10),dimension(ndat) :: rtype,rrtype,drtype
    logical print_verbose

    if (.not. iaerosol_fwd) return

    print_verbose=.false.
    if(verbose)print_verbose=.true.
    
    drtype='nonaod'
    do i=1,ndat
       rtype(i)=dtype(i)!     rtype  - observation types to process
       if (rtype(i) == 'modis_aod' .or. rtype(i) == 'viirs_aod') then
           drtype='aod'
       end if 
    end do
   
    k=0
    k2i=0
    first=.true.
    rrtype=''
    do i=1,ndat
       if (drtype(i) /= 'aod') cycle

       found=.false.
       if (first) then
          k=k+1
          rrtype(k)=rtype(i)
          k2i(k)=i
          first=.false.
       else
          do j=1,k
             if (trim(rtype(i)) == trim(rrtype(j))) then
                found=.true.
                exit
             end if
          end do
          if (.not. found) then
             k=k+1
             rrtype(k)=rtype(i)
             k2i(k)=i
          end if
       end if
    end do
    total_aod_type=k

    if (mype==0) write(6,*) 'radiance_obstype_init: total_aod_type=',k,' types are: ', rrtype(1:total_aod_type)

    if (total_aod_type<=0) return
    allocate(aod_type_info(total_aod_type))

    do k=1, total_aod_type
       aod_type_info(k)%rtype=rrtype(k)
       aod_type_info(k)%cld_sea_only=.false.
       aod_type_info(k)%ex_obserr=' '
       aod_type_info(k)%ex_biascor=.false.
       aod_type_info(k)%cld_effect=.false.
       aod_type_info(k)%lcloud_fwd=.false.
       aod_type_info(k)%lallsky=.false.
       aod_type_info(k)%laerosol_fwd=.true.
       aod_type_info(k)%laerosol=.true.

       ii=k2i(k)
       first=.true.
       nn1=0
       nn2=0
       do j=1,jpch_aero
          if (j==jpch_aero) then
             diffistr = .true.
          else
             diffistr = trim(nusis_aero(j))/=trim(nusis_aero(j+1))
          end if
          if (trim(dsis(ii))==trim(nusis_aero(j))) then
             if (first) then
                nn1=j
                first=.false.
             else
                nn2=j
             end if
             if (diffistr) exit
          end if
       end do
       if (nn1/=0 .and. nn2/=0) then
          aod_type_info(k)%nchannel=nn2-nn1+1
       else
          cycle
       end if

       allocate(aod_type_info(k)%lcloud4crtm(aod_type_info(k)%nchannel))
       allocate(aod_type_info(k)%laerosol4crtm(aod_type_info(k)%nchannel))
       aod_type_info(k)%lcloud4crtm=0
       aod_type_info(k)%laerosol4crtm=0

       if (mype==0)  &
                               write(6,*) 'radiance_obstype_init: type=',aod_type_info(k)%rtype, &
                               ' nch=',aod_type_info(k)%nchannel, &
                               ' lcloud_fwd=',aod_type_info(k)%lcloud_fwd, &
                               ' lallsky=',aod_type_info(k)%lallsky, &
                               ' laerosol_fwd=',aod_type_info(k)%laerosol_fwd, &
                               ' laerosol=',aod_type_info(k)%laerosol

       allocate(aod_type_info(k)%cclr(aod_type_info(k)%nchannel))
       allocate(aod_type_info(k)%ccld(aod_type_info(k)%nchannel))
       aod_type_info(k)%cclr(:)=9999.9_r_kind
       aod_type_info(k)%ccld(:)=zero

    end do ! end total_aod_type

  end subroutine radiance_parameter_aerosol_init


  ! CCH

  subroutine radiance_ex_obserr_evolve_gauss(radmod,chan_number, innovation, clwp_obs ,clw_guess_retrieval, &
                                             ng_predictor, tnoise, tnoise_cld, error0, io_cloud_bc)
  !$$$  subprogram documentation block
  !                .      .    .
  ! subprogram:    radiance_ex_obserr_evolve_gauss
  !
  !   prgrmmr:    Chih-Chi Hu                  date: xxxx-xx-xx
  !
  ! abstract:  This routine includes evolving-Gaussian techniques 
  !            to account for the non-Gaussian errors. Note that 
  !            this is different from radiance_ex_obserr1, which works
  !            for all the channels at once. This subroutine works for
  !            one channel at a time. Need a do loop to wrap this 
  !            subroutine for all the channels.
  !
  !            Note: This code should work universally across different
  !                  sensors (atms, amsua, etc)
  !
  ! program history log:
  !
  !   input argument list:
  !      radmod      = the radiance object ( check: %rtype = sensor name (e.g., amsua) )
  !      chan_number = which channel 
  !      innovation  = O-B innovation
  !      clwp_obs, clw_guess_retrieval = cloud amount
  !      ng_predictor = predictor for the obs error model
  !      tnoise, tnoise_cld = original error stdev for clear sky and cloudy sky
  !
  !   output argument list:
  !      error0      = output error stdev (will overwrite the original values)
  !
  ! attributes:
  !   language: f90
  !   machine: 
  !
  !$$$ end documentation block

    use kinds, only: i_kind,r_kind
    implicit none

    !integer(i_kind),intent(in) :: nchanl   ! don't need this anymore
    integer(i_kind), intent(in)    :: chan_number  ! Instead, need to know which channel are we working on 
    real(r_kind),    intent(inout) :: innovation   ! O-B innovation
                                                   ! if io_cloud_bc == .true. , innovation will be changed

    character(len=*), intent(in) :: ng_predictor ! thepredictor; e.g., 'sym_cld', 'obs_cld', 'no_pred'

    real(r_kind),intent(in) :: clwp_obs, clw_guess_retrieval
    real(r_kind),intent(in) :: tnoise,tnoise_cld
    real(r_kind),intent(inout) :: error0
    type(rad_obs_type),intent(in) :: radmod
    logical, intent(in) :: io_cloud_bc

    integer(i_kind) :: i, unit_number, fstatus
    real(r_kind) :: clwtmp
    real(r_kind) :: cclr,ccld
    character(len=10)  :: chan_number_string
    character(len=256) :: ng_table_name
    character(len=1200) :: line, abs_path

    integer(i_kind) :: num_cloud_cat, num_hist_pdf
    real(r_kind)    :: max_range, dx, bdy_slope
    real(r_kind), dimension(:),   allocatable :: cloud_bin, bias
    real(r_kind), dimension(:,:), allocatable :: stdev

    real(r_kind)    :: cld_predictor, dTB, error1, inv_upp, inv_low, wt_upp, wt_low
    integer(i_kind) :: cld_cat, ndx_low, ndx_upp, ndx_left, ndx_right

    real(r_kind)    :: innov_tmp, inov_left, inov_right ! temporary storage for innovation

    !do i=1,nchanl
    !   cclr(i)=radmod%cclr(i)
    !   ccld(i)=radmod%ccld(i)
    !end do

    ! if (radmod%lcloud4crtm(chan_number)<0) cycle  ==> check the meaning!

    !clwtmp=half*(clwp_amsua+clw_guess_retrieval)

    write(chan_number_string, '(I0)') chan_number

    ! needs to be changed here:
    !abs_path ='/scratch2/GFDL/gfdlscr/Chih-Chi.Hu/SHiELD/shield_workflow/fix/fix_gsi/'
    !ng_table_name = trim(abs_path)//'non_Gaussian_'//trim(radmod%rtype)//'_ch'//trim(chan_number_string)//'_OmF_sym_cld.txt'
    ng_table_name = 'ng_'//trim(radmod%rtype)//'_ch'//trim(chan_number_string)//'.txt'
    

    ! read the non-Gaussian table:
    call read_nongauss_table(ng_table_name, max_range, dx, bdy_slope, &
                             num_cloud_cat, num_hist_pdf, cloud_bin,  &
                             bias, stdev ) 

    ! define the predictor:
    if ( trim(ng_predictor) == 'sym_cld' .or. trim(ng_predictor) == 'no_pred' ) then
       cld_predictor = half*(clwp_obs+clw_guess_retrieval)
    elseif ( trim(ng_predictor) == 'obs_cld' ) then
       cld_predictor = clwp_obs
    else
       write(6,*) 'ng_predictor = ', ng_predictor, ' is not supported!'
       stop
    endif

    if (cld_predictor > 1.0) cld_predictor = 1.0

    ! look for which cloud cat:
    ! when cld_cat = 0  : cld_predictor <= 0 --> make this case cld_cat = 1
    ! when cld_cat = 1  : cld_predictor is between [cloud_bin(1), cloud_bin(2)]
    ! when cld_cat = n  : cld_predictor is between [cloud_bin(n), cloud_bin(n+1)]
    do i=1,num_cloud_cat+1
       if (cld_predictor - cloud_bin(i) <= 0) then
          cld_cat = i-1
          exit
       endif
    enddo

    if ( cld_cat == 0 ) cld_cat = 1

    !write(6,*) 'cloud = ', cld_predictor, 'is between', cloud_bin(cld_cat), 'and ', cloud_bin(cld_cat+1)


    ! move the innovation to the 'mode-relative coordinate'
    ! Note: this is assuming all the biases in O-B pdfs are attributed to 
    !       the background error
    innov_tmp = innovation
    if ( .not. io_cloud_bc ) innov_tmp = innov_tmp + bias(cld_cat)

    !write(6,*) 'bias = ', bias

    ! find the stdev based on the cloud category
    if (innov_tmp >= max_range ) then
       dTB = innov_tmp - max_range
       error1 = stdev(cld_cat,num_hist_pdf) + bdy_slope*dTB ! extrapolation

    elseif (innov_tmp <= -max_range ) then
       dTB = -max_range - innov_tmp
       error1 = stdev(cld_cat, 1) + bdy_slope*dTB ! extrapolation

    else
       ! if innovation is within [-max_range, max_range]
       ! interpolate from the table to get the stdev:

       !ndx_low = floor( (innov_tmp - (-max_range))/dx ) ! ndx = how many dx
       !ndx_upp = ndx_low + 1

       !inv_low = -max_range + dx*(ndx_low)
       !inv_upp = -max_range + dx*(ndx_upp)

       !wt_upp  = innov_tmp - inv_low
       !wt_low  = inv_upp - innov_tmp
       !error0 = (wt_upp*stdev(cld_cat,ndx_upp+1) + wt_low*stdev(cld_cat,ndx_low+1))/dx

      ndx_left  = floor( (innov_tmp - (-max_range))/dx ) ! ndx = how many dx
      ndx_right = ndx_left + 1

      inov_left  = -max_range + dx*(ndx_left )
      inov_right = -max_range + dx*(ndx_right)

      call interpolate(inov_left,  stdev(cld_cat, ndx_left+1), &
                       inov_right, stdev(cld_cat, ndx_right+1), &
                       innov_tmp,  error0 )

    endif

    if (io_cloud_bc) innovation = innovation + bias(cld_cat)

    return

  end subroutine radiance_ex_obserr_evolve_gauss

 ! CCH

  subroutine radiance_ex_obserr_evolve_gauss_interp(radmod,chan_number, &
                                                    innovation,clwp_obs ,clw_guess_retrieval, &
                                                    ng_predictor, tnoise, tnoise_cld, error0, io_cloud_bc)
  !$$$  subprogram documentation block
  !                .      .    .
  ! subprogram:    radiance_ex_obserr_evolve_gauss_interp
  !
  !   prgrmmr:    Chih-Chi Hu                  date: 2024-03-14
  !
  ! abstract:  This routine works similar to radiance_ex_obserr_evolve_gauss, 
  !            except that there is an additional interpolation between
  !            different (either symmetric, obs) cloud bins. 
  !
  ! program history log:
  !
  !   input argument list:
  !      radmod      = the radiance object ( check: %rtype = sensor name (e.g.,
  !      amsua) )
  !      chan_number = which channel 
  !      innovation  = O-B innovation
  !      clwp_obs, clw_guess_retrieval = cloud amount
  !      ng_predictor = predictor for the obs error model
  !      tnoise, tnoise_cld = original error stdev for clear sky and cloudy sky
  !
  !   output argument list:
  !      error0      = output error stdev (will overwrite the original values)
  !
  ! attributes:
  !   language: f90
  !   machine: 
  !
  !$$$ end documentation block

    use kinds, only: i_kind,r_kind
    implicit none

    !integer(i_kind),intent(in) :: nchanl   ! don't need this anymore
    integer(i_kind), intent(in)    :: chan_number ! Instead, need to know which channel are we working on 
    real(r_kind),    intent(inout) :: innovation  ! O-B innovation
                                                  ! if io_cloud_bc == .true. ,
                                                  ! innovation will be changed

    character(len=*), intent(in) :: ng_predictor ! thepredictor; e.g., 'sym_cld', 'obs_cld', 'no_pred'

    real(r_kind),intent(in) :: clwp_obs, clw_guess_retrieval
    real(r_kind),intent(in) :: tnoise,tnoise_cld
    real(r_kind),intent(inout) :: error0
    type(rad_obs_type),intent(in) :: radmod
    logical, intent(in) :: io_cloud_bc

    integer(i_kind) :: i, unit_number, fstatus
    real(r_kind) :: clwtmp
    real(r_kind) :: cclr,ccld
    character(len=10)  :: chan_number_string
    character(len=256) :: ng_table_name
    character(len=1200) :: line, abs_path

    integer(i_kind) :: num_cloud_cat, num_hist_pdf
    real(r_kind)    :: max_range, dx, bdy_slope
    real(r_kind), dimension(:),   allocatable :: cloud_bin, bin_center, bias
    real(r_kind), dimension(:,:), allocatable :: stdev

    integer(i_kind) :: cld_cat_low, cld_cat_upp 
    real(r_kind)    :: innov_tmp_low, innov_tmp_upp, stdev_low, stdev_upp
    integer(i_kind) :: ndx_left, ndx_right
    real(r_kind)    :: inov_left, inov_right
    real(r_kind)    :: cld_predictor, dTB, bias_interp


    write(chan_number_string, '(I0)') chan_number
    ng_table_name ='ng_'//trim(radmod%rtype)//'_ch'//trim(chan_number_string)//'.txt'

    ! read the non-Gaussian table:
    call read_nongauss_table(ng_table_name, max_range, dx, bdy_slope, &
                             num_cloud_cat, num_hist_pdf, cloud_bin,  &
                             bias, stdev )

    ! define the predictor:
    if ( trim(ng_predictor) == 'sym_cld' .or. trim(ng_predictor) == 'no_pred' ) then
       cld_predictor = half*(clwp_obs+clw_guess_retrieval)
    elseif ( trim(ng_predictor) == 'obs_cld' ) then
       cld_predictor = clwp_obs
    else
       write(6,*) 'ng_predictor = ', ng_predictor, ' is not supported!'
       stop
    endif

    if (cld_predictor > 1.0) cld_predictor = 1.0
    if (cld_predictor < 0.0) cld_predictor = 0.0

   ! define the bin center:
   ! Note the estimated discretized cloud-dependent pdfs are assumed
   ! to be defined on the bin center

   allocate(bin_center(num_cloud_cat))

   do i=1,num_cloud_cat
      bin_center(i) = 0.5*( cloud_bin(i) + cloud_bin(i+1) )
   enddo

   !if (mype==150) write(6,*) 'In radiance_mod:: Sensor, channel, bin_center =', radmod%rtype, chan_number_string, bin_center

   ! look for which cloud cat:
   ! i.e., cld_predictor is in [lower_cld_cat, upper_cld_cat]
   if ( cld_predictor <= bin_center(1) ) then ! if in the smallest category
      cld_cat_low = 1
      cld_cat_upp = 1
   elseif ( cld_predictor >= bin_center(num_cloud_cat) ) then ! if in the largest category
      cld_cat_low = num_cloud_cat
      cld_cat_upp = num_cloud_cat
   else
      do i=1,num_cloud_cat
         if ( cld_predictor - bin_center(i) <0 ) then
            cld_cat_low = i-1
            cld_cat_upp = i
            exit
         endif
      enddo
   endif

   !if (mype==150) write(6,*) 'In radiance_mod:: cld_predictor = ', cld_predictor
   !if (mype==150) write(6,*) 'In radiance_mod:: innovation = ', innovation
   !if (mype==150) write(6,*) 'In radiance_mod:: cld_cat_low = ', cld_cat_low, 'cld_cat_upp = ',cld_cat_upp

   ! interpolation:
   ! (1) interpolate the innovation in the lower and upper cloud cat, respectively
   !     so we get two different stdevs
   ! (2) interpolate the cld cat based on the bin center

   ! === lower cat ===

   ! move the innovation to the 'mode-relative coordinate'
   ! Note: this is assuming all the biases in O-B pdfs are attributed to 
   !       the background error
   innov_tmp_low = innovation
   if ( .not. io_cloud_bc ) innov_tmp_low = innov_tmp_low + bias(cld_cat_low)

   ! find the stdev based on the cloud category
   if (innov_tmp_low >= max_range ) then
      dTB = innov_tmp_low - max_range
      stdev_low = stdev(cld_cat_low,num_hist_pdf) + bdy_slope*dTB !extrapolation

   elseif (innov_tmp_low <= -max_range ) then
      dTB = -max_range - innov_tmp_low
      stdev_low = stdev(cld_cat_low, 1) + bdy_slope*dTB ! extrapolation

   else
      ! if innovation is within [-max_range, max_range]
      ! interpolate from the table to get the stdev:
      ndx_left  = floor( (innov_tmp_low - (-max_range))/dx ) ! ndx = how many dx
      ndx_right = ndx_left + 1

      inov_left  = -max_range + dx*(ndx_left )
      inov_right = -max_range + dx*(ndx_right)

      call interpolate(inov_left,  stdev(cld_cat_low, ndx_left+1), &
                       inov_right, stdev(cld_cat_low, ndx_right+1), &
                       innov_tmp_low, stdev_low )
   endif

   !if (mype==150) write(6,*) 'In radiance_mod:: lower cat stdev = ', stdev_low



   ! === upper cat ===

   ! move the innovation to the 'mode-relative coordinate'
   ! Note: this is assuming all the biases in O-B pdfs are attributed to 
   !       the background error
   innov_tmp_upp = innovation
   if ( .not. io_cloud_bc ) innov_tmp_upp = innov_tmp_upp + bias(cld_cat_upp)

   ! find the stdev based on the cloud category
   if (innov_tmp_upp >= max_range ) then
      dTB = innov_tmp_upp - max_range
      stdev_upp = stdev(cld_cat_upp,num_hist_pdf) + bdy_slope*dTB !extrapolation

   elseif (innov_tmp_upp <= -max_range ) then
      dTB = -max_range - innov_tmp_upp
      stdev_upp = stdev(cld_cat_upp, 1) + bdy_slope*dTB ! extrapolation

   else
      ! if innovation is within [-max_range, max_range]
      ! interpolate from the table to get the stdev:
      ndx_left  = floor( (innov_tmp_upp - (-max_range))/dx ) ! ndx = how many dx
      ndx_right = ndx_left + 1

      inov_left  = -max_range + dx*(ndx_left )
      inov_right = -max_range + dx*(ndx_right)

      call interpolate(inov_left,  stdev(cld_cat_upp, ndx_left+1), &
                       inov_right, stdev(cld_cat_upp, ndx_right+1), &
                       innov_tmp_upp, stdev_upp )
   endif

   !if (mype==150) write(6,*) 'In radiance_mod:: upper cat stdev = ', stdev_upp


   ! === interpolate the lower and upper cat ===
   call interpolate(bin_center(cld_cat_low), stdev_low, &
                    bin_center(cld_cat_upp), stdev_upp, &
                    cld_predictor, error0 )

   if (io_cloud_bc) then
       call interpolate(bin_center(cld_cat_low), bias(cld_cat_low), &
                        bin_center(cld_cat_upp), bias(cld_cat_upp), &
                        cld_predictor, bias_interp )
       innovation = innovation - bias_interp
   endif

   !if (mype==150) write(6,*) 'In radiance_mod:: final interpolated stdev = ', error0

   return

  end subroutine radiance_ex_obserr_evolve_gauss_interp



  subroutine interpolate(x1,y1,x2,y2,x,y) 
      use kinds, only: i_kind,r_kind
      implicit none
      real(r_kind), intent(in)  :: x1, y1, x2, y2, x
      real(r_kind), intent(out) :: y

      ! Check if x1 equals x2 to avoid division by zero
      if (x1 == x2) then
          y = y1
      else
          ! Perform linear interpolation
          y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
      end if

  end subroutine interpolate


  subroutine read_nongauss_table(filename, max_range, dx, bdy_slope, &
                                 num_cloud_cat, num_hist_pdf, cloud_bin, &
                                 bias, stdev)
   implicit none

   character(len=*), intent(in) :: filename
   real(r_kind), intent(out)    :: max_range, dx, bdy_slope
   integer(i_kind), intent(out) :: num_cloud_cat, num_hist_pdf
   real(r_kind), dimension(:), allocatable, intent(out)   :: cloud_bin, bias
   real(r_kind), dimension(:,:), allocatable, intent(out) :: stdev


   character(len=1200) :: line
   character(len=256), dimension(2) :: title
   integer(i_kind) :: unit_number, status, i, n, num_chars_read
   real(r_kind)    :: cri_costfn
   character(len=200) :: chan_number_string, ng_table_name

   !write(chan_number_string, '(I0)') 15
   !ng_table_name ='non_Gaussian_table_amsua_ch'//adjustl(trim(chan_number_string))//'.txt'
   !write(*,*) ng_table_name

  ! Open the file for reading
  open(newunit=unit_number, file=filename, status='old',action='read',iostat=status)

  if (status /= 0) then
    write(6,*) 'In radiance_mod.f90/read_nongauss table :: Error opening file ', filename
    stop
  end if

 ! Read the parameter block

  n=1 ! index for the line number of the title

  do

    read(unit_number, '(A)', iostat=status) line

    if (status /= 0) exit

    if(trim(line)=='') cycle ! advance empty line

    if(index(line,trim('nonGaussian observation error table::')) /= 0) then !table

       do ! start read table
          read(unit_number, '(A)', iostat=status) line  ! read the next line

          if (line(1:1)=='!') cycle ! skip the comments
          if (line(1:1)==':') exit  ! end of the table

          if(index(line(1:1),'#') /= 0) then
             title(n) = line
             n=n+1
             cycle
          endif

          if(index(line(1:15),'max_range') /= 0) then
             num_chars_read = index(line, '=') + 1
             read(line(num_chars_read:),*) max_range
             cycle
          elseif(index(line(1:15),'dx') /= 0) then
             num_chars_read = index(line, '=') + 1
             read(line(num_chars_read:),*) dx
             cycle
          elseif(index(line(1:15),'bdy_slope') /= 0) then
             num_chars_read = index(line, '=') + 1
             read(line(num_chars_read:),*) bdy_slope
             cycle
          elseif(index(line(1:15),'cri_costfn') /= 0) then
             num_chars_read = index(line, '=') + 1
             read(line(num_chars_read:),*) cri_costfn
             cycle
          elseif(index(line(1:15),'num_cloud_cat') /= 0) then
             num_chars_read = index(line, '=') + 1
             read(line(num_chars_read:),*) num_cloud_cat
             cycle
          endif
       enddo ! end read table

    endif ! table

    if(index(line,trim('cloud_bin::')) /= 0) then ! cloud bin
       if(.not. allocated(cloud_bin)) allocate(cloud_bin(num_cloud_cat+1))
       read(unit_number, '(A)', iostat=status) line  ! read the nextline
       read(line, *) cloud_bin
    endif ! cloud bin


    if(index(line,trim('bias::')) /= 0) then ! bias
       if(.not. allocated(bias)) allocate(bias(num_cloud_cat)  )
       read(unit_number, '(A)', iostat=status) line  ! read the nextline
       read(line, *) bias
    endif ! bias


    if(index(line,trim('stdev::')) /= 0) then ! stdev
       num_hist_pdf = 2*int(max_range/dx)+1 ! number of discretized histograms for a pdf
       if(.not. allocated(stdev)) allocate(stdev(num_cloud_cat, num_hist_pdf))

       do i=1, num_cloud_cat
          read(unit_number, '(A1200)', iostat=status) line  ! read the nextline
          read(line, *) stdev(i,:)
       enddo

    endif ! stdev

  end do


  ! Close the file
  close(unit_number)

  !write(6,*) '=== start of read_and_output.F90 ==='
  !write(6,*) '#   summary of read info:   '
  !write(6,*) trim(title(1))
  !write(6,*) trim(title(2))
  !write(6,*) ''
  !write(6,'(A, f8.2)') 'max_range     = ', max_range
  !write(6,'(A, f8.2)') 'dx            = ', dx
  !write(6,'(A, f8.2)') 'bdy_slope     = ', bdy_slope
  !write(6,'(A, I4)') 'num_cloud_cat = ', num_cloud_cat
  !write(6,'(A, I4)') 'num_hist_pdf  = ', num_hist_pdf
  !write(6,'(A, *(f8.2))') 'cloud bin     = ', cloud_bin
  !write(6,'(A, *(f8.2))') 'bias          = ', bias
  !do i=1,num_cloud_cat
  !   write(6,'(A, I2, A, *(F8.2))') 'stdev (cloud cat',i,') = ', stdev(i,:)
  !enddo
  !write(6,*) ''
  !write(6,*) '=== end of read_and_output.F90 ==='


  end subroutine read_nongauss_table





  subroutine radiance_ex_obserr_1(radmod,nchanl,clwp_amsua,clw_guess_retrieval, &
                                tnoise,tnoise_cld,error0)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_obserr_1
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-09-10
!
! abstract:  This routine includes extra observation error assignment routines.
!
! program history log:
!   2015-09-10  zhu
!   2016-10-27  zhu - add ATMS
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
    implicit none
    
    integer(i_kind),intent(in) :: nchanl
    real(r_kind),intent(in) :: clwp_amsua,clw_guess_retrieval
    real(r_kind),dimension(nchanl),intent(in):: tnoise,tnoise_cld
    real(r_kind),dimension(nchanl),intent(inout) :: error0
    type(rad_obs_type),intent(in) :: radmod 

    integer(i_kind) :: i
    real(r_kind) :: clwtmp
    real(r_kind),dimension(nchanl) :: cclr,ccld

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
       ccld(i)=radmod%ccld(i)
    end do

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       clwtmp=half*(clwp_amsua+clw_guess_retrieval)
       if(clwtmp <= cclr(i)) then
          error0(i) = tnoise(i)
       else if(clwtmp > cclr(i) .and. clwtmp < ccld(i)) then
          error0(i) = tnoise(i) + (clwtmp-cclr(i))* &
                      (tnoise_cld(i)-tnoise(i))/(ccld(i)-cclr(i))
       else
          error0(i) = tnoise_cld(i)
       endif
    end do
    return

  end subroutine radiance_ex_obserr_1

  subroutine radiance_ex_obserr_2(radmod,nchanl,cldeff1,cldeff2,tnoise,tnoise_cld,error0)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_obserr_1
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-09-10
!
! abstract:  This routine includes extra observation error assignment routines.
!            
!
! program history log:
!   2018-04-04  zhu/bi -- adapted from radiance_ex_obserr_1
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
    implicit none

    integer(i_kind),intent(in) :: nchanl
    real(r_kind),dimension(nchanl),intent(in) :: cldeff1,cldeff2
    real(r_kind),dimension(nchanl),intent(in) :: tnoise,tnoise_cld
    real(r_kind),dimension(nchanl),intent(inout) :: error0
    type(rad_obs_type),intent(in) :: radmod

    integer(i_kind) :: i
    real(r_kind) :: cldeff
    real(r_kind),dimension(nchanl) :: cclr,ccld

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
       ccld(i)=radmod%ccld(i)
    end do

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       cldeff=half*(abs(cldeff1(i))+abs(cldeff2(i)))
       if(cldeff <= cclr(i)) then
          error0(i) = tnoise(i)
       else if(cldeff > cclr(i) .and. cldeff < ccld(i)) then
          error0(i) = tnoise(i) + (cldeff-cclr(i))* &
                      (tnoise_cld(i)-tnoise(i))/(ccld(i)-cclr(i))
       else
          error0(i) = tnoise_cld(i)
       endif
    end do
    return

  end subroutine radiance_ex_obserr_2

  subroutine radiance_ex_biascor_1(radmod,nchanl,tsim_bc,tsavg5,zasat, & 
                       clw_guess_retrieval,clwp_amsua,cld_rbc_idx,ierrret)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_biascor_1
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-09-20
!
! abstract:  This routine include extra radiance bias correction routines.
!
! program history log:
!   2015-09-20  zhu
!   2016-10-27  zhu - add ATMS
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
    use clw_mod, only: ret_amsua
    implicit none

    integer(i_kind)                   ,intent(in   ) :: nchanl
    real(r_kind),dimension(nchanl)    ,intent(in   ) :: tsim_bc
    real(r_kind)                      ,intent(in   ) :: tsavg5,zasat
    real(r_kind),dimension(nchanl)    ,intent(inout) :: cld_rbc_idx
    real(r_kind)                      ,intent(inout) :: clwp_amsua
    real(r_kind)                      ,intent(inout) :: clw_guess_retrieval
    type(rad_obs_type)                ,intent(in)    :: radmod
    integer(i_kind)                   ,intent(  out) :: ierrret

    integer(i_kind) :: i
    real(r_kind),dimension(nchanl) :: cclr

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
    end do

!   call ret_amsua(tb_obs,nchanl,tsavg5,zasat,clwp_amsua,ierrret)
    call ret_amsua(tsim_bc,nchanl,tsavg5,zasat,clw_guess_retrieval,ierrret) 

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       if ((clwp_amsua-cclr(i))*(clw_guess_retrieval-cclr(i))<zero  &
          .and. abs(clwp_amsua-clw_guess_retrieval)>=0.005_r_kind) cld_rbc_idx(i)=zero
    end do
    return

  end subroutine radiance_ex_biascor_1

  subroutine radiance_ex_biascor_2(radmod,nchanl,cldeff1,cldeff2,cld_rbc_idx)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_biascor_1
!
!   prgrmmr:    yanqiu zhu      org: np23                date: 2015-09-20
!
! abstract:  This routine include extra radiance bias correction routines using
!            cloud effect.
!
! program history log:
!   2018-04-04  zhu - adapted from radiance_ex_biascor_1
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
    use clw_mod, only: ret_amsua
    implicit none

    integer(i_kind)                   ,intent(in   ) :: nchanl
    real(r_kind),dimension(nchanl)    ,intent(inout) :: cld_rbc_idx
    real(r_kind),dimension(nchanl)    ,intent(in) :: cldeff1
    real(r_kind),dimension(nchanl)    ,intent(in) :: cldeff2
    type(rad_obs_type)                ,intent(in) :: radmod

    integer(i_kind) :: i
    real(r_kind),dimension(nchanl) :: cclr

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
    end do

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       if ((abs(cldeff1(i))-cclr(i))*(abs(cldeff2(i))-cclr(i))<zero  &
          .and. abs(cldeff1(i)-cldeff2(i))>=0.1_r_kind) cld_rbc_idx(i)=zero
    end do
    return

  end subroutine radiance_ex_biascor_2

  subroutine radiance_ex_obserr_gmi(radmod,nchanl,clw_obs,clw_guess_retrieval,tnoise,tnoise_cld,error0)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_obserr_3
!
!   prgrmmr:    min-jeong kim      org: np23                date: 2018-08-10
!
! abstract:  This routine include extra radiance bias correction routines.
!
! program history log:
!   2015-09-20  zhu
!   2016-10-27  zhu - add ATMS
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block

    use kinds, only: i_kind,r_kind
    implicit none 

    integer(i_kind),intent(in) :: nchanl
    real(r_kind),intent(in) :: clw_obs,clw_guess_retrieval
    real(r_kind),dimension(nchanl),intent(in):: tnoise,tnoise_cld
    real(r_kind),dimension(nchanl),intent(inout) :: error0
    type(rad_obs_type),intent(in) :: radmod

    integer(i_kind) :: i 
    real(r_kind) :: clwavg
    real(r_kind),dimension(nchanl) :: cclr,ccld,tnoise_cldval1

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
       ccld(i)=radmod%ccld(i)
       tnoise_cldval1(i)=radmod%cldval1(i)
    end do

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       clwavg=half*(clw_obs+clw_guess_retrieval)
       if(clwavg < cclr(i)) then 
          error0(i) = tnoise(i)
       else if(clwavg >= cclr(i) .and. clwavg < ccld(i)) then 
          error0(i) = tnoise_cldval1(i) + (tnoise(i) - tnoise_cldval1(i))*(ccld(i)-clwavg)/(ccld(i)-cclr(i))
       else if(clwavg >= ccld(i) .and. clwavg < 0.5_r_kind) then
          error0(i) = tnoise_cld(i) + (tnoise_cldval1(i)-tnoise_cld(i)) * (0.5_r_kind-clwavg)/(0.5_r_kind-ccld(i))
       else
          error0(i) = tnoise_cld(i)
       endif
    end do
    return

!  end subroutine radiance_ex_obserr_3
  end subroutine radiance_ex_obserr_gmi

!  subroutine radiance_ex_biascor_3(radmod,nchanl,tsim_bc,tsavg5,zasat, &
  subroutine radiance_ex_biascor_gmi(radmod,clw_obs,clw_guess_retrieval,nchanl,cld_rbc_idx)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    radiance_ex_biascor_gmi
!
!   prgrmmr:    min-jeong kim      org: np23                date: 2085-08-10
!
! abstract:  This routine include extra radiance bias correction routines.
!
! program history log:
!   2018-08-10  mkim
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origin 2000; Compaq/HP
!
!$$$ end documentation block
    use kinds, only: i_kind,r_kind
   use constants, only: zero,one

    implicit none

    integer(i_kind)                   ,intent(in   ) :: nchanl
    real(r_kind),dimension(nchanl)    ,intent(inout) :: cld_rbc_idx
    real(r_kind)                      ,intent(inout) :: clw_obs
    real(r_kind)                      ,intent(inout) :: clw_guess_retrieval
    type(rad_obs_type)                ,intent(in)    :: radmod

    integer(i_kind) :: i
    real(r_kind),dimension(nchanl) :: cclr

    do i=1,nchanl
       cclr(i)=radmod%cclr(i)
    end do

    do i=1,nchanl
       if (radmod%lcloud4crtm(i)<0) cycle
       if ((clw_obs-cclr(i))*(clw_guess_retrieval-cclr(i))<zero .and. abs(clw_obs-clw_guess_retrieval)>=0.005_r_kind) cld_rbc_idx(i)=zero
    end do
    return

!  end subroutine radiance_ex_biascor_3
  end subroutine radiance_ex_biascor_gmi


end module radiance_mod

