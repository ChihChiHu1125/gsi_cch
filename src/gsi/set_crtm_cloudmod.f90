module set_crtm_cloudmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   set_crtm_cloudmod
!  prgmmr: todling          org: gmao                date: 2011-06-01
!
! abstract: module providing interface to set-crtm-cloud procedures
!
! program history log:
!   2011-06-01  todling
!   2011-11-17  zhu     --- merge set_crtm_cloudmod with crtm_cloud
!   2018-05-19  eliu    --- add precipiation components (related to GFDL physics)
!   2022-12-09  mtong   --- add capability to specify hydrometeor habit
!
! subroutines included:
!   sub Set_CRTM_Cloud
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds, only: i_kind,r_kind
  use constants, only: zero,one,two,five,r0_05,t0c,fv,rd,grav
  use CRTM_Cloud_Define, only: CRTM_Cloud_type
  use CRTM_Cloud_Define, only: WATER_CLOUD,ICE_CLOUD,RAIN_CLOUD, &
      SNOW_CLOUD,GRAUPEL_CLOUD,HAIL_CLOUD,PlateType1,ColumnType1, &
      SixBulletRosette,Perpendicular4_BulletRosette,Flat3_BulletRosette, &
      IconCloudIce,SectorSnowflake,EvansSnowAggregate,EightColumnAggregate, &
      LargePlateAggregate,LargeColumnAggregate,LargeBlockAggregate, &
      IconSnow,IconHail,GemGraupel,GemSnow,GemHail,IceSphere,LiquidSphere 
  use mpeu_util, only: die
  use mpimod, only: mype          
  use radiance_mod, only: cw_cv 
  use ncepnems_io, only: imp_physics 
  implicit none

private
public Set_CRTM_Cloud

CONTAINS

  subroutine Set_CRTM_Cloud ( km, nac, cloud_name, icmask, nc, cloud_cont, cloud_efr,jcloud, &
                              dp, tp, pr, qh, cloud, lprecip, ccw_cont) 

  implicit none

  integer(i_kind) , intent(in)    :: km                ! number of levels
  integer(i_kind) , intent(in)    :: nac               ! number of actual clouds
  character(len=*), intent(in)    :: cloud_name(nac)   ! [nac]   Model cloud names: qi, ql, etc.
  logical,          intent(in)    :: icmask            ! mask determining where to consider clouds
  logical,          intent(in)    :: lprecip           ! mask determining where to consider clouds
  integer(i_kind),  intent(in)    :: nc                ! number of clouds
  integer(i_kind),  intent(in)    :: jcloud(nc)        ! cloud index
  real(r_kind),     intent(in)    :: cloud_cont(km,nc) ! cloud content 
  real(r_kind),     intent(in)    :: cloud_efr (km,nc) ! cloud effective radius
  real(r_kind),     intent(in)    :: dp(km)            ! [km]    
  real(r_kind),     intent(in)    :: tp(km)            ! [km]   atmospheric temperature (K)
  real(r_kind),     intent(in)    :: pr(km)            ! [km]   atmospheric pressure  
  real(r_kind),     intent(in)    :: qh(km)            ! [km]   specific humidity
  real(r_kind),intent(in),optional :: ccw_cont(km)     ! convective cloud content

  type(CRTM_Cloud_type), intent(inout) :: cloud(nc)    ! [nc]   CRTM Cloud object

  if (present(ccw_cont))then
     call setCloud (cloud_name, icmask, cloud_cont, cloud_efr, jcloud, dp, tp, &
                    pr, qh, cloud, lprecip, ccw_cont)
  else
     call setCloud (cloud_name, icmask, cloud_cont, cloud_efr, jcloud, dp, tp, pr, qh, cloud, lprecip)
  end if

  end subroutine Set_CRTM_Cloud

 
  subroutine setCloud (cloud_name, icmask, cloud_cont, cloud_efr,jcloud, dp, tp, pr, qh, &
                       cloud, lprecip, ccw_cont)

  use gridmod, only: regional,wrf_mass_regional
  use wrf_params_mod, only: cold_start
  use radinfo, only: hydrotype
  implicit none

! !ARGUMENTS:

  character(len=*), intent(in)    :: cloud_name(:)     ! [nc]    Model cloud names: Water, Ice, etc.
  logical,          intent(in)    :: icmask            !         mask for where to consider clouds  
  logical,          intent(in)    :: lprecip           !         mask for where to consider clouds  
  integer(i_kind),  intent(in)    :: jcloud(:)         !         cloud order
  real(r_kind),     intent(in)    :: cloud_cont(:,:)   ! [km,nc] cloud contents  (kg/m2)
  real(r_kind),     intent(in)    :: cloud_efr (:,:)   ! [km,nc] cloud effective radius (microns)
  real(r_kind),     intent(in)    :: dp(:)             ! [km]    layer thickness   
  real(r_kind),     intent(in)    :: tp(:)             ! [km]    atmospheric temperature (K)
  real(r_kind),     intent(in)    :: pr(:)             ! [km]    atmospheric pressure (??)
  real(r_kind),     intent(in)    :: qh(:)             ! [km]    atmospheric specific humidity (??)
  real(r_kind),intent(in),optional :: ccw_cont(:)      ! [km]    convective cloud contents  (kg/m2) 

  type(CRTM_Cloud_type), intent(inout) :: cloud(:)     ! [nc]   CRTM Cloud object

! !DESCRIPTION: Set the CRTM Cloud object given Model cloud properties.
!
! !REVISION HISTORY:
!
! 03May2011  Min-Jeong  Initial version.
! 14May2011  Todling    Encapsulate Min-Jeong's code in present module.
! 01July2011 Zhu        Add jcloud and cloud_efr; add codes for the regional 
! 19Feb2013  Zhu        Add cold_start for the regional
!
!EOP
!-----------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'setCloud'
  integer(i_kind) :: na, nc, km, n, k
  real(r_kind)    :: tem1,tem2,tem3,tem4
  km = size(cloud_cont,1)
  nc = size(cloud_cont,2)
  na = size(cloud_name)

! Handle hand-split case as particular case
! -----------------------------------------

  if (lprecip) cold_start=.false.

! if (cold_start .or. (na /= nc .and. (.not. regional))) then 
  if (cold_start .or. cw_cv) then                              
!    Initialize Loop over clouds ...
     do n = 1, nc
        Cloud(n)%Type = CloudType_(cloud_name(jcloud(n)),hydrotype)
        Cloud(n)%water_content(:) = zero
        cloud(n)%Effective_Radius(:) = zero
        cloud(n)%effective_variance(:) = two
     enddo

     if(icmask) then
        Cloud(1)%water_content(:) = cloud_cont(:,1)
        Cloud(2)%water_content(:) = cloud_cont(:,2)
     else
        Cloud(1)%water_content(:) = zero
        Cloud(2)%water_content(:) = zero
     endif

!    Calculate effective radius for each cloud component (wired to 2)
!    ----------------------------------------------------------------
     if(icmask) then
        do k=1,km
           ! liquid water cloud drop size
           tem4=max(zero,(t0c-tp(k))*r0_05)
           if (cloud(1)%water_content(k) > 1.0e-6_r_kind) &    
           cloud(1)%effective_radius(k) = five + five * min(one, tem4)

           ! ice water cloud particle size
           tem2 = tp(k) - t0c
           tem1 = grav/rd
           tem3 = tem1 * cloud(2)%water_content(k) * (pr(k)/dp(k)) &
                 /tp(k) * (one + fv * qh(k))

           if (cloud(2)%water_content(k) > 1.0e-6_r_kind) then     
              if (tem2 < -50.0_r_kind ) then
                 cloud(2)%effective_radius(k) =  (1250._r_kind/9.917_r_kind)*tem3**0.109_r_kind
              elseif (tem2 < -40.0_r_kind ) then
                 cloud(2)%effective_radius(k) =  (1250._r_kind/9.337_r_kind)*tem3**0.08_r_kind
              elseif (tem2 < -30.0_r_kind ) then
                 cloud(2)%effective_radius(k) =  (1250._r_kind/9.208_r_kind)*tem3**0.055_r_kind
              else
                 cloud(2)%effective_radius(k) =  (1250._r_kind/9.387_r_kind)*tem3**0.031_r_kind
              endif
           endif  

           cloud(1)%effective_radius(k)=max(zero, cloud(1)%effective_radius(k))
           cloud(2)%effective_radius(k)=max(zero, cloud(2)%effective_radius(k))

        end do
        cloud(1)%effective_variance(:) = two
        cloud(2)%effective_variance(:) = two

     else
        cloud(1)%effective_radius  (:) = zero
        cloud(2)%effective_radius  (:) = zero
        cloud(1)%effective_variance(:) = two
        cloud(2)%effective_variance(:) = two
     endif
  else ! Handle general case with arbitray number of clouds
       ! --------------------------------------------------
!    Loop over clouds ...
!    --------------------
     do n = 1, nc

!       Map Model cloud names into CRTM Cloud indices
!       ---------------------------------------------
        Cloud(n)%Type = CloudType_(cloud_name(jcloud(n)),hydrotype)

        if(icmask) then
           Cloud(n)%water_content(:) = cloud_cont(:,n)
           if (trim(cloud_name(jcloud(n))) == 'ql' .and. &
               present (ccw_cont)) then
              Cloud(n)%water_content(:) = Cloud(n)%water_content(:) &
                                        + ccw_cont(:) 
           end if
        else
           Cloud(n)%water_content(:) = zero
        endif

!       Calculate effective radius of given cloud type
!       ----------------------------------------------
        if(icmask) then
           if (regional .and. (.not. wrf_mass_regional)) then
              cloud(n)%Effective_Radius(:) = cloud_efr(:,n)
           else
             !cloud(n)%Effective_Radius(:) = EftSize_(cloud_name(jcloud(n))) 
              if ( imp_physics==11 .and. lprecip ) then
                 cloud(n)%Effective_Radius(:) = cloud_efr(:,n)
              else
                 do k = 1, km
                    if (cloud(n)%water_content(k) > 1.0e-6_r_kind) &
                    cloud(n)%Effective_Radius(k) = EftSize_(cloud_name(jcloud(n)))
                 enddo
              end if 
           end if
        else
           cloud(n)%Effective_Radius(:) = zero
        endif
        cloud(n)%effective_variance(:) = two

     enddo

  endif 
  end subroutine setCloud

  function CloudType_(name,hydrotype) Result(ctype)
    character(len=*), parameter :: myname = 'CloudType_'
    character(len=*) :: name  ! Model cloud name
    character(len=*) :: hydrotype(6)
    integer(i_kind)  :: ctype ! CRTM cloud type

    if ( trim(name) == 'ql' ) then
       SELECT CASE (hydrotype(1))
         CASE('WATER_CLOUD'); ctype = WATER_CLOUD
         CASE('LiquidSphere'); ctype = LiquidSphere  
       END SELECT
    else if ( trim(name) == 'qi' ) then
       SELECT CASE (hydrotype(2))
         CASE('ICE_CLOUD'); ctype = ICE_CLOUD
         CASE('IceSphere'); ctype = IceSphere
       END SELECT
    else if ( trim(name) == 'qh' ) then
       SELECT CASE (hydrotype(6))
         CASE('HAIL_CLOUD'); ctype = HAIL_CLOUD
         CASE('GemHail'); ctype = GemHail
         CASE('IconHail'); ctype = IconHail
       END SELECT
    else if ( trim(name) == 'qg' ) then
       SELECT CASE (hydrotype(5))
         CASE('GRAUPEL_CLOUD'); ctype = GRAUPEL_CLOUD
         CASE('GemGraupel'); ctype = GemGraupel
       END SELECT
    else if ( trim(name) == 'qr' ) then
       SELECT CASE (hydrotype(3))
         CASE('RAIN_CLOUD'); ctype = RAIN_CLOUD
         CASE('LiquidSphere'); ctype = LiquidSphere
       END SELECT
    else if ( trim(name) == 'qs' ) then
       SELECT CASE (hydrotype(4))
         CASE('SNOW_CLOUD'); ctype = SNOW_CLOUD
         CASE('PlateType1'); ctype = PlateType1
         CASE('ColumnType1'); ctype = ColumnType1
         CASE('SixBulletRosette'); ctype = SixBulletRosette
         CASE('Perpendicular4_BulletRosette'); ctype = Perpendicular4_BulletRosette
         CASE('Flat3_BulletRosette'); ctype = Flat3_BulletRosette
         CASE('SectorSnowflake'); ctype = SectorSnowflake
         CASE('EvansSnowAggregate'); ctype = EvansSnowAggregate
         CASE('EightColumnAggregate'); ctype = EightColumnAggregate
         CASE('LargePlateAggregate'); ctype = LargePlateAggregate
         CASE('LargeColumnAggregate'); ctype = LargeColumnAggregate
         CASE('LargeBlockAggregate'); ctype = LargeBlockAggregate
         CASE('IconSnow'); ctype = IconSnow
         CASE('GemSnow'); ctype = GemSnow
         CASE('IconCloudIce'); ctype = IconCloudIce
       END SELECT
    else
       call die(myname,"cannot recognize cloud name <"//trim(name)//">")
    end if

  end function CloudType_

  function EftSize_(name) Result(csize)                
    character(len=*), parameter :: myname = 'EftSize_'  
    character(len=*) :: name  ! Model cloud name
    real(r_kind)     :: csize ! CRTM cloud type

! Note: Values below from Tom Auligne
    if ( trim(name) == 'ql' ) then
       csize = 10.0_r_kind  ! in micron
    else if ( trim(name) == 'qi' ) then
       csize = 30.0_r_kind
    else if ( trim(name) == 'qh' ) then
    !  csize = zero ! RT: can somebody fill this in? 
       csize = 1000.0_r_kind
    else if ( trim(name) == 'qg' ) then
       csize = 600.0_r_kind
    else if ( trim(name) == 'qr' ) then
       csize = 300.0_r_kind
    else if ( trim(name) == 'qs' ) then
       csize = 600.0_r_kind

    else
       call die(myname,"cannot recognize cloud name <"//trim(name)//">")
    end if

  end function EftSize_    

end module set_crtm_cloudmod
