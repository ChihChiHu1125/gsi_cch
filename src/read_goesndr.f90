subroutine read_goesndr(mype,val_goes,ithin,rmesh,jsatid,infile,&
     lunout,obstype,nread,ndata,nodata,twind,gstime,sis,&
     mype_root,mype_sub,npe_sub,mpi_comm_sub)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_goesndr                   read goes sounder data
!   prgmmr: yang             org: np23                date: 1998-05-15
!
! abstract:  This routine reads GOES sounder radiance (brightness
!            temperature) files.  Optionally, the data are thinned to 
!            a specified resolution using simple quality control checks.
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   1998-05-15 weiyu yang
!   1999-08-24 derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-05-28 kleist - update subroutine call
!   2004-06-16 treadon - update documentation
!   2004-07-23 derber - make changes to eliminate obs. earlier in thinning
!   2004-07-29 treadon - add only to module use, add intent in/out
!   2005-01-26 derber - land/sea determination and weighting for data selection
!   2005-07-08 derber - clean up, fix bugs, and improve observation selection
!   2005-09-08  derber - modify to use input group time window
!   2005-09-28  derber - modify to produce consistent surface info
!   2005-10-17  treadon - add grid and earth relative obs location to output file
!   2005-10-18  treadon - remove array obs_load and call to sumload
!   2005-11-22  derber  - include mean in bias correction
!   2005-11-29  parrish - modify getsfc to work for different regional options
!   2006-02-01  parrish - remove getsfc (different version called now in read_obs)
!   2006-02-03  derber  - modify for new obs control and obs count
!   2006-03-07  derber  - combine reading of 1x1 and 5x5 (prepbufr) files
!   2006-04-27  derber - clean up code
!   2006-05-19  eliu    - add logic to reset relative weight when all channels not used
!   2006-07-28  derber  - add solar and satellite azimuth angles remove isflg from output
!                       - add ability to read g9,g11,g13
!   2006-09-21  treadon - replace serial bufr i/o with parallel bufr i/o (mpi_io)
!   2007-03-01  tremolet - measure time from beginning of assimilation window
!   2008-04-21  safford - rm unused vars
!   2009-04-18  woollen - improve mpi_io interface with bufrlib routines
!   2009-04-21  derber  - add ithin to call to makegrids
!
!   input argument list:
!     mype     - mpi task id
!     val_goes - weighting factor applied to super obs
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window(hours)
!     gstime   - guess time
!     sis      - sensor/instrument/satellite indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!
!   output argument list:
!     nread    - number of BUFR GOES sounder observations read
!     ndata    - number of BUFR GOES sounder profiles retained for further processing
!     nodata   - number of BUFR GOES sounder observations retained for further processing
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
             checkob,finalcheck,score_crit
  use radinfo, only: cbias,newchn,predx,iuse_rad,jpch_rad,nusis
  use gridmod, only: diagnostic_reg,nlat,nlon,regional,tll2xy,txy2ll,rlats,rlons
  use constants, only: deg2rad,zero,one,izero,ione,rad2deg
  use obsmod, only: iadate,offtime_data
  use gsi_4dvar, only: l4dvar,idmodel,iadatebgn,iadateend,time_4dvar,iwinbgn,winlen

  implicit none

! Declare passed variables
  character(len=*),intent(in):: infile,obstype,jsatid
  character(len=*),intent(in):: sis
  integer(i_kind),intent(in):: mype,lunout,ithin
  integer(i_kind),intent(inout):: ndata,nodata,nread
  real(r_kind),intent(in):: rmesh,twind,gstime
  real(r_kind),intent(inout):: val_goes
  integer(i_kind)  ,intent(in) :: mype_root
  integer(i_kind)  ,intent(in) :: mype_sub
  integer(i_kind)  ,intent(in) :: npe_sub
  integer(i_kind)  ,intent(in) :: mpi_comm_sub


! Declare local parameters
  integer(i_kind),parameter:: maxinfo=33
  integer(i_kind),parameter:: mfov=25   ! maximum number of fovs (currently 5x5)

  real(r_kind),parameter:: r360=360.0_r_kind
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind
  real(r_kind),parameter:: R60=60.0_r_kind
  character(80),parameter:: hdstr = &
               'CLON CLAT ELEV SOEL BEARAZ SOLAZI SAID DINU YEAR MNTH DAYS HOUR MINU SECO ACAV' 
  character(80),parameter:: hdstr5 = &
               'XOB YOB ELEV SOEL BEARAZ SOLAZI SAID TYP ACAV DHR SID '
  character(80),parameter:: rbstr = 'TMBR'

! Declare local variables
  logical outside,iuse,g5x5,assim

  character(8)  subset

  integer(i_kind) kx,levs,ldetect
  integer(i_kind) lnbufr,nchanl,nreal,iret,ksatid,lsatid
  integer(i_kind) idate,iout
  integer(i_kind) ilat,ilon,isflg,idomsfc
  integer(i_kind) itx,k,i,itt,iskip,l,ifov,n
  integer(i_kind) ichan8,ich8
  integer(i_kind) nele,iscan,nmind
  integer(i_kind) ntest,ireadsb,ireadmg,irec,isub,next
  integer(i_kind)::  file_handle,ierror,nblocks
  integer(i_kind),dimension(5):: idate5

  real(r_kind) dlon,dlat,timedif,emiss,sfcr
  real(r_kind) dlon_earth,dlat_earth
  real(r_kind) ch8,sstime
  real(r_kind) pred,crit1,tdiff,dist1,toff,t4dv
  real(r_kind) disterr,disterrmax,dlon00,dlat00,r01

  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10

  real(r_kind),allocatable,dimension(:,:):: data_all

  real(r_double),dimension(15):: hdr
  real(r_double),dimension(18):: grad



!**************************************************************************

! Start routine here.  Set constants.  Initialize variables
  lnbufr = 10
  disterrmax=zero
  ntest  = 0
  nreal  = maxinfo
  ich8   = 8        !channel 8
  ndata  = 0
  nchanl = 18
  ifov = -999
  r01 = 0.01_r_kind

  ilon=3
  ilat=4

  rlndsea(0) = 0._r_kind
  rlndsea(1) = 15._r_kind
  rlndsea(2) = 10._r_kind
  rlndsea(3) = 15._r_kind
  rlndsea(4) = 30._r_kind

! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_goes=zero


! Make thinning grids
  call makegrids(rmesh,ithin)
 
!  check to see if prepbufr file

  g5x5 = jsatid == 'g08_prep' .or. jsatid == 'g09_prep' .or.     &
         jsatid == 'g10_prep' .or. jsatid == 'g11_prep' .or.     &
         jsatid == 'g12_prep' .or. jsatid == 'g13_prep'

  if(g5x5)then
       if(jsatid=='g08_prep')lsatid=252
       if(jsatid=='g09_prep')lsatid=253
       if(jsatid=='g10_prep')lsatid=254
       if(jsatid=='g11_prep')lsatid=255
       if(jsatid=='g12_prep')lsatid=256
       if(jsatid=='g13_prep')lsatid=257
   else
       if(jsatid=='g08')lsatid=252
       if(jsatid=='g09')lsatid=253
       if(jsatid=='g10')lsatid=254
       if(jsatid=='g11')lsatid=255
       if(jsatid=='g12')lsatid=256
       if(jsatid=='g13')lsatid=257
       if(obstype == 'sndrd1')ldetect = 1
       if(obstype == 'sndrd2')ldetect = 2
       if(obstype == 'sndrd3')ldetect = 3
       if(obstype == 'sndrd4')ldetect = 4
  end if

! Set array index for surface-sensing channels
  ichan8  = newchn(sis, ich8)


! Open then read the bufr data
  open(lnbufr,file=infile,form='unformatted')
  call openbf(lnbufr,'IN',lnbufr)
  call datelen(10)

! Time offset
  call time_4dvar(idate,toff)

! Allocate arrays to hold data
  nele=nreal+nchanl
  allocate(data_all(nele,itxmax))

! Big loop to read data file
  next=mype_sub+1
  do while(ireadmg(lnbufr,subset,idate)>=0)
  call ufbcnt(lnbufr,irec,isub)
  if(irec<>next)cycle; next=next+npe_sub
  read_loop: do while (ireadsb(lnbufr)==0)

!    Extract type, date, and location information
     if(g5x5)then
!    Prepbufr file
       call ufbint(lnbufr,hdr,11,1,iret,hdstr5)
       kx = hdr(8)
       if(kx /= 164 .and. kx /= 165 .and. kx /= 174 .and. kx /= 175)cycle read_loop
!      If not goes data over ocean , read next bufr record
!      if(kx /= 174 .and. kx /= 175)cycle read_loop

       ksatid=hdr(7)
!      if not proper satellite read next bufr record
       if (ksatid /= lsatid) cycle read_loop

!      Extract number of averaged FOVS
       ifov = hdr(9) ! number of averaged FOVS 
       if(ifov <= 3) cycle read_loop
!      Extract obs time difference. 
       tdiff=hdr(10)  ! relative obs time in hours
       t4dv=toff+tdiff
     else
!      GOES 1x1 or 5x5 file
       call ufbint(lnbufr,hdr,15,1,iret,hdstr)

       ksatid=hdr(7)   !bufr satellite id
!      if not proper satellite/detector read next bufr record
       if (ksatid /=lsatid) cycle read_loop
       if(obstype /= 'sndr')then
         if(ldetect /= nint(hdr(8)))cycle read_loop
       end if

!!     ifov = hdr(15) ! number of averaged FOVS 
       ifov = nint(hdr(15)) ! number of averaged FOVS

       if(ifov < mfov .and. ifov > 0)then
        if(ifov <= 3) cycle read_loop
       end if

!      Extract obs time.  If not within analysis window, skip obs
!      Extract date information.  If time outside window, skip this obs
       idate5(1) = hdr(9) !year
       idate5(2) = hdr(10) !month
       idate5(3) = hdr(11) !day
       idate5(4) = hdr(12) !hour
       idate5(5) = hdr(13) !minute
       call w3fs21(idate5,nmind)
       sstime=real(nmind,r_kind) + hdr(14)/R60
       tdiff=(sstime-gstime)/R60
       t4dv=(real(nmind-iwinbgn,r_kind) + hdr(14)/r60)/r60
     end if

!    If not within analysis window, skip obs
     if (l4dvar) then
       if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
     else
       if (abs(tdiff)>twind) cycle read_loop
     endif

!       Convert obs location to radians
     if (hdr(1)>=r360) hdr(1)=hdr(1)-r360
     if (hdr(1)< zero) hdr(1)=hdr(1)+r360

     dlon_earth = hdr(1)*deg2rad   !convert degrees to radians
     dlat_earth = hdr(2)*deg2rad

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if(diagnostic_reg) then
           call txy2ll(dlon,dlat,dlon00,dlat00)
           ntest=ntest+1
           disterr=acos(sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00)))*rad2deg
           disterrmax=max(disterrmax,disterr)
        end if
      
!          Check to see if in domain
        if(outside) cycle read_loop

     else
        dlon = dlon_earth 
        dlat = dlat_earth 
        call grdcrd(dlat,1,rlats,nlat,1)
        call grdcrd(dlon,1,rlons,nlon,1)
     endif


!    Set common predictor parameters

     if (l4dvar) then
       timedif = 0.0
     else
       timedif = 6.0_r_kind*abs(tdiff)        ! range:  0 to 18
     endif

     nread=nread+nchanl

     crit1=0.01_r_kind+timedif
     if(ifov < mfov .and. ifov > 0)then
       crit1=crit1+2.0_r_kind*float(mfov-ifov)
     end if

     call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis)
     if(.not. iuse)cycle read_loop

!    Increment goes sounder data counter
!    Extract brightness temperatures
     call ufbint(lnbufr,grad,1,18,levs,rbstr)

     iskip = 0
     do l=1,nchanl

        if( grad(l) < tbmin .or. grad(l) > tbmax )then
           iskip = iskip + 1
           if(l == ich8)iskip = nchanl
        endif
     end do

     if( iskip >= nchanl )cycle read_loop

!    "Score" observation.   We use this information to id "best" obs.

!    Locate the observation on the analysis grid.  Get sst and land/sea/ice
!    mask.  

!     isflg    - surface flag
!                0 sea
!                1 land
!                2 sea ice
!                3 snow
!                4 mixed 

     call deter_sfc_type(dlat_earth,dlon_earth,t4dv,isflg,tsavg)

!      If not goes data over ocean , read next bufr record
     if(isflg /= 0) cycle read_loop

     crit1 = crit1 + rlndsea(isflg)  
     call checkob(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle read_loop

!    Set data quality predictor
     iscan   = nint(hdr(3))+0.001_r_kind   ! "scan" position
     ch8     = grad(ich8) -cbias(iscan,ichan8) - r01*predx(1,ichan8)
     emiss=0.992-0.013*(hdr(3)/65.)**3.5-0.026*(hdr(3)/65.)**7.0
     pred = abs(ch8-tsavg*emiss)

!    Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"
     crit1 = crit1+10.0_r_kind*pred  

     call finalcheck(dist1,crit1,itx,iuse)
     if(.not. iuse)cycle read_loop

!    Transfer observation location and other data to local arrays

     data_all(1,itx) = ksatid                       ! satellite id
     data_all(2,itx) = t4dv                         ! time (hours)
     data_all(3,itx) = dlon                         ! grid relative longitude
     data_all(4,itx) = dlat                         ! grid relative latitude
     data_all(5,itx) = hdr(3)*deg2rad               ! satellite zenith angle
     data_all(6,itx) = hdr(5)                       ! satellite zenith angle
     data_all(7,itx) = zero                         ! local zenith angle
     data_all(8,itx) = iscan                        ! "scan" position
     data_all(9,itx) = hdr(4)                       ! solar zenith angle
     data_all(10,itx)= hdr(6)                       ! solar zenith angle

     data_all(30,itx)= dlon_earth                   ! earth relative longitude (rad)
     data_all(31,itx)= dlat_earth                   ! earth relative latitude (rad)

     data_all(32,itx)= val_goes
     data_all(33,itx)= itt

     do k=1,nchanl
        data_all(k+nreal,itx)=grad(k)
     end do


  end do read_loop
  end do
  call closbf(lnbufr)


! If multiple tasks read input bufr file, allow each tasks to write out
! information it retained and then let single task merge files together

  call combine_radobs(mype,mype_sub,mype_root,npe_sub,mpi_comm_sub,&
       nele,itxmax,nread,ndata,data_all,score_crit)


! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then

!    Identify "bad" observation (unreasonable brightness temperatures).
!    Update superobs sum according to observation location

     do n=1,ndata
        do i=1,nchanl
           if(data_all(i+nreal,n) > tbmin .and. &
                data_all(i+nreal,n) < tbmax)nodata=nodata+1
        end do
        itt=nint(data_all(nreal,n))
        super_val(itt)=super_val(itt)+val_goes
        tdiff = data_all(2,n)                ! time (hours)
        dlon=data_all(3,n)                   ! grid relative longitude
        dlat=data_all(4,n)                   ! grid relative latitude
        dlon_earth = data_all(30,n)          ! earth relative longitude (rad)
        dlat_earth = data_all(31,n)          ! earth relative latitude (rad)

        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,tdiff,isflg,idomsfc,sfcpct, &
                ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
        data_all(11,n) = sfcpct(0)           ! sea percentage of
        data_all(12,n) = sfcpct(1)           ! land percentage
        data_all(13,n) = sfcpct(2)           ! sea ice percentage
        data_all(14,n) = sfcpct(3)           ! snow percentage
        data_all(15,n)= ts(0)                ! ocean skin temperature
        data_all(16,n)= ts(1)                ! land skin temperature
        data_all(17,n)= ts(2)                ! ice skin temperature
        data_all(18,n)= ts(3)                ! snow skin temperature
        data_all(19,n)= tsavg                ! average skin temperature
        data_all(20,n)= vty                  ! vegetation type
        data_all(21,n)= vfr                  ! vegetation fraction
        data_all(22,n)= sty                  ! soil type
        data_all(23,n)= stp                  ! soil temperature
        data_all(24,n)= sm                   ! soil moisture
        data_all(25,n)= sn                   ! snow depth
        data_all(26,n)= zz                   ! surface height
        data_all(27,n)= idomsfc + 0.001      ! dominate surface type
        data_all(28,n)= sfcr                 ! surface roughness
        data_all(29,n)= ff10                 ! ten meter wind factor
        data_all(30,n)= data_all(30,n)*rad2deg  ! earth relative longitude (degrees)
        data_all(31,n)= data_all(31,n)*rad2deg  ! earth relative latitude (degrees)

     end do

!    Write final set of "best" observations to output file
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  
  endif

  deallocate(data_all) ! Deallocate data arrays
  call destroygrids    ! Deallocate satthin arrays

  if(diagnostic_reg .and. ntest>0 .and. mype_sub==mype_root) &
       write(6,*)'READ_GOESNDR:  mype,ntest,disterrmax=',&
       mype,ntest,disterrmax

  return

end subroutine read_goesndr