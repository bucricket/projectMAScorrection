!
module rttov_zutility
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  ! Description:
  !   containing routines to load geomagnetic field Lookup table, compute
  !   geomagnetic field and its angles relative to the wave propagation
  !   direction.
  !
  ! History:
  ! Version   Date     Comment
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0   11/07/2007  created by Y. Han, NOAA/JCSDA
  !
#include "throw.h"
  USE PARKIND1  ,ONLY : JPIM, JPRB, JPLM
  Implicit None

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: load_bfield_lut   ! function to load Geomagnetic field LUT
  PUBLIC :: compute_bfield    ! subroutine to calculate the Geomagnetic field 
                              ! using a LUT and two cosine angles of the field
                              ! relative to the wave propagation direction k.
  PUBLIC :: compute_kb_angles ! subroutine to compute two cosine angles of the 
                              ! geomagnetic field relative to the wave 
                              ! propagation direction k.


  INTERFACE compute_bfield
    MODULE PROCEDURE Compute_bfield_F1
    MODULE PROCEDURE compute_bfield_F2
    MODULE PROCEDURE compute_bfield_F3
  END INTERFACE compute_bfield
  
  ! Array for Earth's magnetic field
  INTEGER,      PARAMETER   :: n_Lat = 91, n_Lon = 181               
  INTEGER, SAVE             :: BField(3, n_lat, n_lon)               

  Real(JPRB), parameter :: DEGREES_TO_RADIANS  = 3.141592653589793238462643_JPRB/180.0_JPRB

CONTAINS

  !--------------------------------------------------------------------------------
  !
  ! NAME:
  !       load_bfield_lut
  !
  ! PURPOSE:
  !       Function to the geomagnetic filed LUT
  !
  ! CALLING SEQUENCE:
  !   errorstatus = load_bfield_lut(filename_LUT)
  !
  ! INPUT ARGUMENT:
  ! 
  !   filename_LUT:       file name for the file containing the LUT of
  !                       the Earth magnetic field.
  !                       UNITS:      N/A
  !                       TYPE:       character string
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !----------------------------------------------------------------------------
  Function load_bfield_lut(filename_LUT) Result(ERR)
    Character(*), Intent(in)  :: filename_LUT

    Integer(jpim)             :: ERR

    ! local
    Character(*), parameter :: NameOfRoutine = 'rttov_load_bfield_lut '
    Integer(jpim)      :: file_id
    Logical(Kind=jplm) :: Existence, Is_Open
    Integer(jpim)      :: i, j, iskip
TRY
    ERR = 0
    ! Find a free logical unit number for file access
    file_id = 9
    ! Start open loop for Lun Search
    Lun_Search: DO
      file_id = file_id + 1
      INQUIRE( UNIT = file_id, EXIST = Existence )
      IF ( .NOT. Existence ) THEN
        file_id = -1
        EXIT Lun_Search
      END IF
      INQUIRE( UNIT = file_id, OPENED = Is_Open )
      IF ( .NOT. Is_Open ) EXIT Lun_Search
    END DO Lun_Search

    Open( file_id, FILE   = filename_LUT, &
          & STATUS = 'OLD', &
          & IOSTAT = ERR ) 
    THROWM(ERR .NE. 0,"Error opening "//TRIM(filename_LUT))

    READ(file_id, *, iostat = ERR)iskip, iskip, &
                      ((BField(:, i, j), j=1,n_Lon), i=1,n_Lat)  
!write(6,*) 'Bfield(1,1,1) = ', Bfield(1,1,1)

    THROWM(ERR .NE. 0,"Error reading "//TRIM(filename_LUT))
              
    Close( unit = file_id )

CATCH
  End Function load_bfield_lut

  !--------------------------------------------------------------------------------
  !
  ! NAME:
  !       compute_bfield
  !
  ! PURPOSE:
  !       Subroutine to calculate the Geomagnetic field using a LUT and 
  !       two cosine angles of the field relative to the wave propagation
  !       direction k. 
  !
  ! CALLING SEQUENCE:
  !
  !        CALL compute_bfield(latitude, longitude, & ! Inputs 
  !                            Bx, By, Bz, Be)        ! Outputs
  !
  !   OR
  !
  !        CALL compute_bfield(latitude,          &   ! input
  !                            longitude,         &   ! input   
  !                            sensor_zenang,     &   ! input   
  !                            sensor_aziang,     &   ! input     
  !                            Be,                &   ! output  
  !                            cos_bkang,         &   ! output  
  !                            cos_baziang)           ! output   
  !  
  !   OR       
  !
  !        CALL compute_bfield(latitude,               &   ! input
  !                            longitude,              &   ! input   
  !                            sensor_zenang,          &   ! input   
  !                            sensor_relative_aziang, &   ! input   
  !                            Julian_day,             &   ! input   
  !                            utc_time,               &   ! input   
  !                            Be,                     &   ! output  
  !                            cos_bkang,              &   ! output  
  !                            cos_baziang)                ! output  
  !
  !
  ! INPUT ARGUMENTS:
  !
  !
  !      latitude :       Latitude, -90 - 90 (90 North Pole; -90 - South Pole)
  !                       UNITS:      Degree
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !      longitude:       longitude, 0 - 360 East
  !                       UNITS:      Degree
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !   sensor_zenang:      sensor zenith angle
  !                       UNITS:      Degree
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !   sensor_aziang:      sensor azimuth angle (starts from North,
  !                                             positive clockwise)
  !                       UNITS:      Degree
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  ! sensor_relative_aziang: sensor azimuth angle relative to the sun azimuth
  !                         angle
  !    Sun_azimuth_angle - from north, East positive
  !    sensor_relative_aziang = sun_azimuth_angle - Sensor_aziang
  !                       UNITS:      degree
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !        Julian_day:    Julian_Day 1=Jan 1, 365=Dec31 (366 leap year) 
  !                       UNITS:      day
  !                       TYPE:       integer(JPIM)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !          utc_time:    Universal_Time 0.00-23.999,(GMT,Z time) 
  !                       UNITS:      hour
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  ! OUTPUT ARGUMENTS:
  !       Bx:             Magetic field East component 
  !                       UNITS:      Gauss
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       By:             Magetic field North component 
  !                       UNITS:      Gauss
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       Bz:             Magetic field zenith component (positive upward) 
  !                       UNITS:      Gauss
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       Be:             Magetic field strength (sqrt(BxBx + ByBy + BzBz)) 
  !                       UNITS:      Gauss
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !     cos_bkang:        cosine of the angle between the magnetic field Be 
  !                       vector and the wave propagation direction k.
  !                       UNITS:      N/A
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !   cos_baziang:        cosine of the azimuth angle of the Be vector in the
  !                       (v, h, k) coordinates system, where v, h and k comprise
  !                       a right-hand orthogonal system, similar to the (x, y, z)
  !                       Catesian coordinates. The h vector is normal to the
  !                       plane containing the k and z vectors, where k points
  !                       to the wave propagation direction and z points 
  !                       to the zenith. h = (z cross k)/|z cross k|. The
  !                       azimuth angle is the angle on the (v, h) plane
  !                       from the positive v axis to the projected line of the
  !                       Be vector on this plane, positive counterclockwise.
  !                       UNITS:      N/A
  !                       TYPE:       real(JPRB)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !--------------------------------------------------------------------------------

  SUBROUTINE compute_bfield_F1(latitude, longitude, & ! Inputs 
                               Bx, By, Bz, Be)        ! Outputs
    REAL(JPRB), INTENT(IN)  :: latitude, longitude
    REAL(JPRB), INTENT(OUT) :: Bx, By, Bz, Be

    ! Local
    REAL(JPRB)                :: lat, lon
    REAL(JPRB)                :: x2, w1_lat, w1_lon
    INTEGER                   :: idx_lat, idx_lon
    REAL(JPRB), PARAMETER     :: dlat = 2.0_JPRB, dlon = 2.0_JPRB  

    lat = 90.0_JPRB - latitude  ! lat: 0 - 180 (0 North Pole; 180 South Pole)
    lon = longitude             ! lon: 0 - 360 East
    IF(lon < 0.0_JPRB)lon = lon + 360.0_JPRB
    
    idx_lat = INT(lat/dlat)+1
    IF(idx_lat >= n_Lat)idx_lat = n_lat-1
    idx_lon = INT(lon/dlat)+1
    IF(idx_lon >= n_Lon)idx_lon = n_lon-1

    x2 = REAL(idx_lat, JPRB)*dlat
    w1_lat = (x2 - lat)/dlat
  
    x2 = REAL(idx_lon, JPRB)*dlat
    w1_lon = (x2 - lon)/dlat

    Bx = BField_Component(1, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)
    By = BField_Component(2, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)
    Bz = BField_Component(3, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)
                    
    Be = SQRT(Bx*Bx+By*By+Bz*Bz)

  CONTAINS

    FUNCTION BField_Component(comp, Bfield, w1_lat, w1_lon, &                    
                              idx_lat, idx_lon) Result(B) 

     INTEGER,      INTENT(IN)  :: comp
     INTEGER,      INTENT(IN)  :: BField(:,:,:)
     REAL(JPRB), INTENT(IN)  :: w1_lat, w1_lon
     INTEGER,      INTENT(IN)  :: idx_lat, idx_lon

     REAL(JPRB) :: B                      

     REAL(JPRB), PARAMETER :: Scale = 0.001_JPRB
     REAL(JPRB)            :: w2_lat, w2_lon              

      w2_lat = 1.0_JPRB - w1_lat                                                      
      w2_lon = 1.0_JPRB - w1_lon                                                      
      B = (w1_lon*(w1_lat*REAL(BField(comp,idx_lat,   idx_lon),  JPRB) + &     
                   w2_lat*REAL(BField(comp,idx_lat+1, idx_lon),  JPRB))  &     
         + w2_lon*(w1_lat*REAL(BField(comp,idx_lat,   idx_lon+1),JPRB) + &     
                   w2_lat*REAL(BField(comp,idx_lat+1, idx_lon+1),JPRB)))*Scale 

     END FUNCTION BField_Component
        
  END SUBROUTINE compute_bfield_F1

  Subroutine compute_bfield_F2(latitude, longitude, sensor_zenang, sensor_aziang, & 
                               Be, cos_bkang, cos_baziang)

    !subroutine arguments:
    Real(JPRB), Intent(in)             :: latitude
    Real(JPRB), Intent(in)             :: longitude
    Real(JPRB), Intent(in)             :: sensor_zenang
    Real(JPRB), Intent(in)             :: sensor_aziang
    Real(JPRB), Intent(out)            :: Be
    Real(JPRB), Intent(out)            :: cos_bkang
    Real(JPRB), Intent(out)            :: cos_baziang

    ! local variable
    Real(JPRB)    :: Bx, By, Bz
  
    ! get Earth magnetic filed from LUT
    Call compute_bfield_F1(latitude, longitude,  & ! inputs
                           Bx, By, Bz, Be)         ! outputs

    ! compute the cosines of the angle between the magnetic field Be and
    ! propagation direction and Be's azimuth angle
    Call Compute_kb_Angles(Bx, By, Bz, sensor_zenang, sensor_aziang,  &  ! Input
                           cos_bkang, cos_baziang)                       ! Output

  End Subroutine compute_bfield_F2

  Subroutine compute_bfield_F3(latitude, longitude, sensor_zenang, &
                               sensor_relative_aziang, julian_day, utc_time, &
                               Be, cos_bkang, cos_baziang)
    !subroutine arguments:
    Real(JPRB), Intent(in)             :: latitude
    Real(JPRB), Intent(in)             :: longitude
    Real(JPRB), Intent(in)             :: sensor_zenang
    Real(JPRB), Intent(in)             :: sensor_relative_aziang
    Integer(jpim), Intent(in)          :: julian_day
    Real(JPRB), Intent(in)             :: utc_time
    Real(JPRB), Intent(out)            :: Be
    Real(JPRB), Intent(out)            :: cos_bkang
    Real(JPRB), Intent(out)            :: cos_baziang

    ! Local     
    Real(JPRB)    :: Bx, By, Bz, Solar_ZA, Solar_AZ, lat, lon, sensor_aziang

    ! compute the cosines of the angle between the magnetic field Be and
    ! propagation direction and Be's azimuth angle
    Call compute_bfield_F1(latitude, longitude,  & ! inputs
                            Bx, By, Bz, Be)         ! outputs

    lat = latitude
    lon = longitude
    If(lon > 180.0_JPRB)lon = lon - 360.0_JPRB
    ! Compute Solar azimuth angle
    Call Solar_ZA_AZ(lat,lon, Real(julian_day, JPRB), utc_time, &                     
                     Solar_ZA, Solar_AZ)
    ! Compute satellite azimuth angle (zero north, east +90)
    ! NB relative azimuth is measured as (sun_aziang - sensor_aziang)
    sensor_aziang = solar_AZ - sensor_relative_aziang

    ! compute for cos_bkangle, cos_baziangle
    Call Compute_kb_Angles(Bx, By, Bz, sensor_zenang, sensor_aziang,  &  ! Input
                           cos_bkang, cos_baziang)                       ! Output


  End Subroutine compute_bfield_F3


!--------------------------------------------------------------------------------
!
! NAME:
!       Compute_kb_Angles
!
! PURPOSE:
!       Subroutine to calculate the cosine of the angle between the Geomagnetic 
!       field (Be) and wave propagation direction (k) and the cosine of the 
!       azimuth angle of the Be vector in the (v, h, k) coordinates system (see
!       more detailed description below)
!         
!
! CALLING SEQUENCE:
!      CALL Compute_kb_Angles(Bx, By, Bz,   &                       ! Input
!                             sensor_zenang, sensor_aziang, &       ! Input   
!                             cos_bk_Angle, cos_baziang)            ! Output  
! INPUT ARGUMENTS:
!
!       Bx:             Magetic field East component 
!                       UNITS:      Gauss
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       By:             Magetic field North component 
!                       UNITS:      Gauss
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Bz:             Magetic field zenith component (positive upward) 
!                       UNITS:      Gauss
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!   sensor_zenang :     sensor zenith angle
!                       UNITS:      Degree
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!    
!   sensor_aziang :     sensor zenith angle defined as the
!                       angle from North, positive clockwise.
!                       UNITS:      Degree
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!     cos_bkang:        cosine of the angle between the magnetic field Be   
!                       vector and the wave propagation direction k.       
!                       UNITS:      N/A                                    
!                       TYPE:       real(JPRB)                             
!                       DIMENSION:  Scalar                                 
!                       ATTRIBUTES: INTENT(OUT)                            
!
!   cos_baziang:        cosine of the azimuth angle of the Be vector in the       
!                       (v, h, k) coordinates system, where v, h and k comprise   
!                       a right-hand orthogonal system, similar to the (x, y, z)  
!                       Catesian coordinates. The h vector is normal to the       
!                       plane containing the k and z vectors, where k points      
!                       to the wave propagation direction and z points            
!                       to the zenith. h = (z cross k)/|z cross k|. The           
!                       azimuth angle is the angle on the (v, h) plane            
!                       from the positive v axis to the projected line of the     
!                       Be vector on this plane, positive counterclockwise.       
!                       UNITS:      N/A
!                       TYPE:       Real(JPRB)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE Compute_kb_Angles(Bx, By, Bz,   &                       ! Input
                               sensor_zenang, sensor_aziang, &       ! Input
                               cos_bkang, cos_baziang)   ! Output
    REAL(JPRB), INTENT(IN)  :: Bx, By, Bz
    REAL(JPRB), INTENT(IN)  :: sensor_zenang, sensor_aziang
    REAL(JPRB), INTENT(OUT) :: cos_bkang, cos_baziang

    ! Local
    REAL(JPRB) :: B, B_v, B_h, B_p, kx, ky, kz, azi, &
                  SIN_SenZA, COS_SenAZ, SIN_SenAZ, COS_SenZA 

    ! Input azimuth is zero north, east +90. The code below expects zero east, north +90
    azi = 90._JPRB - sensor_aziang

    SIN_SenZA = SIN(sensor_zenang*DEGREES_TO_RADIANS)
    COS_SenZA = COS(sensor_zenang*DEGREES_TO_RADIANS)
    SIN_SenAZ = SIN(azi*DEGREES_TO_RADIANS)
    COS_SenAZ = COS(azi*DEGREES_TO_RADIANS)

    ! compute k directional vector from satellite's zenith and azimuth angles 
    kx = SIN_SenZA*COS_SenAZ
    ky = SIN_SenZA*SIN_SenAZ
    kz = COS_SenZA

    ! compute consine of the angle between the magnetic field B and k
    B = SQRT(bx*bx+by*by+bz*bz)
    cos_bkang = (kx*bx + ky*by + kz*bz)/B
             
    ! Project the B vector on the V and H plane: B_v - B component on the V
    ! axis; B_h - B component on the H axis. 
    B_v = bx*kx*kz + by*ky*kz - bz*(ky*ky + kx*kx) ;  
    B_h = -bx*ky + by*kx                            

    ! compute the cosine of the azimuth angle
    B_p = SQRT(B_v*B_v + B_h*B_h)
    If(B_p /= 0.0_JPRB)Then
      cos_baziang = B_v / B_p
    Else
      cos_baziang = 0.0   ! not defined (take an arbitrary value)
    Endif  

  END SUBROUTINE Compute_kb_Angles
 
  Subroutine Solar_ZA_Az(latitude,        &      
                         longitude,       &      
                         julian_day,      &      
                         universal_time,  &      
                         ZA,              &      
                         Az)                     
!
!************************************************************************
!*                                                                      *
!*	Module Name:	Solar_Az_Za	                                *
!*                                                                      *
!*	Language:	Fortran 	   Library:	                *
!*	Version.Rev:	1.0  22 Feb 91	Programmer:	Kleespies       *
!*			1.1  28 Feb 91			Kleespies       *
!*			     Put equation of time into hour angle.      *
!*
!*      updated to f95, 4, September 2007, Y. Han                       *
!*                                                                      *
!*	Calling Seq:		Call Solar_Az_ZA(                       *
!*     &				latitude,                       *
!*     &				longitude,                      *
!*     &				julian_day,                     *
!*     &				universal_time,	                *
!*     &				ZA,	                        *
!*     &				Az)	                        *
!*                                                                      *
!*                                                                      *
!*	Description:	Computes solar azimuth and zenith angle for a   *
!*		given place and time.  Solar azimuth is the angle       *
!*		positive east of north form a point to the center of    *
!*		the solar disc. Solar zenith angle is the angle         *
!*		in degrees from zenith to the center of the solar disk. *
!*		The algorithms are taken from "Introduction to Solar    *
!*		Radiation" by Muhamad Iqbal, Academic Press, 1983.      *
!*		Note that lat=0,lon=0, 21Mar 12z will not give sun      *
!*		overhead, because that is not the way things work.      *
!*                                                                      *
!*	Input Args:	R*4  Latitude, +NH, -SH degrees	                *
!*			R*4  Longitude,-W, +E                           *
!*			R*4  Julian_Day 1=Jan 1, 365=Dec31 (366 leap yr)*
!*			R*4  Universal_Time 0.00-23.99,(GMT,Z time)     *
!*                                                                      *
!*	Output Args:	R*4  ZA	Solar Zenith Angle                      *
!*			R*4  AZ	Solar Azmuth Angle                      *
!*                                                                      *
!*	Common Blks:	none                                            *
!*	Include:	none                                            *
!*	Externals:	none                                            *
!*	Data Files:	none                                            *
!*                                                                      *
!*	Restrictions:	Accurate to within .1 deg.                      *
!*			No checking made to the validity of input       *
!*			parameters.                                     *
!*			Since solar zenith angle is a conic angle,      *
!*			it has no sign.                                 *
!*			No correction made for refraction.              *
!*                                                                      *
!*	Error Codes:	none                                            *
!*                                                                      *
!*	Error Messages:	                                                *
!*                                                                      *
!************************************************************************
!
        implicit none

        real(JPRB), intent(in)  :: latitude,longitude,julian_day,universal_time
        real(JPRB), intent(out) :: ZA, Az

        real(JPRB) :: local_sun_time,solar_elevation,equation_of_time
        real(JPRB) :: cosza,cosaz
        real(JPRB) :: hour_angle,day_angle
        real(JPRB) :: solar_declination
        real(JPRB) :: rlatitude
        real(JPRB), parameter :: DEGREES_TO_RADIANS  = 3.141592653589793238462643_JPRB/180.0_JPRB
        real(JPRB), parameter :: one_eighty_over_pi  = 1.0/DEGREES_TO_RADIANS
        real(JPRB), parameter :: threesixty_over_24 = 15.0
        real(JPRB), parameter ::  threesixty_over_365 = 0.98630137
        real(JPRB), parameter ::  min_declination = -23.433
        real(JPRB), parameter ::  day_offset = 10.0  !  original equation had this nine

!*	Compute day angle

        day_angle = threesixty_over_365*(julian_day-1.0)*DEGREES_TO_RADIANS

        rlatitude  = latitude  * DEGREES_TO_RADIANS

!*	Compute equation of Time

        Equation_of_Time =  &           
               ( 0.000075   &
               + 0.001868*Cos(Day_Angle)  &
               - 0.032077*Sin(Day_Angle)  &
               - 0.014615*Cos(2.*Day_Angle) &
               - 0.040890*Sin(2.*Day_Angle) )*229.18/60 ! in hours

!*	Compute local sun time
        local_sun_time  = universal_time   &
                       + Equation_of_Time &
                       + longitude/threesixty_over_24

!*	Compute solar declination

        solar_declination = &
         (       0.006918  &
               - 0.399912 * cos(day_angle)    &
               + 0.070257 * sin(day_angle)    &
               - 0.006758 * cos(2*day_angle)  &
               + 0.000907 * sin(2*day_angle)  &
               - 0.002697 * cos(3*day_angle)  &
               + 0.001480 * sin(3*day_angle) ) 

!*	Compute hour angle
        hour_angle = threesixty_over_24*mod(local_sun_time+12.0_JPRB , 24.0_JPRB)

!*	Compute solar zenith angle
        cosza = sin(rlatitude)*sin(solar_declination)  &
              + cos(rlatitude)*cos(solar_declination)*cos(hour_angle*DEGREES_TO_RADIANS)

        ZA = acos(cosza)*one_eighty_over_pi

!*	Compute solar azimuth angle
        solar_elevation = 90.0 - ZA

        If(Solar_Elevation .eq. 90.0) Then
          Az = 180.0    ! handle arbitrary case
        Else
          cosaz = (sin(solar_elevation*DEGREES_TO_RADIANS)*sin(rlatitude) - &
                  sin(solar_declination)) / &
                 (cos(solar_elevation*DEGREES_TO_RADIANS)*cos(rlatitude))

          If(cosaz .lt. -1.0) cosaz = -1.0
          If(cosaz .gt.  1.0) cosaz =  1.0

          Az = acos(cosaz)*one_eighty_over_pi

!	  The above formula produces azimuth positive east, zero south.
!	  We want positive east, zero north.

!	
          If (Az .ge. 0.0) Then
              Az = 180.0 - Az
          Else
              Az = -180.0 + Az
          EndIf

          If(hour_angle .lt. 180.0) Az = - Az
          If(Az .lt. 0) Az = 360.0 + Az

        EndIf

   end Subroutine Solar_ZA_Az

End module rttov_zutility
