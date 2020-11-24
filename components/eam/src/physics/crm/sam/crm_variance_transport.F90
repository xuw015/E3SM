module variance_transport_mod
#if defined(MMF_VARIANCE_TRANSPORT)
!-------------------------------------------------------------------------------
! Purpose:
!
! Convective Scale Variance Transport (CSVT)
! This module supports the capability of transporting the CRM's internal 
! variance on the large-scale host model's grid using individual wavenumbers
! from a truncated Fouier transformed version of the main prognostics variables.
!
! Author: Walter Hannah - Lawrence Livermore National Lab
!-------------------------------------------------------------------------------
   use params_kind,  only: crm_rknd
   use crmdims,      only: crm_nvark
   use grid,         only: nzm
   use params,       only: pi
   use openacc_utils

   implicit none

   public allocate_CSVT
   public deallocate_CSVT
   public CSVT_diagnose
   public CSVT_forcing

   real(crm_rknd), allocatable :: t_csvt_amp_tend(:,:,:)  ! input tendency of LSE amplitude per wavenumber
   real(crm_rknd), allocatable :: q_csvt_amp_tend(:,:,:)  ! input tendency of QT  amplitude per wavenumber

   real(crm_rknd), allocatable :: t_csvt_amp(:,:,:)       ! LSE amplitude per wavenumber
   real(crm_rknd), allocatable :: q_csvt_amp(:,:,:)       ! QT  amplitude per wavenumber
   real(crm_rknd), allocatable :: t_csvt_phs(:,:,:)       ! LSE phase per wavenumber
   real(crm_rknd), allocatable :: q_csvt_phs(:,:,:)       ! QT  phase per wavenumber

contains

!===============================================================================
!===============================================================================
subroutine allocate_CSVT(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Allocate and initialize variables
   !----------------------------------------------------------------------------
   implicit none
   integer, intent(in) :: ncrms

   allocate( t_csvt_amp_tend( ncrms, nzm, crm_nvark ) )
   allocate( q_csvt_amp_tend( ncrms, nzm, crm_nvark ) )
   allocate( t_csvt_amp( ncrms, nzm, crm_nvark ) )
   allocate( q_csvt_amp( ncrms, nzm, crm_nvark ) )
   allocate( t_csvt_phs( ncrms, nzm, crm_nvark ) )
   allocate( q_csvt_phs( ncrms, nzm, crm_nvark ) )

   call prefetch( t_csvt_amp_tend )
   call prefetch( q_csvt_amp_tend )
   call prefetch( t_csvt_amp )
   call prefetch( q_csvt_amp )
   call prefetch( t_csvt_phs )
   call prefetch( q_csvt_phs )

   t_csvt_amp_tend(:,:,:) = 0.0_crm_rknd
   q_csvt_amp_tend(:,:,:) = 0.0_crm_rknd
   t_csvt_amp(:,:,:) = 0.0_crm_rknd
   q_csvt_amp(:,:,:) = 0.0_crm_rknd
   t_csvt_phs(:,:,:) = 0.0_crm_rknd
   q_csvt_phs(:,:,:) = 0.0_crm_rknd

end subroutine allocate_CSVT

!===============================================================================
!===============================================================================
subroutine deallocate_CSVT()
   !----------------------------------------------------------------------------
   ! Purpose: Deallocate variables
   !----------------------------------------------------------------------------
   deallocate( t_csvt_amp_tend )
   deallocate( q_csvt_amp_tend )
   deallocate( t_csvt_amp )
   deallocate( q_csvt_amp )
   deallocate( t_csvt_phs )
   deallocate( q_csvt_phs )
end subroutine deallocate_CSVT

!===============================================================================
!===============================================================================
subroutine CSVT_diagnose(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Diagnose amplitude for each wavenumber for variance transport
   !----------------------------------------------------------------------------
   use grid,      only: nx,ny
   use vars,      only: t,qv,qcl,qci
   use crmdims,   only: crm_dx
   use fftpack51D
   implicit none

   ! interface arguments
   integer, intent(in) :: ncrms

   ! local variables
   integer, parameter :: lensav = nx+15 ! must be at least N + INT(LOG(REAL(N))) + 4.
   real(crm_rknd), dimension(nx)    :: fft_out_t   ! for FFT input and output
   real(crm_rknd), dimension(nx)    :: fft_out_q   ! for FFT input and output
   real(crm_rknd), dimension(nx)    :: work        ! work array
   real(crm_rknd), dimension(lensav):: wsave       ! prime factors of N and certain trig values used in rfft1f
   ! real(crm_rknd), dimension(nx)    :: wave_num    ! only for debugging
   integer :: i, j, k, icrm   ! loop iterators
   integer :: ier             ! error return code
   !----------------------------------------------------------------------------

   ! initialization for FFT
   call rfft1i(nx,wsave,lensav,ier)
   if(ier /= 0) write(0,*) 'ERROR: rfftmi(): CSVT_diagnose - FFT initialization error ',ier

   ! diagnose amplitude for each wavenumber 
   !$acc parallel loop collapse(3) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do icrm = 1,ncrms

            ! initialize FFT input
            do i = 1,nx
               fft_out_t(i) = t(icrm,i,j,k)
               fft_out_q(i) = qv(icrm,i,j,k)
               ! fft_out_q(i) = qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
            end do

            ! do the forward transform
            call rfft1f( nx, 1, fft_out_t(:), nx, wsave, lensav, work(:), nx, ier )
            if (ier/=0) write(0,*) 'ERROR: rfftmf(): CSVT_diagnose - fft_out_t ',ier

            call rfft1f( nx, 1, fft_out_q(:), nx, wsave, lensav, work(:), nx, ier )
            if (ier/=0) write(0,*) 'ERROR: rfftmf(): CSVT_diagnose - fft_out_q ',ier

            ! calculate amplitude of each wavenumber
            do i = 1,crm_nvark
               ! wave_num(i) = nx/crm_dx * (real(i,crm_rknd)/real(nx,crm_rknd))
               ! phase
               t_csvt_phs(icrm,k,i) = atan2(fft_out_t(2*i+1),fft_out_t(2*i))
               q_csvt_phs(icrm,k,i) = atan2(fft_out_q(2*i+1),fft_out_q(2*i))
               ! amplitude
               t_csvt_amp(icrm,k,i) = sqrt( fft_out_t(2*i)**2. + fft_out_t(2*i+1)**2. )
               q_csvt_amp(icrm,k,i) = sqrt( fft_out_q(2*i)**2. + fft_out_q(2*i+1)**2. )
            end do


            ! ! print FFT weights
            ! write(*,*)
            ! write(*,200) 'i','wave num','FFT cos wgt','FFT sin wgt','amp','phase'
            ! do i = 1,crm_nvark
            !    write(*,201) i, wave_num(i), fft_out_t(2*i), fft_out_t(2*i+1), t_csvt_amp(icrm,k,i), t_csvt_phs(icrm,k,i)
            ! end do

         end do
      end do
   end do
   
! 200 format(A12,A16,A16,A16,A16,A16)
! 201 format(i12,f16.8,f16.8,f16.8,f16.8,f16.8,f16.8)

end subroutine CSVT_diagnose

!===============================================================================
!===============================================================================
! #ifdef SKIP_ME

subroutine CSVT_forcing(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Calculate forcing for variance injection and apply limiters
   !----------------------------------------------------------------------------
   use grid,         only: nx,ny,dtn
   use crmdims,      only: crm_dx
   use vars,         only: t
   use microphysics, only: micro_field, index_water_vapor
   implicit none

   ! interface arguments
   integer, intent(in) :: ncrms

   ! local variables
   real(crm_rknd), allocatable :: ttend_loc(:,:,:,:)
   real(crm_rknd), allocatable :: qtend_loc(:,:,:,:)
   real(crm_rknd) :: angle,freq
   integer :: i, j, k, icrm, iwn   ! loop iterators
   !----------------------------------------------------------------------------

   allocate( ttend_loc( ncrms, nx, ny, nzm ) )
   allocate( qtend_loc( ncrms, nx, ny, nzm ) )

   !$acc parallel loop collapse(4) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do i = 1,nx
            do icrm = 1,ncrms
               ttend_loc(icrm,i,j,k) = 0
               qtend_loc(icrm,i,j,k) = 0
               angle = (2*pi)*(i-1)/nx
               do iwn = 1,crm_nvark
                  freq = nx/crm_dx * (real(i,crm_rknd)/real(nx,crm_rknd))
                  ttend_loc(icrm,i,j,k) = ttend_loc(icrm,i,j,k) + t_csvt_amp_tend(icrm,k,iwn)*cos(freq*angle-t_csvt_phs(icrm,k,iwn))
                  qtend_loc(icrm,i,j,k) = qtend_loc(icrm,i,j,k) + q_csvt_amp_tend(icrm,k,iwn)*cos(freq*angle-q_csvt_phs(icrm,k,iwn))
               end do
            end do
         end do
      end do
   end do

   !$acc parallel loop collapse(4) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do i = 1,nx
            do icrm = 1,ncrms
               t(icrm,i,j,k) = t(icrm,i,j,k) + ttend_loc(icrm,i,j,k) * dtn
               micro_field(icrm,i,j,k,index_water_vapor) = micro_field(icrm,i,j,k,index_water_vapor) + qtend_loc(icrm,i,j,k) * dtn
            end do
         end do
      end do
   end do

   deallocate( ttend_loc )
   deallocate( qtend_loc )

end subroutine CSVT_forcing

! #endif

!===============================================================================
!===============================================================================
#endif
end module variance_transport_mod
