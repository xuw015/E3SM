module crmdims

    use params_kind, only: crm_rknd

    implicit none

    integer, parameter ::  crm_nx=CRM_NX
    integer, parameter ::  crm_ny=CRM_NY
    integer, parameter ::  crm_nz=CRM_NZ

    integer, parameter ::  crm_nx_rad=CRM_NX_RAD
    integer, parameter ::  crm_ny_rad=CRM_NY_RAD

    real(crm_rknd), parameter :: crm_dx=CRM_DX
    real(crm_rknd), parameter :: crm_dy=crm_dx
    real(crm_rknd), parameter :: crm_dt=CRM_DT

#if defined(MMF_BCVT)
    ! use 1 for bulk convective variance transport 
    integer, parameter :: crm_nvark=1 
#else
    ! specify number of tracers for spectral convective variance transport
    integer, parameter :: crm_nvark=CRM_NTRC 
#endif

end module crmdims
