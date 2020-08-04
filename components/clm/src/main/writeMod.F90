module writeMod
contains
subroutine write_vars()
     use fileio_mod, only : fio_open, fio_close
     use clm_instMod, only : canopystate_vars
     use GridcellType, only : grc_pp
     use TopounitDataType, only : top_as
     use TopounitDataType, only : top_af
     use ColumnType, only : col_pp
     use ColumnDataType, only : col_es
     use ColumnDataType, only : col_ws
     use ColumnDataType, only : col_wf
     use ColumnDataType, only : col_cs
     use ColumnDataType, only : col_nf
     use ColumnDataType, only : col_ns
     use ColumnDataType, only : col_ps
     use ColumnDataType, only : col_cf
     use ColumnDataType, only : col_pf
     use ColumnDataType, only : c14_col_cf
     use ColumnDataType, only : c13_col_cf
     use ColumnDataType, only : c13_col_cs
     use ColumnDataType, only : c14_col_cs
     use VegetationType, only : veg_pp
     use VegetationDataType, only : veg_es
     use VegetationDataType, only : veg_cs
     use VegetationDataType, only : veg_nf
     use VegetationDataType, only : veg_cf
     use VegetationDataType, only : veg_ef
     use VegetationDataType, only : veg_pf
     use VegetationDataType, only : c13_veg_cf
     use VegetationDataType, only : c14_veg_cf
     use VegetationDataType, only : c14_veg_cs
     use VegetationDataType, only : veg_ns
     use VegetationDataType, only : veg_ps
     use VegetationPropertiesType, only : veg_vp
     use clm_instMod, only : photosyns_vars
     use clm_instMod, only : atm2lnd_vars
     use clm_instMod, only : ch4_vars
     use clm_instMod, only : soilstate_vars
     use clm_instMod, only : energyflux_vars
     use clm_instMod, only : cnstate_vars
     use clm_instMod, only : soilhydrology_vars
     use clm_instMod, only : crop_vars
     use CNDecompCascadeConType, only : decomp_cascade_con
     use clm_instMod, only : sedflux_vars

     implicit none
      integer :: fid
     character(len=256) :: ofile = "output_EcosystemDyn_vars.txt"
     fid = 23
     call fio_open(fid,ofile, 2)

     !====================== filter ======================!


     !====================== clumps ======================!


     !====================== procinfo ======================!


     !====================== canopystate_vars ======================!


     !====================== grc_pp ======================!

     write (fid, "(A)") "grc_pp%gindex"
     write (fid, *) grc_pp%gindex
     write (fid, "(A)") "grc_pp%area"
     write (fid, *) grc_pp%area
     write (fid, "(A)") "grc_pp%lat"
     write (fid, *) grc_pp%lat
     write (fid, "(A)") "grc_pp%lon"
     write (fid, *) grc_pp%lon
     write (fid, "(A)") "grc_pp%latdeg"
     write (fid, *) grc_pp%latdeg
     write (fid, "(A)") "grc_pp%londeg"
     write (fid, *) grc_pp%londeg
     write (fid, "(A)") "grc_pp%topi"
     write (fid, *) grc_pp%topi
     write (fid, "(A)") "grc_pp%topf"
     write (fid, *) grc_pp%topf
     write (fid, "(A)") "grc_pp%ntopounits"
     write (fid, *) grc_pp%ntopounits
     write (fid, "(A)") "grc_pp%lndi"
     write (fid, *) grc_pp%lndi
     write (fid, "(A)") "grc_pp%lndf"
     write (fid, *) grc_pp%lndf
     write (fid, "(A)") "grc_pp%nlandunits"
     write (fid, *) grc_pp%nlandunits
     write (fid, "(A)") "grc_pp%coli"
     write (fid, *) grc_pp%coli
     write (fid, "(A)") "grc_pp%colf"
     write (fid, *) grc_pp%colf
     write (fid, "(A)") "grc_pp%ncolumns"
     write (fid, *) grc_pp%ncolumns
     write (fid, "(A)") "grc_pp%pfti"
     write (fid, *) grc_pp%pfti
     write (fid, "(A)") "grc_pp%pftf"
     write (fid, *) grc_pp%pftf
     write (fid, "(A)") "grc_pp%npfts"
     write (fid, *) grc_pp%npfts
     write (fid, "(A)") "grc_pp%max_dayl"
     write (fid, *) grc_pp%max_dayl
     write (fid, "(A)") "grc_pp%dayl"
     write (fid, *) grc_pp%dayl
     write (fid, "(A)") "grc_pp%prev_dayl"
     write (fid, *) grc_pp%prev_dayl
     write (fid, "(A)") "grc_pp%landunit_indices"
     write (fid, *) grc_pp%landunit_indices

     !====================== grc_cs ======================!


     !====================== grc_ns ======================!


     !====================== grc_ps ======================!


     !====================== grc_nf ======================!


     !====================== grc_cf ======================!


     !====================== grc_pf ======================!


     !====================== grc_wf ======================!


     !====================== grc_ws ======================!


     !====================== grc_es ======================!


     !====================== grc_ef ======================!


     !====================== c14_grc_cf ======================!


     !====================== c13_grc_cf ======================!


     !====================== top_as ======================!

     write (fid, "(A)") "top_as%tbot"
     write (fid, *) top_as%tbot
     write (fid, "(A)") "top_as%pbot"
     write (fid, *) top_as%pbot
     write (fid, "(A)") "top_as%rhbot"
     write (fid, *) top_as%rhbot
     write (fid, "(A)") "top_as%windbot"
     write (fid, *) top_as%windbot

     !====================== top_af ======================!

     write (fid, "(A)") "top_af%rain"
     write (fid, *) top_af%rain
     write (fid, "(A)") "top_af%snow"
     write (fid, *) top_af%snow
     write (fid, "(A)") "top_af%prec10d"
     write (fid, *) top_af%prec10d
     write (fid, "(A)") "top_af%prec60d"
     write (fid, *) top_af%prec60d

     !====================== lun_pp ======================!


     !====================== col_pp ======================!

     write (fid, "(A)") "col_pp%gridcell"
     write (fid, *) col_pp%gridcell
     write (fid, "(A)") "col_pp%wtgcell"
     write (fid, *) col_pp%wtgcell
     write (fid, "(A)") "col_pp%topounit"
     write (fid, *) col_pp%topounit
     write (fid, "(A)") "col_pp%wttopounit"
     write (fid, *) col_pp%wttopounit
     write (fid, "(A)") "col_pp%landunit"
     write (fid, *) col_pp%landunit
     write (fid, "(A)") "col_pp%wtlunit"
     write (fid, *) col_pp%wtlunit
     write (fid, "(A)") "col_pp%pfti"
     write (fid, *) col_pp%pfti
     write (fid, "(A)") "col_pp%pftf"
     write (fid, *) col_pp%pftf
     write (fid, "(A)") "col_pp%npfts"
     write (fid, *) col_pp%npfts
     write (fid, "(A)") "col_pp%itype"
     write (fid, *) col_pp%itype
     write (fid, "(A)") "col_pp%active"
     write (fid, *) col_pp%active
     write (fid, "(A)") "col_pp%glc_topo"
     write (fid, *) col_pp%glc_topo
     write (fid, "(A)") "col_pp%micro_sigma"
     write (fid, *) col_pp%micro_sigma
     write (fid, "(A)") "col_pp%n_melt"
     write (fid, *) col_pp%n_melt
     write (fid, "(A)") "col_pp%topo_slope"
     write (fid, *) col_pp%topo_slope
     write (fid, "(A)") "col_pp%topo_std"
     write (fid, *) col_pp%topo_std
     write (fid, "(A)") "col_pp%hslp_p10"
     write (fid, *) col_pp%hslp_p10
     write (fid, "(A)") "col_pp%nlevbed"
     write (fid, *) col_pp%nlevbed
     write (fid, "(A)") "col_pp%zibed"
     write (fid, *) col_pp%zibed
     write (fid, "(A)") "col_pp%snl"
     write (fid, *) col_pp%snl
     write (fid, "(A)") "col_pp%dz"
     write (fid, *) col_pp%dz
     write (fid, "(A)") "col_pp%z"
     write (fid, *) col_pp%z
     write (fid, "(A)") "col_pp%zi"
     write (fid, *) col_pp%zi
     write (fid, "(A)") "col_pp%zii"
     write (fid, *) col_pp%zii
     write (fid, "(A)") "col_pp%dz_lake"
     write (fid, *) col_pp%dz_lake
     write (fid, "(A)") "col_pp%z_lake"
     write (fid, *) col_pp%z_lake
     write (fid, "(A)") "col_pp%lakedepth"
     write (fid, *) col_pp%lakedepth
     write (fid, "(A)") "col_pp%hydrologically_active"
     write (fid, *) col_pp%hydrologically_active

     !====================== col_es ======================!

     write (fid, "(A)") "col_es%t_soisno"
     write (fid, *) col_es%t_soisno
     write (fid, "(A)") "col_es%t_soi10cm"
     write (fid, *) col_es%t_soi10cm
     write (fid, "(A)") "col_es%t_soi17cm"
     write (fid, *) col_es%t_soi17cm

     !====================== col_ef ======================!


     !====================== col_ws ======================!

     write (fid, "(A)") "col_ws%h2osoi_liq"
     write (fid, *) col_ws%h2osoi_liq
     write (fid, "(A)") "col_ws%h2osoi_vol"
     write (fid, *) col_ws%h2osoi_vol
     write (fid, "(A)") "col_ws%snow_depth"
     write (fid, *) col_ws%snow_depth
     write (fid, "(A)") "col_ws%wf"
     write (fid, *) col_ws%wf
     write (fid, "(A)") "col_ws%wf2"
     write (fid, *) col_ws%wf2

     !====================== col_wf ======================!

     write (fid, "(A)") "col_wf%qflx_evap_veg"
     write (fid, *) col_wf%qflx_evap_veg
     write (fid, "(A)") "col_wf%qflx_tran_veg"
     write (fid, *) col_wf%qflx_tran_veg

     !====================== col_cs ======================!

     write (fid, "(A)") "col_cs%totvegc"
     write (fid, *) col_cs%totvegc
     write (fid, "(A)") "col_cs%leafc"
     write (fid, *) col_cs%leafc
     write (fid, "(A)") "col_cs%deadstemc"
     write (fid, *) col_cs%deadstemc
     write (fid, "(A)") "col_cs%decomp_cpools_vr"
     write (fid, *) col_cs%decomp_cpools_vr
     write (fid, "(A)") "col_cs%prod1c"
     write (fid, *) col_cs%prod1c
     write (fid, "(A)") "col_cs%prod10c"
     write (fid, *) col_cs%prod10c
     write (fid, "(A)") "col_cs%prod100c"
     write (fid, *) col_cs%prod100c
     write (fid, "(A)") "col_cs%totlitc"
     write (fid, *) col_cs%totlitc
     write (fid, "(A)") "col_cs%totsomc"
     write (fid, *) col_cs%totsomc

     !====================== col_nf ======================!

     write (fid, "(A)") "col_nf%hrv_deadstemn_to_prod10n"
     write (fid, *) col_nf%hrv_deadstemn_to_prod10n
     write (fid, "(A)") "col_nf%hrv_deadstemn_to_prod100n"
     write (fid, *) col_nf%hrv_deadstemn_to_prod100n
     write (fid, "(A)") "col_nf%m_n_to_litr_met_fire"
     write (fid, *) col_nf%m_n_to_litr_met_fire
     write (fid, "(A)") "col_nf%m_n_to_litr_cel_fire"
     write (fid, *) col_nf%m_n_to_litr_cel_fire
     write (fid, "(A)") "col_nf%m_n_to_litr_lig_fire"
     write (fid, *) col_nf%m_n_to_litr_lig_fire
     write (fid, "(A)") "col_nf%harvest_n_to_litr_met_n"
     write (fid, *) col_nf%harvest_n_to_litr_met_n
     write (fid, "(A)") "col_nf%harvest_n_to_litr_cel_n"
     write (fid, *) col_nf%harvest_n_to_litr_cel_n
     write (fid, "(A)") "col_nf%harvest_n_to_litr_lig_n"
     write (fid, *) col_nf%harvest_n_to_litr_lig_n
     write (fid, "(A)") "col_nf%harvest_n_to_cwdn"
     write (fid, *) col_nf%harvest_n_to_cwdn
     write (fid, "(A)") "col_nf%hrv_cropn_to_prod1n"
     write (fid, *) col_nf%hrv_cropn_to_prod1n
     write (fid, "(A)") "col_nf%fire_mortality_n_to_cwdn"
     write (fid, *) col_nf%fire_mortality_n_to_cwdn
     write (fid, "(A)") "col_nf%ndep_to_sminn"
     write (fid, *) col_nf%ndep_to_sminn
     write (fid, "(A)") "col_nf%nfix_to_sminn"
     write (fid, *) col_nf%nfix_to_sminn
     write (fid, "(A)") "col_nf%fert_to_sminn"
     write (fid, *) col_nf%fert_to_sminn
     write (fid, "(A)") "col_nf%soyfixn_to_sminn"
     write (fid, *) col_nf%soyfixn_to_sminn
     write (fid, "(A)") "col_nf%phenology_n_to_litr_met_n"
     write (fid, *) col_nf%phenology_n_to_litr_met_n
     write (fid, "(A)") "col_nf%phenology_n_to_litr_cel_n"
     write (fid, *) col_nf%phenology_n_to_litr_cel_n
     write (fid, "(A)") "col_nf%phenology_n_to_litr_lig_n"
     write (fid, *) col_nf%phenology_n_to_litr_lig_n
     write (fid, "(A)") "col_nf%gap_mortality_n_to_litr_met_n"
     write (fid, *) col_nf%gap_mortality_n_to_litr_met_n
     write (fid, "(A)") "col_nf%gap_mortality_n_to_litr_cel_n"
     write (fid, *) col_nf%gap_mortality_n_to_litr_cel_n
     write (fid, "(A)") "col_nf%gap_mortality_n_to_litr_lig_n"
     write (fid, *) col_nf%gap_mortality_n_to_litr_lig_n
     write (fid, "(A)") "col_nf%gap_mortality_n_to_cwdn"
     write (fid, *) col_nf%gap_mortality_n_to_cwdn
     write (fid, "(A)") "col_nf%decomp_cascade_ntransfer_vr"
     write (fid, *) col_nf%decomp_cascade_ntransfer_vr
     write (fid, "(A)") "col_nf%decomp_cascade_sminn_flux_vr"
     write (fid, *) col_nf%decomp_cascade_sminn_flux_vr
     write (fid, "(A)") "col_nf%potential_immob_vr"
     write (fid, *) col_nf%potential_immob_vr
     write (fid, "(A)") "col_nf%potential_immob"
     write (fid, *) col_nf%potential_immob
     write (fid, "(A)") "col_nf%actual_immob_vr"
     write (fid, *) col_nf%actual_immob_vr
     write (fid, "(A)") "col_nf%actual_immob"
     write (fid, *) col_nf%actual_immob
     write (fid, "(A)") "col_nf%sminn_to_plant_vr"
     write (fid, *) col_nf%sminn_to_plant_vr
     write (fid, "(A)") "col_nf%sminn_to_plant"
     write (fid, *) col_nf%sminn_to_plant
     write (fid, "(A)") "col_nf%supplement_to_sminn_vr"
     write (fid, *) col_nf%supplement_to_sminn_vr
     write (fid, "(A)") "col_nf%gross_nmin_vr"
     write (fid, *) col_nf%gross_nmin_vr
     write (fid, "(A)") "col_nf%net_nmin_vr"
     write (fid, *) col_nf%net_nmin_vr
     write (fid, "(A)") "col_nf%f_nit_vr"
     write (fid, *) col_nf%f_nit_vr
     write (fid, "(A)") "col_nf%f_denit_vr"
     write (fid, *) col_nf%f_denit_vr
     write (fid, "(A)") "col_nf%pot_f_nit_vr"
     write (fid, *) col_nf%pot_f_nit_vr
     write (fid, "(A)") "col_nf%pot_f_denit_vr"
     write (fid, *) col_nf%pot_f_denit_vr
     write (fid, "(A)") "col_nf%n2_n2o_ratio_denit_vr"
     write (fid, *) col_nf%n2_n2o_ratio_denit_vr
     write (fid, "(A)") "col_nf%actual_immob_no3_vr"
     write (fid, *) col_nf%actual_immob_no3_vr
     write (fid, "(A)") "col_nf%actual_immob_nh4_vr"
     write (fid, *) col_nf%actual_immob_nh4_vr
     write (fid, "(A)") "col_nf%smin_no3_to_plant_vr"
     write (fid, *) col_nf%smin_no3_to_plant_vr
     write (fid, "(A)") "col_nf%smin_nh4_to_plant_vr"
     write (fid, *) col_nf%smin_nh4_to_plant_vr
     write (fid, "(A)") "col_nf%sminn_to_denit_decomp_cascade_vr"
     write (fid, *) col_nf%sminn_to_denit_decomp_cascade_vr
     write (fid, "(A)") "col_nf%sminn_to_denit_excess_vr"
     write (fid, *) col_nf%sminn_to_denit_excess_vr
     write (fid, "(A)") "col_nf%dwt_prod10n_gain"
     write (fid, *) col_nf%dwt_prod10n_gain
     write (fid, "(A)") "col_nf%dwt_prod100n_gain"
     write (fid, *) col_nf%dwt_prod100n_gain
     write (fid, "(A)") "col_nf%decomp_npools_sourcesink"
     write (fid, *) col_nf%decomp_npools_sourcesink
     write (fid, "(A)") "col_nf%plant_ndemand"
     write (fid, *) col_nf%plant_ndemand
     write (fid, "(A)") "col_nf%pmnf_decomp_cascade"
     write (fid, *) col_nf%pmnf_decomp_cascade

     !====================== col_ns ======================!

     write (fid, "(A)") "col_ns%decomp_npools_vr"
     write (fid, *) col_ns%decomp_npools_vr
     write (fid, "(A)") "col_ns%sminn_vr"
     write (fid, *) col_ns%sminn_vr
     write (fid, "(A)") "col_ns%smin_no3_vr"
     write (fid, *) col_ns%smin_no3_vr
     write (fid, "(A)") "col_ns%smin_nh4_vr"
     write (fid, *) col_ns%smin_nh4_vr
     write (fid, "(A)") "col_ns%sminn"
     write (fid, *) col_ns%sminn
     write (fid, "(A)") "col_ns%prod1n"
     write (fid, *) col_ns%prod1n
     write (fid, "(A)") "col_ns%prod10n"
     write (fid, *) col_ns%prod10n
     write (fid, "(A)") "col_ns%prod100n"
     write (fid, *) col_ns%prod100n

     !====================== col_ps ======================!

     write (fid, "(A)") "col_ps%decomp_ppools_vr"
     write (fid, *) col_ps%decomp_ppools_vr
     write (fid, "(A)") "col_ps%solutionp_vr"
     write (fid, *) col_ps%solutionp_vr
     write (fid, "(A)") "col_ps%labilep_vr"
     write (fid, *) col_ps%labilep_vr
     write (fid, "(A)") "col_ps%secondp_vr"
     write (fid, *) col_ps%secondp_vr
     write (fid, "(A)") "col_ps%occlp_vr"
     write (fid, *) col_ps%occlp_vr
     write (fid, "(A)") "col_ps%primp_vr"
     write (fid, *) col_ps%primp_vr
     write (fid, "(A)") "col_ps%prod1p"
     write (fid, *) col_ps%prod1p
     write (fid, "(A)") "col_ps%prod10p"
     write (fid, *) col_ps%prod10p
     write (fid, "(A)") "col_ps%prod100p"
     write (fid, *) col_ps%prod100p

     !====================== col_cf ======================!

     write (fid, "(A)") "col_cf%decomp_cpools_sourcesink"
     write (fid, *) col_cf%decomp_cpools_sourcesink
     write (fid, "(A)") "col_cf%decomp_cascade_hr_vr"
     write (fid, *) col_cf%decomp_cascade_hr_vr
     write (fid, "(A)") "col_cf%decomp_cascade_ctransfer_vr"
     write (fid, *) col_cf%decomp_cascade_ctransfer_vr
     write (fid, "(A)") "col_cf%decomp_k"
     write (fid, *) col_cf%decomp_k
     write (fid, "(A)") "col_cf%hr_vr"
     write (fid, *) col_cf%hr_vr
     write (fid, "(A)") "col_cf%o_scalar"
     write (fid, *) col_cf%o_scalar
     write (fid, "(A)") "col_cf%w_scalar"
     write (fid, *) col_cf%w_scalar
     write (fid, "(A)") "col_cf%t_scalar"
     write (fid, *) col_cf%t_scalar
     write (fid, "(A)") "col_cf%phr_vr"
     write (fid, *) col_cf%phr_vr
     write (fid, "(A)") "col_cf%phenology_c_to_litr_met_c"
     write (fid, *) col_cf%phenology_c_to_litr_met_c
     write (fid, "(A)") "col_cf%phenology_c_to_litr_cel_c"
     write (fid, *) col_cf%phenology_c_to_litr_cel_c
     write (fid, "(A)") "col_cf%phenology_c_to_litr_lig_c"
     write (fid, *) col_cf%phenology_c_to_litr_lig_c
     write (fid, "(A)") "col_cf%gap_mortality_c_to_litr_met_c"
     write (fid, *) col_cf%gap_mortality_c_to_litr_met_c
     write (fid, "(A)") "col_cf%gap_mortality_c_to_litr_cel_c"
     write (fid, *) col_cf%gap_mortality_c_to_litr_cel_c
     write (fid, "(A)") "col_cf%gap_mortality_c_to_litr_lig_c"
     write (fid, *) col_cf%gap_mortality_c_to_litr_lig_c
     write (fid, "(A)") "col_cf%gap_mortality_c_to_cwdc"
     write (fid, *) col_cf%gap_mortality_c_to_cwdc
     write (fid, "(A)") "col_cf%m_decomp_cpools_to_fire_vr"
     write (fid, *) col_cf%m_decomp_cpools_to_fire_vr
     write (fid, "(A)") "col_cf%m_c_to_litr_met_fire"
     write (fid, *) col_cf%m_c_to_litr_met_fire
     write (fid, "(A)") "col_cf%m_c_to_litr_cel_fire"
     write (fid, *) col_cf%m_c_to_litr_cel_fire
     write (fid, "(A)") "col_cf%m_c_to_litr_lig_fire"
     write (fid, *) col_cf%m_c_to_litr_lig_fire
     write (fid, "(A)") "col_cf%fire_mortality_c_to_cwdc"
     write (fid, *) col_cf%fire_mortality_c_to_cwdc
     write (fid, "(A)") "col_cf%decomp_cpools_yield_vr"
     write (fid, *) col_cf%decomp_cpools_yield_vr
     write (fid, "(A)") "col_cf%harvest_c_to_litr_met_c"
     write (fid, *) col_cf%harvest_c_to_litr_met_c
     write (fid, "(A)") "col_cf%harvest_c_to_litr_cel_c"
     write (fid, *) col_cf%harvest_c_to_litr_cel_c
     write (fid, "(A)") "col_cf%harvest_c_to_litr_lig_c"
     write (fid, *) col_cf%harvest_c_to_litr_lig_c
     write (fid, "(A)") "col_cf%harvest_c_to_cwdc"
     write (fid, *) col_cf%harvest_c_to_cwdc
     write (fid, "(A)") "col_cf%hrv_deadstemc_to_prod10c"
     write (fid, *) col_cf%hrv_deadstemc_to_prod10c
     write (fid, "(A)") "col_cf%hrv_deadstemc_to_prod100c"
     write (fid, *) col_cf%hrv_deadstemc_to_prod100c
     write (fid, "(A)") "col_cf%hrv_cropc_to_prod1c"
     write (fid, *) col_cf%hrv_cropc_to_prod1c
     write (fid, "(A)") "col_cf%dwt_prod10c_gain"
     write (fid, *) col_cf%dwt_prod10c_gain
     write (fid, "(A)") "col_cf%dwt_prod100c_gain"
     write (fid, *) col_cf%dwt_prod100c_gain
     write (fid, "(A)") "col_cf%rr"
     write (fid, *) col_cf%rr
     write (fid, "(A)") "col_cf%annsum_npp"
     write (fid, *) col_cf%annsum_npp
     write (fid, "(A)") "col_cf%lag_npp"
     write (fid, *) col_cf%lag_npp

     !====================== col_pf ======================!

     write (fid, "(A)") "col_pf%hrv_deadstemp_to_prod10p"
     write (fid, *) col_pf%hrv_deadstemp_to_prod10p
     write (fid, "(A)") "col_pf%hrv_deadstemp_to_prod100p"
     write (fid, *) col_pf%hrv_deadstemp_to_prod100p
     write (fid, "(A)") "col_pf%m_p_to_litr_met_fire"
     write (fid, *) col_pf%m_p_to_litr_met_fire
     write (fid, "(A)") "col_pf%m_p_to_litr_cel_fire"
     write (fid, *) col_pf%m_p_to_litr_cel_fire
     write (fid, "(A)") "col_pf%m_p_to_litr_lig_fire"
     write (fid, *) col_pf%m_p_to_litr_lig_fire
     write (fid, "(A)") "col_pf%harvest_p_to_litr_met_p"
     write (fid, *) col_pf%harvest_p_to_litr_met_p
     write (fid, "(A)") "col_pf%harvest_p_to_litr_cel_p"
     write (fid, *) col_pf%harvest_p_to_litr_cel_p
     write (fid, "(A)") "col_pf%harvest_p_to_litr_lig_p"
     write (fid, *) col_pf%harvest_p_to_litr_lig_p
     write (fid, "(A)") "col_pf%harvest_p_to_cwdp"
     write (fid, *) col_pf%harvest_p_to_cwdp
     write (fid, "(A)") "col_pf%hrv_cropp_to_prod1p"
     write (fid, *) col_pf%hrv_cropp_to_prod1p
     write (fid, "(A)") "col_pf%fire_mortality_p_to_cwdp"
     write (fid, *) col_pf%fire_mortality_p_to_cwdp
     write (fid, "(A)") "col_pf%pdep_to_sminp"
     write (fid, *) col_pf%pdep_to_sminp
     write (fid, "(A)") "col_pf%phenology_p_to_litr_met_p"
     write (fid, *) col_pf%phenology_p_to_litr_met_p
     write (fid, "(A)") "col_pf%phenology_p_to_litr_cel_p"
     write (fid, *) col_pf%phenology_p_to_litr_cel_p
     write (fid, "(A)") "col_pf%phenology_p_to_litr_lig_p"
     write (fid, *) col_pf%phenology_p_to_litr_lig_p
     write (fid, "(A)") "col_pf%gap_mortality_p_to_litr_met_p"
     write (fid, *) col_pf%gap_mortality_p_to_litr_met_p
     write (fid, "(A)") "col_pf%gap_mortality_p_to_litr_cel_p"
     write (fid, *) col_pf%gap_mortality_p_to_litr_cel_p
     write (fid, "(A)") "col_pf%gap_mortality_p_to_litr_lig_p"
     write (fid, *) col_pf%gap_mortality_p_to_litr_lig_p
     write (fid, "(A)") "col_pf%gap_mortality_p_to_cwdp"
     write (fid, *) col_pf%gap_mortality_p_to_cwdp
     write (fid, "(A)") "col_pf%decomp_cascade_ptransfer_vr"
     write (fid, *) col_pf%decomp_cascade_ptransfer_vr
     write (fid, "(A)") "col_pf%decomp_cascade_sminp_flux_vr"
     write (fid, *) col_pf%decomp_cascade_sminp_flux_vr
     write (fid, "(A)") "col_pf%potential_immob_p_vr"
     write (fid, *) col_pf%potential_immob_p_vr
     write (fid, "(A)") "col_pf%potential_immob_p"
     write (fid, *) col_pf%potential_immob_p
     write (fid, "(A)") "col_pf%actual_immob_p"
     write (fid, *) col_pf%actual_immob_p
     write (fid, "(A)") "col_pf%sminp_to_plant_vr"
     write (fid, *) col_pf%sminp_to_plant_vr
     write (fid, "(A)") "col_pf%sminp_to_plant"
     write (fid, *) col_pf%sminp_to_plant
     write (fid, "(A)") "col_pf%gross_pmin_vr"
     write (fid, *) col_pf%gross_pmin_vr
     write (fid, "(A)") "col_pf%net_pmin_vr"
     write (fid, *) col_pf%net_pmin_vr
     write (fid, "(A)") "col_pf%biochem_pmin_vr"
     write (fid, *) col_pf%biochem_pmin_vr
     write (fid, "(A)") "col_pf%primp_to_labilep_vr"
     write (fid, *) col_pf%primp_to_labilep_vr
     write (fid, "(A)") "col_pf%dwt_prod10p_gain"
     write (fid, *) col_pf%dwt_prod10p_gain
     write (fid, "(A)") "col_pf%dwt_prod100p_gain"
     write (fid, *) col_pf%dwt_prod100p_gain
     write (fid, "(A)") "col_pf%decomp_ppools_sourcesink"
     write (fid, *) col_pf%decomp_ppools_sourcesink
     write (fid, "(A)") "col_pf%plant_pdemand"
     write (fid, *) col_pf%plant_pdemand
     write (fid, "(A)") "col_pf%pmpf_decomp_cascade"
     write (fid, *) col_pf%pmpf_decomp_cascade
     write (fid, "(A)") "col_pf%col_plant_pdemand_vr"
     write (fid, *) col_pf%col_plant_pdemand_vr

     !====================== c14_col_cf ======================!

     write (fid, "(A)") "c14_col_cf%decomp_cpools_sourcesink"
     write (fid, *) c14_col_cf%decomp_cpools_sourcesink
     write (fid, "(A)") "c14_col_cf%hrv_cropc_to_prod1c"
     write (fid, *) c14_col_cf%hrv_cropc_to_prod1c
     write (fid, "(A)") "c14_col_cf%dwt_prod10c_gain"
     write (fid, *) c14_col_cf%dwt_prod10c_gain
     write (fid, "(A)") "c14_col_cf%dwt_prod100c_gain"
     write (fid, *) c14_col_cf%dwt_prod100c_gain

     !====================== c13_col_cf ======================!

     write (fid, "(A)") "c13_col_cf%decomp_cpools_sourcesink"
     write (fid, *) c13_col_cf%decomp_cpools_sourcesink
     write (fid, "(A)") "c13_col_cf%hrv_cropc_to_prod1c"
     write (fid, *) c13_col_cf%hrv_cropc_to_prod1c
     write (fid, "(A)") "c13_col_cf%dwt_prod10c_gain"
     write (fid, *) c13_col_cf%dwt_prod10c_gain
     write (fid, "(A)") "c13_col_cf%dwt_prod100c_gain"
     write (fid, *) c13_col_cf%dwt_prod100c_gain

     !====================== c13_col_cs ======================!

     write (fid, "(A)") "c13_col_cs%decomp_cpools_vr"
     write (fid, *) c13_col_cs%decomp_cpools_vr
     write (fid, "(A)") "c13_col_cs%prod1c"
     write (fid, *) c13_col_cs%prod1c
     write (fid, "(A)") "c13_col_cs%prod10c"
     write (fid, *) c13_col_cs%prod10c
     write (fid, "(A)") "c13_col_cs%prod100c"
     write (fid, *) c13_col_cs%prod100c

     !====================== c14_col_cs ======================!

     write (fid, "(A)") "c14_col_cs%decomp_cpools_vr"
     write (fid, *) c14_col_cs%decomp_cpools_vr
     write (fid, "(A)") "c14_col_cs%seedc"
     write (fid, *) c14_col_cs%seedc
     write (fid, "(A)") "c14_col_cs%prod1c"
     write (fid, *) c14_col_cs%prod1c
     write (fid, "(A)") "c14_col_cs%prod10c"
     write (fid, *) c14_col_cs%prod10c
     write (fid, "(A)") "c14_col_cs%prod100c"
     write (fid, *) c14_col_cs%prod100c

     !====================== veg_pp ======================!

     write (fid, "(A)") "veg_pp%gridcell"
     write (fid, *) veg_pp%gridcell
     write (fid, "(A)") "veg_pp%wtgcell"
     write (fid, *) veg_pp%wtgcell
     write (fid, "(A)") "veg_pp%topounit"
     write (fid, *) veg_pp%topounit
     write (fid, "(A)") "veg_pp%wttopounit"
     write (fid, *) veg_pp%wttopounit
     write (fid, "(A)") "veg_pp%landunit"
     write (fid, *) veg_pp%landunit
     write (fid, "(A)") "veg_pp%wtlunit"
     write (fid, *) veg_pp%wtlunit
     write (fid, "(A)") "veg_pp%column"
     write (fid, *) veg_pp%column
     write (fid, "(A)") "veg_pp%wtcol"
     write (fid, *) veg_pp%wtcol
     write (fid, "(A)") "veg_pp%itype"
     write (fid, *) veg_pp%itype
     write (fid, "(A)") "veg_pp%mxy"
     write (fid, *) veg_pp%mxy
     write (fid, "(A)") "veg_pp%active"
     write (fid, *) veg_pp%active
     write (fid, "(A)") "veg_pp%is_veg"
     write (fid, *) veg_pp%is_veg
     write (fid, "(A)") "veg_pp%is_bareground"
     write (fid, *) veg_pp%is_bareground
     write (fid, "(A)") "veg_pp%wt_ed"
     write (fid, *) veg_pp%wt_ed
     write (fid, "(A)") "veg_pp%is_fates"
     write (fid, *) veg_pp%is_fates

     !====================== veg_es ======================!

     write (fid, "(A)") "veg_es%t_ref2m"
     write (fid, *) veg_es%t_ref2m
     write (fid, "(A)") "veg_es%t_a10"
     write (fid, *) veg_es%t_a10
     write (fid, "(A)") "veg_es%t_a10min"
     write (fid, *) veg_es%t_a10min
     write (fid, "(A)") "veg_es%t_a5min"
     write (fid, *) veg_es%t_a5min
     write (fid, "(A)") "veg_es%t_ref2m_min"
     write (fid, *) veg_es%t_ref2m_min
     write (fid, "(A)") "veg_es%t_ref2m_max"
     write (fid, *) veg_es%t_ref2m_max
     write (fid, "(A)") "veg_es%gdd0"
     write (fid, *) veg_es%gdd0
     write (fid, "(A)") "veg_es%gdd8"
     write (fid, *) veg_es%gdd8
     write (fid, "(A)") "veg_es%gdd10"
     write (fid, *) veg_es%gdd10
     write (fid, "(A)") "veg_es%gdd020"
     write (fid, *) veg_es%gdd020
     write (fid, "(A)") "veg_es%gdd820"
     write (fid, *) veg_es%gdd820
     write (fid, "(A)") "veg_es%gdd1020"
     write (fid, *) veg_es%gdd1020

     !====================== veg_ws ======================!


     !====================== veg_wf ======================!


     !====================== veg_cs ======================!

     write (fid, "(A)") "veg_cs%leafc"
     write (fid, *) veg_cs%leafc
     write (fid, "(A)") "veg_cs%leafc_storage"
     write (fid, *) veg_cs%leafc_storage
     write (fid, "(A)") "veg_cs%leafc_xfer"
     write (fid, *) veg_cs%leafc_xfer
     write (fid, "(A)") "veg_cs%frootc"
     write (fid, *) veg_cs%frootc
     write (fid, "(A)") "veg_cs%frootc_storage"
     write (fid, *) veg_cs%frootc_storage
     write (fid, "(A)") "veg_cs%frootc_xfer"
     write (fid, *) veg_cs%frootc_xfer
     write (fid, "(A)") "veg_cs%livestemc"
     write (fid, *) veg_cs%livestemc
     write (fid, "(A)") "veg_cs%livestemc_storage"
     write (fid, *) veg_cs%livestemc_storage
     write (fid, "(A)") "veg_cs%livestemc_xfer"
     write (fid, *) veg_cs%livestemc_xfer
     write (fid, "(A)") "veg_cs%deadstemc"
     write (fid, *) veg_cs%deadstemc
     write (fid, "(A)") "veg_cs%deadstemc_storage"
     write (fid, *) veg_cs%deadstemc_storage
     write (fid, "(A)") "veg_cs%deadstemc_xfer"
     write (fid, *) veg_cs%deadstemc_xfer
     write (fid, "(A)") "veg_cs%livecrootc"
     write (fid, *) veg_cs%livecrootc
     write (fid, "(A)") "veg_cs%livecrootc_storage"
     write (fid, *) veg_cs%livecrootc_storage
     write (fid, "(A)") "veg_cs%livecrootc_xfer"
     write (fid, *) veg_cs%livecrootc_xfer
     write (fid, "(A)") "veg_cs%deadcrootc"
     write (fid, *) veg_cs%deadcrootc
     write (fid, "(A)") "veg_cs%deadcrootc_storage"
     write (fid, *) veg_cs%deadcrootc_storage
     write (fid, "(A)") "veg_cs%deadcrootc_xfer"
     write (fid, *) veg_cs%deadcrootc_xfer
     write (fid, "(A)") "veg_cs%gresp_storage"
     write (fid, *) veg_cs%gresp_storage
     write (fid, "(A)") "veg_cs%gresp_xfer"
     write (fid, *) veg_cs%gresp_xfer
     write (fid, "(A)") "veg_cs%cpool"
     write (fid, *) veg_cs%cpool
     write (fid, "(A)") "veg_cs%xsmrpool"
     write (fid, *) veg_cs%xsmrpool
     write (fid, "(A)") "veg_cs%grainc"
     write (fid, *) veg_cs%grainc
     write (fid, "(A)") "veg_cs%grainc_storage"
     write (fid, *) veg_cs%grainc_storage
     write (fid, "(A)") "veg_cs%grainc_xfer"
     write (fid, *) veg_cs%grainc_xfer
     write (fid, "(A)") "veg_cs%cropseedc_deficit"
     write (fid, *) veg_cs%cropseedc_deficit
     write (fid, "(A)") "veg_cs%totvegc"
     write (fid, *) veg_cs%totvegc

     !====================== veg_nf ======================!

     write (fid, "(A)") "veg_nf%m_leafn_to_litter"
     write (fid, *) veg_nf%m_leafn_to_litter
     write (fid, "(A)") "veg_nf%m_frootn_to_litter"
     write (fid, *) veg_nf%m_frootn_to_litter
     write (fid, "(A)") "veg_nf%m_leafn_storage_to_litter"
     write (fid, *) veg_nf%m_leafn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_frootn_storage_to_litter"
     write (fid, *) veg_nf%m_frootn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_livestemn_storage_to_litter"
     write (fid, *) veg_nf%m_livestemn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_deadstemn_storage_to_litter"
     write (fid, *) veg_nf%m_deadstemn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_livecrootn_storage_to_litter"
     write (fid, *) veg_nf%m_livecrootn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_deadcrootn_storage_to_litter"
     write (fid, *) veg_nf%m_deadcrootn_storage_to_litter
     write (fid, "(A)") "veg_nf%m_leafn_xfer_to_litter"
     write (fid, *) veg_nf%m_leafn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_frootn_xfer_to_litter"
     write (fid, *) veg_nf%m_frootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_livestemn_xfer_to_litter"
     write (fid, *) veg_nf%m_livestemn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_deadstemn_xfer_to_litter"
     write (fid, *) veg_nf%m_deadstemn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_livecrootn_xfer_to_litter"
     write (fid, *) veg_nf%m_livecrootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_deadcrootn_xfer_to_litter"
     write (fid, *) veg_nf%m_deadcrootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%m_livestemn_to_litter"
     write (fid, *) veg_nf%m_livestemn_to_litter
     write (fid, "(A)") "veg_nf%m_deadstemn_to_litter"
     write (fid, *) veg_nf%m_deadstemn_to_litter
     write (fid, "(A)") "veg_nf%m_livecrootn_to_litter"
     write (fid, *) veg_nf%m_livecrootn_to_litter
     write (fid, "(A)") "veg_nf%m_deadcrootn_to_litter"
     write (fid, *) veg_nf%m_deadcrootn_to_litter
     write (fid, "(A)") "veg_nf%m_retransn_to_litter"
     write (fid, *) veg_nf%m_retransn_to_litter
     write (fid, "(A)") "veg_nf%m_npool_to_litter"
     write (fid, *) veg_nf%m_npool_to_litter
     write (fid, "(A)") "veg_nf%hrv_leafn_to_litter"
     write (fid, *) veg_nf%hrv_leafn_to_litter
     write (fid, "(A)") "veg_nf%hrv_frootn_to_litter"
     write (fid, *) veg_nf%hrv_frootn_to_litter
     write (fid, "(A)") "veg_nf%hrv_leafn_storage_to_litter"
     write (fid, *) veg_nf%hrv_leafn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_frootn_storage_to_litter"
     write (fid, *) veg_nf%hrv_frootn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_livestemn_storage_to_litter"
     write (fid, *) veg_nf%hrv_livestemn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadstemn_storage_to_litter"
     write (fid, *) veg_nf%hrv_deadstemn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_livecrootn_storage_to_litter"
     write (fid, *) veg_nf%hrv_livecrootn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadcrootn_storage_to_litter"
     write (fid, *) veg_nf%hrv_deadcrootn_storage_to_litter
     write (fid, "(A)") "veg_nf%hrv_leafn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_leafn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_frootn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_frootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_livestemn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_livestemn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadstemn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_deadstemn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_livecrootn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_livecrootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadcrootn_xfer_to_litter"
     write (fid, *) veg_nf%hrv_deadcrootn_xfer_to_litter
     write (fid, "(A)") "veg_nf%hrv_livestemn_to_litter"
     write (fid, *) veg_nf%hrv_livestemn_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadstemn_to_prod10n"
     write (fid, *) veg_nf%hrv_deadstemn_to_prod10n
     write (fid, "(A)") "veg_nf%hrv_deadstemn_to_prod100n"
     write (fid, *) veg_nf%hrv_deadstemn_to_prod100n
     write (fid, "(A)") "veg_nf%hrv_livecrootn_to_litter"
     write (fid, *) veg_nf%hrv_livecrootn_to_litter
     write (fid, "(A)") "veg_nf%hrv_deadcrootn_to_litter"
     write (fid, *) veg_nf%hrv_deadcrootn_to_litter
     write (fid, "(A)") "veg_nf%hrv_retransn_to_litter"
     write (fid, *) veg_nf%hrv_retransn_to_litter
     write (fid, "(A)") "veg_nf%hrv_npool_to_litter"
     write (fid, *) veg_nf%hrv_npool_to_litter
     write (fid, "(A)") "veg_nf%hrv_leafn_to_prod1n"
     write (fid, *) veg_nf%hrv_leafn_to_prod1n
     write (fid, "(A)") "veg_nf%hrv_livestemn_to_prod1n"
     write (fid, *) veg_nf%hrv_livestemn_to_prod1n
     write (fid, "(A)") "veg_nf%hrv_grainn_to_prod1n"
     write (fid, *) veg_nf%hrv_grainn_to_prod1n
     write (fid, "(A)") "veg_nf%grainn_xfer_to_grainn"
     write (fid, *) veg_nf%grainn_xfer_to_grainn
     write (fid, "(A)") "veg_nf%leafn_xfer_to_leafn"
     write (fid, *) veg_nf%leafn_xfer_to_leafn
     write (fid, "(A)") "veg_nf%frootn_xfer_to_frootn"
     write (fid, *) veg_nf%frootn_xfer_to_frootn
     write (fid, "(A)") "veg_nf%livestemn_xfer_to_livestemn"
     write (fid, *) veg_nf%livestemn_xfer_to_livestemn
     write (fid, "(A)") "veg_nf%deadstemn_xfer_to_deadstemn"
     write (fid, *) veg_nf%deadstemn_xfer_to_deadstemn
     write (fid, "(A)") "veg_nf%livecrootn_xfer_to_livecrootn"
     write (fid, *) veg_nf%livecrootn_xfer_to_livecrootn
     write (fid, "(A)") "veg_nf%deadcrootn_xfer_to_deadcrootn"
     write (fid, *) veg_nf%deadcrootn_xfer_to_deadcrootn
     write (fid, "(A)") "veg_nf%livestemn_to_litter"
     write (fid, *) veg_nf%livestemn_to_litter
     write (fid, "(A)") "veg_nf%grainn_to_food"
     write (fid, *) veg_nf%grainn_to_food
     write (fid, "(A)") "veg_nf%leafn_to_litter"
     write (fid, *) veg_nf%leafn_to_litter
     write (fid, "(A)") "veg_nf%leafn_to_retransn"
     write (fid, *) veg_nf%leafn_to_retransn
     write (fid, "(A)") "veg_nf%frootn_to_retransn"
     write (fid, *) veg_nf%frootn_to_retransn
     write (fid, "(A)") "veg_nf%frootn_to_litter"
     write (fid, *) veg_nf%frootn_to_litter
     write (fid, "(A)") "veg_nf%retransn_to_npool"
     write (fid, *) veg_nf%retransn_to_npool
     write (fid, "(A)") "veg_nf%sminn_to_npool"
     write (fid, *) veg_nf%sminn_to_npool
     write (fid, "(A)") "veg_nf%npool_to_grainn"
     write (fid, *) veg_nf%npool_to_grainn
     write (fid, "(A)") "veg_nf%npool_to_grainn_storage"
     write (fid, *) veg_nf%npool_to_grainn_storage
     write (fid, "(A)") "veg_nf%npool_to_leafn"
     write (fid, *) veg_nf%npool_to_leafn
     write (fid, "(A)") "veg_nf%npool_to_leafn_storage"
     write (fid, *) veg_nf%npool_to_leafn_storage
     write (fid, "(A)") "veg_nf%npool_to_frootn"
     write (fid, *) veg_nf%npool_to_frootn
     write (fid, "(A)") "veg_nf%npool_to_frootn_storage"
     write (fid, *) veg_nf%npool_to_frootn_storage
     write (fid, "(A)") "veg_nf%npool_to_livestemn"
     write (fid, *) veg_nf%npool_to_livestemn
     write (fid, "(A)") "veg_nf%npool_to_livestemn_storage"
     write (fid, *) veg_nf%npool_to_livestemn_storage
     write (fid, "(A)") "veg_nf%npool_to_deadstemn"
     write (fid, *) veg_nf%npool_to_deadstemn
     write (fid, "(A)") "veg_nf%npool_to_deadstemn_storage"
     write (fid, *) veg_nf%npool_to_deadstemn_storage
     write (fid, "(A)") "veg_nf%npool_to_livecrootn"
     write (fid, *) veg_nf%npool_to_livecrootn
     write (fid, "(A)") "veg_nf%npool_to_livecrootn_storage"
     write (fid, *) veg_nf%npool_to_livecrootn_storage
     write (fid, "(A)") "veg_nf%npool_to_deadcrootn"
     write (fid, *) veg_nf%npool_to_deadcrootn
     write (fid, "(A)") "veg_nf%npool_to_deadcrootn_storage"
     write (fid, *) veg_nf%npool_to_deadcrootn_storage
     write (fid, "(A)") "veg_nf%grainn_storage_to_xfer"
     write (fid, *) veg_nf%grainn_storage_to_xfer
     write (fid, "(A)") "veg_nf%leafn_storage_to_xfer"
     write (fid, *) veg_nf%leafn_storage_to_xfer
     write (fid, "(A)") "veg_nf%frootn_storage_to_xfer"
     write (fid, *) veg_nf%frootn_storage_to_xfer
     write (fid, "(A)") "veg_nf%livestemn_storage_to_xfer"
     write (fid, *) veg_nf%livestemn_storage_to_xfer
     write (fid, "(A)") "veg_nf%deadstemn_storage_to_xfer"
     write (fid, *) veg_nf%deadstemn_storage_to_xfer
     write (fid, "(A)") "veg_nf%livecrootn_storage_to_xfer"
     write (fid, *) veg_nf%livecrootn_storage_to_xfer
     write (fid, "(A)") "veg_nf%deadcrootn_storage_to_xfer"
     write (fid, *) veg_nf%deadcrootn_storage_to_xfer
     write (fid, "(A)") "veg_nf%fert"
     write (fid, *) veg_nf%fert
     write (fid, "(A)") "veg_nf%livestemn_to_deadstemn"
     write (fid, *) veg_nf%livestemn_to_deadstemn
     write (fid, "(A)") "veg_nf%livestemn_to_retransn"
     write (fid, *) veg_nf%livestemn_to_retransn
     write (fid, "(A)") "veg_nf%livecrootn_to_deadcrootn"
     write (fid, *) veg_nf%livecrootn_to_deadcrootn
     write (fid, "(A)") "veg_nf%livecrootn_to_retransn"
     write (fid, *) veg_nf%livecrootn_to_retransn
     write (fid, "(A)") "veg_nf%nfix_to_plantn"
     write (fid, *) veg_nf%nfix_to_plantn
     write (fid, "(A)") "veg_nf%crop_seedn_to_leaf"
     write (fid, *) veg_nf%crop_seedn_to_leaf
     write (fid, "(A)") "veg_nf%plant_ndemand"
     write (fid, *) veg_nf%plant_ndemand
     write (fid, "(A)") "veg_nf%smin_no3_to_plant"
     write (fid, *) veg_nf%smin_no3_to_plant
     write (fid, "(A)") "veg_nf%smin_nh4_to_plant"
     write (fid, *) veg_nf%smin_nh4_to_plant
     write (fid, "(A)") "veg_nf%sminn_to_plant"
     write (fid, *) veg_nf%sminn_to_plant
     write (fid, "(A)") "veg_nf%supplement_to_plantn"
     write (fid, *) veg_nf%supplement_to_plantn

     !====================== veg_cf ======================!

     write (fid, "(A)") "veg_cf%m_leafc_to_litter"
     write (fid, *) veg_cf%m_leafc_to_litter
     write (fid, "(A)") "veg_cf%m_leafc_storage_to_litter"
     write (fid, *) veg_cf%m_leafc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_leafc_xfer_to_litter"
     write (fid, *) veg_cf%m_leafc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_frootc_to_litter"
     write (fid, *) veg_cf%m_frootc_to_litter
     write (fid, "(A)") "veg_cf%m_frootc_storage_to_litter"
     write (fid, *) veg_cf%m_frootc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_frootc_xfer_to_litter"
     write (fid, *) veg_cf%m_frootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_livestemc_to_litter"
     write (fid, *) veg_cf%m_livestemc_to_litter
     write (fid, "(A)") "veg_cf%m_livestemc_storage_to_litter"
     write (fid, *) veg_cf%m_livestemc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_livestemc_xfer_to_litter"
     write (fid, *) veg_cf%m_livestemc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_deadstemc_to_litter"
     write (fid, *) veg_cf%m_deadstemc_to_litter
     write (fid, "(A)") "veg_cf%m_deadstemc_storage_to_litter"
     write (fid, *) veg_cf%m_deadstemc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_deadstemc_xfer_to_litter"
     write (fid, *) veg_cf%m_deadstemc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_livecrootc_to_litter"
     write (fid, *) veg_cf%m_livecrootc_to_litter
     write (fid, "(A)") "veg_cf%m_livecrootc_storage_to_litter"
     write (fid, *) veg_cf%m_livecrootc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_livecrootc_xfer_to_litter"
     write (fid, *) veg_cf%m_livecrootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_deadcrootc_to_litter"
     write (fid, *) veg_cf%m_deadcrootc_to_litter
     write (fid, "(A)") "veg_cf%m_deadcrootc_storage_to_litter"
     write (fid, *) veg_cf%m_deadcrootc_storage_to_litter
     write (fid, "(A)") "veg_cf%m_deadcrootc_xfer_to_litter"
     write (fid, *) veg_cf%m_deadcrootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_gresp_storage_to_litter"
     write (fid, *) veg_cf%m_gresp_storage_to_litter
     write (fid, "(A)") "veg_cf%m_gresp_xfer_to_litter"
     write (fid, *) veg_cf%m_gresp_xfer_to_litter
     write (fid, "(A)") "veg_cf%m_cpool_to_litter"
     write (fid, *) veg_cf%m_cpool_to_litter
     write (fid, "(A)") "veg_cf%hrv_leafc_to_litter"
     write (fid, *) veg_cf%hrv_leafc_to_litter
     write (fid, "(A)") "veg_cf%hrv_leafc_storage_to_litter"
     write (fid, *) veg_cf%hrv_leafc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_leafc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_leafc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_frootc_to_litter"
     write (fid, *) veg_cf%hrv_frootc_to_litter
     write (fid, "(A)") "veg_cf%hrv_frootc_storage_to_litter"
     write (fid, *) veg_cf%hrv_frootc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_frootc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_frootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_livestemc_to_litter"
     write (fid, *) veg_cf%hrv_livestemc_to_litter
     write (fid, "(A)") "veg_cf%hrv_livestemc_storage_to_litter"
     write (fid, *) veg_cf%hrv_livestemc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_livestemc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_livestemc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_deadstemc_to_prod10c"
     write (fid, *) veg_cf%hrv_deadstemc_to_prod10c
     write (fid, "(A)") "veg_cf%hrv_deadstemc_to_prod100c"
     write (fid, *) veg_cf%hrv_deadstemc_to_prod100c
     write (fid, "(A)") "veg_cf%hrv_deadstemc_storage_to_litter"
     write (fid, *) veg_cf%hrv_deadstemc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_deadstemc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_deadstemc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_livecrootc_to_litter"
     write (fid, *) veg_cf%hrv_livecrootc_to_litter
     write (fid, "(A)") "veg_cf%hrv_livecrootc_storage_to_litter"
     write (fid, *) veg_cf%hrv_livecrootc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_livecrootc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_livecrootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_deadcrootc_to_litter"
     write (fid, *) veg_cf%hrv_deadcrootc_to_litter
     write (fid, "(A)") "veg_cf%hrv_deadcrootc_storage_to_litter"
     write (fid, *) veg_cf%hrv_deadcrootc_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_deadcrootc_xfer_to_litter"
     write (fid, *) veg_cf%hrv_deadcrootc_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_gresp_storage_to_litter"
     write (fid, *) veg_cf%hrv_gresp_storage_to_litter
     write (fid, "(A)") "veg_cf%hrv_gresp_xfer_to_litter"
     write (fid, *) veg_cf%hrv_gresp_xfer_to_litter
     write (fid, "(A)") "veg_cf%hrv_xsmrpool_to_atm"
     write (fid, *) veg_cf%hrv_xsmrpool_to_atm
     write (fid, "(A)") "veg_cf%hrv_cpool_to_litter"
     write (fid, *) veg_cf%hrv_cpool_to_litter
     write (fid, "(A)") "veg_cf%hrv_leafc_to_prod1c"
     write (fid, *) veg_cf%hrv_leafc_to_prod1c
     write (fid, "(A)") "veg_cf%hrv_livestemc_to_prod1c"
     write (fid, *) veg_cf%hrv_livestemc_to_prod1c
     write (fid, "(A)") "veg_cf%hrv_grainc_to_prod1c"
     write (fid, *) veg_cf%hrv_grainc_to_prod1c
     write (fid, "(A)") "veg_cf%m_leafc_to_fire"
     write (fid, *) veg_cf%m_leafc_to_fire
     write (fid, "(A)") "veg_cf%m_leafc_storage_to_fire"
     write (fid, *) veg_cf%m_leafc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_leafc_xfer_to_fire"
     write (fid, *) veg_cf%m_leafc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_livestemc_to_fire"
     write (fid, *) veg_cf%m_livestemc_to_fire
     write (fid, "(A)") "veg_cf%m_livestemc_storage_to_fire"
     write (fid, *) veg_cf%m_livestemc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_livestemc_xfer_to_fire"
     write (fid, *) veg_cf%m_livestemc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_to_fire"
     write (fid, *) veg_cf%m_deadstemc_to_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_storage_to_fire"
     write (fid, *) veg_cf%m_deadstemc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_xfer_to_fire"
     write (fid, *) veg_cf%m_deadstemc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_frootc_to_fire"
     write (fid, *) veg_cf%m_frootc_to_fire
     write (fid, "(A)") "veg_cf%m_frootc_storage_to_fire"
     write (fid, *) veg_cf%m_frootc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_frootc_xfer_to_fire"
     write (fid, *) veg_cf%m_frootc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_to_fire"
     write (fid, *) veg_cf%m_livecrootc_to_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_storage_to_fire"
     write (fid, *) veg_cf%m_livecrootc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_xfer_to_fire"
     write (fid, *) veg_cf%m_livecrootc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_to_fire"
     write (fid, *) veg_cf%m_deadcrootc_to_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_storage_to_fire"
     write (fid, *) veg_cf%m_deadcrootc_storage_to_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_xfer_to_fire"
     write (fid, *) veg_cf%m_deadcrootc_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_gresp_storage_to_fire"
     write (fid, *) veg_cf%m_gresp_storage_to_fire
     write (fid, "(A)") "veg_cf%m_gresp_xfer_to_fire"
     write (fid, *) veg_cf%m_gresp_xfer_to_fire
     write (fid, "(A)") "veg_cf%m_cpool_to_fire"
     write (fid, *) veg_cf%m_cpool_to_fire
     write (fid, "(A)") "veg_cf%m_leafc_to_litter_fire"
     write (fid, *) veg_cf%m_leafc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_leafc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_leafc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_leafc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_leafc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livestemc_to_litter_fire"
     write (fid, *) veg_cf%m_livestemc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livestemc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_livestemc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livestemc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_livestemc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_to_litter_fire"
     write (fid, *) veg_cf%m_deadstemc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_deadstemc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadstemc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_deadstemc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_frootc_to_litter_fire"
     write (fid, *) veg_cf%m_frootc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_frootc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_frootc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_frootc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_frootc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_to_litter_fire"
     write (fid, *) veg_cf%m_livecrootc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_livecrootc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_livecrootc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_livecrootc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_to_litter_fire"
     write (fid, *) veg_cf%m_deadcrootc_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_storage_to_litter_fire"
     write (fid, *) veg_cf%m_deadcrootc_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_deadcrootc_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_deadcrootc_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_gresp_storage_to_litter_fire"
     write (fid, *) veg_cf%m_gresp_storage_to_litter_fire
     write (fid, "(A)") "veg_cf%m_gresp_xfer_to_litter_fire"
     write (fid, *) veg_cf%m_gresp_xfer_to_litter_fire
     write (fid, "(A)") "veg_cf%m_cpool_to_litter_fire"
     write (fid, *) veg_cf%m_cpool_to_litter_fire
     write (fid, "(A)") "veg_cf%grainc_xfer_to_grainc"
     write (fid, *) veg_cf%grainc_xfer_to_grainc
     write (fid, "(A)") "veg_cf%leafc_xfer_to_leafc"
     write (fid, *) veg_cf%leafc_xfer_to_leafc
     write (fid, "(A)") "veg_cf%frootc_xfer_to_frootc"
     write (fid, *) veg_cf%frootc_xfer_to_frootc
     write (fid, "(A)") "veg_cf%livestemc_xfer_to_livestemc"
     write (fid, *) veg_cf%livestemc_xfer_to_livestemc
     write (fid, "(A)") "veg_cf%deadstemc_xfer_to_deadstemc"
     write (fid, *) veg_cf%deadstemc_xfer_to_deadstemc
     write (fid, "(A)") "veg_cf%livecrootc_xfer_to_livecrootc"
     write (fid, *) veg_cf%livecrootc_xfer_to_livecrootc
     write (fid, "(A)") "veg_cf%deadcrootc_xfer_to_deadcrootc"
     write (fid, *) veg_cf%deadcrootc_xfer_to_deadcrootc
     write (fid, "(A)") "veg_cf%leafc_to_litter"
     write (fid, *) veg_cf%leafc_to_litter
     write (fid, "(A)") "veg_cf%frootc_to_litter"
     write (fid, *) veg_cf%frootc_to_litter
     write (fid, "(A)") "veg_cf%livestemc_to_litter"
     write (fid, *) veg_cf%livestemc_to_litter
     write (fid, "(A)") "veg_cf%grainc_to_food"
     write (fid, *) veg_cf%grainc_to_food
     write (fid, "(A)") "veg_cf%leaf_mr"
     write (fid, *) veg_cf%leaf_mr
     write (fid, "(A)") "veg_cf%froot_mr"
     write (fid, *) veg_cf%froot_mr
     write (fid, "(A)") "veg_cf%livestem_mr"
     write (fid, *) veg_cf%livestem_mr
     write (fid, "(A)") "veg_cf%livecroot_mr"
     write (fid, *) veg_cf%livecroot_mr
     write (fid, "(A)") "veg_cf%grain_mr"
     write (fid, *) veg_cf%grain_mr
     write (fid, "(A)") "veg_cf%leaf_curmr"
     write (fid, *) veg_cf%leaf_curmr
     write (fid, "(A)") "veg_cf%froot_curmr"
     write (fid, *) veg_cf%froot_curmr
     write (fid, "(A)") "veg_cf%livestem_curmr"
     write (fid, *) veg_cf%livestem_curmr
     write (fid, "(A)") "veg_cf%livecroot_curmr"
     write (fid, *) veg_cf%livecroot_curmr
     write (fid, "(A)") "veg_cf%grain_curmr"
     write (fid, *) veg_cf%grain_curmr
     write (fid, "(A)") "veg_cf%leaf_xsmr"
     write (fid, *) veg_cf%leaf_xsmr
     write (fid, "(A)") "veg_cf%froot_xsmr"
     write (fid, *) veg_cf%froot_xsmr
     write (fid, "(A)") "veg_cf%livestem_xsmr"
     write (fid, *) veg_cf%livestem_xsmr
     write (fid, "(A)") "veg_cf%livecroot_xsmr"
     write (fid, *) veg_cf%livecroot_xsmr
     write (fid, "(A)") "veg_cf%grain_xsmr"
     write (fid, *) veg_cf%grain_xsmr
     write (fid, "(A)") "veg_cf%xr"
     write (fid, *) veg_cf%xr
     write (fid, "(A)") "veg_cf%psnsun_to_cpool"
     write (fid, *) veg_cf%psnsun_to_cpool
     write (fid, "(A)") "veg_cf%psnshade_to_cpool"
     write (fid, *) veg_cf%psnshade_to_cpool
     write (fid, "(A)") "veg_cf%cpool_to_xsmrpool"
     write (fid, *) veg_cf%cpool_to_xsmrpool
     write (fid, "(A)") "veg_cf%cpool_to_grainc"
     write (fid, *) veg_cf%cpool_to_grainc
     write (fid, "(A)") "veg_cf%cpool_to_grainc_storage"
     write (fid, *) veg_cf%cpool_to_grainc_storage
     write (fid, "(A)") "veg_cf%cpool_to_leafc"
     write (fid, *) veg_cf%cpool_to_leafc
     write (fid, "(A)") "veg_cf%cpool_to_leafc_storage"
     write (fid, *) veg_cf%cpool_to_leafc_storage
     write (fid, "(A)") "veg_cf%cpool_to_frootc"
     write (fid, *) veg_cf%cpool_to_frootc
     write (fid, "(A)") "veg_cf%cpool_to_frootc_storage"
     write (fid, *) veg_cf%cpool_to_frootc_storage
     write (fid, "(A)") "veg_cf%cpool_to_livestemc"
     write (fid, *) veg_cf%cpool_to_livestemc
     write (fid, "(A)") "veg_cf%cpool_to_livestemc_storage"
     write (fid, *) veg_cf%cpool_to_livestemc_storage
     write (fid, "(A)") "veg_cf%cpool_to_deadstemc"
     write (fid, *) veg_cf%cpool_to_deadstemc
     write (fid, "(A)") "veg_cf%cpool_to_deadstemc_storage"
     write (fid, *) veg_cf%cpool_to_deadstemc_storage
     write (fid, "(A)") "veg_cf%cpool_to_livecrootc"
     write (fid, *) veg_cf%cpool_to_livecrootc
     write (fid, "(A)") "veg_cf%cpool_to_livecrootc_storage"
     write (fid, *) veg_cf%cpool_to_livecrootc_storage
     write (fid, "(A)") "veg_cf%cpool_to_deadcrootc"
     write (fid, *) veg_cf%cpool_to_deadcrootc
     write (fid, "(A)") "veg_cf%cpool_to_deadcrootc_storage"
     write (fid, *) veg_cf%cpool_to_deadcrootc_storage
     write (fid, "(A)") "veg_cf%cpool_to_gresp_storage"
     write (fid, *) veg_cf%cpool_to_gresp_storage
     write (fid, "(A)") "veg_cf%xsmrpool_to_atm"
     write (fid, *) veg_cf%xsmrpool_to_atm
     write (fid, "(A)") "veg_cf%cpool_leaf_gr"
     write (fid, *) veg_cf%cpool_leaf_gr
     write (fid, "(A)") "veg_cf%cpool_leaf_storage_gr"
     write (fid, *) veg_cf%cpool_leaf_storage_gr
     write (fid, "(A)") "veg_cf%transfer_leaf_gr"
     write (fid, *) veg_cf%transfer_leaf_gr
     write (fid, "(A)") "veg_cf%cpool_froot_gr"
     write (fid, *) veg_cf%cpool_froot_gr
     write (fid, "(A)") "veg_cf%cpool_froot_storage_gr"
     write (fid, *) veg_cf%cpool_froot_storage_gr
     write (fid, "(A)") "veg_cf%transfer_froot_gr"
     write (fid, *) veg_cf%transfer_froot_gr
     write (fid, "(A)") "veg_cf%cpool_livestem_gr"
     write (fid, *) veg_cf%cpool_livestem_gr
     write (fid, "(A)") "veg_cf%cpool_livestem_storage_gr"
     write (fid, *) veg_cf%cpool_livestem_storage_gr
     write (fid, "(A)") "veg_cf%transfer_livestem_gr"
     write (fid, *) veg_cf%transfer_livestem_gr
     write (fid, "(A)") "veg_cf%cpool_deadstem_gr"
     write (fid, *) veg_cf%cpool_deadstem_gr
     write (fid, "(A)") "veg_cf%cpool_deadstem_storage_gr"
     write (fid, *) veg_cf%cpool_deadstem_storage_gr
     write (fid, "(A)") "veg_cf%transfer_deadstem_gr"
     write (fid, *) veg_cf%transfer_deadstem_gr
     write (fid, "(A)") "veg_cf%cpool_livecroot_gr"
     write (fid, *) veg_cf%cpool_livecroot_gr
     write (fid, "(A)") "veg_cf%cpool_livecroot_storage_gr"
     write (fid, *) veg_cf%cpool_livecroot_storage_gr
     write (fid, "(A)") "veg_cf%transfer_livecroot_gr"
     write (fid, *) veg_cf%transfer_livecroot_gr
     write (fid, "(A)") "veg_cf%cpool_deadcroot_gr"
     write (fid, *) veg_cf%cpool_deadcroot_gr
     write (fid, "(A)") "veg_cf%cpool_deadcroot_storage_gr"
     write (fid, *) veg_cf%cpool_deadcroot_storage_gr
     write (fid, "(A)") "veg_cf%transfer_deadcroot_gr"
     write (fid, *) veg_cf%transfer_deadcroot_gr
     write (fid, "(A)") "veg_cf%cpool_grain_gr"
     write (fid, *) veg_cf%cpool_grain_gr
     write (fid, "(A)") "veg_cf%cpool_grain_storage_gr"
     write (fid, *) veg_cf%cpool_grain_storage_gr
     write (fid, "(A)") "veg_cf%transfer_grain_gr"
     write (fid, *) veg_cf%transfer_grain_gr
     write (fid, "(A)") "veg_cf%grainc_storage_to_xfer"
     write (fid, *) veg_cf%grainc_storage_to_xfer
     write (fid, "(A)") "veg_cf%leafc_storage_to_xfer"
     write (fid, *) veg_cf%leafc_storage_to_xfer
     write (fid, "(A)") "veg_cf%frootc_storage_to_xfer"
     write (fid, *) veg_cf%frootc_storage_to_xfer
     write (fid, "(A)") "veg_cf%livestemc_storage_to_xfer"
     write (fid, *) veg_cf%livestemc_storage_to_xfer
     write (fid, "(A)") "veg_cf%deadstemc_storage_to_xfer"
     write (fid, *) veg_cf%deadstemc_storage_to_xfer
     write (fid, "(A)") "veg_cf%livecrootc_storage_to_xfer"
     write (fid, *) veg_cf%livecrootc_storage_to_xfer
     write (fid, "(A)") "veg_cf%deadcrootc_storage_to_xfer"
     write (fid, *) veg_cf%deadcrootc_storage_to_xfer
     write (fid, "(A)") "veg_cf%gresp_storage_to_xfer"
     write (fid, *) veg_cf%gresp_storage_to_xfer
     write (fid, "(A)") "veg_cf%livestemc_to_deadstemc"
     write (fid, *) veg_cf%livestemc_to_deadstemc
     write (fid, "(A)") "veg_cf%livecrootc_to_deadcrootc"
     write (fid, *) veg_cf%livecrootc_to_deadcrootc
     write (fid, "(A)") "veg_cf%gpp_before_downreg"
     write (fid, *) veg_cf%gpp_before_downreg
     write (fid, "(A)") "veg_cf%prev_leafc_to_litter"
     write (fid, *) veg_cf%prev_leafc_to_litter
     write (fid, "(A)") "veg_cf%prev_frootc_to_litter"
     write (fid, *) veg_cf%prev_frootc_to_litter
     write (fid, "(A)") "veg_cf%availc"
     write (fid, *) veg_cf%availc
     write (fid, "(A)") "veg_cf%xsmrpool_turnover"
     write (fid, *) veg_cf%xsmrpool_turnover
     write (fid, "(A)") "veg_cf%crop_seedc_to_leaf"
     write (fid, *) veg_cf%crop_seedc_to_leaf
     write (fid, "(A)") "veg_cf%tempsum_npp"
     write (fid, *) veg_cf%tempsum_npp
     write (fid, "(A)") "veg_cf%annsum_npp"
     write (fid, *) veg_cf%annsum_npp
     write (fid, "(A)") "veg_cf%allocation_leaf"
     write (fid, *) veg_cf%allocation_leaf
     write (fid, "(A)") "veg_cf%allocation_stem"
     write (fid, *) veg_cf%allocation_stem
     write (fid, "(A)") "veg_cf%allocation_froot"
     write (fid, *) veg_cf%allocation_froot

     !====================== veg_ef ======================!

     write (fid, "(A)") "veg_ef%eflx_soil_grnd"
     write (fid, *) veg_ef%eflx_soil_grnd
     write (fid, "(A)") "veg_ef%netrad"
     write (fid, *) veg_ef%netrad

     !====================== veg_pf ======================!

     write (fid, "(A)") "veg_pf%m_leafp_to_litter"
     write (fid, *) veg_pf%m_leafp_to_litter
     write (fid, "(A)") "veg_pf%m_frootp_to_litter"
     write (fid, *) veg_pf%m_frootp_to_litter
     write (fid, "(A)") "veg_pf%m_leafp_storage_to_litter"
     write (fid, *) veg_pf%m_leafp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_frootp_storage_to_litter"
     write (fid, *) veg_pf%m_frootp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_livestemp_storage_to_litter"
     write (fid, *) veg_pf%m_livestemp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_deadstemp_storage_to_litter"
     write (fid, *) veg_pf%m_deadstemp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_livecrootp_storage_to_litter"
     write (fid, *) veg_pf%m_livecrootp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_deadcrootp_storage_to_litter"
     write (fid, *) veg_pf%m_deadcrootp_storage_to_litter
     write (fid, "(A)") "veg_pf%m_leafp_xfer_to_litter"
     write (fid, *) veg_pf%m_leafp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_frootp_xfer_to_litter"
     write (fid, *) veg_pf%m_frootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_livestemp_xfer_to_litter"
     write (fid, *) veg_pf%m_livestemp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_deadstemp_xfer_to_litter"
     write (fid, *) veg_pf%m_deadstemp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_livecrootp_xfer_to_litter"
     write (fid, *) veg_pf%m_livecrootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_deadcrootp_xfer_to_litter"
     write (fid, *) veg_pf%m_deadcrootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%m_livestemp_to_litter"
     write (fid, *) veg_pf%m_livestemp_to_litter
     write (fid, "(A)") "veg_pf%m_deadstemp_to_litter"
     write (fid, *) veg_pf%m_deadstemp_to_litter
     write (fid, "(A)") "veg_pf%m_livecrootp_to_litter"
     write (fid, *) veg_pf%m_livecrootp_to_litter
     write (fid, "(A)") "veg_pf%m_deadcrootp_to_litter"
     write (fid, *) veg_pf%m_deadcrootp_to_litter
     write (fid, "(A)") "veg_pf%m_retransp_to_litter"
     write (fid, *) veg_pf%m_retransp_to_litter
     write (fid, "(A)") "veg_pf%m_ppool_to_litter"
     write (fid, *) veg_pf%m_ppool_to_litter
     write (fid, "(A)") "veg_pf%hrv_leafp_to_litter"
     write (fid, *) veg_pf%hrv_leafp_to_litter
     write (fid, "(A)") "veg_pf%hrv_frootp_to_litter"
     write (fid, *) veg_pf%hrv_frootp_to_litter
     write (fid, "(A)") "veg_pf%hrv_leafp_storage_to_litter"
     write (fid, *) veg_pf%hrv_leafp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_frootp_storage_to_litter"
     write (fid, *) veg_pf%hrv_frootp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_livestemp_storage_to_litter"
     write (fid, *) veg_pf%hrv_livestemp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadstemp_storage_to_litter"
     write (fid, *) veg_pf%hrv_deadstemp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_livecrootp_storage_to_litter"
     write (fid, *) veg_pf%hrv_livecrootp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadcrootp_storage_to_litter"
     write (fid, *) veg_pf%hrv_deadcrootp_storage_to_litter
     write (fid, "(A)") "veg_pf%hrv_leafp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_leafp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_frootp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_frootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_livestemp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_livestemp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadstemp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_deadstemp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_livecrootp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_livecrootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadcrootp_xfer_to_litter"
     write (fid, *) veg_pf%hrv_deadcrootp_xfer_to_litter
     write (fid, "(A)") "veg_pf%hrv_livestemp_to_litter"
     write (fid, *) veg_pf%hrv_livestemp_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadstemp_to_prod10p"
     write (fid, *) veg_pf%hrv_deadstemp_to_prod10p
     write (fid, "(A)") "veg_pf%hrv_deadstemp_to_prod100p"
     write (fid, *) veg_pf%hrv_deadstemp_to_prod100p
     write (fid, "(A)") "veg_pf%hrv_livecrootp_to_litter"
     write (fid, *) veg_pf%hrv_livecrootp_to_litter
     write (fid, "(A)") "veg_pf%hrv_deadcrootp_to_litter"
     write (fid, *) veg_pf%hrv_deadcrootp_to_litter
     write (fid, "(A)") "veg_pf%hrv_retransp_to_litter"
     write (fid, *) veg_pf%hrv_retransp_to_litter
     write (fid, "(A)") "veg_pf%hrv_ppool_to_litter"
     write (fid, *) veg_pf%hrv_ppool_to_litter
     write (fid, "(A)") "veg_pf%hrv_leafp_to_prod1p"
     write (fid, *) veg_pf%hrv_leafp_to_prod1p
     write (fid, "(A)") "veg_pf%hrv_livestemp_to_prod1p"
     write (fid, *) veg_pf%hrv_livestemp_to_prod1p
     write (fid, "(A)") "veg_pf%hrv_grainp_to_prod1p"
     write (fid, *) veg_pf%hrv_grainp_to_prod1p
     write (fid, "(A)") "veg_pf%grainp_xfer_to_grainp"
     write (fid, *) veg_pf%grainp_xfer_to_grainp
     write (fid, "(A)") "veg_pf%leafp_xfer_to_leafp"
     write (fid, *) veg_pf%leafp_xfer_to_leafp
     write (fid, "(A)") "veg_pf%frootp_xfer_to_frootp"
     write (fid, *) veg_pf%frootp_xfer_to_frootp
     write (fid, "(A)") "veg_pf%livestemp_xfer_to_livestemp"
     write (fid, *) veg_pf%livestemp_xfer_to_livestemp
     write (fid, "(A)") "veg_pf%deadstemp_xfer_to_deadstemp"
     write (fid, *) veg_pf%deadstemp_xfer_to_deadstemp
     write (fid, "(A)") "veg_pf%livecrootp_xfer_to_livecrootp"
     write (fid, *) veg_pf%livecrootp_xfer_to_livecrootp
     write (fid, "(A)") "veg_pf%deadcrootp_xfer_to_deadcrootp"
     write (fid, *) veg_pf%deadcrootp_xfer_to_deadcrootp
     write (fid, "(A)") "veg_pf%livestemp_to_litter"
     write (fid, *) veg_pf%livestemp_to_litter
     write (fid, "(A)") "veg_pf%grainp_to_food"
     write (fid, *) veg_pf%grainp_to_food
     write (fid, "(A)") "veg_pf%leafp_to_litter"
     write (fid, *) veg_pf%leafp_to_litter
     write (fid, "(A)") "veg_pf%leafp_to_retransp"
     write (fid, *) veg_pf%leafp_to_retransp
     write (fid, "(A)") "veg_pf%frootp_to_retransp"
     write (fid, *) veg_pf%frootp_to_retransp
     write (fid, "(A)") "veg_pf%frootp_to_litter"
     write (fid, *) veg_pf%frootp_to_litter
     write (fid, "(A)") "veg_pf%retransp_to_ppool"
     write (fid, *) veg_pf%retransp_to_ppool
     write (fid, "(A)") "veg_pf%sminp_to_ppool"
     write (fid, *) veg_pf%sminp_to_ppool
     write (fid, "(A)") "veg_pf%ppool_to_grainp"
     write (fid, *) veg_pf%ppool_to_grainp
     write (fid, "(A)") "veg_pf%ppool_to_grainp_storage"
     write (fid, *) veg_pf%ppool_to_grainp_storage
     write (fid, "(A)") "veg_pf%ppool_to_leafp"
     write (fid, *) veg_pf%ppool_to_leafp
     write (fid, "(A)") "veg_pf%ppool_to_leafp_storage"
     write (fid, *) veg_pf%ppool_to_leafp_storage
     write (fid, "(A)") "veg_pf%ppool_to_frootp"
     write (fid, *) veg_pf%ppool_to_frootp
     write (fid, "(A)") "veg_pf%ppool_to_frootp_storage"
     write (fid, *) veg_pf%ppool_to_frootp_storage
     write (fid, "(A)") "veg_pf%ppool_to_livestemp"
     write (fid, *) veg_pf%ppool_to_livestemp
     write (fid, "(A)") "veg_pf%ppool_to_livestemp_storage"
     write (fid, *) veg_pf%ppool_to_livestemp_storage
     write (fid, "(A)") "veg_pf%ppool_to_deadstemp"
     write (fid, *) veg_pf%ppool_to_deadstemp
     write (fid, "(A)") "veg_pf%ppool_to_deadstemp_storage"
     write (fid, *) veg_pf%ppool_to_deadstemp_storage
     write (fid, "(A)") "veg_pf%ppool_to_livecrootp"
     write (fid, *) veg_pf%ppool_to_livecrootp
     write (fid, "(A)") "veg_pf%ppool_to_livecrootp_storage"
     write (fid, *) veg_pf%ppool_to_livecrootp_storage
     write (fid, "(A)") "veg_pf%ppool_to_deadcrootp"
     write (fid, *) veg_pf%ppool_to_deadcrootp
     write (fid, "(A)") "veg_pf%ppool_to_deadcrootp_storage"
     write (fid, *) veg_pf%ppool_to_deadcrootp_storage
     write (fid, "(A)") "veg_pf%grainp_storage_to_xfer"
     write (fid, *) veg_pf%grainp_storage_to_xfer
     write (fid, "(A)") "veg_pf%leafp_storage_to_xfer"
     write (fid, *) veg_pf%leafp_storage_to_xfer
     write (fid, "(A)") "veg_pf%frootp_storage_to_xfer"
     write (fid, *) veg_pf%frootp_storage_to_xfer
     write (fid, "(A)") "veg_pf%livestemp_storage_to_xfer"
     write (fid, *) veg_pf%livestemp_storage_to_xfer
     write (fid, "(A)") "veg_pf%deadstemp_storage_to_xfer"
     write (fid, *) veg_pf%deadstemp_storage_to_xfer
     write (fid, "(A)") "veg_pf%livecrootp_storage_to_xfer"
     write (fid, *) veg_pf%livecrootp_storage_to_xfer
     write (fid, "(A)") "veg_pf%deadcrootp_storage_to_xfer"
     write (fid, *) veg_pf%deadcrootp_storage_to_xfer
     write (fid, "(A)") "veg_pf%livestemp_to_deadstemp"
     write (fid, *) veg_pf%livestemp_to_deadstemp
     write (fid, "(A)") "veg_pf%livestemp_to_retransp"
     write (fid, *) veg_pf%livestemp_to_retransp
     write (fid, "(A)") "veg_pf%livecrootp_to_deadcrootp"
     write (fid, *) veg_pf%livecrootp_to_deadcrootp
     write (fid, "(A)") "veg_pf%livecrootp_to_retransp"
     write (fid, *) veg_pf%livecrootp_to_retransp
     write (fid, "(A)") "veg_pf%biochem_pmin_to_plant"
     write (fid, *) veg_pf%biochem_pmin_to_plant
     write (fid, "(A)") "veg_pf%crop_seedp_to_leaf"
     write (fid, *) veg_pf%crop_seedp_to_leaf
     write (fid, "(A)") "veg_pf%plant_pdemand"
     write (fid, *) veg_pf%plant_pdemand
     write (fid, "(A)") "veg_pf%sminp_to_plant"
     write (fid, *) veg_pf%sminp_to_plant
     write (fid, "(A)") "veg_pf%supplement_to_plantp"
     write (fid, *) veg_pf%supplement_to_plantp

     !====================== c13_veg_cf ======================!

     write (fid, "(A)") "c13_veg_cf%psnsun_to_cpool"
     write (fid, *) c13_veg_cf%psnsun_to_cpool
     write (fid, "(A)") "c13_veg_cf%psnshade_to_cpool"
     write (fid, *) c13_veg_cf%psnshade_to_cpool

     !====================== c14_veg_cf ======================!

     write (fid, "(A)") "c14_veg_cf%psnsun_to_cpool"
     write (fid, *) c14_veg_cf%psnsun_to_cpool
     write (fid, "(A)") "c14_veg_cf%psnshade_to_cpool"
     write (fid, *) c14_veg_cf%psnshade_to_cpool

     !====================== c13_veg_cs ======================!


     !====================== c14_veg_cs ======================!

     write (fid, "(A)") "c14_veg_cs%leafc"
     write (fid, *) c14_veg_cs%leafc
     write (fid, "(A)") "c14_veg_cs%leafc_storage"
     write (fid, *) c14_veg_cs%leafc_storage
     write (fid, "(A)") "c14_veg_cs%leafc_xfer"
     write (fid, *) c14_veg_cs%leafc_xfer
     write (fid, "(A)") "c14_veg_cs%frootc"
     write (fid, *) c14_veg_cs%frootc
     write (fid, "(A)") "c14_veg_cs%frootc_storage"
     write (fid, *) c14_veg_cs%frootc_storage
     write (fid, "(A)") "c14_veg_cs%frootc_xfer"
     write (fid, *) c14_veg_cs%frootc_xfer
     write (fid, "(A)") "c14_veg_cs%livestemc"
     write (fid, *) c14_veg_cs%livestemc
     write (fid, "(A)") "c14_veg_cs%livestemc_storage"
     write (fid, *) c14_veg_cs%livestemc_storage
     write (fid, "(A)") "c14_veg_cs%livestemc_xfer"
     write (fid, *) c14_veg_cs%livestemc_xfer
     write (fid, "(A)") "c14_veg_cs%deadstemc"
     write (fid, *) c14_veg_cs%deadstemc
     write (fid, "(A)") "c14_veg_cs%deadstemc_storage"
     write (fid, *) c14_veg_cs%deadstemc_storage
     write (fid, "(A)") "c14_veg_cs%deadstemc_xfer"
     write (fid, *) c14_veg_cs%deadstemc_xfer
     write (fid, "(A)") "c14_veg_cs%livecrootc"
     write (fid, *) c14_veg_cs%livecrootc
     write (fid, "(A)") "c14_veg_cs%livecrootc_storage"
     write (fid, *) c14_veg_cs%livecrootc_storage
     write (fid, "(A)") "c14_veg_cs%livecrootc_xfer"
     write (fid, *) c14_veg_cs%livecrootc_xfer
     write (fid, "(A)") "c14_veg_cs%deadcrootc"
     write (fid, *) c14_veg_cs%deadcrootc
     write (fid, "(A)") "c14_veg_cs%deadcrootc_storage"
     write (fid, *) c14_veg_cs%deadcrootc_storage
     write (fid, "(A)") "c14_veg_cs%deadcrootc_xfer"
     write (fid, *) c14_veg_cs%deadcrootc_xfer
     write (fid, "(A)") "c14_veg_cs%gresp_storage"
     write (fid, *) c14_veg_cs%gresp_storage
     write (fid, "(A)") "c14_veg_cs%gresp_xfer"
     write (fid, *) c14_veg_cs%gresp_xfer
     write (fid, "(A)") "c14_veg_cs%cpool"
     write (fid, *) c14_veg_cs%cpool
     write (fid, "(A)") "c14_veg_cs%xsmrpool"
     write (fid, *) c14_veg_cs%xsmrpool
     write (fid, "(A)") "c14_veg_cs%ctrunc"
     write (fid, *) c14_veg_cs%ctrunc

     !====================== veg_ns ======================!

     write (fid, "(A)") "veg_ns%leafn"
     write (fid, *) veg_ns%leafn
     write (fid, "(A)") "veg_ns%leafn_storage"
     write (fid, *) veg_ns%leafn_storage
     write (fid, "(A)") "veg_ns%leafn_xfer"
     write (fid, *) veg_ns%leafn_xfer
     write (fid, "(A)") "veg_ns%frootn"
     write (fid, *) veg_ns%frootn
     write (fid, "(A)") "veg_ns%frootn_storage"
     write (fid, *) veg_ns%frootn_storage
     write (fid, "(A)") "veg_ns%frootn_xfer"
     write (fid, *) veg_ns%frootn_xfer
     write (fid, "(A)") "veg_ns%livestemn"
     write (fid, *) veg_ns%livestemn
     write (fid, "(A)") "veg_ns%livestemn_storage"
     write (fid, *) veg_ns%livestemn_storage
     write (fid, "(A)") "veg_ns%livestemn_xfer"
     write (fid, *) veg_ns%livestemn_xfer
     write (fid, "(A)") "veg_ns%deadstemn"
     write (fid, *) veg_ns%deadstemn
     write (fid, "(A)") "veg_ns%deadstemn_storage"
     write (fid, *) veg_ns%deadstemn_storage
     write (fid, "(A)") "veg_ns%deadstemn_xfer"
     write (fid, *) veg_ns%deadstemn_xfer
     write (fid, "(A)") "veg_ns%livecrootn"
     write (fid, *) veg_ns%livecrootn
     write (fid, "(A)") "veg_ns%livecrootn_storage"
     write (fid, *) veg_ns%livecrootn_storage
     write (fid, "(A)") "veg_ns%livecrootn_xfer"
     write (fid, *) veg_ns%livecrootn_xfer
     write (fid, "(A)") "veg_ns%deadcrootn"
     write (fid, *) veg_ns%deadcrootn
     write (fid, "(A)") "veg_ns%deadcrootn_storage"
     write (fid, *) veg_ns%deadcrootn_storage
     write (fid, "(A)") "veg_ns%deadcrootn_xfer"
     write (fid, *) veg_ns%deadcrootn_xfer
     write (fid, "(A)") "veg_ns%retransn"
     write (fid, *) veg_ns%retransn
     write (fid, "(A)") "veg_ns%npool"
     write (fid, *) veg_ns%npool
     write (fid, "(A)") "veg_ns%grainn"
     write (fid, *) veg_ns%grainn
     write (fid, "(A)") "veg_ns%grainn_storage"
     write (fid, *) veg_ns%grainn_storage
     write (fid, "(A)") "veg_ns%grainn_xfer"
     write (fid, *) veg_ns%grainn_xfer
     write (fid, "(A)") "veg_ns%cropseedn_deficit"
     write (fid, *) veg_ns%cropseedn_deficit
     write (fid, "(A)") "veg_ns%pnup_pfrootc"
     write (fid, *) veg_ns%pnup_pfrootc
     write (fid, "(A)") "veg_ns%benefit_pgpp_pleafc"
     write (fid, *) veg_ns%benefit_pgpp_pleafc

     !====================== veg_ps ======================!

     write (fid, "(A)") "veg_ps%grainp"
     write (fid, *) veg_ps%grainp
     write (fid, "(A)") "veg_ps%grainp_storage"
     write (fid, *) veg_ps%grainp_storage
     write (fid, "(A)") "veg_ps%grainp_xfer"
     write (fid, *) veg_ps%grainp_xfer
     write (fid, "(A)") "veg_ps%leafp"
     write (fid, *) veg_ps%leafp
     write (fid, "(A)") "veg_ps%leafp_storage"
     write (fid, *) veg_ps%leafp_storage
     write (fid, "(A)") "veg_ps%leafp_xfer"
     write (fid, *) veg_ps%leafp_xfer
     write (fid, "(A)") "veg_ps%frootp"
     write (fid, *) veg_ps%frootp
     write (fid, "(A)") "veg_ps%frootp_storage"
     write (fid, *) veg_ps%frootp_storage
     write (fid, "(A)") "veg_ps%frootp_xfer"
     write (fid, *) veg_ps%frootp_xfer
     write (fid, "(A)") "veg_ps%livestemp"
     write (fid, *) veg_ps%livestemp
     write (fid, "(A)") "veg_ps%livestemp_storage"
     write (fid, *) veg_ps%livestemp_storage
     write (fid, "(A)") "veg_ps%livestemp_xfer"
     write (fid, *) veg_ps%livestemp_xfer
     write (fid, "(A)") "veg_ps%deadstemp"
     write (fid, *) veg_ps%deadstemp
     write (fid, "(A)") "veg_ps%deadstemp_storage"
     write (fid, *) veg_ps%deadstemp_storage
     write (fid, "(A)") "veg_ps%deadstemp_xfer"
     write (fid, *) veg_ps%deadstemp_xfer
     write (fid, "(A)") "veg_ps%livecrootp"
     write (fid, *) veg_ps%livecrootp
     write (fid, "(A)") "veg_ps%livecrootp_storage"
     write (fid, *) veg_ps%livecrootp_storage
     write (fid, "(A)") "veg_ps%livecrootp_xfer"
     write (fid, *) veg_ps%livecrootp_xfer
     write (fid, "(A)") "veg_ps%deadcrootp"
     write (fid, *) veg_ps%deadcrootp
     write (fid, "(A)") "veg_ps%deadcrootp_storage"
     write (fid, *) veg_ps%deadcrootp_storage
     write (fid, "(A)") "veg_ps%deadcrootp_xfer"
     write (fid, *) veg_ps%deadcrootp_xfer
     write (fid, "(A)") "veg_ps%retransp"
     write (fid, *) veg_ps%retransp
     write (fid, "(A)") "veg_ps%ppool"
     write (fid, *) veg_ps%ppool
     write (fid, "(A)") "veg_ps%cropseedp_deficit"
     write (fid, *) veg_ps%cropseedp_deficit

     !====================== veg_vp ======================!

     write (fid, "(A)") "veg_vp%leafcn"
     write (fid, *) veg_vp%leafcn
     write (fid, "(A)") "veg_vp%woody"
     write (fid, *) veg_vp%woody
     write (fid, "(A)") "veg_vp%lflitcn"
     write (fid, *) veg_vp%lflitcn
     write (fid, "(A)") "veg_vp%frootcn"
     write (fid, *) veg_vp%frootcn
     write (fid, "(A)") "veg_vp%livewdcn"
     write (fid, *) veg_vp%livewdcn
     write (fid, "(A)") "veg_vp%deadwdcn"
     write (fid, *) veg_vp%deadwdcn
     write (fid, "(A)") "veg_vp%graincn"
     write (fid, *) veg_vp%graincn
     write (fid, "(A)") "veg_vp%froot_leaf"
     write (fid, *) veg_vp%froot_leaf
     write (fid, "(A)") "veg_vp%stem_leaf"
     write (fid, *) veg_vp%stem_leaf
     write (fid, "(A)") "veg_vp%croot_stem"
     write (fid, *) veg_vp%croot_stem
     write (fid, "(A)") "veg_vp%flivewd"
     write (fid, *) veg_vp%flivewd
     write (fid, "(A)") "veg_vp%fcur"
     write (fid, *) veg_vp%fcur
     write (fid, "(A)") "veg_vp%lf_flab"
     write (fid, *) veg_vp%lf_flab
     write (fid, "(A)") "veg_vp%lf_fcel"
     write (fid, *) veg_vp%lf_fcel
     write (fid, "(A)") "veg_vp%lf_flig"
     write (fid, *) veg_vp%lf_flig
     write (fid, "(A)") "veg_vp%fr_flab"
     write (fid, *) veg_vp%fr_flab
     write (fid, "(A)") "veg_vp%fr_fcel"
     write (fid, *) veg_vp%fr_fcel
     write (fid, "(A)") "veg_vp%fr_flig"
     write (fid, *) veg_vp%fr_flig
     write (fid, "(A)") "veg_vp%leaf_long"
     write (fid, *) veg_vp%leaf_long
     write (fid, "(A)") "veg_vp%froot_long"
     write (fid, *) veg_vp%froot_long
     write (fid, "(A)") "veg_vp%evergreen"
     write (fid, *) veg_vp%evergreen
     write (fid, "(A)") "veg_vp%stress_decid"
     write (fid, *) veg_vp%stress_decid
     write (fid, "(A)") "veg_vp%season_decid"
     write (fid, *) veg_vp%season_decid
     write (fid, "(A)") "veg_vp%fertnitro"
     write (fid, *) veg_vp%fertnitro
     write (fid, "(A)") "veg_vp%fleafcn"
     write (fid, *) veg_vp%fleafcn
     write (fid, "(A)") "veg_vp%ffrootcn"
     write (fid, *) veg_vp%ffrootcn
     write (fid, "(A)") "veg_vp%fstemcn"
     write (fid, *) veg_vp%fstemcn
     write (fid, "(A)") "veg_vp%presharv"
     write (fid, *) veg_vp%presharv
     write (fid, "(A)") "veg_vp%convfact"
     write (fid, *) veg_vp%convfact
     write (fid, "(A)") "veg_vp%fyield"
     write (fid, *) veg_vp%fyield
     write (fid, "(A)") "veg_vp%leafcp"
     write (fid, *) veg_vp%leafcp
     write (fid, "(A)") "veg_vp%lflitcp"
     write (fid, *) veg_vp%lflitcp
     write (fid, "(A)") "veg_vp%frootcp"
     write (fid, *) veg_vp%frootcp
     write (fid, "(A)") "veg_vp%livewdcp"
     write (fid, *) veg_vp%livewdcp
     write (fid, "(A)") "veg_vp%deadwdcp"
     write (fid, *) veg_vp%deadwdcp
     write (fid, "(A)") "veg_vp%graincp"
     write (fid, *) veg_vp%graincp
     write (fid, "(A)") "veg_vp%km_decomp_nh4"
     write (fid, *) veg_vp%km_decomp_nh4
     write (fid, "(A)") "veg_vp%km_decomp_no3"
     write (fid, *) veg_vp%km_decomp_no3
     write (fid, "(A)") "veg_vp%km_decomp_p"
     write (fid, *) veg_vp%km_decomp_p
     write (fid, "(A)") "veg_vp%km_nit"
     write (fid, *) veg_vp%km_nit
     write (fid, "(A)") "veg_vp%km_den"
     write (fid, *) veg_vp%km_den
     write (fid, "(A)") "veg_vp%km_ptase"
     write (fid, *) veg_vp%km_ptase
     write (fid, "(A)") "veg_vp%lamda_ptase"
     write (fid, *) veg_vp%lamda_ptase

     !====================== photosyns_vars ======================!

     write (fid, "(A)") "photosyns_vars%psnsun_patch"
     write (fid, *) photosyns_vars%psnsun_patch
     write (fid, "(A)") "photosyns_vars%psnsha_patch"
     write (fid, *) photosyns_vars%psnsha_patch
     write (fid, "(A)") "photosyns_vars%c13_psnsun_patch"
     write (fid, *) photosyns_vars%c13_psnsun_patch
     write (fid, "(A)") "photosyns_vars%c13_psnsha_patch"
     write (fid, *) photosyns_vars%c13_psnsha_patch
     write (fid, "(A)") "photosyns_vars%c14_psnsun_patch"
     write (fid, *) photosyns_vars%c14_psnsun_patch
     write (fid, "(A)") "photosyns_vars%c14_psnsha_patch"
     write (fid, *) photosyns_vars%c14_psnsha_patch
     write (fid, "(A)") "photosyns_vars%lmrsun_patch"
     write (fid, *) photosyns_vars%lmrsun_patch
     write (fid, "(A)") "photosyns_vars%lmrsha_patch"
     write (fid, *) photosyns_vars%lmrsha_patch

     !====================== frictionvel_vars ======================!


     !====================== atm2lnd_vars ======================!

     write (fid, "(A)") "atm2lnd_vars%forc_hdm"
     write (fid, *) atm2lnd_vars%forc_hdm
     write (fid, "(A)") "atm2lnd_vars%forc_lnfm"
     write (fid, *) atm2lnd_vars%forc_lnfm
     write (fid, "(A)") "atm2lnd_vars%forc_ndep_grc"
     write (fid, *) atm2lnd_vars%forc_ndep_grc
     write (fid, "(A)") "atm2lnd_vars%forc_pdep_grc"
     write (fid, *) atm2lnd_vars%forc_pdep_grc

     !====================== aerosol_vars ======================!


     !====================== ch4_vars ======================!

     write (fid, "(A)") "ch4_vars%finundated_col"
     write (fid, *) ch4_vars%finundated_col
     write (fid, "(A)") "ch4_vars%o2stress_unsat_col"
     write (fid, *) ch4_vars%o2stress_unsat_col
     write (fid, "(A)") "ch4_vars%o2stress_sat_col"
     write (fid, *) ch4_vars%o2stress_sat_col
     write (fid, "(A)") "ch4_vars%conc_o2_sat_col"
     write (fid, *) ch4_vars%conc_o2_sat_col
     write (fid, "(A)") "ch4_vars%conc_o2_unsat_col"
     write (fid, *) ch4_vars%conc_o2_unsat_col
     write (fid, "(A)") "ch4_vars%o2_decomp_depth_sat_col"
     write (fid, *) ch4_vars%o2_decomp_depth_sat_col
     write (fid, "(A)") "ch4_vars%o2_decomp_depth_unsat_col"
     write (fid, *) ch4_vars%o2_decomp_depth_unsat_col

     !====================== soilstate_vars ======================!

     write (fid, "(A)") "soilstate_vars%cellorg_col"
     write (fid, *) soilstate_vars%cellorg_col
     write (fid, "(A)") "soilstate_vars%bd_col"
     write (fid, *) soilstate_vars%bd_col
     write (fid, "(A)") "soilstate_vars%bsw_col"
     write (fid, *) soilstate_vars%bsw_col
     write (fid, "(A)") "soilstate_vars%watsat_col"
     write (fid, *) soilstate_vars%watsat_col
     write (fid, "(A)") "soilstate_vars%watfc_col"
     write (fid, *) soilstate_vars%watfc_col
     write (fid, "(A)") "soilstate_vars%sucsat_col"
     write (fid, *) soilstate_vars%sucsat_col
     write (fid, "(A)") "soilstate_vars%soilpsi_col"
     write (fid, *) soilstate_vars%soilpsi_col
     write (fid, "(A)") "soilstate_vars%rootfr_patch"
     write (fid, *) soilstate_vars%rootfr_patch

     !====================== surfalb_vars ======================!


     !====================== solarabs_vars ======================!


     !====================== surfrad_vars ======================!


     !====================== urbanparams_vars ======================!


     !====================== energyflux_vars ======================!

     write (fid, "(A)") "energyflux_vars%btran2_patch"
     write (fid, *) energyflux_vars%btran2_patch
     write (fid, "(A)") "energyflux_vars%rresis_patch"
     write (fid, *) energyflux_vars%rresis_patch

     !====================== top_pp ======================!

     !====================== cnstate_vars ======================!

     write (fid, "(A)") "cnstate_vars%lfpftd_patch"
     write (fid, *) cnstate_vars%lfpftd_patch
     write (fid, "(A)") "cnstate_vars%hdidx_patch"
     write (fid, *) cnstate_vars%hdidx_patch
     write (fid, "(A)") "cnstate_vars%cumvd_patch"
     write (fid, *) cnstate_vars%cumvd_patch
     write (fid, "(A)") "cnstate_vars%gddmaturity_patch"
     write (fid, *) cnstate_vars%gddmaturity_patch
     write (fid, "(A)") "cnstate_vars%huileaf_patch"
     write (fid, *) cnstate_vars%huileaf_patch
     write (fid, "(A)") "cnstate_vars%huigrain_patch"
     write (fid, *) cnstate_vars%huigrain_patch
     write (fid, "(A)") "cnstate_vars%aleaf_patch"
     write (fid, *) cnstate_vars%aleaf_patch
     write (fid, "(A)") "cnstate_vars%astem_patch"
     write (fid, *) cnstate_vars%astem_patch
     write (fid, "(A)") "cnstate_vars%peaklai_patch"
     write (fid, *) cnstate_vars%peaklai_patch
     write (fid, "(A)") "cnstate_vars%leaf_prof_patch"
     write (fid, *) cnstate_vars%leaf_prof_patch
     write (fid, "(A)") "cnstate_vars%froot_prof_patch"
     write (fid, *) cnstate_vars%froot_prof_patch
     write (fid, "(A)") "cnstate_vars%croot_prof_patch"
     write (fid, *) cnstate_vars%croot_prof_patch
     write (fid, "(A)") "cnstate_vars%stem_prof_patch"
     write (fid, *) cnstate_vars%stem_prof_patch
     write (fid, "(A)") "cnstate_vars%gdp_lf_col"
     write (fid, *) cnstate_vars%gdp_lf_col
     write (fid, "(A)") "cnstate_vars%peatf_lf_col"
     write (fid, *) cnstate_vars%peatf_lf_col
     write (fid, "(A)") "cnstate_vars%abm_lf_col"
     write (fid, *) cnstate_vars%abm_lf_col
     write (fid, "(A)") "cnstate_vars%fpi_vr_col"
     write (fid, *) cnstate_vars%fpi_vr_col
     write (fid, "(A)") "cnstate_vars%fpg_col"
     write (fid, *) cnstate_vars%fpg_col
     write (fid, "(A)") "cnstate_vars%isoilorder"
     write (fid, *) cnstate_vars%isoilorder
     write (fid, "(A)") "cnstate_vars%fpi_p_vr_col"
     write (fid, *) cnstate_vars%fpi_p_vr_col
     write (fid, "(A)") "cnstate_vars%fpg_p_col"
     write (fid, *) cnstate_vars%fpg_p_col
     write (fid, "(A)") "cnstate_vars%rf_decomp_cascade_col"
     write (fid, *) cnstate_vars%rf_decomp_cascade_col
     write (fid, "(A)") "cnstate_vars%pathfrac_decomp_cascade_col"
     write (fid, *) cnstate_vars%pathfrac_decomp_cascade_col
     write (fid, "(A)") "cnstate_vars%nfixation_prof_col"
     write (fid, *) cnstate_vars%nfixation_prof_col
     write (fid, "(A)") "cnstate_vars%ndep_prof_col"
     write (fid, *) cnstate_vars%ndep_prof_col
     write (fid, "(A)") "cnstate_vars%som_adv_coef_col"
     write (fid, *) cnstate_vars%som_adv_coef_col
     write (fid, "(A)") "cnstate_vars%som_diffus_coef_col"
     write (fid, *) cnstate_vars%som_diffus_coef_col
     write (fid, "(A)") "cnstate_vars%tempavg_t2m_patch"
     write (fid, *) cnstate_vars%tempavg_t2m_patch
     write (fid, "(A)") "cnstate_vars%annavg_t2m_patch"
     write (fid, *) cnstate_vars%annavg_t2m_patch
     write (fid, "(A)") "cnstate_vars%annavg_t2m_col"
     write (fid, *) cnstate_vars%annavg_t2m_col
     write (fid, "(A)") "cnstate_vars%scalaravg_col"
     write (fid, *) cnstate_vars%scalaravg_col
     write (fid, "(A)") "cnstate_vars%annsum_counter_col"
     write (fid, *) cnstate_vars%annsum_counter_col
     write (fid, "(A)") "cnstate_vars%lfc_col"
     write (fid, *) cnstate_vars%lfc_col
     write (fid, "(A)") "cnstate_vars%dtrotr_col"
     write (fid, *) cnstate_vars%dtrotr_col
     write (fid, "(A)") "cnstate_vars%trotr1_col"
     write (fid, *) cnstate_vars%trotr1_col
     write (fid, "(A)") "cnstate_vars%trotr2_col"
     write (fid, *) cnstate_vars%trotr2_col
     write (fid, "(A)") "cnstate_vars%cropf_col"
     write (fid, *) cnstate_vars%cropf_col
     write (fid, "(A)") "cnstate_vars%baf_crop_col"
     write (fid, *) cnstate_vars%baf_crop_col
     write (fid, "(A)") "cnstate_vars%baf_peatf_col"
     write (fid, *) cnstate_vars%baf_peatf_col
     write (fid, "(A)") "cnstate_vars%fbac_col"
     write (fid, *) cnstate_vars%fbac_col
     write (fid, "(A)") "cnstate_vars%fbac1_col"
     write (fid, *) cnstate_vars%fbac1_col
     write (fid, "(A)") "cnstate_vars%farea_burned_col"
     write (fid, *) cnstate_vars%farea_burned_col
     write (fid, "(A)") "cnstate_vars%onset_flag_patch"
     write (fid, *) cnstate_vars%onset_flag_patch
     write (fid, "(A)") "cnstate_vars%onset_counter_patch"
     write (fid, *) cnstate_vars%onset_counter_patch
     write (fid, "(A)") "cnstate_vars%onset_gddflag_patch"
     write (fid, *) cnstate_vars%onset_gddflag_patch
     write (fid, "(A)") "cnstate_vars%onset_fdd_patch"
     write (fid, *) cnstate_vars%onset_fdd_patch
     write (fid, "(A)") "cnstate_vars%onset_gdd_patch"
     write (fid, *) cnstate_vars%onset_gdd_patch
     write (fid, "(A)") "cnstate_vars%offset_flag_patch"
     write (fid, *) cnstate_vars%offset_flag_patch
     write (fid, "(A)") "cnstate_vars%offset_counter_patch"
     write (fid, *) cnstate_vars%offset_counter_patch
     write (fid, "(A)") "cnstate_vars%offset_fdd_patch"
     write (fid, *) cnstate_vars%offset_fdd_patch
     write (fid, "(A)") "cnstate_vars%offset_swi_patch"
     write (fid, *) cnstate_vars%offset_swi_patch
     write (fid, "(A)") "cnstate_vars%grain_flag_patch"
     write (fid, *) cnstate_vars%grain_flag_patch
     write (fid, "(A)") "cnstate_vars%bglfr_leaf_patch"
     write (fid, *) cnstate_vars%bglfr_leaf_patch
     write (fid, "(A)") "cnstate_vars%bglfr_froot_patch"
     write (fid, *) cnstate_vars%bglfr_froot_patch
     write (fid, "(A)") "cnstate_vars%bgtr_patch"
     write (fid, *) cnstate_vars%bgtr_patch
     write (fid, "(A)") "cnstate_vars%c_allometry_patch"
     write (fid, *) cnstate_vars%c_allometry_patch
     write (fid, "(A)") "cnstate_vars%n_allometry_patch"
     write (fid, *) cnstate_vars%n_allometry_patch
     write (fid, "(A)") "cnstate_vars%tempsum_potential_gpp_patch"
     write (fid, *) cnstate_vars%tempsum_potential_gpp_patch
     write (fid, "(A)") "cnstate_vars%annsum_potential_gpp_patch"
     write (fid, *) cnstate_vars%annsum_potential_gpp_patch
     write (fid, "(A)") "cnstate_vars%tempmax_retransn_patch"
     write (fid, *) cnstate_vars%tempmax_retransn_patch
     write (fid, "(A)") "cnstate_vars%annmax_retransn_patch"
     write (fid, *) cnstate_vars%annmax_retransn_patch
     write (fid, "(A)") "cnstate_vars%p_allometry_patch"
     write (fid, *) cnstate_vars%p_allometry_patch
     write (fid, "(A)") "cnstate_vars%tempmax_retransp_patch"
     write (fid, *) cnstate_vars%tempmax_retransp_patch
     write (fid, "(A)") "cnstate_vars%annmax_retransp_patch"
     write (fid, *) cnstate_vars%annmax_retransp_patch
     write (fid, "(A)") "cnstate_vars%fpg_nh4_vr_col"
     write (fid, *) cnstate_vars%fpg_nh4_vr_col
     write (fid, "(A)") "cnstate_vars%cn_scalar"
     write (fid, *) cnstate_vars%cn_scalar
     write (fid, "(A)") "cnstate_vars%cp_scalar"
     write (fid, *) cnstate_vars%cp_scalar
     write (fid, "(A)") "cnstate_vars%cn_scalar_runmean"
     write (fid, *) cnstate_vars%cn_scalar_runmean
     write (fid, "(A)") "cnstate_vars%cp_scalar_runmean"
     write (fid, *) cnstate_vars%cp_scalar_runmean
     write (fid, "(A)") "cnstate_vars%r_mort_cal_patch"
     write (fid, *) cnstate_vars%r_mort_cal_patch

     !====================== lakestate_vars ======================!


     !====================== dust_vars ======================!


     !====================== lun_es ======================!


     !====================== lun_ef ======================!


     !====================== lun_ws ======================!


     !====================== soilhydrology_vars ======================!

     write (fid, "(A)") "soilhydrology_vars%fsat_col"
     write (fid, *) soilhydrology_vars%fsat_col

     !====================== crop_vars ======================!

     write (fid, "(A)") "crop_vars%nyrs_crop_active_patch"
     write (fid, *) crop_vars%nyrs_crop_active_patch
     write (fid, "(A)") "crop_vars%croplive_patch"
     write (fid, *) crop_vars%croplive_patch
     write (fid, "(A)") "crop_vars%harvdate_patch"
     write (fid, *) crop_vars%harvdate_patch
     write (fid, "(A)") "crop_vars%gddplant_patch"
     write (fid, *) crop_vars%gddplant_patch
     write (fid, "(A)") "crop_vars%gddtsoi_patch"
     write (fid, *) crop_vars%gddtsoi_patch
     write (fid, "(A)") "crop_vars%cvt_patch"
     write (fid, *) crop_vars%cvt_patch
     write (fid, "(A)") "crop_vars%cvp_patch"
     write (fid, *) crop_vars%cvp_patch
     write (fid, "(A)") "crop_vars%plantmonth_patch"
     write (fid, *) crop_vars%plantmonth_patch
     write (fid, "(A)") "crop_vars%xt_bar_patch"
     write (fid, *) crop_vars%xt_bar_patch
     write (fid, "(A)") "crop_vars%xp_bar_patch"
     write (fid, *) crop_vars%xp_bar_patch

     !====================== decomp_cascade_con ======================!

     write (fid, "(A)") "decomp_cascade_con%cascade_step_name"
     write (fid, *) decomp_cascade_con%cascade_step_name
     write (fid, "(A)") "decomp_cascade_con%cascade_donor_pool"
     write (fid, *) decomp_cascade_con%cascade_donor_pool
     write (fid, "(A)") "decomp_cascade_con%cascade_receiver_pool"
     write (fid, *) decomp_cascade_con%cascade_receiver_pool
     write (fid, "(A)") "decomp_cascade_con%floating_cn_ratio_decomp_pools"
     write (fid, *) decomp_cascade_con%floating_cn_ratio_decomp_pools
     write (fid, "(A)") "decomp_cascade_con%floating_cp_ratio_decomp_pools"
     write (fid, *) decomp_cascade_con%floating_cp_ratio_decomp_pools
     write (fid, "(A)") "decomp_cascade_con%decomp_pool_name_restart"
     write (fid, *) decomp_cascade_con%decomp_pool_name_restart
     write (fid, "(A)") "decomp_cascade_con%decomp_pool_name_history"
     write (fid, *) decomp_cascade_con%decomp_pool_name_history
     write (fid, "(A)") "decomp_cascade_con%decomp_pool_name_long"
     write (fid, *) decomp_cascade_con%decomp_pool_name_long
     write (fid, "(A)") "decomp_cascade_con%decomp_pool_name_short"
     write (fid, *) decomp_cascade_con%decomp_pool_name_short
     write (fid, "(A)") "decomp_cascade_con%is_litter"
     write (fid, *) decomp_cascade_con%is_litter
     write (fid, "(A)") "decomp_cascade_con%is_soil"
     write (fid, *) decomp_cascade_con%is_soil
     write (fid, "(A)") "decomp_cascade_con%is_cwd"
     write (fid, *) decomp_cascade_con%is_cwd
     write (fid, "(A)") "decomp_cascade_con%initial_cn_ratio"
     write (fid, *) decomp_cascade_con%initial_cn_ratio
     write (fid, "(A)") "decomp_cascade_con%initial_cp_ratio"
     write (fid, *) decomp_cascade_con%initial_cp_ratio
     write (fid, "(A)") "decomp_cascade_con%initial_stock"
     write (fid, *) decomp_cascade_con%initial_stock
     write (fid, "(A)") "decomp_cascade_con%is_metabolic"
     write (fid, *) decomp_cascade_con%is_metabolic
     write (fid, "(A)") "decomp_cascade_con%is_cellulose"
     write (fid, *) decomp_cascade_con%is_cellulose
     write (fid, "(A)") "decomp_cascade_con%is_lignin"
     write (fid, *) decomp_cascade_con%is_lignin
     write (fid, "(A)") "decomp_cascade_con%spinup_factor"
     write (fid, *) decomp_cascade_con%spinup_factor
     write (fid, "(A)") "decomp_cascade_con%decomp_k_pools"
     write (fid, *) decomp_cascade_con%decomp_k_pools


     !====================== lnd2atm_vars ======================!


     !====================== glc2lnd_vars ======================!


     !====================== ldomain ======================!


     !====================== params_inst ======================!


     !====================== sedflux_vars ======================!


     !====================== drydepvel_vars ======================!


     !====================== prior_weights ======================!


     !====================== patch_state_updater ======================!


     !====================== column_state_updater ======================!

     call fio_close(fid)
end subroutine write_vars
end module writeMod
