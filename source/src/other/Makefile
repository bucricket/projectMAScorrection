
TOPDIR=../..
SUBDIR=other

FFLAGS_ARCH_create_aer_clim_prof          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_aer_clim_prof           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_bpr_calc                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_bpr_dealloc             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_bpr_init                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_bpr_mod                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_calc_weighting_fn       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_coef_info               = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_mie_params_aer          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_mie_params_cld          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_obs_to_pc               = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_print_info              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_print_opts              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_print_profile           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scattering_mod          = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_us76_prof               = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_zutility                = $(FFLAGS_ARCH)


F90_INTF=rttov_aer_clim_prof.F90 rttov_bpr_calc.F90 rttov_bpr_dealloc.F90 rttov_bpr_init.F90 rttov_calc_weighting_fn.F90 rttov_print_info.F90 rttov_print_opts.F90 rttov_print_profile.F90 rttov_us76_prof.F90

INTF=$(INCDIR)/rttov_aer_clim_prof.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_calc_weighting_fn.interface $(INCDIR)/rttov_print_info.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_profile.interface $(INCDIR)/rttov_us76_prof.interface

LIB=$(LIBDIR)/librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_aer_clim_prof$(EXT_OBJ) $(OBJDIR)/rttov_bpr_calc$(EXT_OBJ) $(OBJDIR)/rttov_bpr_dealloc$(EXT_OBJ) $(OBJDIR)/rttov_bpr_init$(EXT_OBJ) $(OBJDIR)/rttov_calc_weighting_fn$(EXT_OBJ) $(OBJDIR)/rttov_print_info$(EXT_OBJ) $(OBJDIR)/rttov_print_opts$(EXT_OBJ) $(OBJDIR)/rttov_print_profile$(EXT_OBJ) $(OBJDIR)/rttov_us76_prof$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/rttov_bpr_mod$(EXT_OBJ) $(OBJDIR)/rttov_scattering_mod$(EXT_OBJ) $(OBJDIR)/rttov_zutility$(EXT_OBJ)

OBJBIN=$(OBJDIR)/create_aer_clim_prof$(EXT_OBJ) $(OBJDIR)/rttov_coef_info$(EXT_OBJ) $(OBJDIR)/rttov_mie_params_aer$(EXT_OBJ) $(OBJDIR)/rttov_mie_params_cld$(EXT_OBJ) $(OBJDIR)/rttov_obs_to_pc$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/create_aer_clim_prof$(EXT_BIN) $(BINDIR)/rttov_coef_info$(EXT_BIN) $(BINDIR)/rttov_mie_params_aer$(EXT_BIN) $(BINDIR)/rttov_mie_params_cld$(EXT_BIN) $(BINDIR)/rttov_obs_to_pc$(EXT_BIN)

MODS=$(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_zutility$(EXT_MOD)

include $(TOPDIR)/build/Makefile.inc

lib: dir intf $(LIB)

bin: lib $(BINS)

dir: 
	@mkdir -p $(TMPDIR)/$(SUBDIR) $(BINDIR) $(OBJDIR) $(INCDIR) $(MODDIR) $(LIBDIR)

intf: 
	@$(MKINTF) $(F90_INTF) $(INTF)
	@$(CPINCH) $(INCDIR)

$(LIB): $(OBJSUB) $(OBJMOD)
	\rm -f $(LIB)
	$(AR) $(LIB) $(OBJSUB) $(OBJMOD)

clean:
	\rm -f $(INTF) $(LIB) $(BINS) $(OBJSUB) $(OBJMOD) $(OBJBIN) $(MODS)


$(BINDIR)/create_aer_clim_prof$(EXT_BIN): $(OBJDIR)/create_aer_clim_prof$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/create_aer_clim_prof$(EXT_BIN) $(OBJDIR)/create_aer_clim_prof$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_coef_info$(EXT_BIN): $(OBJDIR)/rttov_coef_info$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_coef_info$(EXT_BIN) $(OBJDIR)/rttov_coef_info$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_mie_params_aer$(EXT_BIN): $(OBJDIR)/rttov_mie_params_aer$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_mie_params_aer$(EXT_BIN) $(OBJDIR)/rttov_mie_params_aer$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_mie_params_cld$(EXT_BIN): $(OBJDIR)/rttov_mie_params_cld$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_mie_params_cld$(EXT_BIN) $(OBJDIR)/rttov_mie_params_cld$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_obs_to_pc$(EXT_BIN): $(OBJDIR)/rttov_obs_to_pc$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_obs_to_pc$(EXT_BIN) $(OBJDIR)/rttov_obs_to_pc$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov_aer_clim_prof$(EXT_OBJ): rttov_aer_clim_prof.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_aer_clim_prof) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_aer_clim_prof.F90 -o $(OBJDIR)/rttov_aer_clim_prof$(EXT_OBJ)


$(OBJDIR)/rttov_bpr_calc$(EXT_OBJ): rttov_bpr_calc.F90 $(INCDIR)/rttov_bpr_init.interface $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_bpr_calc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_bpr_calc.F90 -o $(OBJDIR)/rttov_bpr_calc$(EXT_OBJ)


$(OBJDIR)/rttov_bpr_dealloc$(EXT_OBJ): rttov_bpr_dealloc.F90  $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_bpr_dealloc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_bpr_dealloc.F90 -o $(OBJDIR)/rttov_bpr_dealloc$(EXT_OBJ)


$(OBJDIR)/rttov_bpr_init$(EXT_OBJ): rttov_bpr_init.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_bpr_init) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_bpr_init.F90 -o $(OBJDIR)/rttov_bpr_init$(EXT_OBJ)


$(OBJDIR)/rttov_calc_weighting_fn$(EXT_OBJ): rttov_calc_weighting_fn.F90  $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_calc_weighting_fn) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_calc_weighting_fn.F90 -o $(OBJDIR)/rttov_calc_weighting_fn$(EXT_OBJ)


$(OBJDIR)/rttov_print_info$(EXT_OBJ): rttov_print_info.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_print_info) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_print_info.F90 -o $(OBJDIR)/rttov_print_info$(EXT_OBJ)


$(OBJDIR)/rttov_print_opts$(EXT_OBJ): rttov_print_opts.F90  $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_print_opts) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_print_opts.F90 -o $(OBJDIR)/rttov_print_opts$(EXT_OBJ)


$(OBJDIR)/rttov_print_profile$(EXT_OBJ): rttov_print_profile.F90  $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_print_profile) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_print_profile.F90 -o $(OBJDIR)/rttov_print_profile$(EXT_OBJ)


$(OBJDIR)/rttov_us76_prof$(EXT_OBJ): rttov_us76_prof.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_us76_prof) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_us76_prof.F90 -o $(OBJDIR)/rttov_us76_prof$(EXT_OBJ)


$(OBJDIR)/create_aer_clim_prof$(EXT_OBJ): create_aer_clim_prof.F90 $(INCDIR)/rttov_aer_clim_prof.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_create_aer_clim_prof) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/create_aer_clim_prof.F90 -o $(OBJDIR)/create_aer_clim_prof$(EXT_OBJ)


$(OBJDIR)/rttov_coef_info$(EXT_OBJ): rttov_coef_info.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_print_info.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_coef_info) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_coef_info.F90 -o $(OBJDIR)/rttov_coef_info$(EXT_OBJ)


$(OBJDIR)/rttov_mie_params_aer$(EXT_OBJ): rttov_mie_params_aer.F90 $(INCDIR)/rttov_read_ascii_coef.interface $(INCDIR)/rttov_read_binary_coef.interface $(INCDIR)/rttov_init_coef.interface $(INCDIR)/rttov_nullify_coef_scatt_ir.interface $(INCDIR)/rttov_write_ascii_scaercoef.interface $(INCDIR)/rttov_write_binary_scaercoef.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_sskin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_mie_params_aer) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_mie_params_aer.F90 -o $(OBJDIR)/rttov_mie_params_aer$(EXT_OBJ)


$(OBJDIR)/rttov_mie_params_cld$(EXT_OBJ): rttov_mie_params_cld.F90 $(INCDIR)/rttov_read_ascii_coef.interface $(INCDIR)/rttov_read_binary_coef.interface $(INCDIR)/rttov_init_coef.interface $(INCDIR)/rttov_nullify_coef_scatt_ir.interface $(INCDIR)/rttov_write_ascii_sccldcoef.interface $(INCDIR)/rttov_write_binary_sccldcoef.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_hdf_load.interface $(INCDIR)/rttov_hdf_save.interface $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance2_io$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_pccomp_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_transmission_io$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_sskin_io$(EXT_MOD) $(MODDIR)/rttov_hdf_radiance_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_mie_params_cld) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_mie_params_cld.F90 -o $(OBJDIR)/rttov_mie_params_cld$(EXT_OBJ)


$(OBJDIR)/rttov_obs_to_pc$(EXT_OBJ): rttov_obs_to_pc.F90 $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_math_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_obs_to_pc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_obs_to_pc.F90 -o $(OBJDIR)/rttov_obs_to_pc$(EXT_OBJ)


$(OBJDIR)/rttov_bpr_mod$(EXT_OBJ): rttov_bpr_mod.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_bpr_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_bpr_mod.F90 -o $(OBJDIR)/rttov_bpr_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD)

$(MODDIR)/rttov_bpr_mod$(EXT_MOD): rttov_bpr_mod.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_bpr_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_bpr_mod.F90 -o $(OBJDIR)/rttov_bpr_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD)


$(OBJDIR)/rttov_scattering_mod$(EXT_OBJ): rttov_scattering_mod.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scattering_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scattering_mod.F90 -o $(OBJDIR)/rttov_scattering_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD)

$(MODDIR)/rttov_scattering_mod$(EXT_MOD): rttov_scattering_mod.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scattering_mod) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scattering_mod.F90 -o $(OBJDIR)/rttov_scattering_mod$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_scattering_mod$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD)


$(OBJDIR)/rttov_zutility$(EXT_OBJ): rttov_zutility.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_zutility) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_zutility.F90 -o $(OBJDIR)/rttov_zutility$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_zutility$(EXT_MOD) $(MODDIR)/rttov_zutility$(EXT_MOD)

$(MODDIR)/rttov_zutility$(EXT_MOD): rttov_zutility.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_zutility) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_zutility.F90 -o $(OBJDIR)/rttov_zutility$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_zutility$(EXT_MOD) $(MODDIR)/rttov_zutility$(EXT_MOD)

