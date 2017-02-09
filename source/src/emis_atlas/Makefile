
TOPDIR=../..
SUBDIR=emis_atlas

FFLAGS_ARCH_mod_cnrm_mw_atlas                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_iratlas                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_mwatlas                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_rttov_emis_atlas              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_cnrmmwatlas_test            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_deallocate_emis_atlas       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_get_emis                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_iratlas_test                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_mwatlas_test                = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_setup_emis_atlas            = $(FFLAGS_ARCH)


F90_INTF=rttov_deallocate_emis_atlas.F90 rttov_get_emis.F90 rttov_setup_emis_atlas.F90

INTF=$(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_setup_emis_atlas.interface

LIB=$(LIBDIR)/librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_deallocate_emis_atlas$(EXT_OBJ) $(OBJDIR)/rttov_get_emis$(EXT_OBJ) $(OBJDIR)/rttov_setup_emis_atlas$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/mod_cnrm_mw_atlas$(EXT_OBJ) $(OBJDIR)/mod_iratlas$(EXT_OBJ) $(OBJDIR)/mod_mwatlas$(EXT_OBJ) $(OBJDIR)/mod_rttov_emis_atlas$(EXT_OBJ)

OBJBIN=$(OBJDIR)/rttov_cnrmmwatlas_test$(EXT_OBJ) $(OBJDIR)/rttov_iratlas_test$(EXT_OBJ) $(OBJDIR)/rttov_mwatlas_test$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov_cnrmmwatlas_test$(EXT_BIN) $(BINDIR)/rttov_iratlas_test$(EXT_BIN) $(BINDIR)/rttov_mwatlas_test$(EXT_BIN)

MODS=$(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD)

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


$(BINDIR)/rttov_cnrmmwatlas_test$(EXT_BIN): $(OBJDIR)/rttov_cnrmmwatlas_test$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_cnrmmwatlas_test$(EXT_BIN) $(OBJDIR)/rttov_cnrmmwatlas_test$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_iratlas_test$(EXT_BIN): $(OBJDIR)/rttov_iratlas_test$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_iratlas_test$(EXT_BIN) $(OBJDIR)/rttov_iratlas_test$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_mwatlas_test$(EXT_BIN): $(OBJDIR)/rttov_mwatlas_test$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_mwatlas_test$(EXT_BIN) $(OBJDIR)/rttov_mwatlas_test$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov_deallocate_emis_atlas$(EXT_OBJ): rttov_deallocate_emis_atlas.F90  $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_deallocate_emis_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_deallocate_emis_atlas.F90 -o $(OBJDIR)/rttov_deallocate_emis_atlas$(EXT_OBJ)


$(OBJDIR)/rttov_get_emis$(EXT_OBJ): rttov_get_emis.F90  $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_get_emis) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_get_emis.F90 -o $(OBJDIR)/rttov_get_emis$(EXT_OBJ)


$(OBJDIR)/rttov_setup_emis_atlas$(EXT_OBJ): rttov_setup_emis_atlas.F90  $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_setup_emis_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_setup_emis_atlas.F90 -o $(OBJDIR)/rttov_setup_emis_atlas$(EXT_OBJ)


$(OBJDIR)/rttov_cnrmmwatlas_test$(EXT_OBJ): rttov_cnrmmwatlas_test.F90 $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_cnrmmwatlas_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_cnrmmwatlas_test.F90 -o $(OBJDIR)/rttov_cnrmmwatlas_test$(EXT_OBJ)


$(OBJDIR)/rttov_iratlas_test$(EXT_OBJ): rttov_iratlas_test.F90 $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_iratlas_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_iratlas_test.F90 -o $(OBJDIR)/rttov_iratlas_test$(EXT_OBJ)


$(OBJDIR)/rttov_mwatlas_test$(EXT_OBJ): rttov_mwatlas_test.F90 $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_mwatlas_test) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_mwatlas_test.F90 -o $(OBJDIR)/rttov_mwatlas_test$(EXT_OBJ)


$(OBJDIR)/mod_cnrm_mw_atlas$(EXT_OBJ): mod_cnrm_mw_atlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_opencoeff.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_cnrm_mw_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_cnrm_mw_atlas.F90 -o $(OBJDIR)/mod_cnrm_mw_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD)

$(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD): mod_cnrm_mw_atlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_opencoeff.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_cnrm_mw_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_cnrm_mw_atlas.F90 -o $(OBJDIR)/mod_cnrm_mw_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD)


$(OBJDIR)/mod_iratlas$(EXT_OBJ): mod_iratlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_hdf_load.interface $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_sskin_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_iratlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_iratlas.F90 -o $(OBJDIR)/mod_iratlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_iratlas$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD)

$(MODDIR)/mod_iratlas$(EXT_MOD): mod_iratlas.F90 $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_hdf_load.interface $(MODDIR)/rttov_hdf_chanprof_io$(EXT_MOD) $(MODDIR)/rttov_hdf_sskin_io$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_hdf_options_interp_io$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_fast_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_profile_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_mw_io$(EXT_MOD) $(MODDIR)/rttov_hdf_reflectance_io$(EXT_MOD) $(MODDIR)/rttov_hdf_s2m_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_pc_io$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_ir_io$(EXT_MOD) $(MODDIR)/rttov_hdf_rttov_coef_io$(EXT_MOD) $(MODDIR)/rttov_hdf_coefs$(EXT_MOD) $(MODDIR)/rttov_hdf_options_io$(EXT_MOD) $(MODDIR)/rttov_hdf_emissivity_io$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/rttov_hdf_profiles$(EXT_MOD) $(MODDIR)/rttov_hdf_options_config_io$(EXT_MOD) $(MODDIR)/rttov_hdf_options_rt_all_io$(EXT_MOD) $(MODDIR)/rttov_hdf_opt_param_io$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_iratlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_iratlas.F90 -o $(OBJDIR)/mod_iratlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_iratlas$(EXT_MOD) $(MODDIR)/mod_iratlas$(EXT_MOD)


$(OBJDIR)/mod_mwatlas$(EXT_OBJ): mod_mwatlas.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_mwatlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_mwatlas.F90 -o $(OBJDIR)/mod_mwatlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD)

$(MODDIR)/mod_mwatlas$(EXT_MOD): mod_mwatlas.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_mwatlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_mwatlas.F90 -o $(OBJDIR)/mod_mwatlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_mwatlas$(EXT_MOD) $(MODDIR)/mod_mwatlas$(EXT_MOD)


$(OBJDIR)/mod_rttov_emis_atlas$(EXT_OBJ): mod_rttov_emis_atlas.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_rttov_emis_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_rttov_emis_atlas.F90 -o $(OBJDIR)/mod_rttov_emis_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD)

$(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD): mod_rttov_emis_atlas.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_rttov_emis_atlas) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_rttov_emis_atlas.F90 -o $(OBJDIR)/mod_rttov_emis_atlas$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD)

