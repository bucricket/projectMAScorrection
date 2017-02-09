
TOPDIR=../..
SUBDIR=coef_io_789

FFLAGS_ARCH_rttov789_coeffname              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_conv_coef              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_conv_coef_11to9        = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_opencoeff              = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_readcoeffs             = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_readcoeffs_ascii       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov789_write_ascii_coef       = $(FFLAGS_ARCH)


F90_INTF=rttov789_coeffname.F90 rttov789_opencoeff.F90 rttov789_readcoeffs.F90 rttov789_readcoeffs_ascii.F90 rttov789_write_ascii_coef.F90

INTF=$(INCDIR)/rttov789_coeffname.interface $(INCDIR)/rttov789_opencoeff.interface $(INCDIR)/rttov789_readcoeffs.interface $(INCDIR)/rttov789_readcoeffs_ascii.interface $(INCDIR)/rttov789_write_ascii_coef.interface

LIB=$(LIBDIR)/librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov789_coeffname$(EXT_OBJ) $(OBJDIR)/rttov789_opencoeff$(EXT_OBJ) $(OBJDIR)/rttov789_readcoeffs$(EXT_OBJ) $(OBJDIR)/rttov789_readcoeffs_ascii$(EXT_OBJ) $(OBJDIR)/rttov789_write_ascii_coef$(EXT_OBJ) 

OBJMOD=

OBJBIN=$(OBJDIR)/rttov789_conv_coef$(EXT_OBJ) $(OBJDIR)/rttov789_conv_coef_11to9$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io_789$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov789_conv_coef$(EXT_BIN) $(BINDIR)/rttov789_conv_coef_11to9$(EXT_BIN)

MODS=

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


$(BINDIR)/rttov789_conv_coef$(EXT_BIN): $(OBJDIR)/rttov789_conv_coef$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov789_conv_coef$(EXT_BIN) $(OBJDIR)/rttov789_conv_coef$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_coef_io_789 -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov789_conv_coef_11to9$(EXT_BIN): $(OBJDIR)/rttov789_conv_coef_11to9$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov789_conv_coef_11to9$(EXT_BIN) $(OBJDIR)/rttov789_conv_coef_11to9$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_coef_io_789 -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/rttov789_coeffname$(EXT_OBJ): rttov789_coeffname.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_coeffname) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_coeffname.F90 -o $(OBJDIR)/rttov789_coeffname$(EXT_OBJ)


$(OBJDIR)/rttov789_opencoeff$(EXT_OBJ): rttov789_opencoeff.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_opencoeff) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_opencoeff.F90 -o $(OBJDIR)/rttov789_opencoeff$(EXT_OBJ)


$(OBJDIR)/rttov789_readcoeffs$(EXT_OBJ): rttov789_readcoeffs.F90 $(INCDIR)/rttov789_coeffname.interface $(INCDIR)/rttov789_opencoeff.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_nullify_coef.interface $(INCDIR)/rttov789_readcoeffs_ascii.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_readcoeffs) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_readcoeffs.F90 -o $(OBJDIR)/rttov789_readcoeffs$(EXT_OBJ)


$(OBJDIR)/rttov789_readcoeffs_ascii$(EXT_OBJ): rttov789_readcoeffs_ascii.F90 $(INCDIR)/rttov789_opencoeff.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_skipcommentline.interface $(INCDIR)/rttov_deletecomment.interface $(INCDIR)/rttov_cmpuc.interface $(INCDIR)/rttov_findnextsection.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_readcoeffs_ascii) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_readcoeffs_ascii.F90 -o $(OBJDIR)/rttov789_readcoeffs_ascii$(EXT_OBJ)


$(OBJDIR)/rttov789_write_ascii_coef$(EXT_OBJ): rttov789_write_ascii_coef.F90 $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_rttov_fastem3_coef$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_write_ascii_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_write_ascii_coef.F90 -o $(OBJDIR)/rttov789_write_ascii_coef$(EXT_OBJ)


$(OBJDIR)/rttov789_conv_coef$(EXT_OBJ): rttov789_conv_coef.F90 $(INCDIR)/rttov789_readcoeffs.interface $(INCDIR)/rttov_init_coef.interface $(INCDIR)/rttov_write_ascii_coef.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_conv_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_conv_coef.F90 -o $(OBJDIR)/rttov789_conv_coef$(EXT_OBJ)


$(OBJDIR)/rttov789_conv_coef_11to9$(EXT_OBJ): rttov789_conv_coef_11to9.F90 $(INCDIR)/rttov_read_ascii_coef.interface $(INCDIR)/rttov789_write_ascii_coef.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/mod_rttov_fastem3_coef$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov789_conv_coef_11to9) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov789_conv_coef_11to9.F90 -o $(OBJDIR)/rttov789_conv_coef_11to9$(EXT_OBJ)

