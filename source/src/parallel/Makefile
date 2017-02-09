
TOPDIR=../..
SUBDIR=parallel

FFLAGS_ARCH_rttov_parallel_ad           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_parallel_direct       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_parallel_k            = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_parallel_tl           = $(FFLAGS_ARCH)


F90_INTF=rttov_parallel_ad.F90 rttov_parallel_direct.F90 rttov_parallel_k.F90 rttov_parallel_tl.F90

INTF=$(INCDIR)/rttov_parallel_ad.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_parallel_tl.interface

LIB=$(LIBDIR)/librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_parallel_ad$(EXT_OBJ) $(OBJDIR)/rttov_parallel_direct$(EXT_OBJ) $(OBJDIR)/rttov_parallel_k$(EXT_OBJ) $(OBJDIR)/rttov_parallel_tl$(EXT_OBJ) 

OBJMOD=

OBJBIN=

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=

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


$(OBJDIR)/rttov_parallel_ad$(EXT_OBJ): rttov_parallel_ad.F90 $(INCDIR)/rttov_ad.interface $(INCDIR)/rttov_alloc_prof.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_add_prof.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_parallel_ad) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_parallel_ad.F90 -o $(OBJDIR)/rttov_parallel_ad$(EXT_OBJ)


$(OBJDIR)/rttov_parallel_direct$(EXT_OBJ): rttov_parallel_direct.F90 $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_parallel_direct) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_parallel_direct.F90 -o $(OBJDIR)/rttov_parallel_direct$(EXT_OBJ)


$(OBJDIR)/rttov_parallel_k$(EXT_OBJ): rttov_parallel_k.F90 $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_parallel_k) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_parallel_k.F90 -o $(OBJDIR)/rttov_parallel_k$(EXT_OBJ)


$(OBJDIR)/rttov_parallel_tl$(EXT_OBJ): rttov_parallel_tl.F90 $(INCDIR)/rttov_tl.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_parallel_tl) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_parallel_tl.F90 -o $(OBJDIR)/rttov_parallel_tl$(EXT_OBJ)

