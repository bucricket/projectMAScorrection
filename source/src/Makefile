TOPDIR=..

include $(TOPDIR)/build/Makefile.inc


default: test coef_io other mw_scatt mw_scatt_coef coef_io_789 wrapper emis_atlas brdf_atlas gui 


all: bin

bin: main/bin hdf/bin parallel/bin coef_io/bin mw_scatt/bin other/bin emis_atlas/bin coef_io_789/bin mw_scatt_coef/bin brdf_atlas/bin test/bin wrapper/bin gui/bin


lib: main/lib hdf/lib parallel/lib coef_io/lib mw_scatt/lib other/lib emis_atlas/lib coef_io_789/lib mw_scatt_coef/lib brdf_atlas/lib test/lib wrapper/lib gui/lib


intf: main/intf hdf/intf parallel/intf coef_io/intf mw_scatt/intf other/intf emis_atlas/intf coef_io_789/intf mw_scatt_coef/intf brdf_atlas/intf test/intf wrapper/intf gui/intf


dir: main/dir hdf/dir parallel/dir coef_io/dir mw_scatt/dir other/dir emis_atlas/dir coef_io_789/dir mw_scatt_coef/dir brdf_atlas/dir test/dir wrapper/dir gui/dir



mw_scatt: mw_scatt/bin

mw_scatt/bin: mw_scatt/lib
	cd mw_scatt && make bin

mw_scatt/lib: mw_scatt/intf coef_io/lib hdf/lib main/lib
	cd mw_scatt && make lib

mw_scatt/intf: mw_scatt/dir
	cd mw_scatt && make intf

mw_scatt/dir:
	cd mw_scatt && make dir

mw_scatt/clean:
	cd mw_scatt && make clean

coef_io: coef_io/bin

coef_io/bin: coef_io/lib
	cd coef_io && make bin

coef_io/lib: coef_io/intf main/lib hdf/lib
	cd coef_io && make lib

coef_io/intf: coef_io/dir
	cd coef_io && make intf

coef_io/dir:
	cd coef_io && make dir

coef_io/clean:
	cd coef_io && make clean

hdf: hdf/bin

hdf/bin: hdf/lib
	cd hdf && make bin

hdf/lib: hdf/intf main/lib
	cd hdf && make lib

hdf/intf: hdf/dir
	cd hdf && make intf

hdf/dir:
	cd hdf && make dir

hdf/clean:
	cd hdf && make clean

wrapper: wrapper/bin

wrapper/bin: wrapper/lib
	cd wrapper && make bin

wrapper/lib: wrapper/intf other/lib parallel/lib emis_atlas/lib brdf_atlas/lib coef_io/lib hdf/lib main/lib
	cd wrapper && make lib

wrapper/intf: wrapper/dir
	cd wrapper && make intf

wrapper/dir:
	cd wrapper && make dir

wrapper/clean:
	cd wrapper && make clean

other: other/bin

other/bin: other/lib
	cd other && make bin

other/lib: other/intf main/lib hdf/lib coef_io/lib
	cd other && make lib

other/intf: other/dir
	cd other && make intf

other/dir:
	cd other && make dir

other/clean:
	cd other && make clean

mw_scatt_coef: mw_scatt_coef/bin

mw_scatt_coef/bin: mw_scatt_coef/lib
	cd mw_scatt_coef && make bin

mw_scatt_coef/lib: mw_scatt_coef/intf mw_scatt/lib main/lib hdf/lib coef_io/lib
	cd mw_scatt_coef && make lib

mw_scatt_coef/intf: mw_scatt_coef/dir
	cd mw_scatt_coef && make intf

mw_scatt_coef/dir:
	cd mw_scatt_coef && make dir

mw_scatt_coef/clean:
	cd mw_scatt_coef && make clean

parallel: parallel/bin

parallel/bin: parallel/lib
	cd parallel && make bin

parallel/lib: parallel/intf main/lib
	cd parallel && make lib

parallel/intf: parallel/dir
	cd parallel && make intf

parallel/dir:
	cd parallel && make dir

parallel/clean:
	cd parallel && make clean

emis_atlas: emis_atlas/bin

emis_atlas/bin: emis_atlas/lib
	cd emis_atlas && make bin

emis_atlas/lib: emis_atlas/intf hdf/lib coef_io/lib main/lib
	cd emis_atlas && make lib

emis_atlas/intf: emis_atlas/dir
	cd emis_atlas && make intf

emis_atlas/dir:
	cd emis_atlas && make dir

emis_atlas/clean:
	cd emis_atlas && make clean

brdf_atlas: brdf_atlas/bin

brdf_atlas/bin: brdf_atlas/lib
	cd brdf_atlas && make bin

brdf_atlas/lib: brdf_atlas/intf main/lib coef_io/lib hdf/lib emis_atlas/lib other/lib parallel/lib
	cd brdf_atlas && make lib

brdf_atlas/intf: brdf_atlas/dir
	cd brdf_atlas && make intf

brdf_atlas/dir:
	cd brdf_atlas && make dir

brdf_atlas/clean:
	cd brdf_atlas && make clean

test: test/bin

test/bin: test/lib
	cd test && make bin

test/lib: test/intf mw_scatt/lib parallel/lib other/lib main/lib hdf/lib coef_io/lib
	cd test && make lib

test/intf: test/dir
	cd test && make intf

test/dir:
	cd test && make dir

test/clean:
	cd test && make clean

coef_io_789: coef_io_789/bin

coef_io_789/bin: coef_io_789/lib
	cd coef_io_789 && make bin

coef_io_789/lib: coef_io_789/intf hdf/lib coef_io/lib main/lib
	cd coef_io_789 && make lib

coef_io_789/intf: coef_io_789/dir
	cd coef_io_789 && make intf

coef_io_789/dir:
	cd coef_io_789 && make dir

coef_io_789/clean:
	cd coef_io_789 && make clean

gui: gui/bin

gui/bin: gui/lib
	cd gui && make bin

gui/lib: gui/intf main/lib hdf/lib coef_io/lib brdf_atlas/lib emis_atlas/lib parallel/lib other/lib
	cd gui && make lib

gui/intf: gui/dir
	cd gui && make intf

gui/dir:
	cd gui && make dir

gui/clean:
	cd gui && make clean

main: main/bin

main/bin: main/lib
	cd main && make bin

main/lib: main/intf 
	cd main && make lib

main/intf: main/dir
	cd main && make intf

main/dir:
	cd main && make dir

main/clean:
	cd main && make clean

clean:
	\rm -f $(OBJDIR)/*$(EXT_OBJ) $(LIBDIR)/*$(EXT_LIB) $(LIBDIR)/*$(EXT_SO) $(MODDIR)/*$(EXT_MOD) $(BINDIR)/*$(EXT_BIN) $(INCDIR)/*
	\rm -rf $(TMPDIR)/*

dist:
	cd .. ; tar cf - readme.txt docs gui src build data brdf_data emis_data rtcoef_rttov11 rttov_test wrapper --exclude src/dev --exclude src/lbl --exclude src/lbl-add --exclude src/baran --exclude "src/*/*.pl" --exclude "src/parallel/*.sh" --exclude "src/parallel/*.txt" --exclude "docs/*lbl*" --exclude "docs/rttov-compile.*" --exclude "rtcoef_rttov11/extract*srf*.*" --exclude "rtcoef_rttov11/plot_lbl_check.py" --exclude "docs/*.odt" --exclude "docs/*.doc" --exclude "docs/rttov_developer_guide.*" --exclude "docs/doxygen_doc_*" --exclude "docs/sharp_1983_*" --exclude gui/rttov_tests.tar.bz2 --exclude rttov_test/rttov_remote.pl --exclude "rttov_test/*_dev*" --exclude "rttov_test/*.tar.gz" --exclude "rtcoef_rttov11/rttov*pred*/*airs*" --exclude "rtcoef_rttov11/rttov*pred*/*iasi*" --exclude "rtcoef_rttov11/rttov*pred*/*cris*" --exclude "rtcoef_rttov11/rttov*pred*/*iris*" --exclude "rtcoef_rttov11/rttov*pred*/*mtg*irs*" --exclude "rtcoef_rttov11/cldaer/*.dat" --exclude "rtcoef_rttov11/mietable/*.dat" --exclude "rtcoef_rttov11/pc/*.dat" --exclude "emis_data/*" --exclude "brdf_data/*" --exclude "wrapper/doxygen_doc_*" --exclude "wrapper/*_dev*" | gzip > rttov_dist.tgz

distdev:
	cd .. ; tar cf - readme.txt docs gui src build data brdf_data emis_data rttov_test wrapper --exclude "docs/doxygen_doc_*" --exclude "rttov_test/*.tar.gz" --exclude "emis_data/*" --exclude "brdf_data/*" --exclude "wrapper/doxygen_doc_*" | gzip > rttov_dist_dev.tgz

distlbl:
	cd .. ; tar cf - docs/*lbl* src/lbl src/lbl-add rttov_lbl-data --exclude rttov_lbl-data/scripts_pascal --exclude rttov_lbl-data/rttov_lbl_instrument.kcarta --exclude rttov_lbl-data/instrument-data.kcarta | gzip > rttov_dist_lbl.tgz
