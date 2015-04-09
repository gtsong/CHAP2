
# By default (if you just run "make"), CHAP is configured to keep its scripts
# and Java programs in the main "package directory" (where the distribution
# was unpacked), while compiled binaries and resource data files are located
# in the "bin" and "resources" subdirectories, respectively.  If you want to
# install it elsewhere, edit the following lines to specify the desired
# locations, and then run "make install".
#
CHAP_SCRIPT_DIR=/usr/local/chap
CHAP_JAVA_DIR=/usr/local/chap
CHAP_BINARY_DIR=/usr/local/chap
CHAP_RESOURCE_DIR=/usr/local/chap

INSTALL = install
TARGETS = build_lastz build_cage build_conversion build_post_process build_ortho_map build_utils \
          build_infer_annot build_piptools

build : clean bindir $(TARGETS)
	@echo '=== Build complete ==='

bindir :
	@echo ---$@---
	mkdir -p bin

build_lastz :
	@echo ---$@---
	cd lastz/src; make lastz
	cp lastz/src/lastz bin/lastz

build_cage :
	@echo ---$@---
	cd cage; make

build_conversion :
	@echo ---$@---
	cd conversion; make

build_post_process :
	@echo ---$@---
	cd post_process; make

build_ortho_map :
	@echo ---$@---
	cd ortho_map; make

build_infer_annot :
	@echo ---$@---
	cd infer_annot; make

build_piptools :
	@echo ---$@---
	cd piptools; make

build_utils :
	@echo ---$@---
	cd utils; make

clean :
	@echo ---$@---
	cd lastz/src; make clean
	rm -rf bin

install : build
	@echo ---$@---
	rm -rf bin/*.dSYM    # debug info on Macs
	umask 022; mkdir -p $(CHAP_SCRIPT_DIR) $(CHAP_JAVA_DIR) $(CHAP_BINARY_DIR) $(CHAP_RESOURCE_DIR)
	$(INSTALL) -m 755 *.sh $(CHAP_SCRIPT_DIR)
	$(INSTALL) -m 644 *.jar $(CHAP_JAVA_DIR)
	$(INSTALL) -m 755 bin/* $(CHAP_BINARY_DIR)
	@echo '====>  One of the resource files is large, so this may take some time.  <===='
	$(INSTALL) -m 644 resources/* $(CHAP_RESOURCE_DIR)
	for f in *.sh; do sed -i'~' -e 's!^SCRIPTS=.*!SCRIPTS=$(CHAP_SCRIPT_DIR)!' $(CHAP_SCRIPT_DIR)/$$f; done
	for f in *.sh; do sed -i'~' -e 's!^JARS=.*!JARS=$(CHAP_JAVA_DIR)!' $(CHAP_SCRIPT_DIR)/$$f; done
	for f in *.sh; do sed -i'~' -e 's!^BIN=.*!BIN=$(CHAP_BINARY_DIR)!' $(CHAP_SCRIPT_DIR)/$$f; done
	for f in *.sh; do sed -i'~' -e 's!^RESOURCES=.*!RESOURCES=$(CHAP_RESOURCE_DIR)!' $(CHAP_SCRIPT_DIR)/$$f; done
	for f in *.sh; do rm -f $(CHAP_SCRIPT_DIR)/$$f~; done
	@echo '=== Install complete ==='

