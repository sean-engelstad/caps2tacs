#	Makefile for udpKulfan
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
BDIR  = $(ESP_ROOT)/bin
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

$(LDIR)/udpTwistFoil.so:	$(ODIR)/udpTwistFoil.o
	touch $(LDIR)/udpTwistFoil.so
	rm $(LDIR)/udpTwistFoil.so
	$(CC) $(SOFLGS) -o $(LDIR)/twistFoil.so $(ODIR)/udpTwistFoil.o \
		-L$(LDIR) -legads -lm $(RPATH)

$(ODIR)/udpTwistFoil.o:	udpTwistFoil.c $(IDIR)/udpUtilities.c $(IDIR)/udpUtilities.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I. \
		udpTwistFoil.c -o $(ODIR)/udpTwistFoil.o

test : $(LDIR)/udpTwistFoil.so
	@date
	@for f in design*.csm; do \
	    date; \
	    sensCSM $$f -geom -outLevel 0; \
	    printf "Finished sensCSM (config) %s\n\n" $$f; \
	done
	@date
	@printf "Complete\n"

clean:
	(cd $(ODIR); rm -f udpTwistFoil.o )

cleanall:	clean
	(cd $(LDIR); rm -f udpTwistFoil.so )
