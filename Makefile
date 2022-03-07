export OPENGL=0
REBSRC = /home/sara/Rebound/rebound/src/
include $(REBSRC)Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I $(REBSRC) -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REBSRC)
	@-rm -f librebound.so
	@ln -s $(REBSRC)librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REBSRC) clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
