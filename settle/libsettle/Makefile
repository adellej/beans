# detect OS (Linux or OSX)
UNAME := $(shell uname)

LIBRARY = libsettle.so

# CSRCS  := $(shell find . -maxdepth 1 -name '*.c')
# CPPSRCS := $(shell find . -maxdepth 1 -name '*.cc')
CSRCS	:= nrutil.c gasdev.c root.c
CPPSRCS := settle.cc odeint.cc eos.cc spline.cc useful.cc
COBJS	:= $(CSRCS:.c=.o)
CPPOBJS := $(CPPSRCS:.cc=.o)
OBJS	:= $(COBJS) $(CPPOBJS)

CFLAGS		:= -g -Wall -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-variable -Ofast -fPIC -c
CPPFLAGS 	:= $(CFLAGS)
LDFLAGS		:= -Wall -Ofast -shared -fPIC

ifeq ($(UNAME), Linux)
	LDFLAGS := ${LDFLAGS} -Wl,-soname,${LIBRARY}
	CC 	= gcc
	CPPC 	= g++
	LD 	= g++
	INSTALL_DIR=/usr/local/lib
else
ifeq ($(UNAME), Darwin)
	LDFLAGS := ${LDFLAGS} -Wl,-install_name,${LIBRARY}
# add suitable install dir for OSX here
	CC 	= clang
	CPPC 	= clang++
	LD 	= clang++
	INSTALL_DIR=/usr/local/lib
else
	@echo "ERROR: unsupported platform, this ${LIBRARY} library can be build only on Linux or Mac!"
	exit 1
endif
endif

settle:	${OBJS}
	@echo "Linking $(LIBRARY) ..."
	$(LD) ${LDFLAGS} -o ${LIBRARY} ${OBJS}

.PHONY: clean
clean:
	rm -fv *.o

cleaner: clean
	rm -fv *.so *.pyc

install: settle
	 sudo rm -f ${INSTALL_DIR}/${LIBRARY}
	 sudo cp -av ${LIBRARY} ${INSTALL_DIR}
ifeq ($(UNAME), Linux)
	sudo chown root:root ${INSTALL_DIR}/${LIBRARY}
	sudo ldconfig
else
ifeq ($(UNAME), Darwin)
	sudo chown root:wheel ${INSTALL_DIR}/${LIBRARY}
endif
endif

.cc.o: 
	@echo Compiling ...
	$(CPPC) -c $(CPPFLAGS) $<

.c.o:
	@echo Compiling ...
	$(CC) -c $(CFLAGS) $<
