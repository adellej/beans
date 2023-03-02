# detect OS (Linux or OSX)
UNAME := $(shell uname)

LIBRARY = libsettle.so

CSRCS	:= $(shell find . -maxdepth 1 -name '*.c')
CPPSRCS := $(shell find . -maxdepth 2 -name '*.cc')
COBJS	:= $(CSRCS:.c=.o)
CPPOBJS := $(CPPSRCS:.cc=.o)
OBJS	:= $(COBJS) $(CPPOBJS)

CFLAGS		:= -g -Wall -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-variable -Ofast -fPIC -c
CPPFLAGS 	:= $(CFLAGS)
LDFLAGS		:= -Wall -Ofast -shared -fPIC

CC 	= gcc
CPPC 	= g++
LD 	= g++

ifeq ($(UNAME), Linux)
	LDFLAGS := ${LDFLAGS} -Wl,-soname,${LIBRARY}
	INSTALL_DIR=/usr/local/lib
else
ifeq ($(UNAME), Darwin)
	LDFLAGS := ${LDFLAGS} -Wl,-install_name,${LIBRARY}
# add suitable install dir for OSX here
	INSTALL_DIR=.
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
	sudo cp -arv --preserve=mode,timestamps ${LIBRARY} ${INSTALL_DIR}
	sudo ldconfig

test:	install
	time python try.py

.cc.o: 
	@echo Compiling ...
	$(CPPC) -c $(CPPFLAGS) $<

.c.o:
	@echo Compiling ...
	$(CC) -c $(CFLAGS) $<
