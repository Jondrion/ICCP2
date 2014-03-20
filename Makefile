FC = gfortran
FFLAGS =
LDFLAGS = 
LIBS = -llapack
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)


COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = 
OBJS += parameters.o
OBJS += polymer.o
OBJS += poly.o

all: poly

poly: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) poly $(OBJS) *.mod
