LIBNAME=gnlm

OBJS=	dist.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(SHLIBLD) $(SHLIBLDFLAGS) $(FLIBS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
