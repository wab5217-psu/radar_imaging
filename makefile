DEBUG=0


ifeq (1, ${DEBUG})
	CFLAGS=-L /usr/local/lib -I$(RSTPATH)/include/base -I$(RSTPATH)/include/superdarn -I$(RSTPATH)/include/general -I$(RSTPATH)/include/analysis -L$(RSTPATH)/lib -L/opt/intel -L/usr/local/lib -L/Users/bristow/lib -O3 -ipo -g \
	$(info ********DEBUG*********)
else

CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -O3 -ipo  -qmkl=parallel -mcmodel medium -fPIC -Wall -D_GNU_SOURCE -D_LINUX -Wno-unused-but-set-variable

endif


LIBS=-lm


OBJS=vel_est.o calculate_beam_matrix.o rd_image.o corr_lag_mtx.o estimate_amp.o estimate_vel.o pinv_c.o select_data.o x_corr_lag_mtx.o isgood.o


.c.o:
	icx $(CFLAGS) -v -c -o $@ $<


vel_est:	$(OBJS)
	icx -o ~/bin/vel_est $(CFLAGS) $(OBJS) $(LIBS)
