# Default format for making a file with only one cpp code which links to other headers infile
ALL_OBJS = $(OBJS0) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7) $(OBJS8)
ALL_PROG = $(PROG0) $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7) $(PROG8)

all: $(ALL_PROG)

$(PROG0): $(OBJS0)
	$(LD) $? $(LDFLAGS) -o $@

$(OBJS0): $(SRCS0)  
	$(CC) $(CFLAGS) -c $?
	
$(PROG1): $(OBJS1)
	$(LD) $? $(LDFLAGS) -o $@

$(OBJS1): $(SRCS1)  
	$(CC) $(CFLAGS) -c $?

$(PROG2): $(OBJS2)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS2): $(SRCS2)  
		$(CC) $(CFLAGS) -c $?

$(PROG3): $(OBJS3)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS3): $(SRCS3)  
		$(CC) $(CFLAGS) -c $?

$(PROG4): $(OBJS4)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS4): $(SRCS4)  
		$(CC) $(CFLAGS) -c $?	
		
$(PROG5): $(OBJS5)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS5): $(SRCS5)  
		$(CC) $(CFLAGS) -c $?
		
$(PROG6): $(OBJS6)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS6): $(SRCS6)  
		$(CC) $(CFLAGS) -c $?	

$(PROG7): $(OBJS7)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS7): $(SRCS7)  
		$(CC) $(CFLAGS) -c $?	

$(PROG8): $(OBJS8)
		$(LD) $? $(LDFLAGS) -o $@

$(OBJS8): $(SRCS8)  
		$(CC) $(CFLAGS) -c $?
		
clean:
	$(RM) $(ALL_OBJS) $(ALL_PROG)
	
info:
	@echo 'Hostname is $(HOSTNAME)'
	@echo 'Username is $(USER)'
	@echo 'Build type is $(BUILD_TYPE)'
	@echo 'operating system is $(MAKE_UNAME)'

