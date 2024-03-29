MV = /bin/mv
RM = /bin/rm -f

# c and lib flags
CFLAGS = # empty
LDFLAGS = # empty
	
ifndef CLUSTER

	#user names and file location, default is in home/Cpp/include
	ifeq ($(USER),jmgambet)
		# options set on my laptop
		CFLAGS+= -I $(HOME)/Dropbox/CppLib-felix/include
	else
		ifeq ($(USER),fmotzoi)
			# options set on my laptop
			CFLAGS+= -I /svndev/CppLib/include
		else 
			ifeq ($(USER),anastasiamctaggart)
				# options set on my laptop
				CFLAGS+= -I /home/anastasiamctaggart/InternshipSummer14/code/CppLib/include

			else
				# default
				CFLAGS+= -I $(HOME)/CppLib/include
			endif
		endif
	endif

	#build options, default is O3 and linux version of blas and Lapack
	#
	# Platform specific tweaks
	#
	OSTYPE := $(strip $(shell uname -s))
	
	ifneq (,$(strip $(findstring $(OSTYPE), Darwin)))
	  	# Mac OSX
		CC      = g++
		LD      = g++
	  	LDFLAGS += -framework Accelerate
	  	ifeq ($(BUILD_TYPE), debug)
			CFLAGS += -g -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -Wno-deprecated
		else
			CFLAGS += -fast -Wno-deprecated
		endif 
	else
		ifneq (,$(strip $(findstring $(OSTYPE), Linux)))
			# Linux
	  		CC      = clang++
			LD      = clang++
			LDFLAGS +=  -lblas -llapack 
	  		ifeq ($(BUILD_TYPE), debug)
				CFLAGS += -g -ansi -w -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -Wno-deprecated -std=c++0x
			else
				CFLAGS += -O3 -Wno-deprecated -std=c++0x
			endif
			
		else 
			# default 
	  		CC      = g++
			LD      = g++
			LDFLAGS += -lgsl -lgslcblas -llapack -lblas -lg2c -lm -m64 
			ifeq ($(BUILD_TYPE), debug)
				CFLAGS += -g -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -m64
			else
				CFLAGS += -O3 -m64 -std=c++0x
			endif	
				LDFLAGS +=  -L /usr/lib -lacml -lpathfortran 
		endif
	endif
endif	

ifdef CLUSTER 		#or: ifneq ($(HOSTNAME),)
	ifdef ACMLPATH
		# we should be on sharcnet: options set on sharcnet cluster for computers using the AMD processos
		CC      = pathCC
	  	LD      = pathCC
		CFLAGS += -I $(HOME)/CppLib/include 
		ifeq ($(BUILD_TYPE), debug)
			CFLAGS += -g -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
		else
			CFLAGS += -O3
		endif
		LDFLAGS += -L /opt/sharcnet/acml/current/pathscale64/lib -lacml -lpathfortran
		# should be but doesnt work:  LDFLAGS += -L $(ACMLPATH) -lacml -lpathfortran 
	endif
endif

