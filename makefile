NAME = run

CLEANNAME = clean

CC = gcc
LFLAGS = -lm

SOURCE = fbgk6.c general.c derivatives.c

OBJECT = fbgk6.o general.o derivatives.o

FILES = output/*

all: $(NAME)
		echo All done

$(CLEANNAME):
		$(RM) $(OBJECT) 
		$(RM) -f $(FILES)
						
$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)
                        
$(OBJECT): defs.h