NAME = run

CC = gcc
LFLAGS = -lm

SOURCE = fbgk6.c general.c derivatives.c

OBJECT = fbgk6.o general.o derivatives.o


all: $(NAME)
		echo All done

clean:
		$(RM) $(OBJECT)
		$(RM) -f output/*
                        
$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)
                        
$(OBJECT): defs.h