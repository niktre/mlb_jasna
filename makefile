NAME = run

CLEANNAME = clean

CC = gcc
LFLAGS = -lm

SOURCE = d2q21.c fbgk6.c general.c derivatives.c

OBJECT = d2q21.o fbgk6.o general.o derivatives.o

FILES = output/*

all: $(NAME)
		echo All done

$(CLEANNAME):
		$(RM) $(OBJECT) 
		$(RM) -f $(FILES)

$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)

$(OBJECT): defs.h
