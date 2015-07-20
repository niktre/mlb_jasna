NAME = run

CLEANNAME = clean

CC = gcc
LFLAGS += -lm

HEADERS = mlb.h d2q21.h defs.h derivFD.h

SOURCE = mlb.c derivFD.c d2q21.c fbgk6.c general.c derivatives.c

OBJECT = mlb.o derivFD.o d2q21.o fbgk6.o general.o derivatives.o

FILES = output/*

all: $(NAME)
		echo All done

$(CLEANNAME):
		$(RM) $(OBJECT) 
		$(RM) -f $(FILES)

$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)

$(OBJECT): $(HEADERS)
