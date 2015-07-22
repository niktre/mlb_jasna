NAME = run

CLEANNAME = clean

CC = gcc
LFLAGS += -lm

HEADERS = mlb.h d2q21.h derivFD.h

SOURCE = mlb.c derivFD.c d2q21.c

OBJECT = mlb.o derivFD.o d2q21.o 

FILES = output/*

all: $(NAME)
		echo All done

$(CLEANNAME):
		$(RM) $(OBJECT) 
		$(RM) -f $(FILES)

$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)

$(OBJECT): $(HEADERS)
