CC      = gcc
CFLAGS  = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRCDIR  = src
OBJDIR  = obj
TESTDIR = tests

SRCS    = $(SRCDIR)/conservation.c $(SRCDIR)/topology.c $(SRCDIR)/spectral.c $(SRCDIR)/cst_unified.c
OBJS    = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

LIB     = libcstcore.a
TEST_BIN = test_cst

.PHONY: all clean test

all: $(LIB)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	ar rcs $@ $^

$(TEST_BIN): $(TESTDIR)/test_cst.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lcstcore $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -rf $(OBJDIR) $(LIB) $(TEST_BIN)
