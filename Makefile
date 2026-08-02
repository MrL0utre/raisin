# SPDX-License-Identifier: GPL-3.0-only

CC ?= cc

CPPFLAGS ?= -Iinclude
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
LDLIBS ?= -lm

SRC_DIR := src
BUILD_DIR := build
BIN_DIR := exe
TARGET := $(BIN_DIR)/raisin
TEST_TARGET := $(BUILD_DIR)/test_invariants
PYTHON ?= python3

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SOURCES))
CORE_OBJECTS := $(filter-out $(BUILD_DIR)/raisin.o,$(OBJECTS))
DEPS := $(OBJECTS:.o=.d)

.PHONY: all raisin check debug sanitize clean

all: raisin

raisin: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_invariants.o: test/test_invariants.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(TEST_TARGET): $(BUILD_DIR)/test_invariants.o $(CORE_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

check: raisin $(TEST_TARGET)
	$(TEST_TARGET) hypergraphs/b14.rzn2
	$(PYTHON) test/test_cli.py $(TARGET) .

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

debug: CFLAGS := -O0 -g3 -std=c11 -Wall -Wextra -Wpedantic -Werror
debug: clean raisin

sanitize: CFLAGS := -O1 -g3 -std=c11 -Wall -Wextra -Wpedantic -Werror \
	-fsanitize=address,undefined -fno-omit-frame-pointer
sanitize: LDFLAGS := -fsanitize=address,undefined
sanitize: clean raisin

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(DEPS)
