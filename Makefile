CC ?= cc

CPPFLAGS ?= -Iinclude
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic
LDLIBS ?= -lm

SRC_DIR := src
BUILD_DIR := build
BIN_DIR := exe
TARGET := $(BIN_DIR)/raisin

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS := $(OBJECTS:.o=.d)

.PHONY: all raisin debug sanitize clean

all: raisin

raisin: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

debug: CFLAGS := -O0 -g3 -std=c11 -Wall -Wextra -Wpedantic
debug: clean raisin

sanitize: CFLAGS := -O1 -g3 -std=c11 -Wall -Wextra -Wpedantic \
	-fsanitize=address,undefined -fno-omit-frame-pointer
sanitize: LDFLAGS := -fsanitize=address,undefined
sanitize: clean raisin

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(DEPS)
