CC = g++
CFLAGS = -std=c++17 -g -O0
INCLUDE = -IDDEINT 

BUILD_DIR = build
PERF_DIR = perf_output
PLOT_DIR = plots

# Create the necessary directories
.PHONY: directories
directories:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(PERF_DIR)
	mkdir -p $(PLOT_DIR)

.PHONY: all BCTest clean

all: directories BCTest

# Compile the Breast Cancer model test and run performance test
BCTest: $(BUILD_DIR)/bc_test

$(BUILD_DIR)/bc_test: directories test/Breast_Cancer_Model/bc_model.cpp test/Breast_Cancer_Model/bc_functions.cpp
	@echo "Compiling Breast Cancer Model Test"
	$(CC) $(CFLAGS) $(INCLUDE) test/Breast_Cancer_Model/bc_model.cpp test/Breast_Cancer_Model/bc_functions.cpp -o $(BUILD_DIR)/bc_test
	@echo "Running Breast Cancer Model Test"
	@perf record -g ./$(BUILD_DIR)/bc_test
	@if [ -f $(PERF_DIR)/perf.data ]; then \
		echo "Generating Flame Graph"; \
		perf script > $(PERF_DIR)/out.perf; \
		./utils/FlameGraph/stackcollapse-perf.pl $(PERF_DIR)/out.perf > $(PERF_DIR)/out.folded; \
		./utils/FlameGraph/flamegraph.pl --color=java --hash --title="Breast Cancer Model" $(PERF_DIR)/out.folded > $(PLOT_DIR)/bc_model.svg; \
		echo "Cleaning up"; \
		rm -f $(PERF_DIR)/*; \
	else \
		echo "Error: perf.data file not found"; \
		exit 1; \
	fi

clean:
	rm -rf $(BUILD_DIR) $(PERF_DIR) $(PLOT_DIR)

