SHELL := /bin/bash
CC = g++
CFLAGS = -std=c++17 -IDDEINT -Icpp_utils
PYTHON = python3
PIP = .venv/bin/pip
PYTHON_RUN = .venv/bin/python
VENV = .venv
CPP_TARGET = test/run_cpp
JL_SCRIPT = test/main.jl
PY_SCRIPT = test/main.py

.PHONY: all run run_cpp run_jl setup_py run_py clean help

all: setup_py 

run: run_cpp run_jl run_py

run_cpp: $(CPP_TARGET)

$(CPP_TARGET): test/main.cpp
	@echo "Compiling C++ code..."
	@$(CC) $(CFLAGS) test/main.cpp -o $(CPP_TARGET)
	@echo "Running C++ code..."
	@./$(CPP_TARGET)

run_jl: $(JL_SCRIPT)
	@echo "Running Julia code..."
	@julia $(JL_SCRIPT)

setup_py:
	@if [ ! -d "$(VENV)" ]; then \
		echo "Setting up Python environment..."; \
		$(PYTHON) -m venv $(VENV); \
		$(PIP) install pandas matplotlib; \
	else \
		echo "Python environment already set up."; \
	fi

run_py: $(PY_SCRIPT) setup_py
	@echo "Running Python code..."
	@$(PYTHON_RUN) $(PY_SCRIPT)

clean:
	@echo "Cleaning up..."
	@rm -rf $(VENV) $(CPP_TARGET)
	@rm -f test/run_cpp