CXX ?= g++

PYTHON = "python"
REQUIRED_PYTHON = "2.7"

PIP = "pip"

TAR = "tar"

CURL = "curl"
CURL_FLAGS = "--remote-name"

WGET = "wget"
WGET_FLAGS = ""

REQUIREMENTS = "requirements.txt"

VIRTUALENV = "virtualenv"
VIRTUALENV_NAME = "AutotuneVE"
VIRTUALENV_URL = "https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.tar.gz"
VIRTUALENV_DIR = `basename $(VIRTUALENV_URL) .tar.gz`
VIRTUALENV_FILE = `basename $(VIRTUALENV_URL)`

all: convert pythondeps

install: all

convert:
	@printf "Compiling required C packages\n";
	$(CXX) convertruns.c -o convertruns;

pythondeps:
	@printf "\nInstalling required Python packages\n";
	@if [ ! -x "`which $(PYTHON) 2> /dev/null`" ]; \
	then \
		printf "ERROR: $(PYTHON) cannot be found in PATH. Please install Python v2.7 to continue.\n"; \
		exit 2; \
	fi
	@if [ `$(PYTHON) -c 'import sys; print "%d.%d" % (sys.version_info[0], sys.version_info[1])'` != $(REQUIRED_PYTHON) ]; \
	then \
		printf "ERROR: Python was detected, but is not version $(REQUIRED_PYTHON). Please install Python version $(REQUIRED_PYTHON).\n"; \
		exit 2; \
	fi
	@if [ ! -x "`which $(PIP) 2> /dev/null`" ]; \
	then \
		printf "ERROR: $(PIP) cannot be found in PATH. Please install pip or run 'make virtualenv' to create a virtual environment to install Python packages into.\n"; \
		exit 2; \
	fi
	@$(PIP) install -r $(REQUIREMENTS);
	@printf "\nThe installation of required Python packages was successful.\n";

virtualenv:
	@if [ ! -x "`which $(VIRTUALENV) 2> /dev/null`" ]; \
	then \
		printf "ERROR: $(VIRTUALENV) cannot be found in PATH. Please run 'make virtualenv_install' download and install it to the current directory.\n"; \
		exit 2; \
	fi
	@printf "Attempting to create a new virtualenv (virtual environment) at: AutotuneVE\n";
	@$(VIRTUALENV) $(VIRTUALENV_NAME);
	@printf "\nAttempting to install dependencies using new virtualenv (virtual environment)\n";
	@source $(VIRTUALENV_NAME)/bin/activate; \
	pip install -r $(REQUIREMENTS);
	@printf "\nThe installation of required Python packages was successful.\n";
	@printf "\n*********************************\n";
	@printf " PLEASE READ: USAGE INSTRUCTIONS \n";
	@printf "*********************************\n";
	@printf "Type to update your PATH to use virtualenv (virtual environment):\n\n";
	@printf "\tcd ./process;\n";
	@printf "\tsource AutotuneVE/bin/activate;\n\n";
	@printf "To stop using this virtualenv.\n\n";
	@printf "\tdeactivate;\n";
	@if [ ! $(PYTHON_HOME) ]; \
	then \
		printf "\nYou can also use Python directly using 'AutotuneVE/bin/python' rather than the usual '$(PYTHON)'.\n"; \
	fi
	@printf "\nPlease see 'process/README.md' for more information.\n";

virtualenv_install:
	@printf "Downloading virtualenv from pypi.python.org\n";
	@if [ -x "`which $(CURL) 2> /dev/null`" ]; \
	then \
		curl $(CURL_FLAGS) "$(VIRTUALENV_URL)"; \
	else \
		if [ -x "`which $(WGET) 2> /dev/null`" ]; \
		then \
			 wget $(WGET_FLAGS) "$(VIRTUALENV_URL)"; \
		else \
			printf "Could not download virtualenv. Please install 'curl' or 'wget' to continue.\n"; \
			exit 2; \
		fi \
	fi
	@if [ ! -f "$(VIRTUALENV_FILE)" ]; \
	then \
		printf "Unable to connect to Pypi. Failed to download virtualenv. Please check your proxy settings.\n"; \
		exit 2; \
	fi
	@printf "\nUnzipping virtualenv\n";
	@if [ ! -x "`which $(TAR) 2> /dev/null`" ]; \
	then \
		printf "Could not find $(TAR) in PATH. Please install tar.\n"; \
		exit 2; \
	fi
	@$(TAR) xvfz "$(VIRTUALENV_FILE)";
	@printf "\nAttempting to create a new virtualenv (virtual environment) at: AutotuneVE\n";
	@python "$(VIRTUALENV_DIR)/virtualenv.py" --setuptools $(VIRTUALENV_NAME);
	@printf "\nAttempting to install dependencies using new virtualenv (virtual environment)\n";
	@source $(VIRTUALENV_NAME)/bin/activate; \
	pip install -r $(REQUIREMENTS);
	@printf "\nThe installation of required Python packages was successful.\n";
	@printf "\n*********************************\n";
	@printf " PLEASE READ: USAGE INSTRUCTIONS \n";
	@printf "*********************************\n";
	@printf "Type to update your PATH to use virtualenv (virtual environment):\n\n";
	@printf "\tcd ./process;\n";
	@printf "\tsource AutotuneVE/bin/activate;\n\n";
	@printf "To stop using this virtualenv.\n\n";
	@printf "\tdeactivate;\n";
	@if [ ! $(PYTHON_HOME) ]; \
	then \
		printf "\nYou can also use Python directly using 'AutotuneVE/bin/python' rather than the usual python.\n"; \
	fi
	@printf "\nPlease see 'process/README.md' for more information.\n";

clean:
	-rm convertruns

.PHONY:
	clean
