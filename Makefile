CURL = "curl"
CURL_FLAGS = "--remote-name"

WGET = "wget"
WGET_FLAGS = ""

TAR = "tar"

BLOSSOMV_URL = "http://pub.ist.ac.at/~vnk/software/blossom5-v2.04.src.tar.gz"
BLOSSOMV_DIR = `basename $(BLOSSOMV_URL) .tar.gz`
BLOSSOMV_FILE = `basename $(BLOSSOMV_URL)`

all: blossomv libautotune polyestimate ex process

install: download all
	@echo "\n#\n# Please read and understand LICENSE.md and blossomv/LICENSE.TXT before continuing.\n#"

autotune: libautotune ex

.PHONY: process
process:
	make -C process all

.PHONY: blossomv
blossomv:
	make -C blossomv

.PHONY: libautotune
libautotune:
	make -C libautotune all

.PHONY: ex
ex:
	make -C ex all

.PHONY: polyestimate
polyestimate:
	make -C tools/polyestimate all

.PHONY: download
download:
	@printf "Downloading BlossomV from pub.ist.ac.at: $(BLOSSOMV_URL)\n";
	@if [ ! -f "$(BLOSSOMV_FILE)" ]; \
	then \
		if [ -x "`which $(CURL) 2> /dev/null`" ]; \
		then \
			$(CURL) $(CURL_FLAGS) "$(BLOSSOMV_URL)"; \
		else \
			if [ -x "`which $(WGET) 2> /dev/null`" ]; \
			then \
				$(WGET) $(WGET_FLAGS) "$(BLOSSOMV_URL)"; \
			else \
				printf "ERROR: Could not download Blossom V. Please install 'curl' or 'wget' to continue, or download to this directory manually.\n"; \
				exit 2; \
			fi \
		fi; \
		if [ ! -f "$(BLOSSOMV_FILE)" ]; \
		then \
			printf "ERROR: Failed to download BlossomV. Please check your proxy settings and that write permissions to this folder are available.\n"; \
			exit 2; \
		fi \
	fi
	@printf "\nUnzipping BlossomV\n";
	@if [ ! -x "`which $(TAR) 2> /dev/null`" ]; \
	then \
		printf "Could not find $(TAR) in PATH. Please install tar.\n"; \
		exit 2; \
	fi;
	@mkdir -p blossomv
	@if [ ! -d "blossomv" ]; \
	then \
		printf "Insufficient privileges to create directory 'blossomv' in this directory.\n"; \
		exit 2; \
	fi
	@$(TAR) -xzf $(BLOSSOMV_FILE) -C blossomv --strip-components 1
	@rm $(BLOSSOMV_FILE)
	@printf "\nDownloading/Unzipping BlossomV Complete.\n"

.PHONY: virtualenv
virtualenv:
	make -C process virtualenv

.PHONY: virtualenv_install
virtualenv_install:
	make -C process virtualenv_install

.PHONY: clean
clean:
	make -C libautotune clean
	make -C ex clean
