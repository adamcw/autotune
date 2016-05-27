CURL = "curl"
CURL_FLAGS = "--remote-name"

WGET = "wget"
WGET_FLAGS = ""

TAR = "tar"

all: libautotune polyestimate ex process

install: all

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
