all: process blossomv libautotune ex

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

.PHONY: download
download:
	wget http://pub.ist.ac.at/~vnk/software/blossom5-v2.04.src.tar.gz
	mkdir -p blossomv
	tar -xzvf blossom5-v2.04.src.tar.gz -C blossomv --strip-components 1
	rm blossom5-v2.04.src.tar.gz
