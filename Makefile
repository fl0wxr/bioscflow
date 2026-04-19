SHELL := /bin/bash

.PHONY: noop clean reset

noop:
	mkdir -p tmp data/raw data/interim data/cluster libexec vendor

clean:
	rm -rf tmp/*

reset: clean
	rm -f data/interim/*.casta
	rm -f data/interim/*.adj
	rm -f data/interim/*.mci
	rm -f data/interim/*.tab
	rm -f data/cluster/*.cl
	rm -f data/cluster/session_report*.json
