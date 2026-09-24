# PFv2 build, test and image targets.
#
# The namespace defaults to the Docker Hub account the images are published under. Override
# it for a private registry:
#
#   make image REGISTRY=ghcr.io/dysplasiadx NAMESPACE=
#
NAMESPACE ?= conidiobolus
REGISTRY  ?=
IMAGE     ?= pfv2
VERSION   ?= $(shell grep -m1 '^readonly VERSION=' PFv2.sh | cut -d'"' -f2)
PLATFORMS ?= linux/amd64,linux/arm64

ifeq ($(strip $(REGISTRY)),)
TAG_BASE := $(NAMESPACE)/$(IMAGE)
else
TAG_BASE := $(REGISTRY)/$(NAMESPACE)/$(IMAGE)
endif

.PHONY: help build test smoke e2e install uninstall image image-local push clean

help:
	@printf 'PFv2 %s\n\n' '$(VERSION)'
	@printf '  make build        compile, run the test suite and package PFv2.jar\n'
	@printf '  make test         run the unit test suite\n'
	@printf '  make smoke        run the end-to-end smoke test\n'
	@printf '  make e2e          full run against the committed fixture (needs STAR and Bowtie2)\n'
	@printf '  make install      symlink the ptesfinder CLI into PREFIX/bin (default /usr/local)\n'
	@printf '  make uninstall    remove those symlinks\n'
	@printf '  make image-local  build the image for this machine only, and self test it\n'
	@printf '  make image        build and push a multi-arch image (%s)\n' '$(PLATFORMS)'
	@printf '  make push         alias for image\n'
	@printf '  make clean        remove build output\n\n'
	@printf '  image tags: %s:%s and %s:latest\n' '$(TAG_BASE)' '$(VERSION)' '$(TAG_BASE)'

build:
	bash setup.sh

test:
	bash setup.sh

smoke: build
	bash test/smoke/run.sh

# Runs the real aligners against the committed fixture, so it needs them on PATH. Unlike
# `smoke`, this exercises candidate nomination and all four Bowtie2 passes.
e2e: build
	bash test/e2e/run.sh

PREFIX ?= /usr/local

# A symlink rather than a copy, so the CLI keeps resolving PFV2_HOME back to this checkout
# and an edit here takes effect immediately.
install: build
	install -d $(PREFIX)/bin
	ln -sf $(CURDIR)/bin/ptesfinder $(PREFIX)/bin/ptesfinder
	ln -sf ptesfinder $(PREFIX)/bin/pfv2
	@printf '>>> installed %s/bin/{ptesfinder,pfv2} -> %s/bin/ptesfinder\n' '$(PREFIX)' '$(CURDIR)'
	@$(PREFIX)/bin/ptesfinder version

uninstall:
	rm -f $(PREFIX)/bin/ptesfinder $(PREFIX)/bin/pfv2

# Single-architecture build that stays in the local daemon, for iterating on the image.
image-local:
	docker build -t $(TAG_BASE):$(VERSION) -t $(TAG_BASE):latest .
	docker run --rm $(TAG_BASE):$(VERSION) version
	docker run --rm $(TAG_BASE):$(VERSION) selftest

# Multi-arch images cannot be loaded into the local daemon, so this pushes. Requires
# `docker login` and a buildx builder: docker buildx create --use
image:
	docker buildx build --platform $(PLATFORMS) \
	  -t $(TAG_BASE):$(VERSION) -t $(TAG_BASE):latest \
	  --push .

push: image

clean:
	rm -rf build
