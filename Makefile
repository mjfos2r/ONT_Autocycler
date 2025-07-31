IMAGE_NAME = mjf_workflow
VERSION := $(shell cat .VERSION)
TAG1 = mjfos2r/$(IMAGE_NAME):$(VERSION)
TAG2 = mjfos2r/$(IMAGE_NAME):latest
CONTAINERS_FILE = containers.txt
CONTAINERS = $(shell cat $(CONTAINERS_FILE) 2>/dev/null | grep -v "^$$" | tr -d ' \t\r')

all: | check tag
build-containers: refresh-modules img-build-all
refresh-modules: clean-containers-file generate-containers-list

check:
	find . -name '.venv' -prune -o -name '.git' -prune -o -regex  '.*/*.wdl' -print0 | xargs -0 miniwdl check
	find . -name '.venv' -prune -o -name '.git' -prune -o -regex  '.*\.\(ya?ml\)' -print0 | xargs -0 yamllint -d relaxed

tag:
	git tag -s v$(VERSION) -m "Workflow version $(VERSION)"
	git push origin tag v$(VERSION)

generate-containers-list:
	find . -maxdepth 1 -type d -name "[a-z]*" -not -name ".git" -printf "%f\n" |\
	       	sort |\
	       	uniq > $(CONTAINERS_FILE)
	@echo "Generated list of containers to build in $(CONTAINERS_FILE)"

img-build-all:
	for dir in $(CONTAINERS); do \
		echo begin $$dir; \
		(cd $$dir; ${MAKE}); \
		echo end $$dir; \
		echo ==========; \
	done

clean-containers-file:
	rm -f $(CONTAINERS_FILE)

