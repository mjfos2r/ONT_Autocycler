IMAGE_NAME = mjf_workflow
VERSION := $(shell cat .VERSION)

all: | check tag

check:
	find . -name '.venv' -prune -o -name '.git' -prune -o -regex  '.*/*.wdl' -print0 | xargs -0 miniwdl check
	find . -name '.venv' -prune -o -name '.git' -prune -o -regex  '.*\.\(ya?ml\)' -print0 | xargs -0 yamllint -d relaxed

tag:
	git tag -s v$(VERSION) -m "Workflow version $(VERSION)"
	git push origin tag v$(VERSION)

