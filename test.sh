#!/bin/sh
docker run -it \
	-v `pwd`/test:/data \
	seglh/dmlpa:latest \
	/data/test.xlsx \
	/data/test.pdf
