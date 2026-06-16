#!/usr/bin/env bash

# builds the image multi-dwpc:latest and runs it on port 8501

IMAGE_NAME=multi-dwpc:latest

docker build -t ${IMAGE_NAME} . && \
docker run --rm -it \
	-p 8501:8501 \
	-v ./data/:/app/data/ \
	${IMAGE_NAME}
