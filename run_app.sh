#!/usr/bin/env bash

IMAGE_NAME=multi-dwpc-streamlit:latest

docker build -t "${IMAGE_NAME}" . && \
docker run --rm -it --name multi-dwpc-streamlit-app \
	-p 8501:8501 \
	-v ./data/:/app/data/ \
	"${IMAGE_NAME}" \
	${@}
