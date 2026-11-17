#!/usr/bin/env bash

# runs the compose stack, which launches the API and the React frontend instead
# of the streamlit app
docker compose build && docker compose up --remove-orphans
